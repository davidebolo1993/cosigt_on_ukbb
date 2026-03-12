#!/usr/bin/env bash

set -euo pipefail

################################################################################
# COSIGT Pipeline - Region Processing (Docker/native version)
################################################################################

# CHANGED: added optional [EXTRA_BED] argument
if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <REGION> <CONFIG> <THREADS> [EXTRA_BED]"
    echo ""
    echo "Example: $0 chr1_100000_200000 config.yaml 12"
    echo "         $0 chr1_100000_200000 config.yaml 12 /path/to/extra.bed"
    exit 1
fi

REGION="$1"
CONFIG="$2"
THREADS="$3"
EXTRA_BED="${4:-}"   # CHANGED: optional 4th argument, empty string if not provided

################################################################################
# Simple YAML Parser
################################################################################

get_yaml_value() {
    local key="$1"
    grep "^${key}:" "$CONFIG" | sed "s/^${key}:[[:space:]]*//" | tr -d "'" | tr -d '"'
}

get_yaml_list() {
    local key="$1"
    awk "/^${key}:/{flag=1;next}/^[a-zA-Z]/{flag=0}flag" "$CONFIG" | \
        sed 's/^[[:space:]]*-[[:space:]]*//' | tr -d "'" | tr -d '"' | grep -v '^[[:space:]]*$'
}

# Get graph path from TSV file
get_graph_path() {
    local region="$1"
    local tsv_file="$2"

    if [ ! -f "$tsv_file" ]; then
        echo "ERROR: Graph map file not found: $tsv_file" >&2
        return 1
    fi

    local graph_path=$(awk -v r="$region" '$1 == r {print $2}' "$tsv_file")
    if [ -z "$graph_path" ]; then
        graph_path=$(awk -v r="$region" '$2 == r {print $1}' "$tsv_file")
    fi

    echo "$graph_path"
}

# Get alignment path from TSV file
get_sample_alignment() {
    local sample="$1"
    local tsv_file="$2"

    if [ ! -f "$tsv_file" ]; then
        echo "ERROR: Alignment map file not found: $tsv_file" >&2
        return 1
    fi

    local alignment=$(awk -F'\t' -v s="$sample" '$2 == s {print $1}' "$tsv_file")
    if [ -z "$alignment" ]; then
        alignment=$(awk -F'\t' -v s="$sample" '$1 == s {print $2}' "$tsv_file")
    fi

    echo "$alignment"
}

################################################################################
# Load Configuration
################################################################################

OUTPUT_DIR=$(get_yaml_value "output")
INPUT_DIR=$(get_yaml_value "input")
REFERENCE=$(get_yaml_value "reference")
ALIGNMENT_MAP=$(get_yaml_value "alignment_map")
GRAPH_MAP=$(get_yaml_value "graph_map")

# Load samples
SAMPLES=()
while IFS= read -r line; do
    [ -n "$line" ] && SAMPLES+=("$line")
done < <(get_yaml_list "samples")

# Extract chromosome from region
CHROM=$(echo "$REGION" | rev | cut -d'_' -f3- | rev)

echo "==================================================================="
echo "COSIGT Pipeline - Region: $REGION (Chromosome: $CHROM)"
echo "==================================================================="
echo "Samples: ${#SAMPLES[@]} | Threads: $THREADS"
echo "Graph map: $GRAPH_MAP"
echo "Alignment map: $ALIGNMENT_MAP"
# CHANGED: print EXTRA_BED info if provided
if [ -n "$EXTRA_BED" ]; then
    echo "Extra BED: $EXTRA_BED"
fi
echo "==================================================================="
echo ""

################################################################################
# Validate Prerequisites
################################################################################

MERYL_REF_DB="${INPUT_DIR}/meryl/reference"
if [ ! -d "$MERYL_REF_DB" ] || [ ! -f "${MERYL_REF_DB}/merylIndex" ]; then
    echo "ERROR: Reference k-mer database not found!"
    echo "  Expected: $MERYL_REF_DB"
    echo ""
    echo "Please run preprocess_reference first"
    exit 1
fi

################################################################################
# Step 1: Extract Assemblies from Graph
################################################################################

echo "[Step 1] Processing graph and extracting assemblies"

ALLELES_DIR="${OUTPUT_DIR}/alleles/${CHROM}/${REGION}"
mkdir -p "$ALLELES_DIR"

INPUT_GRAPH=$(get_graph_path "$REGION" "$GRAPH_MAP")

if [ -z "$INPUT_GRAPH" ]; then
    echo "ERROR: No graph found for region: $REGION"
    echo "Please check $GRAPH_MAP contains an entry for $REGION"
    exit 1
fi

echo "  Graph: $INPUT_GRAPH"

if [ ! -f "$INPUT_GRAPH" ]; then
    echo "ERROR: Graph file not found: $INPUT_GRAPH"
    exit 1
fi

ALLELES_FASTA="${ALLELES_DIR}/${REGION}.fasta.gz"

if [ ! -f "${INPUT_DIR}/alleles/${CHROM}/${REGION}/${REGION}.fasta.gz" ]; then
    echo "  Extracting assemblies from graph..."
    odgi paths -i "$INPUT_GRAPH" -f | bgzip -c > "$ALLELES_FASTA"
    samtools faidx "$ALLELES_FASTA"
    echo "  Building BWA-MEM2 index.."
    bwa-mem2 index "$ALLELES_FASTA"
else
    ALLELES_DIR="${INPUT_DIR}/alleles/${CHROM}/${REGION}"
    ALLELES_FASTA="${ALLELES_DIR}/${REGION}.fasta.gz"
fi

################################################################################
# Step 2: ODGI Graph Utilities
################################################################################

echo "[Step 2] ODGI utilities and graph analysis"

ODGI_DIR="${OUTPUT_DIR}/odgi"
mkdir -p "${ODGI_DIR}/view/${CHROM}/${REGION}"
mkdir -p "${ODGI_DIR}/paths/${CHROM}/${REGION}"

GFA_FILE="${ODGI_DIR}/view/${CHROM}/${REGION}/${REGION}.gfa.gz"
PATHS_FILE="${ODGI_DIR}/paths/${CHROM}/${REGION}/${REGION}.tsv.gz"
LENGTH_FILE="${ODGI_DIR}/view/${CHROM}/${REGION}/${REGION}.node.length.tsv"

if [ ! -f "${INPUT_DIR}/odgi/view/${CHROM}/${REGION}/${REGION}.gfa.gz" ]; then
    echo "  Converting to GFA and extracting paths..."
    odgi view -i "$INPUT_GRAPH" -g --threads "$THREADS" | gzip > "$GFA_FILE"
    zgrep '^S' "$GFA_FILE" | awk '{print "node."$2"\t"length($3)}' > "$LENGTH_FILE"
    odgi paths -i "$INPUT_GRAPH" -H | cut -f 1,4- | gzip > "$PATHS_FILE"
else
    ODGI_DIR="${INPUT_DIR}/odgi"
    GFA_FILE="${ODGI_DIR}/view/${CHROM}/${REGION}/${REGION}.gfa.gz"
    PATHS_FILE="${ODGI_DIR}/paths/${CHROM}/${REGION}/${REGION}.tsv.gz"
    LENGTH_FILE="${ODGI_DIR}/view/${CHROM}/${REGION}/${REGION}.node.length.tsv"
fi

################################################################################
# Step 3: Prepare BED Files
################################################################################
# CHANGED: REF_BED is always just the region coordinates.
# ALIGN_BED is built from EXTRA_BED column 5 (seqname:start-end entries) if a
# matching line is found (col1==CHROM, col2==REGION_START, col3==REGION_END).
# If no match or no EXTRA_BED is provided, ALIGN_BED is a copy of REF_BED.

echo "[Step 3] Preparing BED files"

BEDTOOLS_DIR="${OUTPUT_DIR}/bedtools"
mkdir -p "${BEDTOOLS_DIR}/reference_bed/${CHROM}/${REGION}"
mkdir -p "${BEDTOOLS_DIR}/alignment_bed/${CHROM}/${REGION}"

REF_BED="${BEDTOOLS_DIR}/reference_bed/${CHROM}/${REGION}/${REGION}.bed.gz"
ALIGN_BED="${BEDTOOLS_DIR}/alignment_bed/${CHROM}/${REGION}/${REGION}.bed.gz"

REGION_START=$(echo "$REGION" | rev | cut -d'_' -f2 | rev)
REGION_END=$(echo "$REGION" | rev | cut -d'_' -f1 | rev)

if [ ! -f "${INPUT_DIR}/bedtools/alignment_bed/${CHROM}/${REGION}/${REGION}.bed.gz" ]; then

    echo "  Creating reference BED from region coordinates..."
    echo -e "${CHROM}\t${REGION_START}\t${REGION_END}" | gzip > "$REF_BED"

    if [ -n "$EXTRA_BED" ] && [ -f "$EXTRA_BED" ]; then
        EXTRA_ENTRIES=$(awk \
            -v c="$CHROM" -v s="$REGION_START" -v e="$REGION_END" \
            '$1 == c && $2 == s && $3 == e { print $5 }' \
            "$EXTRA_BED")

        if [ -n "$EXTRA_ENTRIES" ]; then
            echo "  Found extra alignment regions in $EXTRA_BED, extending ALIGN_BED..."
            {
                # include the reference
                echo -e "${CHROM}\t${REGION_START}\t${REGION_END}"

                echo "$EXTRA_ENTRIES" | tr ',' '\n' | while IFS= read -r entry; do
                    [ -z "$entry" ] && continue
                    seq=$(echo "$entry" | rev | cut -d: -f2- | rev)
                    coords=$(echo "$entry" | rev | cut -d: -f1 | rev)
                    start=$(echo "$coords" | cut -d- -f1)
                    end=$(echo "$coords" | cut -d- -f2)
                    echo -e "${seq}\t${start}\t${end}"
                done
            } | gzip > "$ALIGN_BED"
        else
            echo "  No matching entry in $EXTRA_BED for $REGION, ALIGN_BED = REF_BED"
            cp "$REF_BED" "$ALIGN_BED"
        fi
    else
        cp "$REF_BED" "$ALIGN_BED"
    fi

else
    BEDTOOLS_DIR="${INPUT_DIR}/bedtools"
    REF_BED="${BEDTOOLS_DIR}/reference_bed/${CHROM}/${REGION}/${REGION}.bed.gz"
    ALIGN_BED="${BEDTOOLS_DIR}/alignment_bed/${CHROM}/${REGION}/${REGION}.bed.gz"
fi

################################################################################
# Step 4a: Compute Allele-Specific K-mers (Meryl)
################################################################################

echo "[Step 4a] Computing allele-specific k-mers"

MERYL_DIR="${OUTPUT_DIR}/meryl/${CHROM}/${REGION}"
mkdir -p "$MERYL_DIR"

UNIQUE_KMERS="${MERYL_DIR}/${REGION}.unique_kmers.txt"

if [ ! -f "${INPUT_DIR}/meryl/${CHROM}/${REGION}/${REGION}.unique_kmers.txt" ]; then

    MERYL_ALLELES="${MERYL_DIR}/${REGION}"
    MERYL_DIFF="${MERYL_DIR}/${REGION}_unique"

    echo "  Computing allele-specific k-mers..."
    MEM_GB=$((THREADS < 8 ? THREADS * 2 : 16))

    meryl count k=31 threads="$THREADS" memory="$MEM_GB" "$ALLELES_FASTA" output "$MERYL_ALLELES"
    meryl difference "$MERYL_ALLELES" "$MERYL_REF_DB" output "$MERYL_DIFF"
    meryl print "$MERYL_DIFF" > "$UNIQUE_KMERS"

    rm -rf "$MERYL_ALLELES" "$MERYL_DIFF"
else
    MERYL_DIR="${INPUT_DIR}/meryl/${CHROM}/${REGION}"
    UNIQUE_KMERS="${MERYL_DIR}/${REGION}.unique_kmers.txt"
fi

################################################################################
# Step 4b: Build K-mer Filter Index (kfilt)
################################################################################

echo "[Step 4b] Building k-mer filter index"

KFILT_DIR="${OUTPUT_DIR}/kfilt/index/${CHROM}/${REGION}"
mkdir -p "$KFILT_DIR"

KFILT_IDX="${KFILT_DIR}/${REGION}.kfilt.idx"

if [ ! -f "${INPUT_DIR}/kfilt/index/${CHROM}/${REGION}/${REGION}.kfilt.idx" ]; then
    kfilt build -k "$UNIQUE_KMERS" -K 31 -o "$KFILT_IDX"
else
    KFILT_DIR="${INPUT_DIR}/kfilt/index/${CHROM}/${REGION}"
    KFILT_IDX="${KFILT_DIR}/${REGION}.kfilt.idx"
fi
    
################################################################################
# Step 5: Node Filtering
################################################################################

echo "[Step 5] Filtering low-complexity nodes"

PANPLEXITY_DIR="${OUTPUT_DIR}/panplexity/${CHROM}/${REGION}"
mkdir -p "$PANPLEXITY_DIR"

PANPLEXITY_MASK="${PANPLEXITY_DIR}/${REGION}.mask.tsv"

if [ ! -f "${INPUT_DIR}/panplexity/${CHROM}/${REGION}/${REGION}.mask.tsv" ]; then
    panplexity \
        --input-gfa "$GFA_FILE" \
        -t auto -k 16 -w 100 -d 100 \
        --complexity linguistic \
        -m "$PANPLEXITY_MASK" \
        --threads "$THREADS"
else
    PANPLEXITY_DIR="${INPUT_DIR}/panplexity/${CHROM}/${REGION}"
    PANPLEXITY_MASK="${PANPLEXITY_DIR}/${REGION}.mask.tsv"
fi
    
COMBINED_MASK="${OUTPUT_DIR}/odgi/paths/${CHROM}/${REGION}/${REGION}.mask.tsv"

if [ ! -f "${INPUT_DIR}/odgi/paths/${CHROM}/${REGION}/${REGION}.mask.tsv" ]; then
    echo "  Filtering coverage outliers..."
    Rscript /usr/local/bin/coverage_outliers.r \
        "$PATHS_FILE" \
        "${ODGI_DIR}/paths/${CHROM}/${REGION}/${REGION}" \
        "$LENGTH_FILE" \
        "$PANPLEXITY_MASK"
else
    COMBINED_MASK="${INPUT_DIR}/odgi/paths/${CHROM}/${REGION}/${REGION}.mask.tsv"
fi

################################################################################
# Step 6: Clustering
################################################################################

echo "[Step 6] Computing path dissimilarities and clustering"

DISSIM_DIR="${OUTPUT_DIR}/odgi/dissimilarity/${CHROM}/${REGION}"
mkdir -p "$DISSIM_DIR"

DISSIM_FILE="${DISSIM_DIR}/${REGION}.tsv.gz"

if [ ! -f "${INPUT_DIR}/odgi/dissimilarity/${CHROM}/${REGION}/${REGION}.tsv.gz" ]; then
    odgi similarity -i "$INPUT_GRAPH" --all --distances --threads "$THREADS" | gzip > "$DISSIM_FILE"
else
    DISSIM_DIR="${INPUT_DIR}/odgi/dissimilarity/${CHROM}/${REGION}"
    DISSIM_FILE="${DISSIM_DIR}/${REGION}.tsv.gz"
fi
    
CLUSTER_DIR="${OUTPUT_DIR}/cluster/${CHROM}/${REGION}"
mkdir -p "$CLUSTER_DIR"

CLUSTER_JSON="${CLUSTER_DIR}/${REGION}.clusters.json"

if [ ! -f "${INPUT_DIR}/cluster/${CHROM}/${REGION}/${REGION}.clusters.json" ]; then
    Rscript /usr/local/bin/cluster.r \
        "$DISSIM_FILE" "$CLUSTER_JSON" "automatic" "100.0" "1"
else
    CLUSTER_JSON="${INPUT_DIR}/cluster/${CHROM}/${REGION}/${REGION}.clusters.json"
fi

################################################################################
# Step 7: Process Samples
################################################################################

echo "[Step 7] Processing ${#SAMPLES[@]} samples"
echo ""

process_sample() {
    local SAMPLE="$1"
    echo "  [Sample: $SAMPLE]"

    SAMPLE_ALIGNMENT=$(get_sample_alignment "$SAMPLE" "$ALIGNMENT_MAP")

    if [ -z "$SAMPLE_ALIGNMENT" ]; then
        echo "    WARNING: No alignment found for $SAMPLE in $ALIGNMENT_MAP, skipping"
        return
    fi

    echo "    Alignment: $SAMPLE_ALIGNMENT"

    if [ ! -f "$SAMPLE_ALIGNMENT" ]; then
        echo "    ERROR: Alignment file not found: $SAMPLE_ALIGNMENT"
        return 1
    fi

    # 7.1: Extract mapped reads
    SAMTOOLS_DIR="${OUTPUT_DIR}/samtools/fasta/${SAMPLE}/${CHROM}/${REGION}"
    mkdir -p "$SAMTOOLS_DIR"

    MAPPED_FASTA="${SAMTOOLS_DIR}/${REGION}.mapped.fasta.gz"

    if [ ! -f "$MAPPED_FASTA" ]; then
        echo "    Extracting mapped reads..."
        samtools view -T "$REFERENCE" -@ "$THREADS" -L "$ALIGN_BED" -M -b "$SAMPLE_ALIGNMENT" | \
            samtools sort -n -@ "$THREADS" -T "${SAMTOOLS_DIR}/${REGION}" - | \
            samtools fasta -@ "$THREADS" - | gzip > "$MAPPED_FASTA"
    fi

    # 7.2: Use pre-extracted unmapped reads
    UNMAPPED_FASTA="${INPUT_DIR}/samtools/fasta/${SAMPLE}/unmapped.fasta.gz"

    if [ ! -f "$UNMAPPED_FASTA" ]; then
        echo "    ERROR: Unmapped reads not found for $SAMPLE"
        echo "    Expected: $UNMAPPED_FASTA"
        echo "    Please run preprocess_sample for $SAMPLE first"
        return 1
    fi

    # 7.3: Filter unmapped reads
    KFILT_SAMPLE_DIR="${OUTPUT_DIR}/kfilt/${SAMPLE}/${CHROM}/${REGION}"
    mkdir -p "$KFILT_SAMPLE_DIR"

    FILTERED_UNMAPPED="${KFILT_SAMPLE_DIR}/${REGION}.unmapped.fasta.gz"

    if [ ! -f "$FILTERED_UNMAPPED" ]; then
        echo "    Filtering unmapped reads..."
        kfilt filter \
            -I "$UNMAPPED_FASTA" -o "$FILTERED_UNMAPPED" \
            -f fasta -z -i "$KFILT_IDX" -n 1 -m 0 -t "$THREADS"
    fi

    # 7.4: Combine reads
    COMBINE_DIR="${OUTPUT_DIR}/combine/${SAMPLE}/${CHROM}/${REGION}"
    mkdir -p "$COMBINE_DIR"

    COMBINED_FASTA="${COMBINE_DIR}/${REGION}.fasta.gz"

    if [ ! -f "$COMBINED_FASTA" ]; then
        cat "$MAPPED_FASTA" "$FILTERED_UNMAPPED" > "$COMBINED_FASTA"
    fi

    # 7.5: Realign and project
    BWA_DIR="${OUTPUT_DIR}/bwa-mem2/${SAMPLE}/${CHROM}/${REGION}"
    mkdir -p "$BWA_DIR"

    REALIGNED_SAM="${BWA_DIR}/${REGION}.realigned.sam"

    GFAINJECT_DIR="${OUTPUT_DIR}/gfainject/${SAMPLE}/${CHROM}/${REGION}"
    mkdir -p "$GFAINJECT_DIR"

    GAF_FILE="${GFAINJECT_DIR}/${REGION}.gaf.gz"

    if [ ! -f "$GAF_FILE" ]; then
        echo "    Realigning reads..."
        bwa-mem2 mem -t "$THREADS" -p -h 10000 -o "$REALIGNED_SAM" "$ALLELES_FASTA" "$COMBINED_FASTA"

        echo "    Projecting to graph..."
        gfainject --gfa "$GFA_FILE" --sam "$REALIGNED_SAM" --alt-hits 10000 | gzip > "$GAF_FILE"
        rm "$REALIGNED_SAM"
    fi

    # 7.6: Calculate coverage
    GAFPACK_DIR="${OUTPUT_DIR}/gafpack/${SAMPLE}/${CHROM}/${REGION}"
    mkdir -p "$GAFPACK_DIR"

    COVERAGE_FILE="${GAFPACK_DIR}/${REGION}.gafpack.gz"

    if [ ! -f "$COVERAGE_FILE" ]; then
        echo "    Computing node coverage..."
        gafpack \
            --gfa "$GFA_FILE" --gaf "$GAF_FILE" \
            --len-scale --weight-queries | gzip > "$COVERAGE_FILE"
    fi

    # 7.7: Genotype
    COSIGT_DIR="${OUTPUT_DIR}/cosigt/${SAMPLE}/${CHROM}/${REGION}"
    mkdir -p "$COSIGT_DIR"

    GENOTYPE_FILE="${COSIGT_DIR}/${REGION}.cosigt_genotype.tsv"

    if [ ! -f "$GENOTYPE_FILE" ]; then
        echo "    Genotyping..."
        cosigt \
            -p "$PATHS_FILE" \
            -g "$COVERAGE_FILE" \
            -c "$CLUSTER_JSON" \
            -o "$COSIGT_DIR" \
            -i "$SAMPLE" \
            -m "$COMBINED_MASK"
    fi

    echo "    Cleaning up intermediaries..."
    rm -rf \
        "${OUTPUT_DIR}/samtools/fasta/${SAMPLE}/${CHROM}/${REGION}" \
        "${OUTPUT_DIR}/kfilt/${SAMPLE}/${CHROM}/${REGION}" \
        "${OUTPUT_DIR}/combine/${SAMPLE}/${CHROM}/${REGION}" \
        "${OUTPUT_DIR}/bwa-mem2/${SAMPLE}/${CHROM}/${REGION}" \
        "${OUTPUT_DIR}/gfainject/${SAMPLE}/${CHROM}/${REGION}" \
        "${OUTPUT_DIR}/gafpack/${SAMPLE}/${CHROM}/${REGION}"

    echo "    ✓ Completed"
}

TOTAL_SAMPLES=${#SAMPLES[@]}
PROCESSED=0

for SAMPLE in "${SAMPLES[@]}"; do
    process_sample "$SAMPLE"
    PROCESSED=$((PROCESSED + 1))
    PCT=$(awk "BEGIN { printf \"%.1f\", ($PROCESSED / $TOTAL_SAMPLES) * 100 }")
    echo " Progress: ${PROCESSED}/${TOTAL_SAMPLES} samples (${PCT}%)"
    echo ""
done

# Small clean-up because at the moment regenerating those file is 

rm ${OUTPUT_DIR}/kfilt/index/${CHROM}/${REGION}/${REGION}.kfilt.idx
rm ${OUTPUT_DIR}/cluster/${CHROM}/${REGION}/${REGION}.clusters.json
rm ${OUTPUT_DIR}/odgi/paths/${CHROM}/${REGION}/${REGION}.mask.tsv

################################################################################
# Completion
################################################################################

echo ""
echo "==================================================================="
echo "✓ Region $REGION completed successfully!"
echo "Results: ${OUTPUT_DIR}/cosigt/*/${CHROM}/${REGION}/"
echo "==================================================================="
