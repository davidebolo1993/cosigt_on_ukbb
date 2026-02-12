#!/usr/bin/env bash

################################################################################
# COSIGT Pipeline - Region-Parallelized (Starting from .og graphs)
# 
# Usage: ./process_region.sh <region> <config_file> <threads> [cache_dir] [tmp_dir]
#
################################################################################

if [ "$#" -lt 3 ] || [ "$#" -gt 5 ]; then
    echo "Usage: $0 <region> <config_file> <threads> [cache_dir] [tmp_dir]"
    echo ""
    echo "Arguments:"
    echo "  region        Region to process (e.g., chr1_100000_200000)"
    echo "  config_file   Path to config.yaml"
    echo "  threads       Number of threads to use"
    echo "  cache_dir     Singularity cache directory [default: \$PWD/singularity_cache]"
    echo "  tmp_dir       Singularity tmp directory [default: /tmp]"
    echo ""
    echo "Example: $0 chr1_100000_200000 config/config.yaml 12 ./cache /scratch/tmp"
    exit 1
fi

REGION="$1"
CONFIG="$2"
THREADS="$3"
CACHE_DIR="${4:-$PWD/singularity_cache}"
TMP_DIR="${5:-/tmp}"

################################################################################
# Setup Singularity Environment
################################################################################

# Create and set cache directory
CACHE_DIR=$(readlink -f "$CACHE_DIR" 2>/dev/null || echo "$CACHE_DIR")
mkdir -p "$CACHE_DIR"
export SINGULARITY_CACHEDIR="$CACHE_DIR"

# Set tmp directory
TMP_DIR=$(readlink -f "$TMP_DIR" 2>/dev/null || echo "$TMP_DIR")
if [ ! -d "$TMP_DIR" ]; then
    echo "ERROR: tmp directory does not exist: $TMP_DIR"
    exit 1
fi
if [ ! -w "$TMP_DIR" ]; then
    echo "ERROR: tmp directory is not writable: $TMP_DIR"
    exit 1
fi
export SINGULARITY_TMPDIR="$TMP_DIR"

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

# Load configuration
OUTPUT_DIR=$(get_yaml_value "output")
REFERENCE=$(get_yaml_value "reference")
PANSN_PREFIX=$(get_yaml_value "pansn_prefix")

# Load samples
SAMPLES=()
while IFS= read -r line; do
    [ -n "$line" ] && SAMPLES+=("$line")
done < <(get_yaml_list "samples")

# Extract chromosome from region (format: chr_start_end)
CHROM=$(echo "$REGION" | rev | cut -d'_' -f3- | rev)

echo "==================================================================="
echo "COSIGT Pipeline - Region: $REGION (Chromosome: $CHROM)"
echo "==================================================================="
echo "Samples: ${#SAMPLES[@]} | Threads: $THREADS"
echo "Singularity cache: $CACHE_DIR"
echo "Singularity tmp: $TMP_DIR"
echo "==================================================================="
echo ""

################################################################################
# Setup Bind Paths for Singularity
################################################################################

# Collect all directories that need to be bound
BIND_PATHS="$PWD"

# Add reference directory
REFERENCE_DIR=$(dirname "$(readlink -f "$REFERENCE")")
BIND_PATHS="${BIND_PATHS},${REFERENCE_DIR}"

# Add output directory if it's not under PWD
OUTPUT_DIR_REAL=$(readlink -f "$OUTPUT_DIR" 2>/dev/null || echo "$OUTPUT_DIR")
if [[ "$OUTPUT_DIR_REAL" != "$PWD"* ]]; then
    BIND_PATHS="${BIND_PATHS},${OUTPUT_DIR_REAL}"
fi

# Add tmp directory
BIND_PATHS="${BIND_PATHS},${TMP_DIR}"

# Add resources directory (for alignments, graphs, regions, etc.)
if [ -d "resources" ]; then
    RESOURCES_REAL=$(readlink -f "resources")
    BIND_PATHS="${BIND_PATHS},${RESOURCES_REAL}"
    
    # Also add any symlinked directories within resources
    if [ -d "resources/alignments" ]; then
        for item in resources/alignments/*; do
            if [ -L "$item" ]; then
                REAL_PATH=$(readlink -f "$item")
                REAL_DIR=$(dirname "$REAL_PATH")
                BIND_PATHS="${BIND_PATHS},${REAL_DIR}"
            fi
        done
    fi
    
    if [ -d "resources/assemblies" ]; then
        for item in resources/assemblies/*/*; do
            if [ -L "$item" ] || [ -f "$item" ]; then
                REAL_PATH=$(readlink -f "$item" 2>/dev/null || echo "$item")
                REAL_DIR=$(dirname "$REAL_PATH")
                BIND_PATHS="${BIND_PATHS},${REAL_DIR}"
            fi
        done
    fi
    
    if [ -d "resources/reference" ]; then
        for item in resources/reference/*; do
            if [ -L "$item" ]; then
                REAL_PATH=$(readlink -f "$item")
                REAL_DIR=$(dirname "$REAL_PATH")
                BIND_PATHS="${BIND_PATHS},${REAL_DIR}"
            fi
        done
    fi
fi

# Remove duplicates and trailing commas
BIND_PATHS=$(echo "$BIND_PATHS" | tr ',' '\n' | sort -u | tr '\n' ',' | sed 's/,$//')

echo "Singularity bind paths: $BIND_PATHS"
echo ""

################################################################################
# Validate Prerequisites
################################################################################

MERYL_REF_DB="${OUTPUT_DIR}/meryl/reference"

if [ ! -d "$MERYL_REF_DB" ] || [ ! -f "${MERYL_REF_DB}/merylIndex" ]; then
    echo "ERROR: Reference k-mer database not found!"
    echo "       Expected: $MERYL_REF_DB"
    echo ""
    echo "Please run preprocessing first:"
    echo "bash preprocess.sh $CONFIG $THREADS $CACHE_DIR $TMP_DIR"
    exit 1
fi

################################################################################
# Singularity Container Management
################################################################################

CONTAINER_DIR="${OUTPUT_DIR}/containers"
mkdir -p "$CONTAINER_DIR"

get_container() {
    local docker_uri="$1"
    local name=$(echo "$docker_uri" | sed 's|docker://||' | tr '/:' '_').sif
    local path="${CONTAINER_DIR}/${name}"
    
    if [ ! -f "$path" ]; then
        echo "  [Container] WARNING: Container not found: $docker_uri" >&2
        echo "              Expected: $path" >&2
        echo "              Downloading now (this may cause issues if run in parallel)..." >&2
        singularity pull "$path" "$docker_uri" || {
            echo "  [Container] ERROR: Failed to download container!" >&2
            echo "              Please run: bash src/download_containers.sh $CONFIG first" >&2
            exit 1
        }
    fi
    echo "$path"
}

run_container() {
    local container="$1"
    shift
    singularity exec --cleanenv --bind "$BIND_PATHS" "$container" "$@"
}

# Download all containers
echo "[Setup] Preparing containers..."
SAMTOOLS_C=$(get_container "docker://davidebolo1993/samtools:1.22")
BWA_C=$(get_container "docker://davidebolo1993/bwa-mem2:2.2.1")
BEDTOOLS_C=$(get_container "docker://davidebolo1993/bedtools:2.31.0")
KFILT_C=$(get_container "docker://davidebolo1993/kfilt:0.1.1")
ODGI_C=$(get_container "docker://pangenome/odgi:1753347183")
GFAINJECT_C=$(get_container "docker://davidebolo1993/gfainject:0.2.1")
GAFPACK_C=$(get_container "docker://davidebolo1993/gafpack:0.1.3")
PANPLEXITY_C=$(get_container "docker://davidebolo1993/panplexity:0.1.1")
COSIGT_C=$(get_container "docker://davidebolo1993/cosigt:0.1.7")
RENV_C=$(get_container "docker://davidebolo1993/renv:4.3.3")

################################################################################
# Step 1: Extract Assemblies from Graph
################################################################################

echo "[Step 1] Processing graph and extracting assemblies"

ALLELES_DIR="${OUTPUT_DIR}/alleles/${CHROM}/${REGION}"
mkdir -p "$ALLELES_DIR"

# Input: pre-built .og graph
INPUT_GRAPH="resources/assemblies/${CHROM}/${REGION}.og"

if [ ! -f "$INPUT_GRAPH" ]; then
    echo "ERROR: Graph not found: $INPUT_GRAPH"
    exit 1
fi

# Extract assemblies as FASTA using odgi paths
ALLELES_FASTA="${ALLELES_DIR}/${REGION}.fasta.gz"

if [ ! -f "$ALLELES_FASTA" ]; then
    echo "  Extracting assemblies from graph..."
    run_container "$ODGI_C" odgi paths -i "$INPUT_GRAPH" -f | \
        run_container "$SAMTOOLS_C" bgzip -c > "$ALLELES_FASTA"
    run_container "$SAMTOOLS_C" samtools faidx "$ALLELES_FASTA"
fi

# Build BWA-MEM2 index
echo "  Building BWA-MEM2 index..."
if [ ! -f "${ALLELES_FASTA}.bwt.2bit.64" ]; then
    run_container "$BWA_C" bwa-mem2 index "$ALLELES_FASTA"
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

if [ ! -f "$GFA_FILE" ]; then
    echo "  Converting to GFA and extracting paths..."
    run_container "$ODGI_C" odgi view -i "$INPUT_GRAPH" -g --threads "$THREADS" | gzip > "$GFA_FILE"
    zgrep '^S' "$GFA_FILE" | awk '{print "node."$2,length($3)}' OFS='\t' > "$LENGTH_FILE"
    run_container "$ODGI_C" odgi paths -i "$INPUT_GRAPH" -H | cut -f 1,4- | gzip > "$PATHS_FILE"
fi

################################################################################
# Step 3: Prepare BED Files
################################################################################

echo "[Step 3] Preparing BED files"

BEDTOOLS_DIR="${OUTPUT_DIR}/bedtools"
mkdir -p "${BEDTOOLS_DIR}/reference_bed/${CHROM}/${REGION}"
mkdir -p "${BEDTOOLS_DIR}/alignment_bed/${CHROM}/${REGION}"

REF_BED="${BEDTOOLS_DIR}/reference_bed/${CHROM}/${REGION}/${REGION}.bed.gz"
ALIGN_BED="${BEDTOOLS_DIR}/alignment_bed/${CHROM}/${REGION}/${REGION}.bed.gz"

INPUT_BED="resources/regions/${CHROM}/${REGION}.bed"

if [ ! -f "$REF_BED" ]; then
    grep -w "$CHROM" "$INPUT_BED" | gzip > "$REF_BED"
fi

if [ ! -f "$ALIGN_BED" ]; then
    # Create temporary files for bedtools
    TMP_SORTED=$(mktemp -p "$TMP_DIR")
    run_container "$BEDTOOLS_C" bedtools sort -i "$INPUT_BED" > "$TMP_SORTED"
    
    run_container "$BEDTOOLS_C" bedtools intersect -a "$REF_BED" -b "$TMP_SORTED" -nonamecheck -u | gzip > "$ALIGN_BED"
    run_container "$BEDTOOLS_C" bedtools intersect -a "$TMP_SORTED" -b "$REF_BED" -nonamecheck -v | gzip >> "$ALIGN_BED"
    
    rm "$TMP_SORTED"
fi

################################################################################
# Step 4: Build K-mer Filter Index (uses pre-built reference DB)
################################################################################

echo "[Step 4] Building k-mer filter index"

KFILT_DIR="${OUTPUT_DIR}/kfilt/index/${CHROM}/${REGION}"
mkdir -p "$KFILT_DIR"
KFILT_IDX="${KFILT_DIR}/${REGION}.kfilt.idx"

if [ ! -f "$KFILT_IDX" ]; then

    MERYL_ALLELES="${OUTPUT_DIR}/meryl/${CHROM}/${REGION}/${REGION}"
    mkdir -p "$MERYL_ALLELES"
    MERYL_DIFF="${OUTPUT_DIR}/meryl/${CHROM}/${REGION}/${REGION}_unique"
    UNIQUE_KMERS="${OUTPUT_DIR}/meryl/${CHROM}/${REGION}/${REGION}.unique_kmers.txt"
    
    echo "  Computing allele-specific k-mers (using pre-built reference DB)..."
    MEM_GB=$((${THREADS} < 8 ? ${THREADS} * 2 : 16))  # Cap at 16GB for this step
    
    run_container "$KFILT_C" meryl count k=31 threads="$THREADS" memory="$MEM_GB" "$ALLELES_FASTA" output "$MERYL_ALLELES"
    run_container "$KFILT_C" meryl difference "$MERYL_ALLELES" "$MERYL_REF_DB" output "$MERYL_DIFF"
    run_container "$KFILT_C" meryl print "$MERYL_DIFF" > "$UNIQUE_KMERS"
    run_container "$KFILT_C" kfilt build -k "$UNIQUE_KMERS" -K 31 -o "$KFILT_IDX"
    
    rm -rf "$MERYL_ALLELES" "$MERYL_DIFF"
fi

################################################################################
# Step 5: Node Filtering (Panplexity + Coverage)
################################################################################

echo "[Step 5] Filtering low-complexity nodes"

PANPLEXITY_DIR="${OUTPUT_DIR}/panplexity/${CHROM}/${REGION}"
mkdir -p "$PANPLEXITY_DIR"
PANPLEXITY_MASK="${PANPLEXITY_DIR}/${REGION}.mask.tsv"

if [ ! -f "$PANPLEXITY_MASK" ]; then
    run_container "$PANPLEXITY_C" panplexity \
        --input-gfa "$GFA_FILE" \
        -t auto -k 16 -w 100 -d 100 \
        --complexity linguistic \
        -m "$PANPLEXITY_MASK" \
        --threads "$THREADS"
fi

COMBINED_MASK="${ODGI_DIR}/paths/${CHROM}/${REGION}/${REGION}.mask.tsv"

if [ ! -f "$COMBINED_MASK" ]; then
    echo "  Filtering coverage outliers..."
    run_container "$RENV_C" Rscript src/coverage_outliers.r \
        "$PATHS_FILE" \
        "${ODGI_DIR}/paths/${CHROM}/${REGION}/${REGION}" \
        "$LENGTH_FILE" \
        "$PANPLEXITY_MASK"
fi

################################################################################
# Step 6: Clustering
################################################################################

echo "[Step 6] Computing path dissimilarities and clustering"

DISSIM_DIR="${ODGI_DIR}/dissimilarity/${CHROM}/${REGION}"
mkdir -p "$DISSIM_DIR"
DISSIM_FILE="${DISSIM_DIR}/${REGION}.tsv.gz"

if [ ! -f "$DISSIM_FILE" ]; then
    run_container "$ODGI_C" odgi similarity -i "$INPUT_GRAPH" --all --distances --threads "$THREADS" | gzip > "$DISSIM_FILE"
fi

CLUSTER_DIR="${OUTPUT_DIR}/cluster/${CHROM}/${REGION}"
mkdir -p "$CLUSTER_DIR"
CLUSTER_JSON="${CLUSTER_DIR}/${REGION}.clusters.json"

if [ ! -f "$CLUSTER_JSON" ]; then
    run_container "$RENV_C" Rscript src/cluster.r \
        "$DISSIM_FILE" "$CLUSTER_JSON" "automatic" "100.0" "1"
fi

################################################################################
# Step 7: Process Samples
################################################################################

echo "[Step 7] Processing ${#SAMPLES[@]} samples"

process_sample() {
    local SAMPLE="$1"
    
    echo "  [Sample: $SAMPLE]"
    
    SAMPLE_ALIGNMENT=$(find resources/alignments -name "${SAMPLE}.*am" 2>/dev/null | head -1)
    
    if [ -z "$SAMPLE_ALIGNMENT" ]; then
        echo "    WARNING: No alignment found for $SAMPLE, skipping"
        return
    fi
    
    # 7.1: Extract mapped reads
    SAMTOOLS_DIR="${OUTPUT_DIR}/samtools/fasta/${SAMPLE}/${CHROM}/${REGION}"
    mkdir -p "$SAMTOOLS_DIR"
    MAPPED_FASTA="${SAMTOOLS_DIR}/${REGION}.mapped.fasta.gz"
    
    if [ ! -f "$MAPPED_FASTA" ]; then
        echo "    Extracting mapped reads..."
        run_container "$SAMTOOLS_C" samtools view -T "$REFERENCE" -@ "$THREADS" -L "$ALIGN_BED" -M -b "$SAMPLE_ALIGNMENT" | \
        run_container "$SAMTOOLS_C" samtools sort -n -@ "$THREADS" -T "${SAMTOOLS_DIR}/${REGION}" - | \
        run_container "$SAMTOOLS_C" samtools fasta -@ "$THREADS" - | gzip > "$MAPPED_FASTA"
    fi
    
    # 7.2: Use pre-extracted unmapped reads
    UNMAPPED_FASTA="${OUTPUT_DIR}/samtools/fasta/${SAMPLE}/unmapped.fasta.gz"
    
    if [ ! -f "$UNMAPPED_FASTA" ]; then
        echo "    ERROR: Unmapped reads not found for $SAMPLE"
        echo "           Expected: $UNMAPPED_FASTA"
        echo "           Please run: bash src/preprocess_sample.sh $SAMPLE $CONFIG <threads>"
        return 1
    fi
    
    echo "    Using pre-extracted unmapped reads"
    
    # 7.3: Filter unmapped reads with k-mer index
    KFILT_SAMPLE_DIR="${OUTPUT_DIR}/kfilt/${SAMPLE}/${CHROM}/${REGION}"
    mkdir -p "$KFILT_SAMPLE_DIR"
    FILTERED_UNMAPPED="${KFILT_SAMPLE_DIR}/${REGION}.unmapped.fasta.gz"
    
    if [ ! -f "$FILTERED_UNMAPPED" ]; then
        echo "    Filtering unmapped reads..."
        run_container "$KFILT_C" kfilt filter \
            -I "$UNMAPPED_FASTA" -o "$FILTERED_UNMAPPED" \
            -f fasta -z -i "$KFILT_IDX" -n 1 -m 0 -t "$THREADS"
    fi
    
    # 7.4: Combine mapped and filtered unmapped
    COMBINE_DIR="${OUTPUT_DIR}/combine/${SAMPLE}/${CHROM}/${REGION}"
    mkdir -p "$COMBINE_DIR"
    COMBINED_FASTA="${COMBINE_DIR}/${REGION}.fasta.gz"
    
    if [ ! -f "$COMBINED_FASTA" ]; then
        cat "$MAPPED_FASTA" "$FILTERED_UNMAPPED" > "$COMBINED_FASTA"
    fi
    
    # 7.5: Realign to alleles and project to graph
    BWA_DIR="${OUTPUT_DIR}/bwa-mem2/${SAMPLE}/${CHROM}/${REGION}"
    mkdir -p "$BWA_DIR"
    REALIGNED_SAM="${BWA_DIR}/${REGION}.realigned.sam"
    
    GFAINJECT_DIR="${OUTPUT_DIR}/gfainject/${SAMPLE}/${CHROM}/${REGION}"
    mkdir -p "$GFAINJECT_DIR"
    GAF_FILE="${GFAINJECT_DIR}/${REGION}.gaf.gz"
    
    if [ ! -f "$GAF_FILE" ]; then
        echo "    Realigning reads..."
        run_container "$BWA_C" bwa-mem2 mem -t "$THREADS" -p -h 10000 -o "$REALIGNED_SAM" "$ALLELES_FASTA" "$COMBINED_FASTA"
        
        echo "    Projecting to graph..."
        run_container "$GFAINJECT_C" gfainject --gfa "$GFA_FILE" --sam "$REALIGNED_SAM" --alt-hits 10000 | gzip > "$GAF_FILE"
        
        #rm "$REALIGNED_SAM"
    fi
    
    # 7.6: Calculate coverage (gafpack)
    GAFPACK_DIR="${OUTPUT_DIR}/gafpack/${SAMPLE}/${CHROM}/${REGION}"
    mkdir -p "$GAFPACK_DIR"
    COVERAGE_FILE="${GAFPACK_DIR}/${REGION}.gafpack.gz"
    
    if [ ! -f "$COVERAGE_FILE" ]; then
        echo "    Computing node coverage..."
        run_container "$GAFPACK_C" gafpack \
            --gfa "$GFA_FILE" --gaf "$GAF_FILE" \
            --len-scale --weight-queries | gzip > "$COVERAGE_FILE"
    fi
    
    # 7.7: Genotype (COSIGT)
    COSIGT_DIR="${OUTPUT_DIR}/cosigt/${SAMPLE}/${CHROM}/${REGION}"
    mkdir -p "$COSIGT_DIR"
    GENOTYPE_FILE="${COSIGT_DIR}/${REGION}.cosigt_genotype.tsv"
    
    if [ ! -f "$GENOTYPE_FILE" ]; then
        echo "    Genotyping..."
        run_container "$COSIGT_C" cosigt \
            -p "$PATHS_FILE" \
            -g "$COVERAGE_FILE" \
            -c "$CLUSTER_JSON" \
            -o "$COSIGT_DIR" \
            -i "$SAMPLE" \
            -m "$COMBINED_MASK"
    fi
    
    echo "    ✓ Completed"
}

for SAMPLE in "${SAMPLES[@]}"; do
    process_sample "$SAMPLE"
done

################################################################################
# Completion
################################################################################

echo "==================================================================="
echo "✓ Region $REGION completed successfully!"
echo "Results: ${OUTPUT_DIR}/cosigt/*/${CHROM}/${REGION}/"
echo "==================================================================="

