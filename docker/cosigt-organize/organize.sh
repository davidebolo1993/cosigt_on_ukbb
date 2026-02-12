#!/bin/bash
set -e

# organize.sh - Set up COSIGT pipeline directory structure
# Usage: ./organize.sh -g graphs.tsv -r reads.tsv -b regions.bed -f reference.fa -o output_dir

# Parse arguments
while getopts "g:r:b:f:o:p:" opt; do
    case $opt in
        g) GRAPHS_TSV="$OPTARG" ;;
        r) READS_TSV="$OPTARG" ;;
        b) BED_FILE="$OPTARG" ;;
        f) REFERENCE="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        p) PANSN_PREFIX="$OPTARG" ;;
        *) echo "Usage: $0 -g graphs.tsv -r reads.tsv -b bed -f ref.fa -o outdir [-p pansn_prefix]" >&2
           exit 1 ;;
    esac
done

# Set defaults
PANSN_PREFIX="${PANSN_PREFIX:-grch38#1#}"

# Validate required arguments
if [[ -z "$GRAPHS_TSV" || -z "$READS_TSV" || -z "$BED_FILE" || -z "$REFERENCE" || -z "$OUTPUT_DIR" ]]; then
    echo "ERROR: Missing required arguments" >&2
    echo "Usage: $0 -g graphs.tsv -r reads.tsv -b bed -f ref.fa -o outdir [-p pansn_prefix]" >&2
    exit 1
fi

# Function to convert path to absolute, handling container /work prefix
realpath_smart() {
    local path="$1"
    local abs_path

    # First convert to absolute path
    abs_path=$(realpath "$path" 2>/dev/null || echo "$path")

    # If we're in a container (cwd=/work) and path starts with /work/,
    # the realpath is already correct - just return it
    if [[ "$abs_path" == /work/* ]]; then
        echo "$abs_path"
        return
    fi

    # Otherwise return the absolute path
    echo "$abs_path"
}

# Convert all paths to absolute
GRAPHS_TSV=$(realpath_smart "$GRAPHS_TSV")
READS_TSV=$(realpath_smart "$READS_TSV")
BED_FILE=$(realpath_smart "$BED_FILE")
REFERENCE=$(realpath_smart "$REFERENCE")
OUTPUT_DIR=$(realpath_smart "$OUTPUT_DIR")

echo "============================================================"
echo "COSIGT Pipeline - Organize Resources"
echo "============================================================"
echo "Graphs:     $GRAPHS_TSV"
echo "Reads:      $READS_TSV"
echo "Regions:    $BED_FILE"
echo "Reference:  $REFERENCE"
echo "Output:     $OUTPUT_DIR"
echo "PanSN:      $PANSN_PREFIX"
echo "============================================================"

# Create directory structure
RESOURCES="resources"
CONFIG="config"

if [[ -d "$RESOURCES" && -n "$(ls -A $RESOURCES 2>/dev/null)" ]]; then
    echo "ERROR: Directory '$RESOURCES' is not empty. Clean it first!" >&2
    exit 1
fi

mkdir -p "$RESOURCES"/{assemblies,alignments,regions,reference}
mkdir -p "$CONFIG"

echo ""
echo "============================================================"
echo "Processing graphs..."
echo "============================================================"

# Process graphs TSV: region<TAB>graph_path
REGIONS_LIST=""
while IFS=$'\t' read -r region graph_path; do
    [[ -z "$region" ]] && continue

    # Get absolute path
    graph_abs=$(realpath_smart "$graph_path")

    if [[ ! -f "$graph_abs" ]]; then
        echo "ERROR: Graph file not found: $graph_abs" >&2
        exit 1
    fi

    # Extract chromosome from region (e.g., chr1_103304997_103901127 -> chr1)
    chrom=$(echo "$region" | cut -d_ -f1)

    # Create chromosome subdirectory
    mkdir -p "$RESOURCES/assemblies/$chrom"

    # Create symlink
    ln -sf "$graph_abs" "$RESOURCES/assemblies/$chrom/${region}.og"

    echo "  Added: $region -> $graph_abs"
done < "$GRAPHS_TSV"

echo ""
echo "============================================================"
echo "Processing alignments..."
echo "============================================================"

# Process reads TSV: alignment_path<TAB>sample_id
SAMPLES_LIST=""
while IFS=$'\t' read -r aln_path sample_id; do
    [[ -z "$sample_id" ]] && continue

    # Get absolute path
    aln_abs=$(realpath_smart "$aln_path")

    if [[ ! -f "$aln_abs" ]]; then
        echo "ERROR: Alignment file not found: $aln_abs" >&2
        exit 1
    fi

    # Determine file type and index
    if [[ "$aln_abs" == *.bam ]]; then
        ext=".bam"
        if [[ -f "${aln_abs}.bai" ]]; then
            idx_ext=".bai"
        elif [[ -f "${aln_abs}.csi" ]]; then
            idx_ext=".csi"
        else
            echo "ERROR: BAM index not found for $aln_abs" >&2
            exit 1
        fi
    elif [[ "$aln_abs" == *.cram ]]; then
        ext=".cram"
        idx_ext=".crai"
        if [[ ! -f "${aln_abs}.crai" ]]; then
            echo "ERROR: CRAM index not found for $aln_abs" >&2
            exit 1
        fi
    else
        echo "ERROR: Unsupported alignment format: $aln_abs" >&2
        exit 1
    fi

    # Create symlinks
    ln -sf "$aln_abs" "$RESOURCES/alignments/${sample_id}${ext}"
    ln -sf "${aln_abs}${idx_ext}" "$RESOURCES/alignments/${sample_id}${ext}${idx_ext}"

    # Add to samples list
    SAMPLES_LIST="${SAMPLES_LIST}- ${sample_id}\n"

    echo "  Added: $sample_id -> $aln_abs"
done < "$READS_TSV"

echo ""
echo "============================================================"
echo "Processing regions..."
echo "============================================================"

# Process BED file
while IFS=$'\t' read -r chrom start end rest; do
    [[ -z "$chrom" || "$chrom" == "#"* ]] && continue

    region="${chrom}_${start}_${end}"

    # Create chromosome subdirectory
    mkdir -p "$RESOURCES/regions/$chrom"

    # Write BED entry
    echo -e "${chrom}\t${start}\t${end}\t${rest}" > "$RESOURCES/regions/$chrom/${region}.bed"

    # Add to regions list
    REGIONS_LIST="${REGIONS_LIST}- ${region}\n"

    echo "  Added: $region"
done < "$BED_FILE"

echo ""
echo "============================================================"
echo "Processing reference..."
echo "============================================================"

# Create reference symlinks
ref_name=$(basename "$REFERENCE")
ln -sf "$REFERENCE" "$RESOURCES/reference/$ref_name"

if [[ -f "${REFERENCE}.fai" ]]; then
    ln -sf "${REFERENCE}.fai" "$RESOURCES/reference/${ref_name}.fai"
else
    echo "WARNING: Reference index not found: ${REFERENCE}.fai" >&2
fi

echo "  Added: $ref_name -> $REFERENCE"

# Get absolute path to reference symlink
ref_link=$(realpath_smart "$RESOURCES/reference/$ref_name")

echo ""
echo "============================================================"
echo "Writing configuration..."
echo "============================================================"

# Write config.yaml
cat > "$CONFIG/config.yaml" << EOF
output: $OUTPUT_DIR
pansn_prefix: $PANSN_PREFIX
reference: $ref_link
regions:
$(echo -e "$REGIONS_LIST")
samples:
$(echo -e "$SAMPLES_LIST")
EOF

echo "  Config written to: $CONFIG/config.yaml"

echo ""
echo "============================================================"
echo "✓ Setup complete!"
echo "============================================================"
echo "Resources:  $(pwd)/$RESOURCES"
echo "Config:     $(pwd)/$CONFIG"
echo "============================================================"

