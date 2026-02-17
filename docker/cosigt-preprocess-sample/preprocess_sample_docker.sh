#!/usr/bin/env bash

set -euo pipefail

################################################################################
# COSIGT Pipeline - Preprocessing Sample (Docker/native version)
################################################################################

if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <SAMPLE> <CONFIG> <THREADS>"
    exit 1
fi

SAMPLE="$1"
CONFIG="$2"
THREADS="$3"

# Simple YAML Parser
get_yaml_value() {
    local key="$1"
    grep "^${key}:" "$CONFIG" | sed "s/^${key}:[[:space:]]*//" | tr -d "'" | tr -d '"'
}

# Get alignment path from TSV file
get_sample_alignment() {
    local sample="$1"
    local tsv_file="$2"
    
    if [ ! -f "$tsv_file" ]; then
        echo "ERROR: Alignment map file not found: $tsv_file" >&2
        return 1
    fi
    
    local alignment=$(awk -v s="$sample" '$2 == s {print $1}' "$tsv_file")
    if [ -z "$alignment" ]; then
        alignment=$(awk -v s="$sample" '$1 == s {print $2}' "$tsv_file")
    fi
    
    echo "$alignment"
}

OUTPUT_DIR=$(get_yaml_value "output")
REFERENCE=$(get_yaml_value "reference")
ALIGNMENT_MAP=$(get_yaml_value "alignment_map")

echo "==================================================================="
echo "COSIGT Pipeline - Preprocessing Sample: $SAMPLE"
echo "==================================================================="

# Get alignment path from TSV
SAMPLE_ALIGNMENT=$(get_sample_alignment "$SAMPLE" "$ALIGNMENT_MAP")

if [ -z "$SAMPLE_ALIGNMENT" ]; then
    echo "ERROR: No alignment found for sample: $SAMPLE"
    echo "Please check $ALIGNMENT_MAP contains an entry for $SAMPLE"
    exit 1
fi

echo "Sample: $SAMPLE"
echo "Alignment: $SAMPLE_ALIGNMENT"
echo "Reference: $REFERENCE"
echo "Threads: $THREADS"
echo "==================================================================="

# Check if alignment file exists (via FUSE mount)
if [ ! -f "$SAMPLE_ALIGNMENT" ]; then
    echo "ERROR: Alignment file not found: $SAMPLE_ALIGNMENT"
    exit 1
fi

# Check for index file
if [[ "$SAMPLE_ALIGNMENT" == *.cram ]]; then
    INDEX_FILE="${SAMPLE_ALIGNMENT}.crai"
elif [[ "$SAMPLE_ALIGNMENT" == *.bam ]]; then
    if [ -f "${SAMPLE_ALIGNMENT}.bai" ]; then
        INDEX_FILE="${SAMPLE_ALIGNMENT}.bai"
    else
        INDEX_FILE="${SAMPLE_ALIGNMENT}.csi"
    fi
fi

if [ -n "${INDEX_FILE:-}" ] && [ ! -f "$INDEX_FILE" ]; then
    echo "WARNING: Index file not found: $INDEX_FILE"
    echo "Processing may be slow without an index"
fi

# Extract unmapped reads
UNMAPPED_DIR="${OUTPUT_DIR}/samtools/fasta/${SAMPLE}"
mkdir -p "$UNMAPPED_DIR"

UNMAPPED_FASTA="${UNMAPPED_DIR}/unmapped.fasta.gz"

if [ -f "$UNMAPPED_FASTA" ]; then
    echo "✓ Unmapped reads already exist: $UNMAPPED_FASTA"
    echo "Skipping extraction..."
else
    echo "Extracting unmapped reads..."
    samtools view -u -f 4 -@ "$THREADS" -T "$REFERENCE" "$SAMPLE_ALIGNMENT" | \
        samtools sort -n -@ "$THREADS" -T "${UNMAPPED_DIR}/unmapped_tmp" - | \
        samtools fasta -0 /dev/null -@ "$THREADS" - | gzip > "$UNMAPPED_FASTA"
    
    echo "✓ Unmapped reads extracted successfully!"
fi

echo "==================================================================="
echo "✓ Sample preprocessing complete"
echo "Output: $UNMAPPED_FASTA"
echo "==================================================================="

