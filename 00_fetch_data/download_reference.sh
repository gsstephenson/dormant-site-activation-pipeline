#!/bin/bash
# Download GRCh38 reference genome
# Part of Dormant Site Activation Pipeline - Module 00

set -euo pipefail

# Get script directory and pipeline root
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PIPELINE_ROOT="$(dirname "$SCRIPT_DIR")"

# Configuration
REFERENCE_DIR="${PIPELINE_ROOT}/data/reference"
ENSEMBL_URL="ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna"
GENOME_FILE="Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
OUTPUT_FA="${REFERENCE_DIR}/GRCh38.fa"

echo "============================================"
echo "Downloading GRCh38 Reference Genome"
echo "============================================"

# Create directory
mkdir -p "${REFERENCE_DIR}"

# Download if not already present
if [ -f "${OUTPUT_FA}" ]; then
    echo "✓ Reference genome already exists: ${OUTPUT_FA}"
else
    echo "Downloading from Ensembl..."
    wget -O "${REFERENCE_DIR}/${GENOME_FILE}" "${ENSEMBL_URL}/${GENOME_FILE}"
    
    echo "Decompressing..."
    gunzip "${REFERENCE_DIR}/${GENOME_FILE}"
    
    echo "✓ Reference genome downloaded: ${OUTPUT_FA}"
fi

# Index with samtools
if command -v samtools &> /dev/null; then
    if [ -f "${OUTPUT_FA}.fai" ]; then
        echo "✓ FASTA index already exists"
    else
        echo "Creating FASTA index..."
        samtools faidx "${OUTPUT_FA}"
        echo "✓ FASTA index created: ${OUTPUT_FA}.fai"
    fi
else
    echo "⚠ samtools not found - skipping FASTA indexing"
    echo "  Install samtools and run: samtools faidx ${OUTPUT_FA}"
fi

# Print summary
echo ""
echo "============================================"
echo "Download Complete"
echo "============================================"
echo "Reference genome: ${OUTPUT_FA}"

if [ -f "${OUTPUT_FA}.fai" ]; then
    n_chromosomes=$(wc -l < "${OUTPUT_FA}.fai")
    echo "Number of sequences: ${n_chromosomes}"
fi

echo ""
echo "You can now proceed with motif scanning."
