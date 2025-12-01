#!/bin/bash
#
# Run Disease Overlap Analysis for Dormant Site Activation Pipeline
#
# This script intersects AP1 activating variants with ClinVar and GWAS databases
#
# Usage: ./run_disease_overlap.sh [TF_NAME]
#

set -euo pipefail

# Configuration
TF_NAME="${1:-AP1}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(dirname "$SCRIPT_DIR")"

# Paths
AP1_LANDSCAPE="${BASE_DIR}/results/landscape/${TF_NAME}/${TF_NAME}_activation_landscape.tsv"
CLINVAR_VCF="${BASE_DIR}/data/clinvar/raw/clinvar.vcf.gz"
GWAS_CATALOG="${BASE_DIR}/data/gwas/raw/gwas_catalog_associations.tsv"
OUTPUT_DIR="${BASE_DIR}/results/disease_overlap/${TF_NAME}"

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "=============================================="
echo "Disease Overlap Analysis - ${TF_NAME}"
echo "=============================================="
echo ""
echo "Input files:"
echo "  - AP1 landscape: ${AP1_LANDSCAPE}"
echo "  - ClinVar VCF: ${CLINVAR_VCF}"
echo "  - GWAS Catalog: ${GWAS_CATALOG}"
echo ""
echo "Output directory: ${OUTPUT_DIR}"
echo ""

# Verify inputs exist
if [[ ! -f "$AP1_LANDSCAPE" ]]; then
    echo "ERROR: AP1 landscape file not found: $AP1_LANDSCAPE"
    exit 1
fi

if [[ ! -f "$CLINVAR_VCF" ]]; then
    echo "ERROR: ClinVar VCF not found: $CLINVAR_VCF"
    echo "Please ensure data symlinks are set up correctly"
    exit 1
fi

if [[ ! -f "$GWAS_CATALOG" ]]; then
    echo "ERROR: GWAS Catalog not found: $GWAS_CATALOG"
    echo "Please ensure data symlinks are set up correctly"
    exit 1
fi

# Run the analysis
echo "Starting disease overlap analysis..."
echo ""

python "${SCRIPT_DIR}/disease_overlap.py" \
    --ap1-landscape "$AP1_LANDSCAPE" \
    --clinvar-vcf "$CLINVAR_VCF" \
    --gwas-catalog "$GWAS_CATALOG" \
    --output-dir "$OUTPUT_DIR" \
    --gwas-window 1000

echo ""
echo "=============================================="
echo "Analysis complete!"
echo "=============================================="
echo ""
echo "Results saved to: ${OUTPUT_DIR}"
echo ""
echo "Key outputs:"
echo "  - disease_overlap_report.txt  : Summary statistics"
echo "  - clinvar_overlaps.tsv        : ClinVar pathogenic overlaps"
echo "  - gwas_overlaps.tsv           : GWAS associations"
echo "  - gwas_unique_variants.tsv    : Unique GWAS-linked variants"
echo "  - figures/                    : Visualization plots"
