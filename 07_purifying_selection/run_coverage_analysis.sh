#!/bin/bash
#
# Purifying Selection Analysis - Coverage Constraint
# 
# Run this script in your alphagenome-env:
#   conda activate alphagenome-env
#   cd /mnt/work_1/gest9386/CU_Boulder/rotations/LAYER/dormant_site_activation_pipeline/07_purifying_selection
#   bash run_coverage_analysis.sh
#
# What it does:
#   1. Reformats the gnomAD coverage file for tabix (one-time, ~10 min)
#   2. Creates tabix index (one-time, ~5 min)
#   3. Queries coverage for all H=1, H=2, H=3 dormant site positions (~1 min)
#   4. Computes constraint statistics and generates figures
#
# Hardware: Optimized for 32 cores, 188GB RAM
#

set -e  # Exit on error

# Paths
COVERAGE_DIR="/mnt/data_1/gnomAD_data/raw/gnomad_v4.1/coverage"
WORK_DIR="/mnt/work_1/gest9386/CU_Boulder/rotations/LAYER/dormant_site_activation_pipeline"
ORIGINAL_FILE="$COVERAGE_DIR/gnomad.genomes.v4.1.allele_number_all_sites.tsv.bgz"
REFORMATTED_FILE="$COVERAGE_DIR/gnomad_AN_tabix.tsv.bgz"
INDEX_FILE="$REFORMATTED_FILE.tbi"

PATHS_FILE="$WORK_DIR/results/mutation_paths/AP1/paths.tsv"
OUTPUT_DIR="$WORK_DIR/results/purifying_selection/AP1"
FIGURE_DIR="$WORK_DIR/figures/purifying_selection"

THREADS=$(nproc)  # Use all available cores (32 on this machine)

echo "============================================================"
echo "Purifying Selection Analysis - Coverage Constraint"
echo "============================================================"
echo ""

# Check we're in the right conda env
if [[ -z "$CONDA_PREFIX" ]] || [[ ! "$CONDA_PREFIX" == *"alphagenome"* ]]; then
    echo "ERROR: Please activate alphagenome-env first:"
    echo "  conda activate alphagenome-env"
    exit 1
fi

# Check required tools
for tool in bgzip tabix awk python; do
    if ! command -v $tool &> /dev/null; then
        echo "ERROR: $tool not found"
        exit 1
    fi
done
echo "✓ All required tools found"

# Create output directories
mkdir -p "$OUTPUT_DIR" "$FIGURE_DIR"

# ============================================================
# Step 1: Reformat coverage file for tabix (if needed)
# ============================================================
if [[ -f "$REFORMATTED_FILE" ]] && [[ -f "$INDEX_FILE" ]]; then
    echo ""
    echo "✓ Reformatted file and index already exist, skipping Step 1-2"
else
    echo ""
    echo "Step 1: Reformatting coverage file for tabix..."
    echo "  Input:  $ORIGINAL_FILE (12 GB)"
    echo "  Output: $REFORMATTED_FILE"
    echo "  This will take ~10-15 minutes..."
    echo ""
    
    # Remove partial file if exists
    rm -f "$REFORMATTED_FILE"
    
    # Reformat: split chr:pos into separate columns
    # Original: chr1:10001\t16
    # New:      chr1\t10001\t16
    # Using bgzip -d for BGZF-aware parallel decompression (better than pigz for .bgz files)
    time bgzip -dc -@ $THREADS "$ORIGINAL_FILE" | \
        awk -F'\t' 'NR==1 {print "#chr\tpos\tAN"; next} {split($1,a,":"); print a[1]"\t"a[2]"\t"$2}' | \
        bgzip -@ $THREADS > "$REFORMATTED_FILE"
    
    echo ""
    echo "✓ Reformatting complete"
    ls -lh "$REFORMATTED_FILE"
    
    # ============================================================
    # Step 2: Create tabix index
    # ============================================================
    echo ""
    echo "Step 2: Creating tabix index..."
    echo "  This will take ~5 minutes..."
    echo ""
    
    time tabix -s 1 -b 2 -e 2 -S 1 "$REFORMATTED_FILE"
    
    echo ""
    echo "✓ Index created"
    ls -lh "$INDEX_FILE"
fi

# ============================================================
# Step 3: Extract positions and query coverage
# ============================================================
echo ""
echo "Step 3: Querying coverage for dormant site positions..."

# Run the Python script
python3 "$(dirname "$0")/query_coverage.py"

echo ""
echo "============================================================"
echo "Analysis Complete!"
echo "============================================================"
echo ""
echo "Results saved to:"
echo "  $OUTPUT_DIR/constraint_by_hamming.tsv"
echo "  $OUTPUT_DIR/purifying_selection_summary.txt"
echo ""
echo "To generate figures, run:"
echo "  python plot_constraint_evidence.py --input $OUTPUT_DIR --output $FIGURE_DIR"
