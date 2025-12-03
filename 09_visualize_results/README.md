# Module 09: Visualization Suite

## Overview

Generate publication-quality figures summarizing the dormant site activation pipeline findings, with focus on the ~30,000 "forbidden variants" and purifying selection evidence.

## Figures Generated

### Core Figures (3)

#### 1. Selection Gradient (`plot_selection_gradient.py`)
Multi-panel figure showing purifying selection evidence:
- **Panel A**: Survival rate by Hamming distance with inset zoom for H=1 visibility
- **Panel B**: AP1 score distributions with all pairwise comparisons (H=1 vs H=2, H=2 vs H=3, H=1 vs H=3)
- **Panel C**: Fold depletion with log scale and p-values
- **Panel D**: Key findings in scannable bullet format

#### 2. Forbidden Variants Characterization (`plot_forbidden_characterization.py`)
Multi-panel figure characterizing forbidden variants:
- **Panel A**: Mutation spectrum with transitions/transversions grouped and proper y-axis (starts at 0)
- **Panel B**: Top TFs by frequency with family color coding
- **Panel C**: Top tissues/cell types with category color coding
- **Panel D**: AP1 vs Enhancer hexbin plot (fixes overplotting with 30K points)

#### 3. TF-Tissue Heatmap (`plot_tf_tissue_heatmap.py`)
Two-panel heatmap with filtering and clustering:
- **Panel A**: Median AP1 score by TF-tissue (cells with n≥10 only, hatching for low-n)
- **Panel B**: Variant counts by TF-tissue with log color scale
- Hierarchical clustering applied to rows and columns
- Marginal totals shown

### Additional Figures (3)

#### 4. Genomic Distribution (`plot_genomic_distribution.py`)
Multi-panel figure showing WHERE forbidden variants are located:
- **Panel A**: Chromosome distribution (observed vs expected by size)
- **Panel B**: Manhattan-style genomic distribution plot
- **Panel C**: AP1 score by chromosome (box plots)
- **Panel D**: AP1 quantile distribution with percentile annotations

#### 5. Selection Signal Strength (`plot_selection_signal.py`)
Multi-panel figure characterizing selection INTENSITY:
- **Panel A**: Proportion with extreme AP1 scores (>0.9, >0.95, etc.)
- **Panel B**: Allele number (AN) distribution showing coverage quality
- **Panel C**: AP1 raw score distribution (log scale)
- **Panel D**: Top 10 highest-impact candidates table

#### 6. Functional Context (`plot_functional_context.py`)
Multi-panel figure showing functional properties:
- **Panel A**: Enhancer activity (H3K27ac) distribution
- **Panel B**: Chromatin accessibility (ATAC/DNase) distribution
- **Panel C**: Multi-feature scatter (accessibility vs enhancer, colored by AP1)
- **Panel D**: High-priority multi-track candidates table

## Usage

```bash
# Activate environment
conda activate alphagenome-env
cd /path/to/dormant_site_activation_pipeline

# Generate core figures
python 09_visualize_results/plot_selection_gradient.py --results-dir results --output-dir figures/selection_gradient
python 09_visualize_results/plot_forbidden_characterization.py --results-dir results --output-dir figures/forbidden_variants
python 09_visualize_results/plot_tf_tissue_heatmap.py --results-dir results --output-dir figures/forbidden_variants

# Generate additional figures
python 09_visualize_results/plot_genomic_distribution.py --results-dir results --output-dir figures/forbidden_variants
python 09_visualize_results/plot_selection_signal.py --results-dir results --output-dir figures/forbidden_variants
python 09_visualize_results/plot_functional_context.py --results-dir results --output-dir figures/forbidden_variants
```

## Input Data

| File | Source | Description |
|------|--------|-------------|
| `results/purifying_selection/AP1/constraint_by_hamming.tsv` | Module 07 | Selection statistics |
| `results/forbidden_variants/AP1/predictions_summary_ap1.tsv` | Module 08 | AP1-specific scores |
| `results/forbidden_variants/AP1/forbidden_variants.tsv` | Module 08 | Variant metadata |
| `results/landscape/AP1/AP1_activation_landscape.tsv` | Module 05 | Observed variants |

## Output

Figures are saved to:
- `figures/selection_gradient/` — Selection gradient figure
- `figures/forbidden_variants/` — All other figures

Both PNG (300 DPI) and PDF versions are generated for each figure.

## Improvements Over Original Scripts

### High Priority Fixes
- ✅ Fixed Panel A negative y-axis in mutation spectrum (now starts at 0)
- ✅ Fixed survival rate scale (H=1 visible via inset zoom)
- ✅ Fixed violin plot y-axis truncation (auto-scaling)
- ✅ Fixed overplotting with hexbin (30K points now visible)

### Medium Priority Improvements
- ✅ Added statistical comparisons for all H=1 vs H=2 vs H=3 pairs
- ✅ Grouped tissues by category (Cancer, Immune, Stem, etc.)
- ✅ Filtered low-count cells in heatmap (n≥10)
- ✅ Added hierarchical clustering to heatmap

### New Figures
- ✅ Genomic distribution figure
- ✅ Selection signal strength figure
- ✅ Functional context figure

## Key Findings Visualized

1. **76.5× depletion of H=1 variants** — Single mutations that would complete AP1 sites are essentially forbidden
2. **Selection targets completion, not strength** — Forbidden variants have similar predicted impact to observed variants
3. **T↔A transversions dominate** — 58% of forbidden variants are T>A or A>T mutations
4. **Tissue specificity** — HepG2 and K562 show highest enrichment
5. **TF specificity** — FOS, ATF3, and ATF7 are most frequently the "best responders"
6. **Extreme scores** — ~47% of forbidden variants exceed 99.9th percentile of all variants

## Dependencies

- pandas
- numpy
- matplotlib
- seaborn
- scipy
