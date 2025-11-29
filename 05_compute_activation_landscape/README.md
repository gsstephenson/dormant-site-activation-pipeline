# Module 05: Compute Activation Landscape

## Overview

This module combines population genetics data (gnomAD, Module 03) with functional impact predictions (AlphaGenome, Module 04) to create a two-dimensional "activation landscape" that maps dormant AP1 site accessibility versus functional impact.

## Scientific Rationale

### The Activation Landscape Concept

The activation landscape answers the question: **Which dormant TF binding sites are both accessible through human variation AND functionally impactful?**

We construct a 2D coordinate system:
- **X-axis (Accessibility):** How easy is it to activate this site through existing human variants?
- **Y-axis (Impact):** How much would activation affect AP1-family TF binding?

### Biologically-Specific Y-Axis Design

**Critical Innovation:** Unlike generic approaches that use `max()` across all AlphaGenome tracks, we use **AP1-family-specific TF predictions** with **raw scores**:

```
Y = log₁₀(max(raw_score)) for TF ∈ {JUND, JUN, JUNB, FOS, FOSL1, FOSL2, ATF3, ATF2, ATF7, BATF, BATF2, MAFK, MAFF, MAFG}
```

**Why raw scores instead of quantile scores:**
- Quantile scores have ceiling effects (90% of variants are >0.9)
- Raw scores (range: 3-34,000) preserve full dynamic range of effect magnitudes
- Raw scores reveal the purifying selection signal (Spearman r=0.096, p<10⁻¹⁴)

**Why AP1-specific:**
- Using global max would conflate unrelated signals (CTCF, splice sites, etc.)
- AP1-specific scoring directly measures what we hypothesize changes: AP1 binding
- Strong correlation (r=0.579) with enhancer marks validates this approach

### X-Axis: Population Accessibility Score

```
X = -log₁₀(AF) × Hamming_distance
```

Where:
- **AF:** Allele frequency from gnomAD (higher AF = more accessible)
- **Hamming_distance:** Number of mutations needed (1, 2, or 3)

**Interpretation:**
- Low X = accessible (common variant, few mutations)
- High X = hard to access (rare variant, many mutations)

## Scripts

### `compute_activation_landscape.py`

Main script that:
1. Loads AlphaGenome predictions (parquet file, ~186M rows)
2. Computes AP1-specific impact scores using vectorized groupby operations
3. Computes secondary metrics (enhancer marks, accessibility)
4. Joins with gnomAD data for X-axis calculation
5. Creates quadrant classifications
6. Outputs landscape TSV files

**Usage:**
```bash
conda run -n alphagenome-env python compute_activation_landscape.py \
  --predictions ../results/alphagenome/AP1/predictions.parquet \
  --gnomad ../results/gnomad_intersection/AP1/all_observed_variants.tsv \
  --output ../results/landscape/AP1
```

**Performance:** ~30 seconds for 6,810 variants (vectorized implementation)

### `plot_activation_landscape.py`

Generates publication-quality figures:
1. **Main landscape:** 2D scatter with quadrant annotations
2. **Comparison panels:** AP1-specific vs global, enhancer validation
3. **TF breakdown:** Which AP1-family TFs drive the signal
4. **High-priority detail:** Zoom on accessible + impactful variants

**Usage:**
```bash
conda run -n alphagenome-env python plot_activation_landscape.py \
  --input ../results/landscape/AP1/AP1_activation_landscape.tsv \
  --output-dir ../figures/landscape
```

## Outputs

### Data Files (`results/landscape/AP1/`)

| File | Description | Size |
|------|-------------|------|
| `AP1_activation_landscape.tsv` | Complete landscape with all metrics | ~1 MB |
| `AP1_high_priority_candidates.tsv` | High-priority quadrant only | ~250 KB |
| `AP1_landscape_summary.txt` | Summary statistics | ~2 KB |

### Figures (`figures/landscape/`)

| File | Description |
|------|-------------|
| `AP1_activation_landscape_main.png` | Primary landscape visualization |
| `AP1_activation_landscape_comparison.png` | Multi-panel comparison |
| `AP1_tf_breakdown.png` | TF contribution analysis |
| `AP1_high_priority_candidates.png` | High-priority candidate details |

### Scientific Report

See `SCIENTIFIC_REPORT.md` for a complete publication-format analysis including:
- Abstract, Introduction, Methods
- Results with all statistics
- Discussion and conclusions
- Embedded figure references

## Key Results (AP1)

| Metric | Value |
|--------|-------|
| Total variants | 6,810 |
| High-priority candidates | 1,624 (23.8%) |
| Strong AP1 gain (>90th %ile) | 6,170 (90.6%) |
| Ultra-rare variants (AF<0.01%) | 5,798 (85.1%) |
| AP1 vs Enhancer correlation | r = 0.545, p < 10⁻¹⁶ |

## Quadrant Definitions

| Quadrant | X | Y | Interpretation |
|----------|---|---|----------------|
| HIGH PRIORITY | < median | > median | Accessible + High Impact |
| High Impact, Hard to Access | ≥ median | > median | Impactful but rare |
| Accessible, Low Impact | < median | ≤ median | Common but weak effect |
| Low Priority | ≥ median | ≤ median | Rare and weak effect |

## Dependencies

- pandas
- numpy
- matplotlib
- scipy

## Configuration

AP1-family TFs are defined in `compute_activation_landscape.py`:

```python
AP1_FAMILY_TFS = [
    'JUND', 'JUN', 'JUNB',      # Jun family
    'FOS', 'FOSL1', 'FOSL2',    # Fos family
    'ATF3', 'ATF2', 'ATF7',     # ATF family
    'BATF', 'BATF2',            # BATF family
    'MAFK', 'MAFF', 'MAFG',     # Small Maf
]
```

To analyze a different TF, modify this list to include the relevant TF family members.

## Execution Time

| Step | Time |
|------|------|
| Load predictions (186M rows) | ~27 seconds |
| Compute AP1 scores (vectorized) | ~5 seconds |
| Compute X-axis and merge | ~1 second |
| Generate all plots | ~3 seconds |
| **Total** | **~35 seconds** |
