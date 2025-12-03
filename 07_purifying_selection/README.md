# Module 07: Purifying Selection Analysis

## Overview

This module uses gnomAD v4.1 **allele number (AN) data for all sites** to rigorously test the hypothesis that dormant AP1 site activation is under purifying selection. Unlike previous modules that only used VCF files (which contain only observed variants), this module leverages coverage data to distinguish between:

- **Truly constrained sites:** High coverage (AN ≥ 100K), variant absent → Strong purifying selection
- **Uncertain sites:** Low coverage (AN < 50K), variant absent → Cannot conclude constraint
- **Missing sites:** Not in coverage file → Uncallable region, exclude from analysis

## Scientific Rationale

### The Coverage Problem

In Modules 00-06, we identified 31,174 possible single-nucleotide mutations that would create functional AP1 binding sites (Hamming distance = 1). Only 1 was observed in gnomAD. This suggests 36× depletion compared to baseline.

**But there's a problem:** VCF files only tell us what variants EXIST. They don't tell us whether absent variants are:
1. **Truly absent** (site was sequenced in 800K individuals, no one has it)
2. **Uncertain** (site was poorly covered, we can't tell)

### The Solution: AN (Allele Number) for All Sites

From gnomAD v4.1 release notes (April 2024):
> "We have calculated AN information for loci that do not yet have variation 
> observed within gnomAD. This means that even for loci where zero gnomAD 
> samples have a non-reference genotype call, we report the exact AN based 
> on the total number of defined sample genotype calls at that site."

**AN = 2 × number of individuals with a genotype call at that position**

- AN = 1,614,006 → All 807,003 individuals were called → Maximum confidence
- AN = 100,000 → ~50,000 individuals called → High confidence
- AN = 10,000 → ~5,000 individuals called → Low confidence

### Refined Constraint Analysis

With coverage data, we can now calculate:

```
Constraint Score = (Expected variants) / (Observed variants)

Where:
  Expected = Σ (AN_i / max_AN) for all positions  # Coverage-adjusted expectation
  Observed = Number actually seen in gnomAD
```

This gives us a **coverage-corrected** measure of purifying selection.

## Data Requirements

### Input Files

| File | Source | Description |
|------|--------|-------------|
| `data/gnomad_coverage/gnomad.genomes.v4.1.allele_number_all_sites.tsv.bgz` | Module 00 | AN for all callable sites (~12 GB) |
| `results/mutation_paths/AP1/paths.tsv` | Module 02 | All possible mutation paths |
| `results/landscape/AP1/AP1_activation_landscape.tsv` | Module 05 | Observed variants with AlphaGenome scores |

### Output Files

| File | Description |
|------|-------------|
| `results/purifying_selection/AP1/h1_coverage_analysis.tsv` | Coverage data for all H=1 positions |
| `results/purifying_selection/AP1/constraint_by_hamming.tsv` | Constraint scores by Hamming distance |
| `results/purifying_selection/AP1/purifying_selection_summary.txt` | Statistical summary |
| `figures/purifying_selection/constraint_evidence.png` | Main visualization |

## Scripts

### `analyze_coverage_constraint.py`

Main analysis script that:
1. Loads all possible H=1, H=2, H=3 mutation positions from paths.tsv
2. Queries gnomAD coverage file for AN at each position
3. Computes coverage-adjusted expected vs observed variants
4. Performs binomial tests for significance
5. Outputs constraint analysis tables

**Usage:**
```bash
python analyze_coverage_constraint.py \
    --paths ../results/mutation_paths/AP1/paths.tsv \
    --coverage ../data/gnomad_coverage/gnomad.genomes.v4.1.allele_number_all_sites.tsv.bgz \
    --observed ../results/landscape/AP1/AP1_activation_landscape.tsv \
    --output ../results/purifying_selection/AP1 \
    --threads 16
```

**Performance Notes:**
- Uses streaming approach to process 12GB coverage file
- Loads query positions into memory first (set lookup = O(1))
- With 32 cores and 188GB RAM, completes in ~2-3 minutes

### `plot_constraint_evidence.py`

Generates publication-quality figures:
1. **Coverage distribution:** AN values for dormant sites
2. **Constraint by Hamming:** Observed vs expected by mutation distance
3. **Fold depletion:** Selection intensity visualization
4. **Statistical summary:** P-values and confidence intervals

**Usage:**
```bash
python plot_constraint_evidence.py \
    --input ../results/purifying_selection/AP1 \
    --output ../figures/purifying_selection
```

## Expected Results

Based on preliminary analysis:

| Hamming | Possible Sites | High-Conf Sites (AN≥100K) | Observed | Fold Depletion |
|---------|----------------|---------------------------|----------|----------------|
| 1       | 31,174         | TBD                       | 1        | TBD            |
| 2       | 774,398        | TBD                       | 677      | TBD            |
| 3       | 5,519,454      | TBD                       | 6,359    | TBD            |

**Key question this module answers:** Of the ~31,000 possible H=1 mutations, how many are at well-covered positions where we can confidently say the variant is truly absent?

## Confidence Thresholds

| Threshold | AN Value | Samples Called | Interpretation |
|-----------|----------|----------------|----------------|
| High      | ≥ 100,000 | ~50,000+      | Strong confidence in constraint call |
| Medium    | 50,000-100,000 | 25,000-50,000 | Moderate confidence |
| Low       | < 50,000  | < 25,000      | Cannot reliably conclude constraint |
| Missing   | Not in file | 0           | Uncallable region, exclude |

## Dependencies

```
pandas>=1.5.0
numpy>=1.20.0
scipy>=1.7.0
matplotlib>=3.5.0
```

## References

1. gnomAD v4.1 Release Notes: https://gnomad.broadinstitute.org/news/2024-04-gnomad-v4-1/
2. Karczewski et al. (2020). The mutational constraint spectrum quantified from variation in 141,456 humans. Nature.
