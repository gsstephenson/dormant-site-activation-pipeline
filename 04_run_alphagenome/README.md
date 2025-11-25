# Module 04: AlphaGenome Variant Scoring

Score functional impact of activating mutations using AlphaGenome with 1MB context windows.

## Overview

This module uses the AlphaGenome API to predict chromatin accessibility and TF binding changes for each variant matched to gnomAD population data.

**Implementation follows the proven methodology from `alphagenome-qtl-validation` repository.**

## Key Features

- **1MB Context Windows**: Uses `dna_client.SEQUENCE_LENGTH_1MB` (1,048,576 bp) as specified in AlphaGenome paper
- **CHIP_HISTONE Scorer**: Predicts TF binding and chromatin accessibility across multiple cell types/tissues
- **Direct API Usage**: Uses `client.score_variant()` - no manual sequence extraction needed
- **Batch Processing**: Scores all variants from Module 03 gnomAD intersection

## Input

From Module 03:
```
results/gnomad_intersection/AP1/paths_with_gnomad.tsv  # All mutation steps
results/gnomad_intersection/AP1/unique_variants.tsv    # Deduplicated (recommended)
```

**Current dataset:**
- 18,138,332 total mutation steps
- 38,961 steps matched to gnomAD (all chromosomes)
- 6,921 unique genomic variants after deduplication
- Coverage: ALL 24 chromosomes (chr1-22, chrX, chrY) ✅
- Note: Deduplication removes redundant measurements where the same genomic variant appears in multiple mutation paths

Expected columns:
- `chr`, `genomic_position`, `ref_base`, `alt_base`
- `AF`, `AC`, `AN` (gnomAD allele frequency data)
- `is_missing`, `coverage_confidence` (Module 03 coverage quality metrics)
- `path_id`, `step_num`, `site_id`, `tier`

**Note on Module 03 updates:** The input now includes AN-based coverage confidence metrics:
- `is_missing`: Flag for variants not found in gnomAD
- `coverage_confidence`: high/medium/low/missing based on AN (allele number)
- These metrics help distinguish true constraint (AF=0 with high AN) from low coverage

## Output

```
results/alphagenome/AP1/
├── predictions.parquet          # Raw scores (variant-track level)
└── predictions_summary.tsv      # Mean scores per variant
```

### Output Schema

**predictions.parquet** (variant-track pairs):
- `variant_name`: chr_pos_ref_alt
- `track_name`: Cell type/tissue track ID
- `quantile_score`: AlphaGenome functional score
- `gnomad_AF`, `gnomad_AC`, `gnomad_AN`: Population frequency
- `path_id`, `mutation_step`: Traceability to Module 02

**predictions_summary.tsv** (per-variant aggregates):
- `variant_name`: Unique variant identifier
- `quantile_score`: Mean across all tracks
- `n_tracks`: Number of tracks contributing
- Population and path metadata

## Usage

### Prerequisites

1. **AlphaGenome Environment**:
   ```bash
   conda activate alphagenome-env
   ```

2. **API Key**: Set environment variable or pass via flag:
   ```bash
   export ALPHAGENOME_API_KEY="your_api_key_here"
   ```

### Step 1: Prepare Unique Variants (REQUIRED)

**Why deduplicate?** The same genomic variant (chr:pos:ref>alt) can appear in multiple mutation paths to different motif sites. AlphaGenome scores genomic positions, not paths - scoring the same variant 6 times gives identical results but wastes API calls.

**Deduplication is scientifically correct:**
- Preserves all unique genomic variants
- One functional measurement per variant (ΔAlphaGenome is position-specific)
- Scores can be merged back to all paths for downstream analysis
- Reduces 38,961 steps → 6,921 unique variants (saves ~6.6 hours)

**Critical: AF=0 variants are INCLUDED**
- **KEEP:** Variants observed in gnomAD (is_missing=0), including AF=0 with known AN
- **SKIP:** Variants missing from gnomAD (is_missing=1, unknown coverage)
- **Rationale:** AF=0 with high AN (≥50K) represents true purifying selection - these are the most interesting constraint candidates and MUST be scored by AlphaGenome

```bash
python 04_run_alphagenome/prepare_unique_variants.py
```

This creates `results/gnomad_intersection/AP1/unique_variants.tsv` with 6,921 unique variants.

**Output includes:**
- Coverage breakdown (observed vs missing variants)
- Coverage confidence distribution (high/medium/low based on AN)
- AN (allele number) statistics for unique variants
- Identifies high-confidence variants (AN≥50K) suitable for downstream analysis

**Expected runtime:** ~84 minutes (~1.4 hours) at 1.37 variants/sec

### Step 2: Test Mode

Test on a small batch first:

```bash
python 04_run_alphagenome/run_alphagenome_scoring.py \
    --config pipeline_config.yaml \
    --limit 50
```

### Step 3: Full Run

Score all unique variants (~84 minutes / 1.4 hours):

```bash
python 04_run_alphagenome/run_alphagenome_scoring.py \
    --config pipeline_config.yaml
```

The script automatically uses `unique_variants.tsv` if available, otherwise falls back to `paths_with_gnomad.tsv`.

## Performance Estimates

Based on testing:
- **Speed**: ~1.37 variants/second
- **2,701 unique variants**: ~33 minutes
- **Memory**: <4GB typical
- **Output size**: ~3MB parquet + summary TSV per 1,000 variants

## Implementation Details

### Code Pattern (from alphagenome-qtl-validation)

```python
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers

# Initialize
client = dna_client.create(api_key=api_key)
scorers = [variant_scorers.RECOMMENDED_VARIANT_SCORERS['CHIP_HISTONE']]

# For each variant:
variant = genome.Variant(chromosome, position, ref, alt, name)
interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
scores = client.score_variant(interval, variant, scorers, organism)

# Aggregate results
df = variant_scorers.tidy_scores(results)
```

### Why 1MB Context?

From AlphaGenome paper:
- Long-range regulatory interactions can span hundreds of kb
- 1MB captures distal enhancers and topological domains
- Validated optimal window size for functional predictions

## Quality Control

Monitor for:
1. **API Errors**: Network issues, rate limits
2. **Reference Mismatches**: Variant position issues
3. **Missing Scores**: Failed predictions
4. **Score Distribution**: Expect quantile_score range [0, 1]

Check logs:
```bash
tail -f logs/04_alphagenome_scoring_AP1.log
```

## Next Step

Module 05: Compute 2D activation landscape (population frequency × functional impact)

## References

- AlphaGenome Paper: Nature 2024
- Working Implementation: `alphagenome-qtl-validation/scripts/02_predict.py`
- QTL Validation: r=0.40 correlation with GTEx caQTLs
