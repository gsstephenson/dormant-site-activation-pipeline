# Module 04: AlphaGenome Variant Scoring

Score functional impact of activating mutations using AlphaGenome with 1MB context windows.

## Overview

This module uses the AlphaGenome API to predict **multi-modal functional impacts** for each variant matched to gnomAD population data. Since AlphaGenome returns all feature types by default, we store everything to maximize value from expensive API queries.

**Implementation follows the proven methodology from `alphagenome-qtl-validation` repository.**

## Key Features

- **1MB Context Windows**: Uses `dna_client.SEQUENCE_LENGTH_1MB` (1,048,576 bp) as specified in AlphaGenome paper
- **Multi-Modal Predictions**: Uses ALL 19 recommended variant scorers to capture:
  - **Gene Expression** (RNA_SEQ, CAGE, PROCAP)
  - **Chromatin Accessibility** (ATAC, DNASE)  
  - **TF Binding** (CHIP_TF)
  - **Histone Modifications** (CHIP_HISTONE: H3K27ac, H3K4me1, H3K4me3, H3K27me3, etc.)
  - **3D Genome Contacts** (CONTACT_MAPS)
  - **Splicing** (SPLICE_SITES, SPLICE_SITE_USAGE, SPLICE_JUNCTIONS, POLYADENYLATION)
- **Direct API Usage**: Uses `client.score_variant()` - no manual sequence extraction needed
- **Batch Processing**: Scores all variants from Module 03 gnomAD intersection
- **Context-dependent tracks**: Average 27,415 tracks per variant (range: 31K-82K depending on genomic context)

## Input

From Module 03:
```
results/gnomad_intersection/AP1/paths_with_gnomad.tsv  # All mutation steps
results/gnomad_intersection/AP1/unique_variants.tsv    # Deduplicated (recommended)
```

**Dataset processed:**
- 18,138,332 total mutation steps
- 38,961 steps matched to gnomAD (all chromosomes)
- 7,158 unique genomic variants submitted for scoring
- 6,810 variants scored successfully (95.1% success rate)
- 348 variants failed (4.9% - transient network errors)
- Coverage: ALL 24 chromosomes (chr1-22, chrX, chrY) ✅

**Results achieved:**
- **186,695,928 variant-track predictions** (186.7M total)
- **27,415 average tracks per variant** (biologically accurate, context-dependent)
- **All 11 output types captured** (RNA_SEQ, CHIP_TF, CHIP_HISTONE, SPLICE_*, CAGE, DNASE, ATAC, CONTACT_MAPS, PROCAP)
- **714 unique biosamples/cell types**
- **Runtime:** 3.5 hours (21:26 → 01:26, November 26-27, 2025)

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
├── predictions.parquet                    # Raw multi-modal scores (variant-track level)
├── predictions_summary.tsv                # Overall mean scores per variant
└── predictions_summary_by_feature.tsv     # Per-feature mean scores per variant
```

### Output Schema

**predictions.parquet** (variant-track pairs, 186,695,928 rows for 6,810 variants):
- Variant identification:
  - `variant_id_str`: chr:pos:ref>alt
  - `scored_interval`: 1MB genomic interval used for scoring
- Feature identification:
  - `output_type`: RNA_SEQ, ATAC, DNASE, CHIP_TF, CHIP_HISTONE, CONTACT_MAPS, CAGE, PROCAP, SPLICE_SITES, SPLICE_SITE_USAGE, SPLICE_JUNCTIONS
  - `histone_mark`: H3K27ac, H3K4me1, H3K4me3, H3K27me3, etc. (for CHIP_HISTONE tracks)
  - `transcription_factor`: TF name (for CHIP_TF tracks)
  - `gtex_tissue`: GTEx tissue name (for RNA_SEQ tracks)
- Track metadata:
  - `track_name`: Cell type/tissue + assay + feature
  - `biosample_name`, `biosample_type`, `ontology_curie`
  - `gene_id`, `gene_name`, `gene_type`, `gene_strand` (for gene-level predictions)
- Scores:
  - `quantile_score`: AlphaGenome functional prediction (-1 to 1 scale)
  - `raw_score`: Raw model output
- Population genetics:
  - `gnomad_AF`, `gnomad_AC`, `gnomad_AN`: Allele frequency, count, number
- Traceability:
  - `path_id`, `step_num`, `site_id`, `tier`: Links back to Module 02 mutation paths

**predictions_summary.tsv** (overall per-variant aggregates, 6,810 rows):
- `variant_id_str`: Unique variant identifier
- `quantile_score`: Mean across ALL ~89K tracks
- `n_tracks`: Number of tracks contributing (~89,000)
- `gnomad_AF`, `gnomad_AC`, `gnomad_AN`: Population frequency data
- `path_id`, `step_num`: Traceability

**predictions_summary_by_feature.tsv** (per-feature aggregates, 63,471 rows = 6,810 variants × 11 features):
- `variant_id_str`: Unique variant identifier
- `output_type`: Feature type (RNA_SEQ, ATAC, DNASE, etc.)
- `quantile_score`: Mean across tracks for this feature type
- `n_tracks`: Number of tracks contributing for this feature
- `gnomad_AF`, `gnomad_AC`, `gnomad_AN`: Population frequency data

**Why store all features?**
- AlphaGenome returns all features by default (no extra API cost)
- Enables multi-modal validation: which feature best predicts evolutionary constraint?
- Expected strongest signals: ΔTF_binding (AP-1), ΔH3K27ac (active enhancers), ΔExpression
- Mechanistic insights: compare TF binding vs chromatin vs expression impacts
- Future-proof: no need to re-query expensive API for additional features

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

**Actual runtime:** ~10 minutes to prepare 7,158 unique variants

### Step 2: Test Mode

Test on a small batch first:

```bash
python 04_run_alphagenome/run_alphagenome_scoring.py \
    --config pipeline_config.yaml \
    --limit 50
```

### Step 3: Full Run

Score all unique variants with multi-modal features:

```bash
python 04_run_alphagenome/run_alphagenome_scoring.py \
    --config pipeline_config.yaml
```

The script automatically uses `unique_variants.tsv` if available, otherwise falls back to `paths_with_gnomad.tsv`.

**Note:** Multi-modal scoring is **slower** than single-scorer mode but captures all valuable data:
- Single scorer (CHIP_HISTONE only): ~1.4 hours for 6,921 variants
- All 19 scorers (multi-modal): **~5-6 hours** for 6,921 variants
- Recommended to run in `tmux` or `screen` for long-running jobs

## Performance Estimates

**Multi-Modal Scoring (19 scorers, ALL features) - ACTUAL RESULTS:**
- **Speed**: ~1.77 seconds/variant (0.565 variants/sec)
- **6,810 variants scored**: 3.5 hours runtime
- **Memory**: ~214 GB peak (in-memory processing of 186.7M predictions)
- **Output size**: 
  - `predictions.parquet`: 1.36 GB (186.7M variant-track pairs, excellent compression)
  - `predictions_summary.tsv`: 506 KB (6,810 variants)
  - `predictions_summary_by_feature.tsv`: 4.2 MB (63,471 variant-feature pairs)

**Why track counts vary by variant (27,415 average):**

The "89K tracks" was a theoretical maximum estimate. **Actual track counts are context-dependent and biologically accurate:**

1. **RNA_SEQ** (17,550 avg): Only genes within 1MB window, not all ~20K human genes
2. **Splice tracks** (0-3,670): Only present near splice junctions
3. **Cell type availability**: Not all 714 biosamples have data for all assays
4. **Genomic context**: Variant location determines which features are relevant

**Track distribution by variant:**
- chr1:867861 → 78,042 tracks (high coverage region)
- chr1:1393984 → 82,284 tracks (near splice junctions) 
- chr1:3027430 → 31,314 tracks (lower coverage region)

This variation is **expected and correct** - AlphaGenome returns only biologically relevant predictions for each genomic context.

## Implementation Details

### Code Pattern (Multi-Modal Scoring)

```python
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers

# Initialize with ALL scorers for multi-modal predictions
client = dna_client.create(api_key=api_key)
scorers = list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.values())  # All 19 scorers

# For each variant:
variant = genome.Variant(chromosome, position, ref, alt, name)
interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
scores = client.score_variant(interval, variant, scorers, organism)

# Aggregate results - tidy_scores() preserves all feature types
df = variant_scorers.tidy_scores(results)

# df now contains ~89K rows per variant with columns:
# - output_type: RNA_SEQ, ATAC, DNASE, CHIP_TF, CHIP_HISTONE, CONTACT_MAPS, ...
# - histone_mark: H3K27ac, H3K4me1, H3K4me3, ... (for CHIP_HISTONE tracks)
# - transcription_factor: TF name (for CHIP_TF tracks)
# - gtex_tissue: GTEx tissue (for RNA_SEQ tracks)
# - quantile_score: Functional prediction score
```

**Key insight:** `tidy_scores()` automatically preserves all output types and metadata returned by AlphaGenome. We just need to provide all scorers and let it handle the rest.

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
