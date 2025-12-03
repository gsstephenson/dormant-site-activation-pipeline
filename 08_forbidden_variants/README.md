# Module 08: Forbidden Variants Scoring

## Overview

Score the ~30,000 "forbidden variants" — single-nucleotide mutations that would create AP1 binding sites but are completely absent from 807,000 humans despite high sequencing coverage.

## Memory-Safe Design

**Problem:** Full predictions for 30K variants = 842M rows = 966 GB RAM (impossible with 188 GB)

**Solution:** Streaming/chunked approach:
- Process 500 variants per batch (~30 GB RAM)
- Write each batch to parquet immediately
- Keep only summary statistics in memory
- Final output: partitioned parquet directory + summary TSV

## Scientific Rationale

From Module 07, we identified:
- **31,174** possible H=1 mutations (single change → AP1 consensus)
- **29,957** at high-confidence callable sites (AN ≥ 100K)
- **Only 1** variant observed in gnomAD
- **76.5× depletion** (p = 4.2×10⁻³²)

The ~30,000 absent variants are "forbidden" by natural selection. This module scores them with AlphaGenome to predict:
1. Which would create the strongest AP1 binding?
2. Which would have the greatest enhancer activity changes?
3. Which are in regulatory regions of disease-relevant genes?

## Input

- `results/mutation_paths/AP1/paths.tsv` — All mutation paths (filter H=1)
- `results/gnomad_intersection/AP1/unique_variants.tsv` — Observed gnomAD variants (to exclude)
- Coverage file with AN data — To filter for high-confidence positions

## Output

- `results/forbidden_variants/AP1/forbidden_variants.tsv` — List of forbidden variants
- `results/forbidden_variants/AP1/predictions_summary.tsv` — Summary stats per variant
- `results/forbidden_variants/AP1/top_candidates.tsv` — Top 1000 by predicted impact  
- `results/forbidden_variants/AP1/partitions/` — Full predictions (partitioned parquet)

## Resource Requirements

| Metric | Value |
|--------|-------|
| Variants | ~30,000 |
| Batch size | 500 |
| RAM per batch | ~30 GB |
| Peak RAM | ~40 GB |
| Disk (partitions) | ~6 GB |
| Time estimate | ~11 hours |

## Usage

```bash
# Step 1: Extract forbidden variants
python 08_forbidden_variants/extract_forbidden_variants.py \
  --paths results/mutation_paths/AP1/paths.tsv \
  --observed results/gnomad_intersection/AP1/unique_variants.tsv \
  --coverage /mnt/data_1/gnomAD_data/raw/gnomad_v4.1/coverage/gnomad_AN_tabix.tsv.bgz \
  --output results/forbidden_variants/AP1/forbidden_variants.tsv

# Step 2: Score with AlphaGenome (streaming mode)
python 08_forbidden_variants/score_forbidden_variants.py \
  --input results/forbidden_variants/AP1/forbidden_variants.tsv \
  --output results/forbidden_variants/AP1/ \
  --batch-size 500

# Resume from crash (e.g., resume from batch 20)
python 08_forbidden_variants/score_forbidden_variants.py \
  --input results/forbidden_variants/AP1/forbidden_variants.tsv \
  --output results/forbidden_variants/AP1/ \
  --batch-size 500 \
  --resume-from 20
```

## Results (Run: Dec 2-3, 2025)

### Extraction
- **29,957** forbidden variants extracted (H=1, absent from gnomAD, AN ≥ 100K)
- Median coverage: AN = 152,254 (high confidence)

### AlphaGenome Scoring
- **Runtime:** 11.4 hours
- **Total predictions:** 804,584,291 track-level scores
- **Disk usage:** 5.85 GB (60 parquet partitions)
- **Failures:** 0

### Top Candidates by Maximum Raw Score
| Variant | Raw Max | AN | Description |
|---------|---------|-----|-------------|
| chr1:204410195:A>C | 260,224 | 152,272 | Highest impact |
| chr11:60843382:T>A | 243,904 | 152,352 | Histone ChIP-seq |
| chr5:69369671:A>T | 220,576 | 152,364 | Multi-tissue |
| chr18:63366657:T>A | 209,504 | 152,354 | Keratinocyte |
| chr2:186591142:T>A | 202,416 | 152,332 | Broad impact |

### Key Findings
- **Mutation bias:** T↔A transversions dominate (58% of forbidden variants)
- **Most affected assay:** Histone ChIP-seq (raw_max up to 243,904)
- **Quantile distribution:** 47% have quantile_max > 0.999 (extreme predicted impact)
- **Broad impact:** Mean of ~26,858 tracks affected per variant

## Data Dependencies

### gnomAD Coverage Data
The coverage file with Allele Number (AN) values is critical for determining which genomic positions have sufficient sequencing depth:

```bash
# Location (pre-downloaded)
/mnt/data_1/gnomAD_data/raw/gnomad_v4.1/coverage/gnomad_AN_tabix.tsv.bgz
/mnt/data_1/gnomAD_data/raw/gnomad_v4.1/coverage/gnomad_AN_tabix.tsv.bgz.tbi

# Download script available at:
# data/gnomad_coverage/download_gnomad_coverage.sh
```

AN ≥ 100,000 means the position was successfully sequenced in ≥50,000 individuals (diploid), providing high confidence that variant absence is biological rather than technical.
