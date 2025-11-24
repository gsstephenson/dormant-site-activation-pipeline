# Module 03: gnomAD Intersection

This module queries gnomAD population variation data for variants matching mutation path positions.

## Script

### `query_gnomad_vcfs.py`
Complete workflow that:
1. Creates BED file of query positions from mutation paths
2. Uses bcftools to query gnomAD VCF for variants at those positions (parallel processing)
3. Matches mutation steps to gnomAD variants (chr, pos, ref, alt)
4. Summarizes allele frequencies per mutation path

**Performance Optimizations:**
- **Parallel Processing:** 30 chromosomes processed simultaneously (default)
- **Timeout:** 6 hours per chromosome (extended from 2 hours)
- **Total CPU Utilization:** ~30 cores (parallel chromosome processing)
- **Resource Requirements:** Best performance with 32+ CPU cores, 64+ GB RAM

**Optimization Results (Verified Nov 2025):**
- **Success Rate:** 24/24 chromosomes completed (100%)
  - Previous 2-hour timeout: 14/24 success (58%), 10 timeouts on large chromosomes
  - Current 6-hour timeout: 24/24 success, no timeouts
- **Data Retrieved:** 867,406 variants (2.66x improvement)
- **Matched Mutations:** 40,195 mutation steps matched to gnomAD
- **Actual Runtime:** ~3h 13min on 32-core system with ~3.2M query positions
- **Largest Chromosomes:** chr1 (67K variants), chr2 (75K variants) completed successfully

**Outputs:**
- `results/gnomad_intersection/{TF}/query_regions.bed` - BED file of query positions
- `results/gnomad_intersection/{TF}/gnomad_query_results.tsv` - Raw bcftools query output
- `results/gnomad_intersection/{TF}/paths_with_gnomad.tsv` - Mutation steps matched to gnomAD (includes `is_af_zero` flag)
- `results/gnomad_intersection/{TF}/path_af_summary.tsv` - Per-path AF summary (max_AF, mean_AF, sum_AF, total_AC, total_nhomalt, **num_af_zero**)

## AF=0 Handling (Critical for Constraint Analysis)

**Biological Significance:**  
Variants with AF=0 (not observed in gnomAD's 807K individuals) represent the **most constrained** mutations. These are likely under strong purifying selection and will appear at the far right of the X-axis landscape.

**Implementation:**
1. **Missing variants** (no match in gnomAD) are assigned **AF=0.0** (literal zero)
2. **Storage:** AF=0.0 is stored as a true zero in all output files
3. **Tracking:** The `is_af_zero` flag (0/1) and `num_af_zero` count enable constraint analysis
4. **X-axis calculation:** Epsilon substitution (1e-12) is applied **ONLY** during landscape X-axis computation to avoid log(0) errors - NOT during data storage

**Output Columns:**
- `paths_with_gnomad.tsv`: 
  * `is_af_zero` (1 if AF=0, 0 otherwise)
  * `is_missing` (1 if no gnomAD record, 0 if matched)
  * `coverage_confidence` (missing/high/medium/low/zero_an)
  * `AN` (allele number = 2× sample size at this position)

- `path_af_summary.tsv`: 
  * `num_af_zero` (count of AF=0 steps per path)
  * `num_missing` (count of missing variants per path)
  * `min_AN`, `mean_AN` (coverage metrics across steps)
  * `all_af_zero_high_conf` (1 if all steps AF=0 with AN≥50K)

**Constraint Metrics:**
- **High-confidence constraint**: ALL steps AF=0 **AND** min_AN ≥ 50,000 (well-covered in >25K individuals)
- **Uncertain constraint**: ALL steps AF=0 **BUT** low AN or missing gnomAD records
- Paths with SOME steps AF=0: Partially constrained activation route
- Mean num_af_zero per path: Overall constraint level across the dataset

**Coverage Confidence Categories (based on AN):**
- `high`: AN ≥ 50,000 (>25K individuals, strong signal)
- `medium`: AN ≥ 10,000 (>5K individuals, moderate signal)
- `low`: AN > 0 but < 10,000 (poor coverage, weak signal)
- `missing`: No gnomAD record (unknown coverage)
- `zero_an`: AN = 0 (edge case, shouldn't occur)

**Interpreting Results:**
- Use `all_af_zero_high_conf=1` to identify **true constrained sites** with high statistical power
- Use AN percentiles to determine appropriate coverage cutoffs for your analysis
- Filter by `coverage_confidence='high'` for most reliable constraint inferences

## Usage

```bash
python 03_intersect_gnomad/query_gnomad_vcfs.py --config pipeline_config.yaml
```

**Optional Arguments:**
- `--n-jobs N` - Override number of parallel jobs (default: 30)
- `--skip-query` - Skip bcftools query, use existing results

**Expected Runtime:**
- Small chromosomes (19-22, X, Y): 15-45 minutes each
- Medium chromosomes (9-18): 1-3 hours each  
- Large chromosomes (1-8): 3-4.5 hours each
- **Total: ~3-4 hours** with 24-30 parallel jobs on 32-core system

## Requirements

- bcftools (for VCF querying with threading support)
- gnomAD v4.1 VCF files (per-chromosome, configured in pipeline_config.yaml)
- 32+ CPU cores recommended for optimal performance
- 64+ GB RAM recommended
