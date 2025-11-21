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
- `results/gnomad_intersection/{TF}/paths_with_gnomad.tsv` - Mutation steps matched to gnomAD
- `results/gnomad_intersection/{TF}/path_af_summary.tsv` - Per-path AF summary (max_AF, mean_AF, accessibility_score)

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
