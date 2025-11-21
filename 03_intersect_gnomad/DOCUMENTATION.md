# Module 03: gnomAD Intersection

This module queries gnomAD population variation data for variants matching mutation path positions.

## Script

### `query_gnomad_vcfs.py`
Complete workflow that:
1. Creates BED file of query positions from mutation paths
2. Uses bcftools to query gnomAD VCF for variants at those positions  
3. Matches mutation steps to gnomAD variants (chr, pos, ref, alt)
4. Summarizes allele frequencies per mutation path

**Outputs:**
- `results/gnomad_intersection/{TF}/query_regions.bed` - BED file of query positions
- `results/gnomad_intersection/{TF}/gnomad_query_results.tsv` - Raw bcftools query output
- `results/gnomad_intersection/{TF}/paths_with_gnomad.tsv` - Mutation steps matched to gnomAD
- `results/gnomad_intersection/{TF}/path_af_summary.tsv` - Per-path AF summary (max_AF, mean_AF, accessibility_score)

## Usage

```bash
python 03_intersect_gnomad/query_gnomad_vcfs.py --config pipeline_config.yaml
```

## Requirements

- bcftools (for VCF querying)
- gnomAD v4.1 VCF file (configured in pipeline_config.yaml)
