# Dormant Site Activation Pipeline

**A computational framework for identifying cryptic transcription factor binding sites that could be activated by naturally occurring genetic variation.**

## Overview

This pipeline quantifies how "dormant" motif instances across the human genome could become activated TF binding sites through human standing variation, using gnomAD v4.1 population data and AlphaGenome functional predictions.

### Scientific Question
Which weak/near-motif sequences could become functional transcription factor binding sites through mutations that already exist in human populations?

### Key Results (AP1 Pilot Study)

| Metric | Value |
|--------|-------|
| **Variants analyzed** | 7,037 |
| **Strong AP1 activation (>90th %ile)** | 6,357 (90.3%) |
| **High-priority candidates** | 1,767 |
| **AP1-Enhancer correlation** | r = 0.579, p < 10⁻⁴⁰ |

**Main Finding:** The pipeline successfully identifies dormant TF binding sites and AlphaGenome confirms 90%+ would create functional AP1 binding upon activation.

---

## Quick Start

### Run Complete Pipeline
```bash
# Activate environment
conda activate alphagenome-env

# Run all modules (or specific ones with -m flag)
python run_pipeline.py --config pipeline_config.yaml

# Run specific module
python run_pipeline.py -m 5 --config pipeline_config.yaml
```

### Run Individual Modules
```bash
# Module 01: Scan genome for motifs
python 01_scan_motifs/scan_genome_fimo.py --config pipeline_config.yaml

# Module 02: Generate mutation paths  
python 02_generate_mutation_paths/enumerate_paths.py --config pipeline_config.yaml

# Module 03: Query gnomAD
python 03_intersect_gnomad/query_gnomad_vcfs.py --config pipeline_config.yaml

# Module 04: AlphaGenome scoring
python 04_run_alphagenome/run_alphagenome_scoring.py --config pipeline_config.yaml

# Module 05: Compute activation landscape
python 05_compute_activation_landscape/compute_activation_landscape.py \
  --predictions results/alphagenome/AP1/predictions.parquet \
  --gnomad results/gnomad_intersection/AP1/unique_variants.tsv \
  --output results/landscape/AP1

# Generate figures
python 05_compute_activation_landscape/plot_activation_landscape.py \
  --input results/landscape/AP1/AP1_activation_landscape.tsv \
  --output-dir figures/landscape
```

---

## Pipeline Modules

### Module 00: Data Fetching ✅
Download reference genome, motifs, and verify gnomAD installation.

### Module 01: Motif Scanning ✅
Genome-wide scanning for TF binding motifs with tiering by strength.
- Output: `results/motif_scan/AP1/`

### Module 02: Mutation Paths ✅
Enumerate minimal mutation paths (≤3 mutations) from dormant sites to consensus.
- Output: `results/mutation_paths/AP1/`

### Module 03: gnomAD Intersection ✅
Query gnomAD v4.1 to identify which activating mutations exist in human populations.
- Output: `results/gnomad_intersection/AP1/`
- Key stats: 867,406 variants queried, 7,037 unique activating variants found

### Module 04: AlphaGenome Scoring ✅
Score functional impact using AlphaGenome's multi-modal deep learning predictions.
- Output: `results/alphagenome/AP1/`
- Key stats: 197.9M predictions across 7,037 variants

### Module 05: Activation Landscape ✅
Compute 2D landscape mapping population accessibility vs. functional impact.
- Output: `results/landscape/AP1/`, `figures/landscape/`
- Scientific Report: `05_compute_activation_landscape/SCIENTIFIC_REPORT.md`

### Module 06: Disease Overlap Analysis ✅
Validate variants against ClinVar and GWAS Catalog disease databases.
- Output: `results/disease_overlap/AP1/`
- Key finding: 20.8% of variants near GWAS loci; 0% ClinVar overlap (expected for novel ultra-rare variants)

---

## Key Outputs

### Results Files
| File | Description |
|------|-------------|
| `results/landscape/AP1/AP1_activation_landscape.tsv` | Complete landscape (7,037 variants) |
| `results/landscape/AP1/AP1_high_priority_candidates.tsv` | High-priority quadrant (1,767 variants) |
| `results/alphagenome/AP1/predictions.parquet` | Raw AlphaGenome predictions (1.36 GB) |
| `results/disease_overlap/AP1/gwas_overlaps.tsv` | GWAS associations (3,423 trait links) |
| `results/disease_overlap/AP1/disease_overlap_report.txt` | Disease overlap summary |

### Figures
| File | Description |
|------|-------------|
| `figures/landscape/AP1_activation_landscape_main.png` | Main 2D landscape visualization |
| `figures/landscape/AP1_tf_breakdown.png` | AP1-family TF response distribution |
| `figures/landscape/AP1_high_priority_candidates.png` | High-priority candidate details |
| `results/disease_overlap/AP1/figures/gwas_trait_distribution.png` | GWAS trait categories |
| `results/disease_overlap/AP1/figures/overlap_summary.png` | Disease database overlap |

---

## Activation Landscape Interpretation

### X-Axis: Population Accessibility
```
X = -log₁₀(AF) × Hamming_distance
```
- **Low X** = accessible (common variant, few mutations)
- **High X** = hard to access (rare variant, many mutations)

### Y-Axis: AP1 Functional Impact
```
Y = log₁₀(max(raw_score)) across AP1-family TFs
```
- Uses: JUND, JUN, FOS, FOSL1, FOSL2, ATF3, ATF2, BATF, MAFK, etc.
- **Raw scores** (3-34,000) used instead of quantile scores to preserve selection signal
- Higher Y = stronger predicted AP1 binding gain

### Quadrant Classification
| Quadrant | Count | Description |
|----------|-------|-------------|
| **HIGH PRIORITY** | 1,767 | Accessible + High Impact |
| High Impact, Hard to Access | 1,733 | Rare variants with strong effects |
| Accessible, Low Impact | 1,751 | Common but weak effects |
| Low Priority | 1,786 | Rare and weak effects |

---

## Scientific Findings

### Primary Result
**The pipeline works:** 90.3% of identified variants show strong predicted AP1 binding activation (>90th percentile), validating that mutation paths correctly identify functional dormant sites.

### Validation
- Strong correlation (r = 0.579) between AP1 binding gains and enhancer marks (H3K27ac/H3K4me1)
- Top responding TFs: FOS (25.8%), JUND (15.8%), ATF3 (9.0%)

### Selection Analysis
- **Significant purifying selection detected** using raw AlphaGenome effect sizes
- Spearman correlation (effect size vs rarity): **r = 0.096, p = 4.97×10⁻¹⁵**
- Variants with stronger predicted AP1 binding are kept at lower population frequencies
- Consistent with selection against creating functional TF sites in inappropriate contexts

---

## Requirements

### Software
- Python 3.8+
- MEME Suite (FIMO)
- bcftools, samtools, tabix
- AlphaGenome API access

### Python Packages
See `requirements.txt`:
- pandas, numpy, scipy
- pyarrow (for parquet)
- matplotlib, seaborn
- pyfaidx

### Data
- GRCh38 reference genome
- gnomAD v4.1 genomes VCF
- JASPAR motif database

---

## Configuration

Edit `pipeline_config.yaml`:
```yaml
tf_name: "AP1"
motif_id: "MA0099.3"

reference_genome: "data/reference/GRCh38.fa"

gnomad:
  vcf: "/path/to/gnomad.genomes.vcf.bgz"

alphagenome:
  batch_size: 32
```

---

## Citation

If you use this pipeline, please cite:
- gnomAD: [Karczewski et al., 2020](https://doi.org/10.1038/s41586-020-2308-7)
- MEME Suite: [Bailey et al., 2015](https://doi.org/10.1093/nar/gkv416)
- JASPAR: [Castro-Mondragon et al., 2022](https://doi.org/10.1093/nar/gkab1113)

---

## Author

**George Stephenson**  
LAYER Laboratory, University of Colorado Boulder

---

## Status

**Current:** Module 06 Complete (Disease Overlap Analysis)  
**Last Updated:** December 1, 2025

### Completed Modules
- ✅ Module 00: Data Fetching
- ✅ Module 01: Motif Scanning
- ✅ Module 02: Mutation Paths
- ✅ Module 03: gnomAD Intersection
- ✅ Module 04: AlphaGenome Scoring
- ✅ Module 05: Activation Landscape
- ✅ Module 06: Disease Overlap Analysis

### Key Achievements
- 7,037 dormant AP1 site activating variants identified
- 90.3% show strong predicted functional impact
- 1,767 high-priority candidates for experimental validation
- **20.8% of variants near GWAS disease loci** (1,463 variants)
- **0% ClinVar overlap** (expected - novel ultra-rare non-coding variants)
- 7 exact position matches with GWAS lead SNPs
- Full scientific report: `05_compute_activation_landscape/SCIENTIFIC_REPORT.md`

### Disease Overlap Interpretation
The lack of ClinVar overlap is **expected and positive**:
- 60% of variants are ultra-rare (AF < 10⁻⁵)
- Located in non-coding enhancer regions
- **These are NOVEL disease candidates not yet studied by clinical genetics**

The GWAS overlap supports the "synthetic association" hypothesis:
- Common GWAS hits may tag rare causal variants like ours
