# Dormant Site Activation Pipeline

**A computational framework for identifying cryptic transcription factor binding sites that could be activated by naturally occurring genetic variation.**

## Overview

This pipeline quantifies how "dormant" motif instances across the human genome could become activated TF binding sites through human standing variation, using gnomAD population data and AlphaGenome functional predictions.

### Scientific Question
Which weak/near-motif sequences could become functional transcription factor binding sites through mutations that already exist in human populations?

### Approach
1. **Scan genome** for all TF-like sequences (strong + weak + near-motif)
2. **Enumerate mutation paths** from weak sites to consensus motif
3. **Query gnomAD** to check if activating mutations exist in populations
4. **Score with AlphaGenome** to predict functional impact
5. **Map activation landscape**: population accessibility vs. functional impact

---

##  Quick Start

### 1. Setup Environment

```bash
# Clone repository
git clone <your-repo>
cd dormant_site_activation_pipeline

# Install Python dependencies
pip install -r requirements.txt

# Install system dependencies (via conda)
conda install -c bioconda meme samtools bcftools tabix
```

### 2. Download Data

```bash
# Download reference genome
bash 00_fetch_data/download_reference.sh

# Download AP1 motif from JASPAR
python 00_fetch_data/download_motifs.py --motif-id MA0099.3

# Verify gnomAD installation (assumes data already downloaded)
python 00_fetch_data/download_gnomad_index_info.py \
    --config pipeline_config.yaml \
    --create-symlink
```

### 3. Run Pipeline (Module by Module)

```bash
# Module 01: Scan genome for motifs
python 01_scan_motifs/scan_genome_fimo.py --config pipeline_config.yaml
python 01_scan_motifs/tier_sites.py --config pipeline_config.yaml

# Module 02: Generate mutation paths
python 02_generate_mutation_paths/enumerate_paths.py --config pipeline_config.yaml

# Module 03: Query gnomAD (WARNING: large chromosomes may timeout)
python 03_intersect_gnomad/query_parallel.py --config pipeline_config.yaml

# Module 04: AlphaGenome scoring (next)
python 04_run_alphagenome/make_variant_seqs.py --config pipeline_config.yaml
```

---

## Directory Structure

```
dormant_site_activation_pipeline/
├── README.md                    # This file
├── directive.md                 # Full specification
├── pipeline_config.yaml         # Central configuration
├── requirements.txt             # Python dependencies
│
├── data/                        # Input data
│   ├── reference/               # GRCh38 genome
│   ├── motifs/JASPAR/          # TF motif models
│   └── gnomad/                 # Symlink to population variants
│
├── 00_fetch_data/              # Data acquisition
├── 01_scan_motifs/             # Genome-wide motif scanning
├── 02_generate_mutation_paths/ # Enumerate activating mutations
├── 03_intersect_gnomad/        # Query population variants
├── 04_run_alphagenome/         # Functional predictions
├── 05_compute_activation_landscape/  # 2D landscape
├── 06_visualization/           # Plots and figures
│
├── utils/                      # Shared utilities
│   ├── constants.py
│   ├── io.py
│   ├── logging_utils.py
│   └── DOCUMENTATION.md
│
├── results/                    # Pipeline outputs
│   ├── motif_scan/
│   ├── mutation_paths/
│   ├── gnomad_intersection/
│   ├── alphagenome/
│   └── landscape/
│
├── figures/                    # Publication-ready plots
└── logs/                       # Execution logs
```

---

## Pipeline Modules

### Module 00: Data Fetching [Complete]
Download reference genome, motifs, and verify gnomAD installation.

**Status:** Complete  
**Documentation:** `00_fetch_data/DOCUMENTATION.md`

---

### Module 01: Motif Scanning [Complete]
Genome-wide scanning for TF binding motifs with tiering by strength.

**Status:** Complete  
**Documentation:** `01_scan_motifs/DOCUMENTATION.md`

**Key outputs:**
- All motif matches across genome
- Tiered sites (strong → very weak)
- Statistics and reports

**Usage:**
```bash
python 01_scan_motifs/scan_genome_fimo.py --config pipeline_config.yaml
python 01_scan_motifs/tier_sites.py --config pipeline_config.yaml
```

---

### Module 02: Mutation Paths [Complete]
Enumerate minimal mutation paths to consensus motif.

---

### Module 03: gnomAD Intersection [Complete]
Query population variants that match activating mutations with AN-based coverage confidence.

**Status:** ✅ Complete (ALL 24 chromosomes) with coverage quality metrics  
**Documentation:** `03_intersect_gnomad/DOCUMENTATION.md`

**Key features:**
- **AF=0 handling:** Missing variants assigned AF=0 (literal zero, epsilon only for X-axis)
- **Coverage confidence:** Uses AN (allele number) to distinguish true constraint from low coverage
  - High confidence: AN≥50K (>25K individuals)
  - Medium: AN≥10K
  - Low: AN>0 but <10K
  - Missing: No gnomAD record
- **Constraint metrics:** `num_af_zero`, `num_missing`, `all_af_zero_high_conf` flag

**Key outputs:**
- 867,406 gnomAD variants retrieved (ALL 24 chromosomes)
- 38,961 mutation steps matched to population data (0.21% of total)
- Allele frequency distribution: 39.6% rare, 44.7% low, 9.3% moderate, 6.3% common
- Coverage quality tracked via `min_AN`, `mean_AN`, `coverage_confidence` fields

**Coverage:**
- ✅ Successfully queried: ALL chromosomes (chr1-22, X, Y) - 100% success
- Extended timeout: 6 hours per chromosome (was 2 hours)
- Increased parallelism: 30 jobs (was 24)
- Runtime: ~3 hours 13 minutes on 32-core system

**Usage:**
```bash
python 03_intersect_gnomad/query_gnomad_vcfs.py --config pipeline_config.yaml
```

**Improvement:** 2.66x more variants than partial run (867K vs 326K). Complete chromosome coverage achieved by extending timeout and optimizing parallelization.

---

### Module 04: AlphaGenome Scoring [In Progress]
**Score functional impact using AlphaGenome with multi-modal feature capture (expression, ATAC, DNase, TF binding, histones, Hi-C).**

**Status:** 🔄 Running full dataset (7,158 unique variants, ~2.5 hours estimated)

**Implementation:** Follows proven methodology from `alphagenome-qtl-validation` repository
- Uses `dna_client.SEQUENCE_LENGTH_1MB` (1,048,576 bp) context windows
- Direct API usage via `client.score_variant()` 
- **ALL 19 recommended variant scorers** for comprehensive feature coverage
- Captures: RNA_SEQ, ATAC, DNASE, CHIP_TF, CHIP_HISTONE, CONTACT_MAPS, CAGE, PROCAP, SPLICE_*
- ~89,000 tracks per variant across all cell types and feature types

**Why store all features?**
- AlphaGenome returns all features by default (no extra API cost)
- Enables multi-modal validation: which feature best predicts evolutionary constraint?
- Expected strongest signals: ΔTF_binding (AP-1), ΔH3K27ac, ΔExpression
- Future-proof: no need to re-query expensive API

**Documentation:** `04_run_alphagenome/README.md`

**Prerequisites:**
```bash
conda activate alphagenome-env
export ALPHA_GENOME_KEY="your_api_key"
```

**Current dataset:**
- 38,961 mutation steps matched to gnomAD (ALL chromosomes)
- 6,921 unique genomic variants (deduplicated)
- Coverage: ALL 24 chromosomes ✅
- Complete genome-wide analysis

**Deduplication rationale:** Same genomic variant (chr:pos:ref>alt) can appear in multiple mutation paths. Since AlphaGenome scores genomic positions (not paths), deduplication avoids redundant API calls while preserving all information.

**Test setup:**
```bash
python 04_run_alphagenome/test_alphagenome_setup.py
```

**Usage:**
```bash
# Step 1: Prepare unique variants (deduplicates mutation paths)
python 04_run_alphagenome/prepare_unique_variants.py

# Step 2: Test on small batch
python 04_run_alphagenome/run_alphagenome_scoring.py \
    --config pipeline_config.yaml \
    --limit 2

# Step 3: Full run (6,921 unique variants, ~5-6 hours with multi-modal)
python 04_run_alphagenome/run_alphagenome_scoring.py \
    --config pipeline_config.yaml
```

**Expected outputs:**
- `results/alphagenome/AP1/predictions.parquet` - Raw multi-modal scores (~199M variant-track pairs, ~1.9 GB)
- `results/alphagenome/AP1/predictions_summary.tsv` - Overall mean scores per variant (7,158 variants)
- `results/alphagenome/AP1/predictions_summary_by_feature.tsv` - Per-feature means (~70K variant-feature pairs)
- Output size: 7,158 variants × ~27,850 tracks/variant = 199M predictions

**Performance:** 
- Multi-modal (19 scorers): ~1.27 sec/variant (much faster than expected)
- Full run: 7,158 variants ≈ **2.5 hours**
- **Recommendation:** Run in tmux/screen for long jobs

**Note:** Final output includes checkpoint/recovery mechanism for large dataset stability

### Module 05: Activation Landscape [Planned]
**Combine gnomAD constraint with AlphaGenome predictions to create 2D activation landscape.**

**Status:** 📋 Not yet implemented  
**Prerequisites:** ✅ Module 03 (gnomAD), 🔄 Module 04 (AlphaGenome)

**What it does:**
- Joins `paths_with_gnomad.tsv` with `predictions.parquet`
- Selects "best path" per motif site (shortest path → lowest AF tiebreaker)
- Computes landscape coordinates:
  - **X-axis (constraint):** `-log10(max(AF, 1e-12)) × hamming_distance`
  - **Y-axis (impact):** `max(ΔAlphaGenome)` [signed, preserves negatives]

**Key outputs:**
- `results/landscape/AP1/activation_landscape.tsv` - Final site-level table
  - Columns: site_id, best_path_id, coordinates, max_AF_step, X_constraint, max_delta, Y_impact
  - ~7,158 rows (one per unique variant with predictions)

**Usage:**
```bash
python 05_compute_activation_landscape/combine_population_and_impact.py \
    --config pipeline_config.yaml
python 05_compute_activation_landscape/classify_quadrants.py \
    --config pipeline_config.yaml
```

---

### Module 06: Visualization [Planned]
**Generate publication-ready plots of the activation landscape.**

**Status:** 📋 Not yet implemented  
**Prerequisites:** ✅ Module 05 (activation_landscape.tsv)

**What it does:**
- Creates 2D scatter/hexbin plots
- Highlights AF=0 sites, GWAS/ClinVar overlaps
- Generates genome browser-style tracks for top candidates

**Key outputs:**
- `figures/AP1/activation_landscape.png` - Main 2D plot
- `figures/AP1/activation_landscape_by_feature.png` - Separate plots per output_type
- `figures/AP1/top_candidates_tracks.png` - IGV-style genome tracks

**Usage:**
```bash
python 06_visualization/plot_landscape.py --config pipeline_config.yaml
python 06_visualization/plot_genome_tracks.py --config pipeline_config.yaml
```

---

### Module 07: Population Statistics [Planned]
**Describe global allele frequency distribution across all AP-1 motif sites.**

**Status:** 📋 Not yet implemented  
**Prerequisites:** ✅ Module 03 (gnomAD data only - does NOT require AlphaGenome)

**What it does:**
- Quantifies constraint landscape-wide
- Generates AF histogram with special handling for AF=0 sites
- Computes summary statistics on population accessibility

**Key outputs:**
- `results/population_stats/AP1/max_af_histogram.tsv`
- `figures/AP1/max_af_histogram.png`
- Summary: fraction with AF=0, ultra-rare sites, AN coverage distribution

**Usage:**
```bash
python 07_population_statistics/compute_af_distribution.py \
    --config pipeline_config.yaml
```

**Note:** This module can run independently of Module 04 - it only analyzes gnomAD data.

---

### Module 08: GWAS & ClinVar Integration [Planned]
**Annotate high-impact sites with disease associations for biological validation.**

**Status:** 📋 Not yet implemented  
**Prerequisites:** ✅ Module 05 (activation_landscape.tsv), ✅ ClinVar + GWAS data (downloaded)

**What it does:**
- Annotates motif sites with ClinVar pathogenic variants (exact position matches)
- Checks proximity to GWAS associations (±10kb window)
- Performs enrichment analysis: are high-constraint + high-impact sites enriched near disease variants?

**Key outputs:**
- `results/landscape/AP1/activation_landscape_annotated.tsv` - Adds columns:
  - `near_gwas`, `gwas_traits`, `gwas_distance`
  - `near_clinvar`, `clinvar_significance`, `clinvar_variant_id`
- `results/validation/AP1/gwas_enrichment.tsv` - Statistical tests
- Updated visualizations with disease overlays

**Usage:**
```bash
python 08_gwas_clinvar_integration/annotate_with_gwas_clinvar.py \
    --config pipeline_config.yaml
```

**Data sources:**
- ClinVar GRCh38: 4.1M variants (downloaded to `/mnt/data_1/clinvar_data/`)
- GWAS Catalog v1.0: 1M associations (downloaded to `/mnt/data_1/gwas_data/`)

---erate plots and ranked candidate lists.

---

## Configuration

Edit `pipeline_config.yaml` to customize:

```yaml
# Transcription factor
tf_name: "AP1"
motif_id: "MA0099.3"

# Data paths
reference_genome: "data/reference/GRCh38.fa"
gnomad:
  vcf: "/mnt/data_1/gnomAD_data/raw/gnomad_v4.1/genomes_vcf/genomes.vcf.bgz"

# Motif scanning
motif_scan:
  pvalue_threshold: 1e-3
  tier_cutoffs:
    tier0: 0.95  # strong
    tier1: 0.85  # medium
    tier2: 0.70  # weak

# AlphaGenome
alphagenome:
  model_path: "/mnt/models/alphagenome/"
**Latest Update:** 
- ✅ Module 03: Complete with 40,195 observed variants (ALL 24 chromosomes), 99.1% high coverage (AN≥50K)
- 🔄 Module 04: Running full AlphaGenome scoring (7,158 unique variants, ~2.5 hours)
- 📋 Module 05-08: Specifications complete, ready for implementation once Module 04 finishes

**Module Dependencies:**
- **Module 07** (Population Statistics): Only needs Module 03 ✅ - can run now
- **Module 05** (Activation Landscape): Needs Module 03 ✅ + Module 04 🔄
- **Module 06** (Visualization): Needs Module 05
- **Module 08** (GWAS/ClinVar): Needs Module 05 + downloaded data ✅
```

---

## Expected Outputs

### Activation Landscape
A 2D plot showing:
- **X-axis**: Population accessibility/constraint  
  Formula: `-log10(max_AF + 1e-12) × hamming_distance`
- **Y-axis**: Functional impact  
  Formula: `max(ΔAlphaGenome_feature)`

### Ranked Candidates
List of dormant sites with:
- High functional impact potential
- Activating mutations present in populations
- Low selection constraint

---

## Use Cases

### 1. AP1 Prototype (Default)
Identify dormant AP1 binding sites across the genome.

### 2. Any Transcription Factor
Generalize to other TFs by changing config:
```yaml
tf_name: "CTCF"
motif_id: "MA0139.1"
```

### 3. Disease Association
Cross-reference with GWAS hits to find:
- Disease variants that activate cryptic sites
- Regulatory explanations for non-coding variants

### 4. Evolution Studies
Examine:
- Constraint on weak binding sites
- Selection against activating mutations
- TF binding site gain/loss dynamics

---

## Documentation

Each module has detailed documentation:
- `00_fetch_data/DOCUMENTATION.md` - Data acquisition
- `01_scan_motifs/DOCUMENTATION.md` - Motif scanning
- `utils/DOCUMENTATION.md` - Utility functions

---

## Troubleshooting

### FIMO not found
```bash
conda install -c bioconda meme
```

### gnomAD files missing
```bash
python 00_fetch_data/download_gnomad_index_info.py --show-download-instructions
```

### Out of memory
- Scan specific chromosomes
- Increase system memory
- Use stricter P-value threshold

---

## Citation

If you use this pipeline, please cite:
- gnomAD: [Karczewski et al., 2020](https://doi.org/10.1038/s41586-020-2308-7)
- MEME Suite: [Bailey et al., 2015](https://doi.org/10.1093/nar/gkv416)
- JASPAR: [Castro-Mondragon et al., 2022](https://doi.org/10.1093/nar/gkab1113)
- AlphaGenome: [your AlphaGenome citation]

---

## Author

**George Stephenson**  
TF Activation Landscape Project  
CU Boulder LAYER Lab

---

## Version

**Current Status:** Module 03 Complete (partial gnomAD coverage)  
**Next:** Module 04 - AlphaGenome Functional Scoring

**Latest Update:** Module 03 completed with 326K variants from 14 chromosomes. Large chromosomes (chr1-8, chr10, chr12) timed out and would require chunked queries for complete coverage.
