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

### Module 04: AlphaGenome Scoring [Complete]
**Score functional impact using AlphaGenome with multi-modal feature capture (expression, ATAC, DNase, TF binding, histones, Hi-C).**

**Status:** ✅ Complete (6,810 variants scored successfully, 95.1% success rate)

**Implementation:** Follows proven methodology from `alphagenome-qtl-validation` repository
- Uses `dna_client.SEQUENCE_LENGTH_1MB` (1,048,576 bp) context windows
- Direct API usage via `client.score_variant()` 
- **ALL 19 recommended variant scorers** for comprehensive feature coverage
- Captures: RNA_SEQ, ATAC, DNASE, CHIP_TF, CHIP_HISTONE, CONTACT_MAPS, CAGE, PROCAP, SPLICE_*
- Average 27,415 tracks per variant (context-dependent, biologically accurate)

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

**Results:**
- **186,695,928 variant-track predictions** (186.7M)
- **6,810 unique variants scored** (95.1% success rate)
- **348 variants failed** (4.9% - transient network errors during scoring)
- **27,415 average tracks per variant** (range: 31K-82K depending on genomic context)
- **All 11 output types present:** RNA_SEQ, CHIP_TF, CHIP_HISTONE, SPLICE_JUNCTIONS, CAGE, DNASE, ATAC, SPLICE_SITE_USAGE, CONTACT_MAPS, PROCAP, SPLICE_SITES
- **714 unique biosamples/cell types** captured
- **Runtime:** 3.5 hours (21:26 → 01:26)

**Data Quality:**
- ✅ 100% metadata completeness (gnomAD AF, path_id, step_num)
- ✅ No null scores (186.7M valid predictions)
- ✅ Score distribution normal (mean=0.27, std=0.50, range=-1 to 1)
- ✅ Context-appropriate track counts (RNA_SEQ only includes genes within 1MB window)

**Why track counts vary by variant:**
- RNA_SEQ predictions are gene-specific (only genes within 1MB window)
- Splice tracks only appear near splice junctions
- Cell type availability varies by assay and genomic region
- This variation is **expected and biologically correct**

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

**Output files:**
- `results/alphagenome/AP1/predictions.parquet` - Raw multi-modal scores (186.7M variant-track pairs, 1.36 GB)
- `results/alphagenome/AP1/predictions_summary.tsv` - Overall mean scores per variant (6,810 variants, 506 KB)
- `results/alphagenome/AP1/predictions_summary_by_feature.tsv` - Per-feature means (63,471 variant-feature pairs, 4.2 MB)

**Performance achieved:** 
- Multi-modal (19 scorers): ~1.77 sec/variant actual
- Runtime: 3.5 hours for 6,810 variants
- Memory usage: ~214 GB peak (in-memory processing)
- File size: 1.36 GB parquet (excellent compression)

**Failed variants (348/7,158 = 4.9%):**
- Distribution: chr18 (64), chr19 (138), chr2 (146)
- Cause: Transient network errors (timeouts, connection failures) during 1-hour window
- Impact: Minimal - 95.1% success rate provides robust coverage
- Can be re-scored later if needed (failed positions saved)

**Next Step:** Module 05 - Compute Activation Landscape

---

### Module 05: Activation Landscape [Complete]
**Combine gnomAD constraint with AlphaGenome predictions to create 2D activation landscape.**

**Status:** ✅ Complete  
**Prerequisites:** ✅ Module 03 (gnomAD), ✅ Module 04 (AlphaGenome)  
**Documentation:** `05_compute_activation_landscape/README.md`  
**Scientific Report:** `05_compute_activation_landscape/SCIENTIFIC_REPORT.md`

**Key Innovation - Biologically-Specific Y-Axis:**
Instead of using `max()` across all AlphaGenome tracks (which conflates unrelated signals), we use **AP1-family-specific TF predictions**:
```
Y = max(quantile_score) for TF ∈ {JUND, JUN, JUNB, FOS, FOSL1, FOSL2, ATF3, ATF2, ATF7, BATF, MAFK}
```

**Landscape Coordinates:**
- **X-axis (Accessibility):** `-log10(AF) × Hamming_distance`
  - Low X = accessible (common variant, few mutations)
  - High X = hard to access (rare variant, many mutations)
- **Y-axis (Impact):** Max AP1-family TF binding quantile score
  - Measures direct biological effect on AP1 binding

**Key Results:**

| Metric | Value |
|--------|-------|
| **Total variants** | 6,810 |
| **High-priority candidates** | 1,624 (23.8%) |
| **Strong AP1 gain (>90th %ile)** | 6,170 (90.6%) |
| **Ultra-rare (AF < 0.01%)** | 5,798 (85.1%) |
| **AP1 vs Enhancer correlation** | r = 0.545, p < 10⁻¹⁶ |

**Quadrant Distribution:**
- HIGH PRIORITY (Accessible + High Impact): 1,624 (23.8%)
- High Impact, Hard to Access: 1,779 (26.1%)
- Accessible, Low Impact: 1,781 (26.2%)
- Low Priority: 1,626 (23.9%)

**Top AP1-Family TFs:**
- FOS: 25.8%
- JUND: 16.8%
- ATF3: 8.8%
- FOSL2: 8.5%
- ATF2: 8.1%

**Validation:** Strong correlation (r=0.545) between AP1 binding impact and enhancer marks (H3K27ac/H3K4me1) confirms biological relevance.

**Usage:**
```bash
# Compute landscape (vectorized, ~30 seconds)
conda run -n alphagenome-env python 05_compute_activation_landscape/compute_activation_landscape.py \
  --predictions results/alphagenome/AP1/predictions.parquet \
  --gnomad results/gnomad_intersection/AP1/all_observed_variants.tsv \
  --output results/landscape/AP1

# Generate figures
conda run -n alphagenome-env python 05_compute_activation_landscape/plot_activation_landscape.py \
  --input results/landscape/AP1/AP1_activation_landscape.tsv \
  --output-dir figures/landscape
```

**Output Files:**
- `results/landscape/AP1/AP1_activation_landscape.tsv` - Complete landscape (6,810 variants)
- `results/landscape/AP1/AP1_high_priority_candidates.tsv` - High-priority quadrant (1,624 variants)
- `results/landscape/AP1/AP1_landscape_summary.txt` - Summary statistics
- `figures/landscape/AP1_activation_landscape_main.png` - Main visualization
- `figures/landscape/AP1_activation_landscape_comparison.png` - Multi-panel comparison
- `figures/landscape/AP1_tf_breakdown.png` - TF contribution analysis
- `figures/landscape/AP1_high_priority_candidates.png` - High-priority details

**Performance:** ~35 seconds total (vectorized implementation processes 186M predictions efficiently)

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
  batch_size: 32
```

---

## Expected Outputs

### Activation Landscape
A 2D plot showing:
- **X-axis**: Population accessibility × constraint  
  Formula: `-log10(AF + 1e-12) × Hamming_distance`
  - Combines allele frequency rarity with mutational distance from reference
  - Higher values = rarer variants requiring more mutations
- **Y-axis**: AP1-specific functional impact  
  Formula: `max(quantile_score)` across 11 AP1-family TFs
  - JUND, JUN, JUNB, FOS, FOSL1, FOSL2, ATF3, ATF2, ATF7, BATF, MAFK
  - Uses quantile normalization for cross-track comparability
  - Values 0-1, with >0.9 indicating top 10% predicted binding

### High-Priority Candidates
Variants in the upper-right quadrant:
- Strong AP1 binding gain (>90th percentile)
- Ultra-rare or absent in population (AF < 0.0001)
- Represent most promising dormant sites for experimental validation

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

**Current Status:** Module 05 Complete (Activation Landscape)  
**Next:** Module 07 - Population Statistics (can run independently)

**Latest Update (Module 05 - November 27, 2025):**
- ✅ Activation landscape computed for 6,810 variants
- ✅ 1,624 high-priority candidates identified (23.8%)
- ✅ 90.6% variants show strong AP1 binding gain (>90th percentile)
- ✅ 4 publication-quality figures generated
- ✅ Full scientific report available: `05_compute_activation_landscape/SCIENTIFIC_REPORT.md`

**Key Results:**
- X-axis formula: `-log10(AF) × Hamming_distance` (population accessibility)
- Y-axis formula: `max(AP1 quantile)` across 11 AP1-family TFs
- Top contributing TFs: FOS (25.8%), JUND (16.8%), ATF3 (8.8%)
- Strong correlation between AP1 binding and enhancer marks (r=0.545)
