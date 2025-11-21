# Dormant Site Activation Pipeline

**A computational framework for identifying cryptic transcription factor binding sites that could be activated by naturally occurring genetic variation.**

## 🎯 Overview

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

## 📋 Quick Start

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

# Module 02: Generate mutation paths (coming next)
# python 02_generate_mutation_paths/enumerate_paths.py --config pipeline_config.yaml

# ... (additional modules to be implemented)
```

---

## 📁 Directory Structure

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
├── 01_scan_motifs/             # ✅ Genome-wide motif scanning
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

## 🧬 Pipeline Modules

### ✅ Module 00: Data Fetching
Download reference genome, motifs, and verify gnomAD installation.

**Status:** Complete  
**Documentation:** `00_fetch_data/DOCUMENTATION.md`

---

### ✅ Module 01: Motif Scanning
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

### 🚧 Module 02: Mutation Paths
*Coming next - Enumerate minimal mutation paths to consensus motif*

---

### 🚧 Module 03: gnomAD Intersection
*Query population variants that match activating mutations*

---

### 🚧 Module 04: AlphaGenome Scoring
*Predict functional impact of activating variants*

---

### 🚧 Module 05: Activation Landscape
*Compute 2D landscape coordinates*

---

### 🚧 Module 06: Visualization
*Generate plots and ranked candidate lists*

---

## ⚙️ Configuration

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
  context_bp: 1000
```

---

## 📊 Expected Outputs

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

## 🔬 Use Cases

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

## 📖 Documentation

Each module has detailed documentation:
- `00_fetch_data/DOCUMENTATION.md` - Data acquisition
- `01_scan_motifs/DOCUMENTATION.md` - Motif scanning
- `utils/DOCUMENTATION.md` - Utility functions

---

## 🐛 Troubleshooting

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

## 📝 Citation

If you use this pipeline, please cite:
- gnomAD: [Karczewski et al., 2020](https://doi.org/10.1038/s41586-020-2308-7)
- MEME Suite: [Bailey et al., 2015](https://doi.org/10.1093/nar/gkv416)
- JASPAR: [Castro-Mondragon et al., 2022](https://doi.org/10.1093/nar/gkab1113)
- AlphaGenome: [your AlphaGenome citation]

---

## 👤 Author

**George Stephenson**  
TF Activation Landscape Project  
CU Boulder LAYER Lab

---

## 🔄 Version

**Current Status:** Module 01 Complete  
**Next:** Module 02 - Mutation Path Enumeration
