# Module 00: Data Fetching

This module handles downloading and verification of all required data sources for the Dormant Site Activation Pipeline.

## Overview

Module 00 ensures that all necessary data files are present and properly configured:
- GRCh38 human reference genome
- Transcription factor motif models from JASPAR
- gnomAD v4.1 population variant data

## Scripts

### `download_reference.sh`
Downloads the GRCh38 reference genome from Ensembl.

**Usage:**
```bash
cd 00_fetch_data
bash download_reference.sh
```

**What it does:**
1. Downloads `Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz` from Ensembl
2. Decompresses the genome
3. Creates FASTA index with `samtools faidx`
4. Places files in `data/reference/`

**Output:**
- `data/reference/GRCh38.fa` - Reference genome
- `data/reference/GRCh38.fa.fai` - FASTA index

**Requirements:**
- `wget` - Download utility
- `samtools` (optional) - For FASTA indexing

---

### `download_motifs.py`
Downloads transcription factor motif matrices from JASPAR database.

**Usage:**
```bash
# Download AP1 motif (default)
python download_motifs.py

# Download different motif
python download_motifs.py --motif-id MA0106.3  # JUN

# Specify output directory
python download_motifs.py --motif-id MA0099.3 --output-dir ../data/motifs/JASPAR
```

**Arguments:**
- `--motif-id` - JASPAR motif ID (default: MA0099.3 for AP1)
- `--output-dir` - Output directory (default: ../data/motifs/JASPAR)
- `--log-file` - Log file path

**What it does:**
1. Queries JASPAR API for motif information
2. Downloads Position Frequency Matrix (PFM)
3. Converts to MEME format (compatible with FIMO)
4. Saves raw JSON metadata
5. Creates `motif_metadata.yaml` for pipeline

**Output:**
- `data/motifs/JASPAR/MA0099.3.meme` - Motif in MEME format
- `data/motifs/JASPAR/MA0099.3.json` - Raw motif data
- `data/motifs/JASPAR/motif_metadata.yaml` - Pipeline metadata

**Requirements:**
- `requests` - HTTP library for API calls

---

### `download_gnomad_index_info.py`
Verifies gnomAD installation and creates symlinks.

**Usage:**
```bash
# Show download instructions
python download_gnomad_index_info.py --show-download-instructions

# Verify installation
python download_gnomad_index_info.py --config ../pipeline_config.yaml

# Create symlink to gnomAD data
python download_gnomad_index_info.py --config ../pipeline_config.yaml --create-symlink
```

**Arguments:**
- `--config` - Path to pipeline configuration (default: ../pipeline_config.yaml)
- `--create-symlink` - Create symlink from data/gnomad to gnomAD root
- `--show-download-instructions` - Display download instructions
- `--log-file` - Log file path

**What it does:**
1. Reads gnomAD paths from configuration
2. Verifies that required files exist:
   - `genomes.vcf.bgz` - Main VCF file (~500 GB)
   - `genomes.vcf.bgz.tbi` - Tabix index
3. Checks optional files:
   - Coverage data
   - Constraint scores
4. Creates symbolic link from `data/gnomad/` to `/mnt/data_1/gnomAD_data/raw/gnomad_v4.1`

**Requirements:**
- gnomAD data already downloaded to `/mnt/data_1/gnomAD_data/`

---

## gnomAD Download Instructions

The gnomAD v4.1 dataset is large (~750 GB total). Follow these steps:

### 1. Create directory structure
```bash
sudo mkdir -p /mnt/data_1/gnomAD_data/raw/gnomad_v4.1
sudo chown $USER:$USER /mnt/data_1/gnomAD_data
cd /mnt/data_1/gnomAD_data/raw/gnomad_v4.1
```

### 2. Download required files (~500 GB)
```bash
mkdir -p genomes_vcf
gsutil cp gs://gcp-public-data--gnomad/release/4.1/genomes_vcf/genomes.vcf.bgz genomes_vcf/
gsutil cp gs://gcp-public-data--gnomad/release/4.1/genomes_vcf/genomes.vcf.bgz.tbi genomes_vcf/
```

### 3. Download optional coverage data (~200 GB, recommended)
```bash
mkdir -p coverage
gsutil -m cp -r gs://gcp-public-data--gnomad/release/4.1/coverage/genomes/* coverage/
```

### 4. Download optional constraint scores (~50 GB)
```bash
mkdir -p constraint
gsutil cp gs://gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint.json.bgz constraint/
```

### 5. Verify installation
```bash
cd /path/to/pipeline
python 00_fetch_data/download_gnomad_index_info.py --config pipeline_config.yaml
```

---

## Complete Setup Workflow

Run all data fetching steps in order:

```bash
# 1. Download reference genome
cd 00_fetch_data
bash download_reference.sh

# 2. Download motifs
python download_motifs.py --motif-id MA0099.3

# 3. Verify gnomAD (assumes already downloaded)
python download_gnomad_index_info.py --show-download-instructions
# ... follow instructions to download gnomAD ...

# 4. Once downloaded, verify and create symlink
python download_gnomad_index_info.py --config ../pipeline_config.yaml --create-symlink
```

---

## Directory Structure After Setup

```
data/
├── reference/
│   ├── GRCh38.fa              # Reference genome (3.1 GB)
│   └── GRCh38.fa.fai          # FASTA index
├── motifs/
│   └── JASPAR/
│       ├── MA0099.3.meme      # AP1 motif in MEME format
│       ├── MA0099.3.json      # Raw motif data
│       └── motif_metadata.yaml # Pipeline metadata
└── gnomad/                     # Symlink to /mnt/data_1/gnomAD_data/raw/gnomad_v4.1
```

---

## Troubleshooting

### Reference genome download fails
- Check internet connection
- Try alternative mirror: `rsync -av rsync://ftp.ensembl.org/ensembl/pub/release-110/fasta/homo_sapiens/dna/`

### JASPAR API timeout
- Check internet connection
- Try again later (occasional API issues)
- Manually download from https://jaspar.elixir.no

### gnomAD download is slow
- Use `gsutil -m` for parallel downloads
- Download during off-peak hours
- Consider using a cloud instance in the same region as the data

### Disk space issues
- gnomAD requires ~750 GB
- Check available space: `df -h /mnt/data_1`
- Consider downloading only required files (skip optional coverage/constraint)

---

## Dependencies

Python packages:
- `pyyaml` - Configuration file handling
- `requests` - HTTP requests for JASPAR API

System tools:
- `wget` - Download utility
- `samtools` - FASTA indexing
- `gsutil` - Google Cloud Storage download (for gnomAD)

Install with:
```bash
# Python packages
pip install pyyaml requests

# System tools (Ubuntu/Debian)
sudo apt-get install wget samtools

# gsutil (Google Cloud SDK)
# See: https://cloud.google.com/storage/docs/gsutil_install
```

---

## Next Steps

After completing Module 00, proceed to:
- **Module 01**: Genome-wide motif scanning
