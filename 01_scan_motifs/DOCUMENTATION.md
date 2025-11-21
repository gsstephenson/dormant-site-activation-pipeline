# Module 01: Motif Scanning

This module performs genome-wide scanning for transcription factor binding motifs and classifies matches into tiers based on PWM score strength.

## Overview

Module 01 identifies all potential binding sites for a transcription factor across the human genome, including:
- **Strong matches** (tier 0): Near-perfect motif matches
- **Medium matches** (tier 1): Good quality matches
- **Weak matches** (tier 2): Lower quality but significant matches  
- **Very weak matches** (tier 3): Borderline matches that could become active with mutations

## Scripts

### `scan_genome_fimo.py`

Scans the entire genome for motif matches using FIMO (Find Individual Motif Occurrences) from the MEME Suite.

**Usage:**
```bash
# Basic usage
python 01_scan_motifs/scan_genome_fimo.py --config pipeline_config.yaml

# Override P-value threshold
python 01_scan_motifs/scan_genome_fimo.py --config pipeline_config.yaml --pvalue 1e-4

# Use custom motif file
python 01_scan_motifs/scan_genome_fimo.py \
    --config pipeline_config.yaml \
    --motif-file data/motifs/custom_motif.meme
```

**Arguments:**
- `--config` - Path to pipeline configuration file (required)
- `--motif-file` - Override motif file from config
- `--output-dir` - Override output directory
- `--pvalue` - P-value threshold (default from config: 1e-3)
- `--log-file` - Custom log file path

**What it does:**
1. Loads motif PWM from MEME format file
2. Runs FIMO against GRCh38 reference genome
3. Filters matches by P-value threshold
4. Converts output to BED format
5. Saves full results as TSV

**Output files:**
- `results/motif_scan/{TF}/fimo_output.tsv` - Raw FIMO output
- `results/motif_scan/{TF}/motif_hits.bed` - BED format (8 columns)
- `results/motif_scan/{TF}/motif_hits_full.tsv` - Full results with all annotations

**BED format columns:**
1. chr - Chromosome
2. start - Start position (0-based)
3. end - End position
4. name - Site identifier
5. score - PWM score
6. strand - Strand (+/-)
7. sequence - Matched sequence
8. pvalue - P-value

---

### `tier_sites.py`

Classifies motif sites into tiers based on PWM score percentiles.

**Usage:**
```bash
# Basic usage
python 01_scan_motifs/tier_sites.py --config pipeline_config.yaml

# Use custom input file
python 01_scan_motifs/tier_sites.py \
    --config pipeline_config.yaml \
    --input-file results/motif_scan/AP1/motif_hits_full.tsv
```

**Arguments:**
- `--config` - Path to pipeline configuration file (required)
- `--input-file` - Override input file (default: from Module 01 scan output)
- `--output-dir` - Override output directory
- `--log-file` - Custom log file path

**What it does:**
1. Loads motif hits from scanning step
2. Calculates score percentiles for all sites
3. Assigns tier based on configurable cutoffs:
   - Tier 0: ≥95th percentile (strong)
   - Tier 1: 85-95th percentile (medium)
   - Tier 2: 70-85th percentile (weak)
   - Tier 3: <70th percentile (very weak)
4. Generates tier-specific BED files
5. Creates summary statistics and report

**Output files:**
- `results/motif_scan/{TF}/motif_hits_tiered.tsv` - All sites with tier assignments
- `results/motif_scan/{TF}/motif_hits_tier0_strong.bed` - Tier 0 sites only
- `results/motif_scan/{TF}/motif_hits_tier1_medium.bed` - Tier 1 sites only
- `results/motif_scan/{TF}/motif_hits_tier2_weak.bed` - Tier 2 sites only
- `results/motif_scan/{TF}/motif_hits_tier3_very_weak.bed` - Tier 3 sites only
- `results/motif_scan/{TF}/tier_summary.json` - Summary statistics (JSON)
- `results/motif_scan/{TF}/tier_report.txt` - Human-readable report

---

## Configuration

Relevant section in `pipeline_config.yaml`:

```yaml
motif_scan:
  tool: "fimo"
  pvalue_threshold: 1e-3
  background: "uniform"
  tier_cutoffs:
    tier0: 0.95  # strong matches (≥95th percentile)
    tier1: 0.85  # medium matches
    tier2: 0.70  # weak matches
    tier3: 0.50  # very weak (not used, implicit <70th)
```

---

## Complete Module 01 Workflow

Run both scripts in sequence:

```bash
# Step 1: Scan genome for motifs
python 01_scan_motifs/scan_genome_fimo.py --config pipeline_config.yaml

# Step 2: Tier the motif sites
python 01_scan_motifs/tier_sites.py --config pipeline_config.yaml
```

**Expected output structure:**
```
results/motif_scan/AP1/
├── fimo_output.tsv              # Raw FIMO output
├── motif_hits.bed               # All sites in BED format
├── motif_hits_full.tsv          # Full annotations
├── motif_hits_tiered.tsv        # All sites with tiers
├── motif_hits_tier0_strong.bed  # Tier-specific BED files
├── motif_hits_tier1_medium.bed
├── motif_hits_tier2_weak.bed
├── motif_hits_tier3_very_weak.bed
├── tier_summary.json            # Statistics
└── tier_report.txt              # Human-readable report
```

---

## Understanding Tiers

The tiering system categorizes motif matches by strength:

**Tier 0 (Strong)**: Top 5% of matches
- These are canonical binding sites
- Likely already functional
- High PWM scores, close to consensus
- Example: `TGAGTCA` (perfect AP1 consensus)

**Tier 1 (Medium)**: 85-95th percentile
- Strong matches with minor deviations
- Likely functional but with reduced affinity
- Example: `TGAGTCC` (1 mismatch from consensus)

**Tier 2 (Weak)**: 70-85th percentile  
- Weaker matches, may be functional in some contexts
- Could be cryptic/dormant sites
- Example: `TGAATCA` (2 mismatches)

**Tier 3 (Very Weak)**: Below 70th percentile
- **Primary focus for dormant site analysis**
- Multiple mismatches from consensus
- Likely non-functional but could be activated by mutations
- Example: `TGAATCC` (3 mismatches)

---

## Performance Notes

**Runtime:**
- Genome scanning: ~2-4 hours for whole genome (depends on motif)
- Tiering: ~5-10 minutes

**Memory:**
- FIMO: ~4-8 GB RAM
- Tiering: ~2 GB RAM

**Disk space:**
- FIMO output: ~100-500 MB (depends on P-value threshold)
- Tiered files: ~200-1000 MB total

**Optimization tips:**
- Use stricter P-value for faster scanning (but may miss weak sites)
- Run on specific chromosomes for testing
- Use background model for more accurate P-values

---

## Troubleshooting

### FIMO not found
```
Error: FIMO not found
```
**Solution:** Install MEME Suite:
```bash
conda install -c bioconda meme
# or
sudo apt-get install meme-suite
```

### Out of memory during scanning
**Solution:** 
- Scan chromosomes individually
- Increase memory allocation
- Use stricter P-value threshold

### No motif hits found
**Possible causes:**
1. P-value too strict → Relax threshold
2. Wrong motif file format → Check MEME format
3. Wrong reference genome → Verify GRCh38

### Very few tier 3 sites
**Explanation:** This is expected! Tier 3 represents very weak matches. If you need more weak sites:
- Relax P-value threshold in scanning
- Adjust tier cutoffs in config

---

## Dependencies

**Required:**
- MEME Suite (for FIMO)
- Python packages: pandas, numpy

**Installation:**
```bash
# MEME Suite
conda install -c bioconda meme

# Python packages (if not already installed)
pip install pandas numpy
```

**Verify installation:**
```bash
fimo --version  # Should show MEME Suite version
```

---

## Next Steps

After completing Module 01:
1. Review tier distribution in `tier_report.txt`
2. Examine tier 3 sites (dormant candidates)
3. Proceed to **Module 02** to generate mutation paths
4. Focus on tier 2-3 sites for activation analysis

**Next command:**
```bash
python 02_generate_mutation_paths/enumerate_paths.py --config pipeline_config.yaml
```
