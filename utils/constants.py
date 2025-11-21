"""
Constants for Dormant Site Activation Pipeline
"""

import os

# Base directories
PIPELINE_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(PIPELINE_ROOT, "data")
RESULTS_DIR = os.path.join(PIPELINE_ROOT, "results")
FIGURES_DIR = os.path.join(PIPELINE_ROOT, "figures")
LOGS_DIR = os.path.join(PIPELINE_ROOT, "logs")

# Data subdirectories
REFERENCE_DIR = os.path.join(DATA_DIR, "reference")
MOTIFS_DIR = os.path.join(DATA_DIR, "motifs")
GNOMAD_DIR = os.path.join(DATA_DIR, "gnomad")

# Standard chromosomes (no patches, alts, or unplaced contigs)
STANDARD_CHROMOSOMES = [
    f"chr{i}" for i in range(1, 23)
] + ["chrX", "chrY"]

# Alternative chromosome naming (without 'chr' prefix)
STANDARD_CHROMOSOMES_NO_CHR = [
    str(i) for i in range(1, 23)
] + ["X", "Y"]

# DNA bases
BASES = ['A', 'C', 'G', 'T']
COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}

# Motif tiers
MOTIF_TIERS = {
    0: "strong",
    1: "medium", 
    2: "weak",
    3: "very_weak"
}

# gnomAD VCF field names
GNOMAD_FIELDS = {
    'AF': 'Allele Frequency',
    'AC': 'Allele Count',
    'AN': 'Allele Number',
    'nhomalt': 'Number of Homozygous Alternate',
    'AF_afr': 'African/African American AF',
    'AF_amr': 'Latino/Admixed American AF',
    'AF_asj': 'Ashkenazi Jewish AF',
    'AF_eas': 'East Asian AF',
    'AF_fin': 'Finnish AF',
    'AF_nfe': 'Non-Finnish European AF',
    'AF_oth': 'Other AF',
    'AF_sas': 'South Asian AF'
}

# AlphaGenome tracks of interest
ALPHAGENOME_TRACKS = [
    'ATAC',
    'H3K27ac',
    'H3K4me3',
    'H3K4me1',
    'H3K9me3',
    'H3K27me3',
    'H3K36me3',
    'CTCF',
    'TF_binding'
]

# File extensions
FASTA_EXTS = ['.fa', '.fasta', '.fna']
BED_EXTS = ['.bed', '.bed.gz']
VCF_EXTS = ['.vcf', '.vcf.gz', '.vcf.bgz']
TSV_EXTS = ['.tsv', '.tsv.gz']

# Default parameters
DEFAULT_CONTEXT_BP = 1000
DEFAULT_PVALUE_THRESHOLD = 1e-3
DEFAULT_MAX_HAMMING_DIST = 3
DEFAULT_MIN_COVERAGE = 20

# JASPAR database URLs
JASPAR_BASE_URL = "https://jaspar.elixir.no/api/v1"
JASPAR_MATRIX_URL = f"{JASPAR_BASE_URL}/matrix"

# HOCOMOCO database URLs
HOCOMOCO_BASE_URL = "https://hocomoco11.autosome.org"

# Color schemes for visualization
TIER_COLORS = {
    0: "#1f77b4",  # strong - blue
    1: "#ff7f0e",  # medium - orange
    2: "#2ca02c",  # weak - green
    3: "#d62728"   # very weak - red
}

QUADRANT_COLORS = {
    'high_impact_accessible': "#e74c3c",      # red
    'high_impact_constrained': "#f39c12",     # orange
    'low_impact_accessible': "#3498db",       # blue
    'low_impact_constrained': "#95a5a6"       # gray
}

# Memory limits
MAX_MEMORY_GB = 64
MAX_SEQUENCES_IN_MEMORY = 100000
