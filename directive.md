🌋 directive.md – Dormant Site Activation Pipeline (Full Specification)

Author: George Stephenson
Project: TF Activation Landscape (AP1 prototype, generalizable to all TFs)
Purpose: Build a fully automated, reproducible pipeline that quantifies how “dormant” motif instances across the human genome could become activated TF binding sites through human-standing variation and predicted molecular consequences (AlphaGenome).
Primary Data Sources: gnomAD v4.1 (GRCh38), JASPAR/HOCOMOCO motifs, GRCh38 genome, AlphaGenome model.

1. 📌 Overall Objective

Implement a generalizable computational pipeline that:

Takes as input any transcription factor (TF) (starting with AP1).

Retrieves its motif model (PWM/PFM).

Scans the entire human genome for all motif-like sequences (strong + weak + near-motif).

Computes mutation paths from each motif instance to the consensus motif.

Intersects each mutation step with gnomAD v4.1 variants to determine whether each activating mutation exists in the population and at what frequency.

Scores the functional impact of activating variants using AlphaGenome (REF vs ALT allele sequences).

Places each genomic site into a 2D landscape of:

X-axis: Population accessibility / selection constraint

Y-axis: Functional impact (ΔAlphaGenome)

Outputs ranked lists of “dormant but potentially functional” motif sites.

Produces publication-ready visualizations and summary reports.

This pipeline must run on ODYSSEUS or any HPC, be modular, and support batch-mode for large models like AlphaGenome.

2. 📁 Directory Structure (Copilot: Create These Exactly)
alphamotif-activation-pipeline/
│
├── directive.md                      # THIS FILE
├── pipeline_config.yaml              # Central config file (to be auto-generated)
│
├── data/
│   ├── reference/                    # GRCh38 genome + indexes
│   ├── motifs/                       # PWMs/PFMs from JASPAR, HOCOMOCO
│   └── gnomad/                       # symlink or pointer to /mnt/data_1/gnomAD_data
│
├── 00_fetch_data/
│   ├── download_reference.sh
│   ├── download_motifs.py
│   ├── download_gnomad_index_info.py
│   └── DOCUMENTATION.md
│
├── 01_scan_motifs/
│   ├── scan_genome_fimo.py
│   ├── tier_sites.py
│   └── DOCUMENTATION.md
│
├── 02_generate_mutation_paths/
│   ├── consensus_from_pwm.py
│   ├── compute_distance.py
│   ├── enumerate_paths.py
│   └── DOCUMENTATION.md
│
├── 03_intersect_gnomad/
│   ├── make_motif_bed.py
│   ├── query_gnomad_vcfs.py
│   ├── summarize_path_AF.py
│   └── DOCUMENTATION.md
│
├── 04_run_alphagenome/
│   ├── make_variant_seqs.py
│   ├── run_alphagenome_batch.py
│   ├── compute_functional_deltas.py
│   └── DOCUMENTATION.md
│
├── 05_compute_activation_landscape/
│   ├── combine_population_and_impact.py
│   ├── classify_quadrants.py
│   └── DOCUMENTATION.md
│
├── 06_visualization/
│   ├── plot_landscape.py
│   ├── plot_genome_tracks.py
│   └── DOCUMENTATION.md
│
└── utils/
    ├── io.py
    ├── logging_utils.py
    ├── constants.py
    └── DOCUMENTATION.md

3. 🌍 Data Requirements (Copilot: Ensure these download commands are used)

All gnomAD data must reside at:

/mnt/data_1/gnomAD_data/raw/gnomad_v4.1/


Copilot must create a symlink inside data/gnomad:

ln -s /mnt/data_1/gnomAD_data/raw/gnomad_v4.1 data/gnomad

3.1 Download gnomAD v4.1 (GRCh38)

Required

mkdir -p /mnt/data_1/gnomAD_data/raw/gnomad_v4.1/genomes_vcf
gsutil cp gs://gcp-public-data--gnomad/release/4.1/genomes_vcf/genomes.vcf.bgz \
    /mnt/data_1/gnomAD_data/raw/gnomad_v4.1/genomes_vcf/
gsutil cp gs://gcp-public-data--gnomad/release/4.1/genomes_vcf/genomes.vcf.bgz.tbi \
    /mnt/data_1/gnomAD_data/raw/gnomad_v4.1/genomes_vcf/


Recommended (for coverage masking)

mkdir -p /mnt/data_1/gnomAD_data/raw/gnomad_v4.1/coverage
gsutil -m cp -r gs://gcp-public-data--gnomad/release/4.1/coverage/genomes \
    /mnt/data_1/gnomAD_data/raw/gnomad_v4.1/coverage/


Optional (constraint scores)

mkdir -p /mnt/data_1/gnomAD_data/raw/gnomad_v4.1/constraint
gsutil cp gs://gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint.json.bgz \
    /mnt/data_1/gnomAD_data/raw/gnomad_v4.1/constraint/

3.2 GRCh38 Reference

Copilot should implement:

wget ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa


Place in:

data/reference/GRCh38.fa

3.3 Motifs

Copilot should pull:

JASPAR 2024 PFMs automatically

specifically AP1 (JASPAR ID MA0099.3)

Place PFMs in:

data/motifs/JASPAR/

4. 🧬 Pipeline Modules (HIGH-LEVEL SPECIFICATION FOR COPILOT)
Module 01 – Genome-Wide Motif Scan

Input:

PWM/PFM

GRCh38 genome

Copilot must:

Run FIMO or MOODS to detect ALL motif-like sites.

Store results in:

results/motif_scan/{TF}/motif_hits.bed


Tier sites by PWM score into tiers 0–3.

Module 02 – Mutation Path Enumeration

For each motif hit:

Infer consensus motif from PWM.

Compute mismatch positions.

Enumerate all minimal mutation paths (distance ≤ 3).

Output:

results/mutation_paths/{TF}/paths.tsv


Columns:

| chr | start | end | site_seq | consensus | hamming_dist | mutation_path | ref_base | alt_base | pos_in_motif |

Module 03 – gnomAD Integration

Copilot must:

Convert motif paths to genomic coordinates.

Query gnomAD VCF by region using:

bcftools view -R motif_windows.bed /mnt/data_1/gnomAD_data/raw/gnomad_v4.1/genomes_vcf/genomes.vcf.bgz


Extract:

AF

AC

AN

nhomalt

Match only those variants where ALT matches the mutation path.

Output:

results/gnomad_intersection/{TF}/gnomad_hits.tsv

Module 04 – AlphaGenome Functional Scoring

Copilot must:

For each motif hit + variant:

Generate reference sequence window (±1kb)

Generate mutated sequence window

Run AlphaGenome on each:

alphagenome predict --model <PATH> \
    --fasta seqs.fa --output predictions.json


Compute:

Δ ATAC

Δ TF binding

Δ enhancer marks

Output:

results/alphagenome/{TF}/functional_deltas.tsv

Module 05 – Activation Landscape

Copilot must compute:

X-axis (population accessibility / constraint)

Formula (conceptual):

X = -log10(max(AF_step_i) + 1e-12) * (hamming_dist)

Y-axis (functional impact)
Y = max(ΔAlphaGenome_feature)


Output:

results/landscape/{TF}/activation_landscape.tsv

Module 06 – Visualization

Copilot must generate:

2D scatter/hexbin of X vs Y

Genome browser-style tracks

Ranked candidate list of top dormant → high-impact sites

Output:

figures/{TF}/activation_landscape.png

5. 🎯 Pipeline Configuration Template (pipeline_config.yaml)

Copilot must generate:

tf_name: "AP1"
motif_id: "MA0099.3"

reference_genome: "data/reference/GRCh38.fa"

gnomad:
  vcf: "/mnt/data_1/gnomAD_data/raw/gnomad_v4.1/genomes_vcf/genomes.vcf.bgz"
  vcf_tbi: "/mnt/data_1/gnomAD_data/raw/gnomad_v4.1/genomes_vcf/genomes.vcf.bgz.tbi"
  coverage: "/mnt/data_1/gnomAD_data/raw/gnomad_v4.1/coverage/genomes"
  constraint: "/mnt/data_1/gnomAD_data/raw/gnomad_v4.1/constraint/gnomad.v4.1.constraint.json.bgz"

alphagenome:
  model_path: "/mnt/models/alphagenome/"
  context_bp: 1000

motif_scan:
  pvalue_threshold: 1e-3
  tier_cutoffs:
    tier0: 0.95
    tier1: 0.85
    tier2: 0.70

6. ⚡ Copilot Implementation Directives

Copilot must:

Implement each module as standalone Python scripts.

Include documentation files in each module folder.

Allow the entire pipeline to run via:

python run_pipeline.py --config pipeline_config.yaml


Support Snakemake or Nextflow auto-generated workflows.

Create unit tests for:

motif scanning

mutation path enumeration

gnomAD querying

AlphaGenome sequencing I/O

All intermediate outputs must be stored under results/.

7. ✔ Success Criteria

Pipeline must:

Successfully process AP1 end-to-end

Produce activation landscape for AP1

Correctly rank dormant → high-impact motif instances

Generalize to any TF simply by changing TF name and motif ID