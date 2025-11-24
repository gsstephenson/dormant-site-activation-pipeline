Here’s an updated, **Copilot-ready** `directive.md` you can drop straight into your repo (replace the old one). It’s written assuming the existing structure, but tells Copilot to *modify* files where they already exist rather than naively rewrite everything.

You can copy–paste this whole thing into `directive.md` at the repo root.

---

````markdown
# Dormant Site Activation Pipeline – Updated Directive

**Author:** George Stephenson  
**Project:** TF Activation Landscape (AP1 prototype, generalizable to all TFs)  

**Purpose:**  
Build a fully automated, reproducible pipeline that quantifies how “dormant” motif instances across the human genome could become activated TF binding sites through human-standing variation and predicted molecular consequences (AlphaGenome). The pipeline must:

- Integrate **gnomAD v4.1** allele frequencies (including **AF = 0** sites).
- Use **AlphaGenome** to compute signed functional impacts (Δ, including negatives).
- Construct a **2D selection–impact landscape**:
  - **X-axis:** population accessibility / constraint, combining allele frequency and mutational steps.
  - **Y-axis:** functional impact (signed ΔAlphaGenome).

This directive **updates and simplifies** the original X/Y definitions and adds **population-level** and **disease-annotation** analyses (GWAS, ClinVar).

---

## 1. Overall Objective

Implement a generalizable computational pipeline that:

1. Takes as input any transcription factor (TF) (starting with **AP1**).
2. Retrieves its motif model (PWM/PFM) from a motif database.
3. Scans the entire human genome for all motif-like sequences (strong + weak + near-motif).
4. Computes mutation paths from each motif instance to the consensus motif.
5. Intersects each mutation step with **gnomAD v4.1** variants to determine:
   - Whether each activating mutation exists in the population.
   - At what **allele frequency (AF)** and **homozygote count**.
6. Scores the functional impact of activating variants using **AlphaGenome** (REF vs ALT sequences).
7. Places each genomic site into a **2D landscape**:
   - **X-axis:** `X = −log10(max_AF_step) × hamming_distance` (population accessibility / constraint).
   - **Y-axis:** `Y = max(ΔAlphaGenome)` (signed functional impact).
8. Produces:
   - Ranked lists of “dormant but potentially functional” motif sites.
   - Population-level AF histograms.
   - Overlaps with **GWAS** and **ClinVar** variants.
   - Publication-ready visualizations and summary reports.

The pipeline must run on **ODYSSEUS or any HPC**, be modular, and support batch-mode for large models like AlphaGenome.

---

## 2. Directory Structure (Copilot: Use / Update These)

The repo is already mostly structured, but Copilot must ensure the following layout (modifying existing files where necessary, not duplicating):

```text
alphamotif-activation-pipeline/
│
├── directive.md                      # THIS FILE (updated)
├── pipeline_config.yaml              # Central config file
│
├── data/
│   ├── reference/                    # GRCh38 genome + indexes
│   ├── motifs/                       # PWMs/PFMs from JASPAR, HOCOMOCO
│   └── gnomad/                       # symlink to /mnt/data_1/gnomAD_data/raw/gnomad_v4.1
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
│   ├── summarize_path_AF.py      # (existing functionality; may be adjusted)
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
├── 07_population_statistics/        # NEW MODULE
│   ├── compute_af_histogram.py
│   └── DOCUMENTATION.md
│
├── 08_gwas_clinvar_integration/     # NEW MODULE
│   ├── annotate_with_gwas_clinvar.py
│   └── DOCUMENTATION.md
│
└── utils/
    ├── io.py
    ├── logging_utils.py
    ├── constants.py
    └── DOCUMENTATION.md
````

---

## 3. Data Requirements (Copilot: Use These Paths/Commands)

All gnomAD data must reside at:

```text
/mnt/data_1/gnomAD_data/raw/gnomad_v4.1/
```

Copilot must ensure there is a symlink:

```bash
ln -s /mnt/data_1/gnomAD_data/raw/gnomad_v4.1 data/gnomad
```

### 3.1 gnomAD v4.1 (GRCh38)

**Required genomes VCF:**

```bash
mkdir -p /mnt/data_1/gnomAD_data/raw/gnomad_v4.1/genomes_vcf
gsutil cp gs://gcp-public-data--gnomad/release/4.1/genomes_vcf/genomes.vcf.bgz \
    /mnt/data_1/gnomAD_data/raw/gnomad_v4.1/genomes_vcf/
gsutil cp gs://gcp-public-data--gnomad/release/4.1/genomes_vcf/genomes.vcf.bgz.tbi \
    /mnt/data_1/gnomAD_data/raw/gnomad_v4.1/genomes_vcf/
```

**Recommended (coverage):**

```bash
mkdir -p /mnt/data_1/gnomAD_data/raw/gnomad_v4.1/coverage
gsutil -m cp -r gs://gcp-public-data--gnomad/release/4.1/coverage/genomes \
    /mnt/data_1/gnomAD_data/raw/gnomad_v4.1/coverage/
```

**Optional (constraint scores):**

```bash
mkdir -p /mnt/data_1/gnomAD_data/raw/gnomad_v4.1/constraint
gsutil cp gs://gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint.json.bgz \
    /mnt/data_1/gnomAD_data/raw/gnomad_v4.1/constraint/
```

### 3.2 GRCh38 Reference

Copilot should implement:

```bash
wget ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
```

Place at:

```text
data/reference/GRCh38.fa
```

### 3.3 Motifs

Copilot should pull:

* JASPAR 2024 PFMs, specifically:

  * AP1 (JASPAR ID `MA0099.3`)

Place PFMs at:

```text
data/motifs/JASPAR/
```

---

## 4. Updated Scientific Design: X and Y Definitions

### 4.1 X-axis – Population Accessibility / Constraint

For each **motif site** and its mutation paths:

* Let `max_AF_step` be the **maximum AF** across all mutation steps along a chosen path.
* Let `hamming_distance` be the number of mutations required to reach consensus (path length).

The X-axis combines **frequency constraint** and **mutational accessibility**:

```text
X = −log10(max_AF_step) × hamming_distance
```

**Biological rationale:**
* **Frequency constraint** (`-log10(max_AF_step)`): Rare variants are harder to access
* **Mutational steps** (`hamming_distance`): Each additional mutation multiplies the constraint
* **Joint probability**: A 3-step path with AF=0.001 per step is ~10⁶× less accessible than a 1-step path at the same AF

Implementation detail (for log(0)):

* If `max_AF_step == 0`, treat it as `max_AF_step = 1e-12` (or another very small epsilon) **only inside the log**, but still store and report AF as 0.
* Conceptually, AF=0 corresponds to **extremely constrained** variants and appears on the far right of the X-axis.

Copilot must:

* **Include Hamming distance as a multiplier in X calculation.**
* Use epsilon (1e-12) only internally to avoid NaNs in log10.
* When summarizing, treat absence of a gnomAD record at (chr, pos) as AF=0 for that REF→ALT mutation.

### 4.2 Y-axis – Functional Impact (AlphaGenome)

For each motif site, using AlphaGenome predictions across REF/ALT alleles:

* Compute signed differences in AlphaGenome scores:

```text
ΔAlphaGenome = AlphaGenome(ALT) − AlphaGenome(REF)
```

The updated Y-axis is:

```text
Y = max(ΔAlphaGenome)
```

Important constraints:

* **Preserve negative values**; **do not** take absolute values.
* **Do not** transform ΔAlphaGenome into quantiles for the primary X–Y landscape (you may create optional secondary plots with quantiles if helpful).

Copilot must:

* Ensure any previous use of `abs(ΔAlphaGenome)` for the landscape is removed.
* Keep the sign for all downstream analyses and visualizations.

---

## 5. Module Specifications (By Number)

### Module 01 – Genome-Wide Motif Scan

**Input:**

* PWM/PFM for TF
* GRCh38 genome

**Tasks (Copilot):**

1. Run FIMO or MOODS to detect **all motif-like sites**.

2. Store results:

   ```text
   results/motif_scan/{TF}/motif_hits.bed
   ```

3. Tier sites by PWM score into tiers 0–3 based on config:

   ```yaml
   motif_scan:
     pvalue_threshold: 1e-3
     tier_cutoffs:
       tier0: 0.95
       tier1: 0.85
       tier2: 0.70
   ```

---

### Module 02 – Mutation Path Enumeration

For each motif hit:

1. Infer the consensus motif from PWM.
2. Compute mismatch positions between site sequence and consensus.
3. Enumerate all **minimal mutation paths** up to a maximum Hamming distance (e.g., ≤ 3).

**Output:**

```text
results/mutation_paths/{TF}/paths.tsv
```

**Columns (Copilot: ensure at least):**

* `site_id`, `path_id`, `path_num`, `step_num`, `total_steps`
* `chr`, `motif_start`, `motif_end`, `strand`, `tier`, `pwm_score`, `hamming_distance`
* `position_in_motif`, `genomic_position`
* `ref_base`, `alt_base`
* `seq_before`, `seq_after`

These already exist in `02_generate_mutation_paths/enumerate_paths.py` – modify rather than recreate.

---

### Module 03 – gnomAD Integration

**Key requirement:** Correct handling of AF=0 (missing variants).

**Biological significance:** Variants with AF=0 (not found in gnomAD's 807K individuals) represent the **most constrained** mutations and will appear at the far right of the X-axis landscape. These are often the most interesting candidates as they suggest strong purifying selection against activation.

**Tasks (Copilot):**

1. Convert mutation paths to genomic coordinates and build a **BED** file of all unique positions.

2. Use **bcftools** to query gnomAD:

   ```bash
   bcftools query \
     -R query_regions.bed \
     -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/AC\t%INFO/AN\t%INFO/nhomalt\n' \
     /mnt/data_1/gnomAD_data/raw/gnomad_v4.1/genomes_vcf/genomes.vcf.bgz
   ```

3. Match queried variants back to mutation paths:

   * Match on `(chr, genomic_position, ref_base, alt_base)`.
   * If no record is found for a mutation step, treat AF=0 for that REF→ALT.

4. Fill NaNs:

   ```text
   AF = 0.0  (for missing records - store as true zero)
   AC = 0
   AN = 0
   nhomalt = 0
   ```
   
   **Critical:** Store AF=0 as literal 0.0 in output files. The epsilon substitution (1e-12) is applied **only** during X-axis calculation to avoid log(0) errors, not during data storage.

5. Summarize by **path_id** (and optionally by site_id):

   * `max_AF` (true value, including 0.0 for missing variants)
   * `mean_AF`
   * `sum_AF`
   * `total_steps`
   * `total_AC`
   * `total_nhomalt`
   * `num_af_zero` (count of steps with AF=0 - important metric for constraint analysis)

**Outputs:**

* `results/gnomad_intersection/{TF}/gnomad_query_results.tsv`
* `results/gnomad_intersection/{TF}/paths_with_gnomad.tsv`
* `results/gnomad_intersection/{TF}/path_af_summary.tsv`

Paths and functions currently exist in `03_intersect_gnomad/query_gnomad_vcfs.py`. Copilot must **update**, not replace, to conform to AF=0 and later X-axis usage.

---

### Module 04 – AlphaGenome Functional Scoring

**Tasks (Copilot):**

1. For each mutation (REF→ALT):

   * Generate **reference** and **mutated** sequence windows (±1 kb or as configured). **Always deduplicate identical windows** prior to scoring by hashing the full `(chr, window_start, window_end, ref_sequence, alt_sequence)` tuple. Treat windows as identical only when both REF and ALT sequences (and their coordinates) match exactly, including ALT base. This deduplication is a fixed default (no CLI/config toggle) and should be documented but not exposed as a parameter.

2. Run AlphaGenome in batch (after deduplication):

   ```bash
   alphagenome predict --model <PATH> \
       --fasta seqs.fa \
       --output predictions.json
   ```

3. Compute **signed** Δ scores per feature and per motif:

   ```text
   ΔFeature = Feature(ALT) − Feature(REF)
   ```

4. Summarize per motif site:

   * `max_delta` = max Δ across features and positions (retain sign).
   * Optionally: per-feature maximums.

**Output:**

```text
results/alphagenome/{TF}/functional_deltas.tsv
```

**Critical:** Copilot must confirm that **no absolute values** are used when computing Δ for the main landscape. Deduplicated predictions must be re-used for every mutation step that maps to the same window so AlphaGenome is never re-run for an identical sequence pair.

---

### Module 05 – Activation Landscape

**Goal:** Combine population and functional impact into X/Y coordinates.

**Inputs:**

* `path_af_summary.tsv` (from Module 03)
* `functional_deltas.tsv` (from Module 04)

**Tasks (Copilot):**

1. For each motif site (via site_id):

   * Select a “best” path using **Option A (default and only supported heuristic)**:
     1. Choose the path with the **fewest mutation steps** (smallest Hamming distance / shortest path length).
     2. If multiple paths tie for shortest length, pick the one whose **maximum allele frequency** (`max_AF`) is **lowest** (i.e., most constrained). There is no fallback to alternative heuristics—do not implement Option B or expose a configuration flag.
   * Compute:

     ```text
     max_AF_step = max_AF for that path
     hamming_distance = number of mutation steps in that path
     X_constraint = −log10(max_AF_step) × hamming_distance    # if max_AF_step == 0, use epsilon internally
     Y_impact     = max_delta (signed ΔAlphaGenome)
     ```

2. Produce a final **site-level** table:

   ```text
   results/landscape/{TF}/activation_landscape.tsv
   ```

   Columns (minimum):

   * `site_id`
   * `best_path_id`
   * `chr`, `motif_start`, `motif_end`, `strand`, `tier`
   * `max_AF_step`, `X_constraint`
   * `max_delta`, `Y_impact`
   * `total_steps`, `hamming_distance`

3. `classify_quadrants.py` may label each site as:

   * High‐constraint / high‐impact
   * High‐constraint / low‐impact
   * Low‐constraint / high‐impact
   * Low‐constraint / low‐impact

Use the updated definition of X and Y.

---

### Module 06 – Visualization

**Tasks (Copilot):**

1. `plot_landscape.py`:

   * Generate 2D scatter or hexbin of:

     * X-axis: `X_constraint = −log10(max_AF_step) × hamming_distance`
     * Y-axis: `Y_impact = max_delta` (signed)

   * Highlight:

     * AF=0 sites (e.g., different color/shape).
     * Optionally, sites with GWAS/ClinVar annotations (after Modules 07–08).

2. `plot_genome_tracks.py`:

   * For selected top candidates, prepare IGV-style or track-style plots using AlphaGenome predictions, motif positions, and variant info.

---

### Module 07 – Population Statistics (NEW)

**Goal:** Describe the global AF distribution across all motif sites.

**Inputs:**

* `path_af_summary.tsv` from Module 03
* Or a derived per‐site max_AF table

**Tasks (Copilot):**

1. Compute per‐site or per‐path:

   * `max_AF_step` = maximum AF observed among its steps.

2. Generate a histogram:

   * X-axis: `max_AF_step` (e.g., log-spaced bins from 0 to 0.5)
   * Y-axis: count of sites in each bin.

3. Special handling:

   * Explicitly count and visualize sites with `max_AF_step = 0` (no observed variation).
   * Report the fraction of sites with AF=0, ultra-low AF, etc.

**Outputs:**

* `results/population_stats/{TF}/max_af_histogram.tsv`
* `figures/{TF}/max_af_histogram.png`

---

### Module 08 – GWAS & ClinVar Integration (NEW)

**Goal:** Determine whether high-impact, low-frequency sites overlap or sit near known disease-associated variants.

**Inputs:**

* `activation_landscape.tsv` (Module 05)
* Public GWAS catalog (e.g., downloaded summary stats)
* ClinVar VCF

**Tasks (Copilot):**

1. `annotate_with_gwas_clinvar.py`:

   * Load GWAS catalog and ClinVar.
   * For each motif site and/or mutation step:

     * Check if there is a ClinVar variant at the same (chr, pos).
     * Check if there is a GWAS variant:

       * Within a configurable window (e.g., ±10 kb).
       * Optionally, in linkage disequilibrium (if LD data available).

2. Add boolean and categorical annotations to sites:

   * `near_gwas` (True/False)
   * `near_clinvar` (True/False)
   * `clinvar_significance` (e.g., pathogenic/benign/VUS)
   * `gwas_traits` (comma-separated list if multiple)

3. Write a final annotated table:

   ```text
   results/landscape/{TF}/activation_landscape_annotated.tsv
   ```

4. Modify visualization (Module 06) to optionally color points by GWAS/ClinVar status.

---

## 6. Pipeline Configuration Template (`pipeline_config.yaml`)

Copilot must ensure `pipeline_config.yaml` includes at least:

```yaml
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
   # AlphaGenome scoring always deduplicates identical REF/ALT windows by genomic position.
   # This behavior is fixed (no user override) and ensures we never rescore the same window twice.

motif_scan:
  pvalue_threshold: 1e-3
  tier_cutoffs:
    tier0: 0.95
    tier1: 0.85
    tier2: 0.70

mutation_paths:
  max_hamming_distance: 3
  tier_filter: [1, 2]

numerical_stability:
   # AF=0 values are stored as true 0.0 in all data files.
   # For X-axis calculation only: X = -log10(max(AF, epsilon)) × hamming_distance
   # where epsilon prevents log(0) errors. This value represents the limit of detection
   # for gnomAD v4.1 (1 / ~800,000 individuals ≈ 1.25e-6, but we use 1e-12 for numerical safety).
   max_af_epsilon: 1e-12  # Applied only in log10 calculations, never stored in data

resources:
  threads: 32
```

---

## 7. Copilot Implementation Directives

Copilot must:

* **Modify existing scripts** to conform to this updated directive rather than creating redundant new ones.

* Implement or adjust each module as standalone Python scripts with clear CLIs.

* Ensure that:

  * AF=0 logic is correctly handled in Module 03 and used downstream.
  * X-axis uses `−log10(max_AF_step) × hamming_distance` with epsilon only for numerical safety in log10.
  * Y-axis uses **signed** ΔAlphaGenome, no absolute values.

* Allow the entire pipeline to be run via:

  ```bash
  python run_pipeline.py --config pipeline_config.yaml
  ```

* Optionally generate a Snakemake or Nextflow workflow.

* Create or maintain unit tests for:

  * motif scanning
  * mutation path enumeration
  * gnomAD querying and AF=0 behavior
  * AlphaGenome I/O and Δ calculations
  * activation landscape combination

All intermediate outputs must be stored under `results/`.

---

## 8. Success Criteria

The pipeline is considered successful when:

1. AP1 (MA0099.3) can be processed **end-to-end**.
2. An activation landscape is produced with:

   * X-axis = `−log10(max_AF_step) × hamming_distance`.
   * Y-axis = signed `max(ΔAlphaGenome)`.
3. The AF histogram (Module 07) shows the distribution of max_AF across all sites, including AF=0.
4. High-impact, low-frequency candidate sites are annotated with GWAS/ClinVar where applicable.
5. The codebase can be generalized to any TF by changing `tf_name` and `motif_id` in the config.

```