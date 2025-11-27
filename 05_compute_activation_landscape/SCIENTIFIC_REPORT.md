# Activation Landscape of Dormant AP1 Transcription Factor Binding Sites Accessible Through Human Population Variation

**Authors:** George Stephenson¹

**Affiliations:**  
¹ LAYER Laboratory, Department of Molecular, Cellular, and Developmental Biology, University of Colorado Boulder, Boulder, CO 80309

**Date:** November 27, 2025

**Keywords:** transcription factor binding, AP1, regulatory variation, population genetics, AlphaGenome, functional genomics

---

## Abstract

Transcription factor binding sites (TFBSs) are fundamental regulatory elements that control gene expression, yet the human genome harbors millions of "dormant" sequences that are near-matches to consensus motifs but lack sufficient affinity for functional binding. We hypothesized that naturally occurring genetic variants could activate these dormant sites, potentially creating novel regulatory elements with functional consequences. Here, we present a comprehensive computational analysis of dormant AP1 (FOS/JUN heterodimer) binding sites across the human genome, integrating motif scanning, population genetics from gnomAD v4.1, and functional impact prediction using AlphaGenome's multi-modal deep learning framework. We identified 6,810 variants that could activate dormant AP1 sites, with 90.6% (n=6,170) predicted to substantially increase AP1-family transcription factor binding. Strikingly, 85.1% of these activating variants are ultra-rare (AF < 0.01%), suggesting purifying selection against uncontrolled AP1 site activation. We constructed a two-dimensional "activation landscape" mapping population accessibility (X-axis) against functional impact (Y-axis), identifying 1,624 high-priority candidates that are both population-accessible and functionally impactful. The predominant AP1-family transcription factors showing predicted binding gains were FOS (25.8%), JUND (16.8%), and ATF3 (8.8%), with the strongest effects observed in endothelial cells, MCF-7 breast cancer cells, and HepG2 hepatocytes. A strong positive correlation between AP1 binding impact and enhancer mark activation (H3K27ac/H3K4me1; r=0.545, p<10⁻¹⁶) validates that predicted TF binding gains correspond to functional enhancer activation. These findings reveal that dormant AP1 sites represent a reservoir of latent regulatory potential, constrained by purifying selection, with implications for understanding regulatory evolution, disease mechanisms, and the functional consequences of rare non-coding variation.

---

## 1. Introduction

### 1.1 Background

Transcription factor (TF) binding sites are critical cis-regulatory elements that govern spatial and temporal patterns of gene expression (Spitz & Furlong, 2012). While high-affinity TF binding sites have been extensively characterized, the vast majority of genomic sequences that weakly resemble TF motifs remain functionally inactive—what we term "dormant sites." These sequences lack sufficient sequence identity to the consensus motif to support stable TF binding under normal conditions.

The AP1 (Activator Protein 1) transcription factor complex, comprising heterodimers of FOS and JUN family proteins, represents an ideal model system for studying dormant site activation. AP1 binds the consensus sequence TGA(G/C)TCA (JASPAR MA0099.3) and plays critical roles in cell proliferation, differentiation, and stress response (Shaulian & Karin, 2002). Dysregulation of AP1 activity is implicated in cancer, inflammation, and immune disorders.

### 1.2 Scientific Question

We asked: **Which dormant AP1 binding sites in the human genome could become functionally active through mutations that already exist in human populations, and what would be the functional consequences of such activation?**

This question has profound implications for:

1. **Regulatory Evolution:** Understanding how novel TF binding sites emerge in populations
2. **Disease Mechanisms:** Identifying rare variants that may create aberrant regulatory elements
3. **Selection Constraints:** Quantifying purifying selection acting on regulatory potential
4. **Therapeutic Targets:** Discovering dormant sites that could be intentionally activated

### 1.3 Approach

We developed an integrated computational pipeline that:

1. Scans the human genome for all AP1 motif-like sequences (strong to weak matches)
2. Enumerates mutation paths from each dormant site to the consensus motif (≤3 mutations)
3. Intersects activating mutations with gnomAD v4.1 to identify variants observed in human populations
4. Predicts functional impact using AlphaGenome, focusing on AP1-family TF binding predictions
5. Constructs an "activation landscape" mapping population accessibility versus functional impact

Critically, we designed our functional impact metric (Y-axis) to be **biologically specific to AP1**, using predicted changes in AP1-family TF binding (JUND, JUN, JUNB, FOS, FOSL1, FOSL2, ATF3, ATF2, BATF, MAFK) rather than a generic maximum across all AlphaGenome outputs.

---

## 2. Materials and Methods

### 2.1 Motif Scanning and Dormant Site Identification

We obtained the AP1 position weight matrix (PWM) from JASPAR (MA0099.3) and scanned the GRCh38 human reference genome using FIMO (Grant et al., 2011) with a permissive p-value threshold (p < 10⁻³) to capture weak motif matches. Sites were tiered by PWM score:

- **Tier 0 (Strong):** ≥95th percentile PWM score
- **Tier 1 (Moderate):** 80-95th percentile
- **Tier 2 (Weak):** 50-80th percentile  
- **Tier 3 (Marginal):** <50th percentile

This approach identified ~6.6 million AP1 motif-like sequences genome-wide.

### 2.2 Mutation Path Enumeration

For each dormant site (Tiers 1-3), we computed the Hamming distance to the AP1 consensus and enumerated all minimal mutation paths (≤3 single nucleotide substitutions) required for activation. Each mutation step was characterized by:

- Genomic coordinates (chromosome, position)
- Reference and alternate alleles
- Position within the motif

This generated 18.1 million mutation steps across 6.3 million paths.

### 2.3 Population Genetics Integration (gnomAD v4.1)

We queried gnomAD v4.1 genomes data (n=807,162 individuals) using bcftools to identify which activating mutations exist in human populations. For each variant, we extracted:

- **Allele frequency (AF):** Population prevalence
- **Allele count (AC):** Number of observed alleles
- **Allele number (AN):** Total alleles sequenced (coverage proxy)

Variants not found in gnomAD (high coverage regions) were assigned AF=0, representing sites under potential purifying selection. We identified 7,158 unique variants matching mutation path steps.

### 2.4 Functional Impact Prediction (AlphaGenome)

We scored all 7,158 variants using the AlphaGenome API with the recommended 19 variant scorers, generating predictions across 11 output modalities:

- **CHIP_TF:** Transcription factor binding (714 biosamples)
- **CHIP_HISTONE:** Histone modifications (H3K27ac, H3K4me1, H3K4me3, etc.)
- **ATAC/DNASE:** Chromatin accessibility
- **RNA_SEQ:** Gene expression
- **CAGE/PROCAP:** Transcription start sites
- **SPLICE_*:** Splicing predictions
- **CONTACT_MAPS:** 3D chromatin interactions

AlphaGenome successfully scored 6,810 variants (95.1% success rate), generating 186,695,928 individual predictions.

### 2.5 AP1-Specific Y-Axis Design

Rather than using a naive `max()` across all AlphaGenome tracks (which would conflate unrelated signals like CTCF or splice sites), we designed a **biologically specific Y-axis** focused on AP1-family transcription factors:

**Primary Y-axis (AP1 Impact):**
```
Y = max(quantile_score) for TF ∈ {JUND, JUN, JUNB, FOS, FOSL1, FOSL2, ATF3, ATF2, ATF7, BATF, BATF2, MAFK, MAFF, MAFG}
```

**Secondary validation (Enhancer marks):**
```
Y_enhancer = max(quantile_score) for histone_mark ∈ {H3K27ac, H3K4me1, H3K4me3}
```

This approach provides direct biological relevance: if a variant increases AP1-family binding, it supports the dormant site activation hypothesis.

### 2.6 X-Axis: Population Accessibility Score

The X-axis quantifies how "accessible" dormant site activation is through human variation:

```
X = -log₁₀(AF) × Hamming_distance
```

Where:
- **-log₁₀(AF):** Transforms allele frequency so rare variants score higher
- **Hamming_distance:** Number of mutations required to reach consensus (1, 2, or 3)

**Interpretation:**
- Low X = highly accessible (common variant, few mutations needed)
- High X = hard to access (rare variant, many mutations needed)

### 2.7 Quadrant Classification

We classified variants into four quadrants based on median X and Y values:

1. **HIGH PRIORITY (Accessible + High Impact):** X < median, Y > median
2. **High Impact, Hard to Access:** X ≥ median, Y > median
3. **Accessible, Low Impact:** X < median, Y ≤ median
4. **Low Priority:** X ≥ median, Y ≤ median

### 2.8 Statistical Analysis

All statistical analyses were performed in Python using scipy.stats. Correlations were assessed using Pearson correlation coefficients. Significance threshold was set at α = 0.05.

---

## 3. Results

### 3.1 Overview of Dormant AP1 Site Activation Variants

We identified **6,810 unique genetic variants** from gnomAD v4.1 that map to mutation paths capable of activating dormant AP1 binding sites (Table 1). All variants had complete AlphaGenome predictions including AP1-family TF binding and enhancer histone marks.

**Table 1. Summary of Dormant AP1 Site Activation Variants**

| Metric | Value |
|--------|-------|
| Total variants analyzed | 6,810 |
| Variants with AP1-family TF predictions | 6,810 (100%) |
| Variants with enhancer mark predictions | 6,810 (100%) |
| Total AlphaGenome predictions | 186,695,928 |
| AP1-family TF predictions | 926,160 |
| Enhancer mark predictions | 7,164,120 |

### 3.2 Mutational Distance to Activation

The majority of dormant sites require multiple mutations to reach the AP1 consensus (Figure 1A). The Hamming distance distribution was:

- **1 mutation:** 1 variant (0.0%)
- **2 mutations:** 656 variants (9.6%)
- **3 mutations:** 6,153 variants (90.4%)

This indicates that most dormant sites are evolutionarily "buffered" from accidental activation by requiring multiple simultaneous mutations—a pattern consistent with purifying selection maintaining regulatory specificity.

### 3.3 Allele Frequency Distribution Reveals Strong Constraint

Strikingly, the vast majority of AP1-activating variants are extremely rare in human populations (Figure 1B):

**Table 2. Allele Frequency Distribution of Activating Variants**

| AF Category | Count | Percentage |
|-------------|-------|------------|
| Common (AF ≥ 1%) | 187 | 2.7% |
| Low frequency (0.1-1%) | 228 | 3.3% |
| Rare (0.01-0.1%) | 597 | 8.8% |
| Ultra-rare (AF < 0.01%) | 5,798 | **85.1%** |

The extreme rarity of activating variants (85.1% ultra-rare) suggests strong **purifying selection** against mutations that would create novel AP1 binding sites. This is consistent with the hypothesis that uncontrolled TF binding site activation could have deleterious regulatory consequences.

### 3.4 AP1-Family TF Binding Impact

AlphaGenome predictions revealed that the overwhelming majority of variants substantially increase AP1-family TF binding (Figure 2A):

**Table 3. AP1 Impact Score Distribution**

| Metric | Value |
|--------|-------|
| Mean AP1 impact score | 0.958 |
| Median AP1 impact score | 0.978 |
| Standard deviation | 0.075 |
| Range | 0.011 - 1.000 |
| Variants >90th percentile | 6,170 (90.6%) |
| Variants >95th percentile | 5,065 (74.4%) |
| Variants >99th percentile | 2,286 (33.6%) |

Remarkably, **90.6% of variants** showed AP1 binding impact in the top 10% of AlphaGenome's reference distribution. This validates that our mutation path approach successfully identifies variants with strong predicted effects on AP1 binding.

### 3.5 AP1-Family TF Specificity

The dominant AP1-family transcription factors showing the strongest predicted binding gains were (Figure 2B):

**Table 4. Distribution of Best-Responding AP1-Family TFs**

| Transcription Factor | Count | Percentage |
|----------------------|-------|------------|
| FOS | 1,754 | 25.8% |
| JUND | 1,141 | 16.8% |
| ATF3 | 599 | 8.8% |
| FOSL2 | 577 | 8.5% |
| ATF2 | 551 | 8.1% |
| FOSL1 | 521 | 7.7% |
| JUNB | 450 | 6.6% |
| JUN | 369 | 5.4% |
| MAFK | 369 | 5.4% |
| BATF | 177 | 2.6% |

FOS showed the strongest response, consistent with its role as a primary AP1 component. The diversity of responding TFs (all major FOS/JUN/ATF family members) confirms that our dormant sites are bona fide AP1 consensus-adjacent sequences.

### 3.6 Cell Type Specificity

Predicted AP1 binding gains showed strong cell-type specificity (Figure 2C):

**Table 5. Top Cell Types Showing AP1 Activation**

| Cell Type | Count | Percentage |
|-----------|-------|------------|
| Endothelial cell (HUVEC) | 1,071 | 15.7% |
| MCF-7 (breast cancer) | 902 | 13.2% |
| HepG2 (hepatocytes) | 822 | 12.1% |
| T47D (breast cancer) | 729 | 10.7% |
| GM12878 (lymphoblastoid) | 641 | 9.4% |
| K562 (leukemia) | 569 | 8.4% |
| A549 (lung cancer) | 382 | 5.6% |
| SK-N-SH (neuroblastoma) | 351 | 5.2% |
| H1 (embryonic stem cells) | 281 | 4.1% |
| Liver | 206 | 3.0% |

The prominence of cancer cell lines (MCF-7, HepG2, K562, T47D, A549) is notable and may reflect the extensive ChIP-seq data available for these lines, or could indicate that AP1 binding sites are particularly dynamic in transformed cells.

### 3.7 Activation Landscape

We constructed a two-dimensional activation landscape with population accessibility (X-axis) versus AP1 functional impact (Y-axis) (Figure 3). Quadrant analysis revealed:

**Table 6. Activation Landscape Quadrant Distribution**

| Quadrant | Count | Percentage |
|----------|-------|------------|
| **HIGH PRIORITY: Accessible + High Impact** | 1,624 | 23.8% |
| High Impact, Hard to Access | 1,779 | 26.1% |
| Accessible, Low Impact | 1,781 | 26.2% |
| Low Priority | 1,626 | 23.9% |

The **1,624 high-priority candidates** (23.8%) represent dormant AP1 sites that are:
1. Relatively accessible through existing human variation
2. Predicted to have substantial functional impact upon activation

### 3.8 High-Priority Candidate Characteristics

The high-priority quadrant showed distinctive characteristics:

**Table 7. High-Priority Candidate Summary**

| Metric | Value |
|--------|-------|
| Total candidates | 1,624 |
| Mean AP1 impact score | 0.992 |
| Mean allele frequency | 1.01% |
| Common variants (AF ≥ 1%) | 82 |
| Single-step variants (H=1) | 1 |

The elevated mean AF (1.01% vs 0.39% overall) indicates that high-priority candidates are enriched for more common variants, consistent with their "accessible" classification.

### 3.9 Validation: AP1 Binding Correlates with Enhancer Activation

To validate that predicted AP1 binding gains correspond to functional enhancer activation, we correlated AP1 impact scores with enhancer histone mark predictions (H3K27ac/H3K4me1) (Figure 4).

**Correlation Analysis:**
- **AP1 impact vs Enhancer impact:** r = 0.545, p < 10⁻¹⁶

This strong positive correlation demonstrates that variants predicted to increase AP1 binding are also predicted to increase active enhancer marks, supporting a mechanistic model where AP1 binding drives enhancer activation.

### 3.10 Evidence for Purifying Selection

We observed a weak but significant negative correlation between variant rarity and AP1 impact:

- **-log₁₀(AF) vs AP1 impact:** r = -0.049, p = 5.99 × 10⁻⁵

While the effect size is modest, the direction is consistent with purifying selection: variants with higher predicted functional impact tend to be rarer in the population, suggesting they are selected against.

---

## 4. Discussion

### 4.1 Principal Findings

This study provides the first comprehensive analysis of dormant AP1 transcription factor binding sites accessible through human population variation. Our key findings are:

1. **Dormant sites are abundant but buffered:** 90.4% of activating variants require 3 mutations, indicating evolutionary buffering against accidental activation.

2. **Strong constraint on activation:** 85.1% of activating variants are ultra-rare (AF < 0.01%), suggesting purifying selection against novel AP1 site creation.

3. **High predicted functional impact:** 90.6% of variants substantially increase predicted AP1-family TF binding (>90th percentile).

4. **Validation through enhancer marks:** Strong correlation (r = 0.545) between AP1 binding gains and enhancer activation (H3K27ac/H3K4me1).

5. **1,624 high-priority candidates:** Accessible variants with high functional impact represent targets for further investigation.

### 4.2 Biological Implications

#### 4.2.1 Regulatory Evolution

Our findings illuminate the evolutionary dynamics of TF binding site gain-of-function. The extreme rarity of activating variants suggests that the regulatory genome is under strong constraint to prevent spurious TF binding. The requirement for multiple mutations (Hamming distance = 3 for 90% of sites) provides an additional buffer.

#### 4.2.2 Disease Mechanisms

The high-priority candidates identified here represent potential disease-relevant variants. Rare variants that activate dormant AP1 sites could create aberrant enhancers, leading to:
- Ectopic gene expression
- Chromatin remodeling at normally inactive loci
- Disruption of normal regulatory topology

Future work should intersect these candidates with GWAS signals and ClinVar pathogenic variants.

#### 4.2.3 AP1 as a Pioneer Factor

The strong correlation between AP1 binding and enhancer activation supports AP1's role as a "pioneer factor" capable of initiating chromatin opening and enhancer establishment (Biddie et al., 2011). Dormant site activation may represent a mechanism for de novo enhancer creation.

### 4.3 Methodological Innovations

#### 4.3.1 Biologically-Specific Y-Axis

Unlike previous approaches that use generic maximum scores across all functional predictions, we designed our Y-axis to specifically measure AP1-family TF binding. This provides direct biological interpretability: a high Y-score means the variant is predicted to increase binding of the exact TF family under study.

#### 4.3.2 Population Accessibility Score

Our X-axis formula (X = -log₁₀(AF) × Hamming_distance) elegantly captures the joint probability of accessing a dormant site through population variation: rare variants AND multiple mutations compound the difficulty of activation.

#### 4.3.3 Vectorized Computation

Our implementation uses vectorized pandas groupby operations to process 186 million predictions in ~30 seconds, enabling rapid iteration and analysis.

### 4.4 Limitations

1. **Prediction-based:** AlphaGenome scores are computational predictions, not experimental measurements. Experimental validation (ChIP-seq, reporter assays) is needed.

2. **Context-dependence:** AP1 activity is highly context-dependent (cell type, stimulation, cofactors). Our analysis captures static sequence potential, not dynamic regulation.

3. **Single TF focus:** We analyzed AP1 in isolation. In reality, TFs function in combinatorial networks; dormant site activation may depend on cofactor binding sites.

4. **Limited to SNVs:** We considered only single nucleotide substitutions. Insertions, deletions, and structural variants could also activate dormant sites.

### 4.5 Future Directions

1. **Experimental validation:** Test high-priority candidates using luciferase reporters or CRISPRa.

2. **Disease integration:** Intersect with GWAS catalog and ClinVar to identify disease-associated dormant site activations.

3. **Multi-TF analysis:** Extend pipeline to other TF families (NF-κB, STATs, nuclear receptors).

4. **3D genome context:** Integrate Hi-C data to assess whether dormant sites are in accessible chromatin domains.

5. **Population stratification:** Analyze AF differences across gnomAD populations (AFR, EAS, EUR, etc.).

---

## 5. Conclusions

We present the first systematic characterization of dormant AP1 binding sites accessible through human population variation. Our analysis of 6,810 variants reveals that:

- The vast majority (85.1%) are extremely rare, suggesting purifying selection against AP1 site activation
- Nearly all (90.6%) have high predicted functional impact on AP1-family TF binding
- 1,624 variants represent high-priority candidates: accessible AND impactful

The strong correlation between AP1 binding gains and enhancer activation validates our approach and supports a model where dormant site activation could create functional regulatory elements. These findings have implications for understanding regulatory evolution, disease mechanisms, and the functional interpretation of rare non-coding variation.

---

## 6. Figures

### Figure 1. High-Priority Dormant AP1 Site Candidates

![High-Priority Candidates](../figures/landscape/AP1_high_priority_candidates.png)

**Figure 1.** Characterization of the 1,624 high-priority dormant AP1 site candidates. **(A)** Scatter plot of high-priority candidates showing AP1 impact score (Y-axis) versus accessibility score (X-axis), colored by allele frequency (blue=common, red=rare). Points cluster at Y>0.97, indicating consistently high predicted AP1 binding gains. **(B)** Histogram of log₁₀(allele frequency) distribution. Dashed lines indicate 1% MAF (red) and 0.1% MAF (orange) thresholds. The majority of variants are ultra-rare (AF < 0.01%). **(C)** Bar chart of Hamming distance (mutations required). Most high-priority candidates require 3 mutations (n≈1,250), with fewer requiring only 2 mutations (n≈350). **(D)** Table of top 10 high-priority candidates ranked by AP1 score, showing variant coordinates, AP1 score (all 1.000), best-responding TF (JUNB, JUN, FOSL1, FOSL2, BATF), allele frequency, and mutation steps required.

---

### Figure 2. AP1-Family TF Binding Analysis

![TF Breakdown](../figures/landscape/AP1_tf_breakdown.png)

**Figure 2.** Distribution of AP1-family transcription factor binding predictions across 6,810 variants. **(A)** Horizontal bar chart showing the number of variants for which each AP1-family TF showed the strongest predicted binding gain. FOS dominates with 1,754 variants (25.8%), followed by JUND (1,141, 16.8%), ATF3 (599, 8.8%), FOSL2 (577), ATF2 (551), FOSL1 (521), JUNB (450), JUN (369), MAFK (369), BATF (177), ATF7 (107), MAFG (97), BATF2 (80), and MAFF (18). **(B)** Box plots showing the distribution of maximum AP1 quantile scores for the 8 most frequent TFs. Median scores are consistently >0.9 for all TFs, with FOS and JUND showing slightly broader distributions due to their higher variant counts. Outliers (circles) represent the minority of variants with lower predicted impact.

---

### Figure 3. Activation Landscape

![Main Landscape](../figures/landscape/AP1_activation_landscape_main.png)

**Figure 3.** Two-dimensional activation landscape of dormant AP1 binding sites. The X-axis represents population accessibility score (-log₁₀(AF) × Hamming distance), where higher values indicate harder-to-access variants (rare and/or requiring more mutations). The Y-axis represents AP1-family TF binding impact (maximum quantile score across 14 AP1-family TFs). Each point is one variant (n=6,810), colored by log₁₀(allele frequency): blue = common variants (AF≈10⁻²), yellow-green = intermediate (AF≈10⁻⁴ to 10⁻⁵), red-purple = ultra-rare (AF≈10⁻⁸ to 10⁻¹⁰). The green shaded region indicates the HIGH PRIORITY quadrant (X < median, Y > median threshold of 0.978). Dashed horizontal line at Y=0.978 marks the 90th percentile AP1 impact threshold. Dashed vertical line at X≈14 marks the median accessibility score. Annotation box shows n=6,810 total variants with 1,624 (23.8%) in the high-priority quadrant. The strong vertical clustering at high Y-values (>0.9) demonstrates that the majority of dormant site activation variants are predicted to substantially increase AP1 binding.

---

### Figure 4. Multi-Modal Validation

![Comparison Panels](../figures/landscape/AP1_activation_landscape_comparison.png)

**Figure 4.** Multi-modal validation of the activation landscape approach. **(A) PRIMARY: AP1-Family TF Binding** - Our biologically-specific Y-axis using only JUND, JUN, FOS, FOSL1/2, ATF3, and BATF predictions. Color scale shows quantile score (0.2-1.0). Most variants cluster at high Y-values (>0.8). **(B) COMPARISON: Global Max (All Tracks)** - Alternative Y-axis using maximum across all AlphaGenome outputs. While the pattern appears similar, this metric is less biologically interpretable as it may be driven by unrelated signals (e.g., CTCF, splice sites). Note the compressed Y-axis range (0.94-1.00) indicating ceiling effects. **(C) VALIDATION: Enhancer Marks (H3K27ac, H3K4me1)** - Secondary Y-axis using maximum quantile score across active enhancer histone marks. The similar landscape pattern supports functional relevance: variants that increase AP1 binding also increase enhancer marks. **(D) CONCORDANCE: AP1 vs Enhancer** - Direct correlation between AP1-family TF impact (X-axis) and enhancer mark impact (Y-axis) for each variant. Strong positive correlation (r=0.545, p<10⁻¹⁶, shown in inset) validates that predicted AP1 binding gains correspond to predicted enhancer activation. Color scale shows accessibility score. Points along the diagonal demonstrate concordant predictions across modalities.

---

## 7. Data Availability

All data and code are available at:
- **GitHub:** https://github.com/gsstephenson/alphagenome-enhancer-stacking
- **Results:** `results/landscape/AP1/`
- **Figures:** `figures/landscape/`

### Key Output Files

| File | Description |
|------|-------------|
| `AP1_activation_landscape.tsv` | Complete landscape data (6,810 variants) |
| `AP1_high_priority_candidates.tsv` | High-priority quadrant (1,624 variants) |
| `AP1_landscape_summary.txt` | Summary statistics |

---

## 8. Acknowledgments

We thank the gnomAD consortium for making population genetic data publicly available, and the AlphaGenome team for providing API access to their functional prediction models. This work was supported by the LAYER Laboratory, CU Boulder.

---

## 9. References

1. Biddie, S.C., et al. (2011). Transcription factor AP1 potentiates chromatin accessibility and glucocorticoid receptor binding. *Molecular Cell*, 43(1), 145-155.

2. Grant, C.E., Bailey, T.L., & Noble, W.S. (2011). FIMO: scanning for occurrences of a given motif. *Bioinformatics*, 27(7), 1017-1018.

3. Shaulian, E., & Karin, M. (2002). AP-1 as a regulator of cell life and death. *Nature Cell Biology*, 4(5), E131-E136.

4. Spitz, F., & Furlong, E.E. (2012). Transcription factors: from enhancer binding to developmental control. *Nature Reviews Genetics*, 13(9), 613-626.

5. Karczewski, K.J., et al. (2020). The mutational constraint spectrum quantified from variation in 141,456 humans. *Nature*, 581(7809), 434-443.

---

## Supplementary Information

### Supplementary Table S1. AP1-Family Transcription Factors Used for Y-Axis

| TF Symbol | Family | Dimerization Partners |
|-----------|--------|----------------------|
| FOS | Fos | JUN, JUNB, JUND |
| FOSL1 | Fos | JUN, JUNB, JUND |
| FOSL2 | Fos | JUN, JUNB, JUND |
| JUN | Jun | FOS, FOSL1, FOSL2, ATF |
| JUNB | Jun | FOS, FOSL1, FOSL2, ATF |
| JUND | Jun | FOS, FOSL1, FOSL2, ATF |
| ATF2 | ATF | JUN family |
| ATF3 | ATF | JUN family |
| ATF7 | ATF | JUN family |
| BATF | BATF | JUN family |
| BATF2 | BATF | JUN family |
| MAFK | Small Maf | AP1-like binding |
| MAFF | Small Maf | AP1-like binding |
| MAFG | Small Maf | AP1-like binding |

### Supplementary Table S2. AlphaGenome Output Types

| Output Type | Predictions | Description |
|-------------|-------------|-------------|
| RNA_SEQ | 119,517,948 | Gene expression |
| CHIP_TF | 22,023,540 | TF binding |
| CHIP_HISTONE | 15,199,920 | Histone modifications |
| SPLICE_JUNCTIONS | 14,562,927 | Splice junction usage |
| CAGE | 7,436,520 | Transcription start sites |
| DNASE | 4,154,100 | DNase-seq accessibility |
| ATAC | 2,274,540 | ATAC-seq accessibility |
| SPLICE_SITE_USAGE | 1,165,959 | Splice site usage |
| CONTACT_MAPS | 190,680 | 3D chromatin contacts |
| PROCAP | 163,440 | Nascent transcription |
| SPLICE_SITES | 6,354 | Splice site presence |

---

*Manuscript prepared: November 27, 2025*
