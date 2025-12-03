# Module 09: Scientific Report — Visualization of Forbidden Variants

## Executive Summary

This module provides comprehensive visualization of the **29,957 "forbidden variants"** identified by the dormant site activation pipeline. These are single-nucleotide mutations that would complete functional AP1 transcription factor binding sites but are essentially absent from 807,000 human genomes sequenced by gnomAD v4.1.

**Key Finding**: Only 1 of 31,174 possible H=1 (single-mutation-away) variants is observed, representing a **76.5-fold depletion** (p < 10⁻³⁰). This provides strong evidence for purifying selection actively preventing the creation of new AP1 binding sites in the human genome.

---

## Figures Generated

| Figure | File | Description |
|--------|------|-------------|
| 1 | `selection_gradient.png` | Selection intensity by Hamming distance |
| 2 | `forbidden_variants_characterization.png` | Multi-panel characterization |
| 3 | `tf_tissue_heatmap.png` | TF-tissue specificity |
| 4 | `genomic_distribution.png` | Chromosomal distribution |
| 5 | `selection_signal_strength.png` | Evidence for selection intensity |
| 6 | `functional_context.png` | Functional properties |

---

## Figure 1: Selection Gradient

**File**: `figures/selection_gradient/selection_gradient.png`

### What It Shows
A 4-panel figure demonstrating the gradient of purifying selection by Hamming distance from canonical AP1 motifs.

### Panels

**A. Survival Rate by Hamming Distance**
- H=1 variants: 0.0032% survival (1 of 31,174 observed)
- H=2 variants: 0.093% survival (677 observed)
- H=3 variants: 0.246% survival (6,359 observed)
- Inset zoom shows the near-zero H=1 bar

**B. AP1 Score Distribution**
- Compares predicted AP1 binding impact between forbidden (H=1) and observed (H=2, H=3) variants
- Key insight: Distributions are similar (median ~500-600), indicating selection targets site completion, not binding strength

**C. Fold Depletion (log scale)**
- H=1: 76.5× depletion (p < 10⁻³⁰)
- H=2: 2.6× depletion (p < 10⁻¹⁹⁹)
- H=3: 1.0× (baseline)

**D. Key Findings Summary**
- Bullet-point summary of the selection paradox

### Scientific Interpretation
The dramatic 76.5-fold depletion of H=1 variants, combined with similar predicted binding impact to observed variants, reveals that **selection purges site COMPLETION, not binding STRENGTH**. Human genomes are actively maintained one mutation away from functional AP1 sites.

---

## Figure 2: Forbidden Variants Characterization

**File**: `figures/forbidden_variants/forbidden_variants_characterization.png`

### What It Shows
A 4-panel characterization of the molecular and genomic properties of forbidden variants.

### Panels

**A. Mutation Spectrum**
- Transversions dominate: **71.8%** (21,519 variants)
- Transitions: 28.2% (8,438 variants)
- A>T is the most common mutation type (29.3%)
- This reflects the core AP1 motif structure (TGA[CG]TCA), where T↔A changes at critical positions complete the binding site

**B. Top Transcription Factors by Frequency**
| TF | Count | Percentage | Family |
|----|-------|------------|--------|
| FOS | 5,362 | 17.9% | FOS/JUN |
| ATF3 | 3,789 | 12.6% | ATF/CREB |
| ATF7 | 2,523 | 8.4% | ATF/CREB |
| ATF2 | 2,425 | 8.1% | ATF/CREB |
| FOSL1 | 2,315 | 7.7% | FOS/JUN |

The FOS/JUN core family accounts for ~40% of best-responding TFs, consistent with classical AP1 biology.

**C. Top Tissues/Cell Types**
| Tissue | Count | Percentage | Category |
|--------|-------|------------|----------|
| HepG2 | 6,312 | 21.1% | Cancer |
| K562 | 6,072 | 20.3% | Cancer |
| HUVEC | 3,224 | 10.8% | Endothelial |
| GM12878 | 2,440 | 8.1% | Immune |
| MCF-7 | 2,211 | 7.4% | Cancer |

Cancer cell lines (HepG2, K562) show highest enrichment, potentially reflecting their high AP1 activity and extensive ChIP-seq profiling.

**D. Functional Potential: AP1 × Enhancer Correlation**
- 2D density heatmap of AP1 quantile vs Enhancer quantile scores
- **Spearman ρ = 0.489** (p < 10⁻³⁰)
- White dashed rectangle highlights the high-AP1/high-enhancer quadrant
- **99.6% of variants** fall in the upper-right quadrant (both scores > 0.5)
- This indicates forbidden variants would create functional binding sites in active regulatory regions

### Scientific Interpretation
The mutation spectrum, TF enrichment, and tissue distribution paint a coherent picture: forbidden variants represent mutations that would complete canonical AP1 sites in functionally active chromatin. The strong AP1-enhancer correlation (ρ = 0.49) confirms these aren't random mutations but would create sites in regulatory hotspots.

---

## Figure 3: TF-Tissue Heatmap

**File**: `figures/forbidden_variants/tf_tissue_heatmap.png`

### What It Shows
Two-panel heatmap showing the joint distribution of forbidden variants across TF-tissue combinations.

### Panels

**A. Median AP1 Score by TF-Tissue**
- Shows median predicted binding score for each TF-tissue pair
- Filtered to cells with n≥10 variants for statistical reliability
- Hatched cells indicate insufficient data
- Hierarchical clustering groups similar TFs and tissues

**B. Variant Counts by TF-Tissue**
- Log-scale color showing variant counts per combination
- Marginal totals shown for rows and columns

### Key Observations
- FOS shows consistently high scores across multiple tissues
- HepG2 and K562 show the broadest TF coverage
- Some TF-tissue pairs show specific enrichment patterns

### Scientific Interpretation
The specificity patterns suggest tissue-specific regulatory programs. The clustering reveals TF families with similar binding preferences, consistent with their sequence specificity.

---

## Figure 4: Genomic Distribution

**File**: `figures/forbidden_variants/genomic_distribution.png`

### What It Shows
The chromosomal and genomic location of forbidden variants.

### Panels

**A. Chromosome Distribution**
- Observed vs expected (by chromosome size) variant counts
- Chi-square test for uniformity
- Most chromosomes show proportional distribution

**B. Manhattan-Style Plot**
- Genomic position vs AP1 score across all chromosomes
- Alternating colors by chromosome
- No obvious clustering at specific loci

**C. AP1 Score by Chromosome**
- Box plots showing score distributions per chromosome (top 10)
- Relatively uniform distributions across chromosomes

**D. AP1 Quantile Distribution**
- Histogram of quantile scores with percentile annotations
- Shows the extreme nature of predicted binding (many > 99th percentile)

### Scientific Interpretation
Forbidden variants are distributed relatively uniformly across the genome, suggesting this is a genome-wide phenomenon rather than localized to specific regions. This supports the model of pervasive purifying selection against AP1 site creation.

---

## Figure 5: Selection Signal Strength

**File**: `figures/forbidden_variants/selection_signal_strength.png`

### What It Shows
Evidence for the intensity of purifying selection acting on forbidden variants.

### Panels

**A. Proportion with Extreme AP1 Scores**
| Threshold | Count | Percentage |
|-----------|-------|------------|
| > 0.90 | 28,234 | 94.2% |
| > 0.95 | 24,789 | 82.7% |
| > 0.99 | 7,578 | 25.3% |
| > 0.999 | 1,106 | 3.7% |
| > 0.9999 | 156 | 0.5% |

**B. Allele Number (AN) Distribution**
- Median AN: ~152,000 (indicating coverage of ~76,000 individuals)
- High coverage ensures reliable absence calls
- >95% of variants have AN > 140,000

**C. AP1 Raw Score Distribution (log scale)**
- Mean: 852
- Median: 560
- 99th percentile: 5,712

**D. Top 10 Highest-Impact Candidates**
- Table of variants with highest predicted AP1 binding
- Top variant: chr19:43596274:A>T with AP1 score of 93,952

### Scientific Interpretation
The extreme quantile scores (94% above 90th percentile) combined with high sequencing coverage (>76K individuals) provides strong confidence in both the predicted functional impact and the observed absence of these variants.

---

## Figure 6: Functional Context

**File**: `figures/forbidden_variants/functional_context.png`

### What It Shows
The functional chromatin context of forbidden variants.

### Panels

**A. Enhancer Activity (H3K27ac) Distribution**
- Histogram of enhancer quantile scores
- Shows distribution of predicted enhancer activity

**B. Chromatin Accessibility (ATAC/DNase) Distribution**
- Histogram of accessibility quantile scores
- Indicates open chromatin status

**C. Multi-Feature Scatter**
- Accessibility vs Enhancer, colored by AP1 score
- High-impact quadrant identification
- Shows clustering in accessible, active enhancer regions

**D. High-Priority Multi-Track Candidates**
- Table of variants exceeding all three thresholds:
  - AP1 quantile > 0.99
  - Enhancer quantile > 0.95
  - Accessibility quantile > 0.95

### Scientific Interpretation
Forbidden variants are concentrated in functionally active chromatin—regions with high enhancer marks and open accessibility. This supports the model that selection specifically removes mutations that would create AP1 sites in regulatory active regions, not just anywhere in the genome.

---

## Summary Statistics

### Core Numbers
| Metric | Value |
|--------|-------|
| Total forbidden variants | 29,957 |
| Fold depletion (H=1 vs expected) | 76.5× |
| P-value (binomial) | 4.2 × 10⁻³² |
| Observed H=1 variants | 1 |
| Expected H=1 variants | 76.5 |

### Molecular Properties
| Property | Value |
|----------|-------|
| Transversions | 71.8% |
| Transitions | 28.2% |
| Dominant mutation type | A>T (29.3%) |
| AP1-Enhancer correlation | ρ = 0.489 |

### AP1 Score Distribution
| Metric | Value |
|--------|-------|
| Mean quantile | 0.959 |
| Median quantile | 0.975 |
| % above 99th percentile | 25.3% |
| Mean raw score | 852 |
| Max raw score | 93,952 |

### Top TF Families
1. **FOS/JUN** (FOS, FOSL1, FOSL2, JUN, JUNB, JUND): ~40%
2. **ATF/CREB** (ATF2, ATF3, ATF7, CREB1): ~35%
3. **MAF** (MAFK, MAFF, MAFG): ~8%
4. **NRF/BACH** (NFE2, NRF1, BACH1): ~5%

### Top Tissues
1. HepG2 (hepatocellular carcinoma): 21.1%
2. K562 (chronic myelogenous leukemia): 20.3%
3. HUVEC (endothelial): 10.8%
4. GM12878 (B-lymphoblast): 8.1%
5. MCF-7 (breast cancer): 7.4%

---

## Scientific Conclusions

### 1. Strong Evidence for Purifying Selection
The 76.5-fold depletion of single-mutation-away variants (p < 10⁻³⁰) provides unambiguous evidence that human genomes are under strong purifying selection to prevent the creation of new AP1 binding sites.

### 2. Selection Targets Completion, Not Strength
The similar predicted binding impact between forbidden (H=1) and observed (H=2/H=3) variants reveals that selection removes mutations based on whether they COMPLETE a functional site, not based on predicted binding affinity.

### 3. Functional Specificity
The strong correlation with enhancer activity (ρ = 0.49) and concentration in accessible chromatin indicates selection specifically targets mutations that would create sites in functionally active regulatory regions.

### 4. Conserved Genome Architecture
The genome-wide distribution of forbidden variants suggests this represents a fundamental constraint on human genome evolution—the architecture of potential AP1 sites is actively maintained across all chromosomes.

### 5. AP1 as a "Dangerous" Motif
The extreme selection against AP1 site creation suggests that uncontrolled AP1 binding poses significant fitness costs, likely through dysregulation of downstream target genes involved in cell proliferation, differentiation, and stress response.

---

## Implications

### For Human Genetics
- Rare variants that create novel AP1 sites should be considered potentially pathogenic
- De novo mutations completing AP1 sites may be enriched in disease cases

### For Cancer Biology  
- AP1 is a proto-oncogenic transcription factor
- The strong selection against site creation is consistent with AP1's role in promoting cell proliferation

### For Gene Regulation
- AP1 sites appear to be under tight evolutionary control
- The regulatory grammar of AP1 may be more constrained than previously appreciated

---

## Methods Summary

### Data Sources
- gnomAD v4.1 exomes/genomes (807,162 individuals)
- AlphaGenome predictions for transcription factor binding
- ENCODE/Roadmap epigenomic annotations

### Analysis Pipeline
1. **Module 01-03**: Motif scanning and mutation path generation
2. **Module 04**: gnomAD intersection
3. **Module 05**: AlphaGenome scoring of observed variants
4. **Module 07**: Purifying selection analysis
5. **Module 08**: AlphaGenome scoring of forbidden variants
6. **Module 09**: Visualization and reporting

### Statistical Tests
- Binomial test for variant depletion
- Spearman correlation for score relationships
- Chi-square test for independence
- Mann-Whitney U test for distribution comparisons

---

## Files Generated

```
figures/
├── selection_gradient/
│   ├── selection_gradient.png
│   └── selection_gradient.pdf
└── forbidden_variants/
    ├── forbidden_variants_characterization.png
    ├── forbidden_variants_characterization.pdf
    ├── tf_tissue_heatmap.png
    ├── tf_tissue_heatmap.pdf
    ├── genomic_distribution.png
    ├── genomic_distribution.pdf
    ├── selection_signal_strength.png
    ├── selection_signal_strength.pdf
    ├── functional_context.png
    └── functional_context.pdf
```

---

*Report generated: December 3, 2025*
*Module 09: Visualization Suite*
*Dormant Site Activation Pipeline*
