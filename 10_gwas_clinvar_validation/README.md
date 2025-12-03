# Module 10: GWAS & ClinVar Validation

## Overview

Orthogonal validation of forbidden variants against disease-associated variant databases. This module tests whether positions harboring "forbidden" AP1-activating mutations show enrichment near known disease loci.

## Scientific Rationale

### Core Hypothesis
If forbidden variants (mutations that would create AP1 sites but are absent from 807K humans) are under strong purifying selection due to disease-causing potential, we might expect:

1. **Proximity enrichment**: Forbidden variant positions may be enriched near GWAS lead SNPs (LD-based)
2. **Trait enrichment**: Certain disease categories (cancer, immune, metabolic) may show stronger enrichment
3. **ClinVar depletion**: Forbidden variants should NOT appear in ClinVar (they don't exist in humans)
4. **Regulatory overlap**: GWAS loci near forbidden sites may be enriched for regulatory/enhancer annotations

### Key Questions

| Question | Analysis | Expected if True |
|----------|----------|------------------|
| Are forbidden sites near GWAS loci? | Proximity analysis (±50kb, ±100kb, ±500kb) | OR > 1, enrichment p < 0.05 |
| Which traits are enriched? | Category-stratified analysis | Cancer/immune/metabolic enriched |
| Do forbidden sites overlap ClinVar? | Exact match analysis | ~0 overlaps (expected) |
| Are nearby GWAS loci regulatory? | Annotation overlap | Enhancer/promoter enriched |

## Analysis Strategy

### Step 1: Data Preprocessing
```
GWAS Catalog → Filter to genome-wide significant (p < 5e-8)
             → Extract chr:pos:ref:alt format
             → Build BED file for interval analysis
             
ClinVar VCF  → Extract pathogenic/likely pathogenic
             → Convert to standard coordinates
             → Build BED file for overlap
```

### Step 2: Proximity Analysis (GWAS)
```
For each window size (50kb, 100kb, 500kb):
  - Count forbidden variants within window of each GWAS hit
  - Generate matched null distribution (genome-wide callable sites)
  - Calculate enrichment odds ratio and p-value
  - Stratify by trait category
```

### Step 3: Exact Overlap Analysis
```
Forbidden variants → Intersect with GWAS lead SNPs (exact match)
                   → Intersect with ClinVar (exact match)
                   → Expected: ~0 (these variants don't exist in humans)
```

### Step 4: LD-Aware Analysis
```
GWAS lead SNPs → Expand to LD proxies (r² > 0.8 in 1000G EUR)
              → Intersect with forbidden variant positions
              → Count forbidden sites that would create AP1 at GWAS-tagged loci
```

### Step 5: Trait Enrichment
```
Group GWAS by EFO category:
  - Cancer (neoplasm)
  - Immune/inflammatory
  - Metabolic/cardiovascular
  - Neurological
  - Anthropometric

For each category:
  - Calculate enrichment vs genome-wide background
  - Test for significant departure from null
```

## Input Data

| File | Location | Description |
|------|----------|-------------|
| Forbidden variants | `results/forbidden_variants/AP1/forbidden_variants.tsv` | 29,957 H=1 absent variants |
| Forbidden scores | `results/forbidden_variants/AP1/predictions_summary_ap1.tsv` | AlphaGenome AP1 scores |
| GWAS Catalog | `/mnt/data_1/gwas_data/raw/gwas_catalog_associations.tsv` | ~1M associations |
| ClinVar VCF | `/mnt/data_1/clinvar_data/raw/clinvar.vcf.gz` | ~4.1M variants |

## Output

| File | Description |
|------|-------------|
| `results/gwas_clinvar/AP1/gwas_proximity_enrichment.tsv` | Proximity enrichment by window size |
| `results/gwas_clinvar/AP1/trait_enrichment.tsv` | Enrichment by disease category |
| `results/gwas_clinvar/AP1/gwas_nearby_forbidden.tsv` | Forbidden variants near GWAS hits |
| `results/gwas_clinvar/AP1/clinvar_overlaps.tsv` | ClinVar overlaps (expected ~0) |
| `results/gwas_clinvar/AP1/validation_report.txt` | Summary statistics |
| `figures/gwas_clinvar/` | Visualization figures |

## Scripts

### Main Analysis Pipeline
```bash
# Step 1: Preprocess databases
python 10_gwas_clinvar_validation/preprocess_gwas_catalog.py \
  --input /mnt/data_1/gwas_data/raw/gwas_catalog_associations.tsv \
  --output results/gwas_clinvar/processed/

python 10_gwas_clinvar_validation/preprocess_clinvar.py \
  --input /mnt/data_1/clinvar_data/raw/clinvar.vcf.gz \
  --output results/gwas_clinvar/processed/

# Step 2: Run proximity and overlap analysis
python 10_gwas_clinvar_validation/analyze_forbidden_enrichment.py \
  --forbidden results/forbidden_variants/AP1/forbidden_variants.tsv \
  --gwas results/gwas_clinvar/processed/gwas_significant.bed \
  --clinvar results/gwas_clinvar/processed/clinvar_pathogenic.bed \
  --output results/gwas_clinvar/AP1/

# Step 3: Stratified trait analysis
python 10_gwas_clinvar_validation/stratify_by_trait.py \
  --forbidden results/forbidden_variants/AP1/forbidden_variants.tsv \
  --gwas results/gwas_clinvar/processed/gwas_significant.bed \
  --efo-mapping results/gwas_clinvar/processed/efo_categories.tsv \
  --output results/gwas_clinvar/AP1/

# Step 4: Generate figures
python 10_gwas_clinvar_validation/plot_validation_results.py \
  --input results/gwas_clinvar/AP1/ \
  --output figures/gwas_clinvar/
```

## Statistical Methods

### Enrichment Testing
We use a permutation-based approach with matched null variants:

1. **Null construction**: Sample random positions from gnomAD-callable sites (AN ≥ 100K) matching the forbidden variant chromosome distribution
2. **Enrichment OR**: `(Observed near GWAS) / (Expected near GWAS)`
3. **P-value**: Empirical from 10,000 permutations

### Multiple Testing Correction
- Bonferroni correction for window sizes (3 tests)
- FDR correction (Benjamini-Hochberg) for trait categories

## Expected Results

### If selection hypothesis is correct:
- Forbidden sites should be **depleted** near GWAS loci (positions under selection for constraint)
- OR: Forbidden sites could be **enriched** near GWAS loci (if these loci tag regulatory regions where new TF sites would be pathogenic)

### Controls:
- ClinVar overlap should be ~0 (forbidden variants don't exist in population)
- Observed AP1-activating variants (H=2/H=3) should show different enrichment pattern

## Dependencies

- bedtools (for interval operations)
- pandas, numpy, scipy
- matplotlib, seaborn

## Resource Requirements

| Step | Time | Memory |
|------|------|--------|
| Preprocessing | ~10 min | 8 GB |
| Proximity analysis | ~30 min | 16 GB |
| Trait stratification | ~15 min | 8 GB |
| Visualization | ~5 min | 4 GB |

## References

- GWAS Catalog: [Sollis et al., 2023](https://doi.org/10.1093/nar/gkac1010)
- ClinVar: [Landrum et al., 2020](https://doi.org/10.1093/nar/gkz972)
