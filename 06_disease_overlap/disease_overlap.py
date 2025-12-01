#!/usr/bin/env python3
"""
Disease Variant Overlap Analysis for Dormant Site Activation Pipeline

This module intersects dormant AP1 site-activating variants with:
1. ClinVar pathogenic/likely pathogenic variants
2. GWAS Catalog disease-associated variants

This validates whether variants that could activate dormant TF binding sites
are implicated in human disease.
"""

import argparse
import gzip
import logging
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def parse_clinvar_vcf(clinvar_path: str) -> pd.DataFrame:
    """
    Parse ClinVar VCF file and extract relevant variant information.
    
    Args:
        clinvar_path: Path to ClinVar VCF file (can be gzipped)
        
    Returns:
        DataFrame with ClinVar variants
    """
    logger.info(f"Parsing ClinVar VCF: {clinvar_path}")
    
    variants = []
    open_func = gzip.open if clinvar_path.endswith('.gz') else open
    
    with open_func(clinvar_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue
            
            chrom, pos, var_id, ref, alt, qual, filt, info = fields[:8]
            
            # Parse INFO field
            info_dict = {}
            for item in info.split(';'):
                if '=' in item:
                    key, value = item.split('=', 1)
                    info_dict[key] = value
            
            # Get clinical significance
            clnsig = info_dict.get('CLNSIG', '')
            clndn = info_dict.get('CLNDN', 'not_specified')
            geneinfo = info_dict.get('GENEINFO', '')
            rs_id = info_dict.get('RS', '')
            clnrevstat = info_dict.get('CLNREVSTAT', '')
            
            # Normalize chromosome
            if not chrom.startswith('chr'):
                chrom = f'chr{chrom}'
            
            # Create variant ID
            variant_id = f"{chrom}:{pos}:{ref}:{alt}"
            
            # Parse gene info
            genes = []
            if geneinfo:
                for gene_entry in geneinfo.split('|'):
                    if ':' in gene_entry:
                        gene_symbol = gene_entry.split(':')[0]
                        genes.append(gene_symbol)
            
            variants.append({
                'variant_id': variant_id,
                'chr': chrom,
                'pos': int(pos),
                'ref': ref,
                'alt': alt,
                'clinvar_id': var_id,
                'rs_id': f'rs{rs_id}' if rs_id else '',
                'clinsig': clnsig,
                'disease': clndn.replace('_', ' '),
                'genes': '|'.join(genes),
                'review_status': clnrevstat.replace('_', ' ')
            })
    
    df = pd.DataFrame(variants)
    logger.info(f"Parsed {len(df):,} ClinVar variants")
    return df


def parse_gwas_catalog(gwas_path: str) -> pd.DataFrame:
    """
    Parse GWAS Catalog TSV file.
    
    Args:
        gwas_path: Path to GWAS Catalog associations file
        
    Returns:
        DataFrame with GWAS variants
    """
    logger.info(f"Parsing GWAS Catalog: {gwas_path}")
    
    # Read GWAS catalog
    df = pd.read_csv(gwas_path, sep='\t', low_memory=False, on_bad_lines='skip')
    
    # Relevant columns - note column names may vary
    # Common format: CHR_ID, CHR_POS, SNPS, DISEASE/TRAIT, P-VALUE, etc.
    
    variants = []
    
    for _, row in df.iterrows():
        try:
            # Get chromosome and position
            chrom = str(row.get('CHR_ID', '')).strip()
            pos = row.get('CHR_POS', '')
            
            if not chrom or chrom == 'nan' or not pos or str(pos) == 'nan':
                continue
            
            # Clean up chromosome
            chrom = str(chrom).replace('chr', '')
            if chrom in ['X', 'Y', 'MT'] or chrom.isdigit():
                chrom = f'chr{chrom}'
            else:
                continue
            
            # Get position
            try:
                pos = int(float(pos))
            except (ValueError, TypeError):
                continue
            
            # Get SNP ID
            snp_id = str(row.get('SNPS', '')).strip()
            
            # Get risk allele
            risk_allele = str(row.get('STRONGEST SNP-RISK ALLELE', '')).strip()
            
            # Get disease/trait
            disease = str(row.get('DISEASE/TRAIT', '')).strip()
            
            # Get p-value
            pvalue = row.get('P-VALUE', '')
            try:
                pvalue = float(pvalue)
            except (ValueError, TypeError):
                pvalue = np.nan
            
            # Get mapped genes
            mapped_genes = str(row.get('MAPPED_GENE', '')).strip()
            reported_genes = str(row.get('REPORTED GENE(S)', '')).strip()
            
            # Get study info
            pubmed_id = str(row.get('PUBMEDID', '')).strip()
            first_author = str(row.get('FIRST AUTHOR', '')).strip()
            journal = str(row.get('JOURNAL', '')).strip()
            
            variants.append({
                'chr': chrom,
                'pos': pos,
                'snp_id': snp_id,
                'risk_allele': risk_allele,
                'disease_trait': disease,
                'pvalue': pvalue,
                'mapped_genes': mapped_genes,
                'reported_genes': reported_genes,
                'pubmed_id': pubmed_id,
                'first_author': first_author,
                'journal': journal
            })
            
        except Exception as e:
            continue
    
    df_out = pd.DataFrame(variants)
    
    # Remove duplicates
    df_out = df_out.drop_duplicates(subset=['chr', 'pos', 'disease_trait'])
    
    logger.info(f"Parsed {len(df_out):,} GWAS associations")
    return df_out


def load_ap1_variants(landscape_path: str) -> pd.DataFrame:
    """
    Load AP1 activation landscape variants.
    
    Args:
        landscape_path: Path to AP1_activation_landscape.tsv
        
    Returns:
        DataFrame with AP1 variants
    """
    logger.info(f"Loading AP1 activation landscape: {landscape_path}")
    
    df = pd.read_csv(landscape_path, sep='\t')
    
    # Parse variant_id to get chr, pos, ref, alt
    def parse_variant_id(vid):
        parts = vid.replace('>', ':').split(':')
        if len(parts) >= 4:
            return parts[0], int(parts[1]), parts[2], parts[3]
        return None, None, None, None
    
    parsed = df['variant_id_str'].apply(parse_variant_id)
    df['chr'] = [p[0] for p in parsed]
    df['pos'] = [p[1] for p in parsed]
    df['ref'] = [p[2] for p in parsed]
    df['alt'] = [p[3] for p in parsed]
    
    logger.info(f"Loaded {len(df):,} AP1 variants")
    return df


def intersect_with_clinvar(
    ap1_df: pd.DataFrame,
    clinvar_df: pd.DataFrame,
    pathogenic_only: bool = True
) -> pd.DataFrame:
    """
    Intersect AP1 variants with ClinVar.
    
    Args:
        ap1_df: AP1 activation landscape DataFrame
        clinvar_df: ClinVar variants DataFrame
        pathogenic_only: If True, only keep pathogenic/likely pathogenic
        
    Returns:
        DataFrame with overlapping variants
    """
    logger.info("Intersecting AP1 variants with ClinVar...")
    
    # Filter ClinVar for pathogenic variants if requested
    if pathogenic_only:
        pathogenic_terms = [
            'Pathogenic',
            'Likely_pathogenic',
            'Pathogenic/Likely_pathogenic'
        ]
        clinvar_filtered = clinvar_df[
            clinvar_df['clinsig'].apply(
                lambda x: any(term in str(x) for term in pathogenic_terms)
            )
        ].copy()
        logger.info(f"Filtered to {len(clinvar_filtered):,} pathogenic/likely pathogenic ClinVar variants")
    else:
        clinvar_filtered = clinvar_df.copy()
    
    # Create variant key for matching
    ap1_df['variant_key'] = ap1_df['chr'] + ':' + ap1_df['pos'].astype(str)
    clinvar_filtered['variant_key'] = clinvar_filtered['chr'] + ':' + clinvar_filtered['pos'].astype(str)
    
    # Also try exact match with ref/alt
    ap1_df['variant_exact'] = ap1_df['variant_id_str']
    clinvar_filtered['variant_exact'] = clinvar_filtered['variant_id']
    
    # Merge on position first (some alleles might differ)
    overlap_pos = pd.merge(
        ap1_df,
        clinvar_filtered,
        on='variant_key',
        how='inner',
        suffixes=('_ap1', '_clinvar')
    )
    
    # Exact matches
    overlap_exact = pd.merge(
        ap1_df,
        clinvar_filtered,
        left_on='variant_exact',
        right_on='variant_exact',
        how='inner',
        suffixes=('_ap1', '_clinvar')
    )
    
    # Combine and deduplicate
    if len(overlap_exact) > 0:
        all_overlaps = pd.concat([overlap_pos, overlap_exact]).drop_duplicates(
            subset=['variant_id_str']
        )
    else:
        all_overlaps = overlap_pos
    
    logger.info(f"Found {len(all_overlaps):,} AP1 variants overlapping ClinVar")
    
    return all_overlaps


def intersect_with_gwas(
    ap1_df: pd.DataFrame,
    gwas_df: pd.DataFrame,
    window_size: int = 1000
) -> pd.DataFrame:
    """
    Intersect AP1 variants with GWAS Catalog.
    
    Uses a window-based approach since GWAS variants may not have exact
    allele information.
    
    Args:
        ap1_df: AP1 activation landscape DataFrame
        gwas_df: GWAS catalog DataFrame
        window_size: Window size for position matching (bp)
        
    Returns:
        DataFrame with overlapping variants
    """
    logger.info(f"Intersecting AP1 variants with GWAS (window={window_size}bp)...")
    
    overlaps = []
    
    # Group by chromosome for efficiency
    ap1_by_chr = ap1_df.groupby('chr')
    gwas_by_chr = gwas_df.groupby('chr')
    
    for chrom in ap1_by_chr.groups.keys():
        if chrom not in gwas_by_chr.groups:
            continue
        
        ap1_chr = ap1_by_chr.get_group(chrom)
        gwas_chr = gwas_by_chr.get_group(chrom).copy()
        
        # Create position index for GWAS
        gwas_positions = gwas_chr['pos'].values
        
        for _, ap1_row in ap1_chr.iterrows():
            ap1_pos = ap1_row['pos']
            
            # Find GWAS variants within window
            close_idx = np.where(
                np.abs(gwas_positions - ap1_pos) <= window_size
            )[0]
            
            for idx in close_idx:
                gwas_row = gwas_chr.iloc[idx]
                distance = abs(ap1_pos - gwas_row['pos'])
                
                overlaps.append({
                    'variant_id_str': ap1_row['variant_id_str'],
                    'chr': chrom,
                    'ap1_pos': ap1_pos,
                    'gwas_pos': gwas_row['pos'],
                    'distance': distance,
                    'snp_id': gwas_row['snp_id'],
                    'disease_trait': gwas_row['disease_trait'],
                    'pvalue': gwas_row['pvalue'],
                    'mapped_genes': gwas_row['mapped_genes'],
                    'pubmed_id': gwas_row['pubmed_id'],
                    'ap1_max_score': ap1_row['ap1_max_score'],
                    'ap1_max_raw_score': ap1_row['ap1_max_raw_score'],
                    'ap1_best_tf': ap1_row['ap1_best_tf'],
                    'gnomad_AF': ap1_row['gnomad_AF'],
                    'quadrant': ap1_row['quadrant'],
                    'hamming_distance': ap1_row['hamming_distance']
                })
    
    overlap_df = pd.DataFrame(overlaps)
    
    if len(overlap_df) > 0:
        # Sort by distance and p-value
        overlap_df = overlap_df.sort_values(['distance', 'pvalue'])
        
        # Get unique AP1 variants with their best GWAS match
        unique_variants = overlap_df.groupby('variant_id_str').first().reset_index()
        
        logger.info(f"Found {len(unique_variants):,} AP1 variants near GWAS loci")
        logger.info(f"Total {len(overlap_df):,} variant-trait associations")
    else:
        unique_variants = overlap_df
        logger.info("No GWAS overlaps found")
    
    return overlap_df, unique_variants


def categorize_diseases(diseases: list) -> dict:
    """
    Categorize diseases into major groups for visualization.
    """
    categories = {
        'Cancer': ['cancer', 'carcinoma', 'tumor', 'melanoma', 'leukemia', 
                   'lymphoma', 'sarcoma', 'neoplasm', 'oncolog'],
        'Cardiovascular': ['heart', 'cardiac', 'cardio', 'coronary', 'artery',
                          'blood pressure', 'hypertension', 'stroke', 'atrial'],
        'Metabolic': ['diabetes', 'metabolic', 'obesity', 'lipid', 'cholesterol',
                     'glucose', 'insulin', 'thyroid'],
        'Neurological': ['alzheimer', 'parkinson', 'schizophrenia', 'autism',
                        'epilepsy', 'neurodegenerative', 'dementia', 'brain'],
        'Autoimmune': ['autoimmune', 'lupus', 'arthritis', 'crohn', 'colitis',
                      'multiple sclerosis', 'psoriasis', 'inflammatory'],
        'Respiratory': ['asthma', 'copd', 'lung', 'pulmonary', 'respiratory'],
        'Other': []
    }
    
    disease_counts = defaultdict(int)
    disease_mapping = {}
    
    for disease in diseases:
        disease_lower = str(disease).lower()
        categorized = False
        
        for category, keywords in categories.items():
            if category == 'Other':
                continue
            for keyword in keywords:
                if keyword in disease_lower:
                    disease_counts[category] += 1
                    disease_mapping[disease] = category
                    categorized = True
                    break
            if categorized:
                break
        
        if not categorized:
            disease_counts['Other'] += 1
            disease_mapping[disease] = 'Other'
    
    return dict(disease_counts), disease_mapping


def generate_summary_report(
    ap1_df: pd.DataFrame,
    clinvar_overlaps: pd.DataFrame,
    gwas_overlaps: pd.DataFrame,
    gwas_unique: pd.DataFrame,
    output_dir: str
) -> str:
    """
    Generate a summary report of disease overlap analysis.
    """
    report = []
    report.append("=" * 80)
    report.append("DISEASE VARIANT OVERLAP ANALYSIS - SUMMARY REPORT")
    report.append("=" * 80)
    report.append("")
    
    # Overall statistics
    report.append("## INPUT DATA")
    report.append(f"Total AP1 activating variants analyzed: {len(ap1_df):,}")
    report.append(f"  - High priority variants: {len(ap1_df[ap1_df['quadrant'].str.contains('HIGH PRIORITY', na=False)]):,}")
    report.append("")
    
    # ClinVar results
    report.append("## CLINVAR OVERLAP")
    report.append(f"AP1 variants with ClinVar pathogenic annotations: {len(clinvar_overlaps):,}")
    
    if len(clinvar_overlaps) > 0:
        # Disease distribution
        diseases = clinvar_overlaps['disease'].dropna().unique()
        report.append(f"Unique diseases represented: {len(diseases)}")
        
        # Top diseases
        disease_counts = clinvar_overlaps['disease'].value_counts().head(10)
        report.append("\nTop 10 diseases in ClinVar overlaps:")
        for disease, count in disease_counts.items():
            report.append(f"  - {disease}: {count}")
        
        # Clinical significance breakdown
        report.append("\nClinical significance breakdown:")
        sig_counts = clinvar_overlaps['clinsig'].value_counts()
        for sig, count in sig_counts.items():
            report.append(f"  - {sig}: {count}")
    
    report.append("")
    
    # GWAS results
    report.append("## GWAS CATALOG OVERLAP")
    report.append(f"AP1 variants near GWAS loci: {len(gwas_unique):,}")
    report.append(f"Total variant-trait associations: {len(gwas_overlaps):,}")
    
    if len(gwas_overlaps) > 0:
        # Trait distribution
        traits = gwas_overlaps['disease_trait'].dropna().unique()
        report.append(f"Unique traits represented: {len(traits)}")
        
        # Categorize diseases
        trait_cats, _ = categorize_diseases(gwas_overlaps['disease_trait'].dropna())
        report.append("\nTrait categories:")
        for cat, count in sorted(trait_cats.items(), key=lambda x: -x[1]):
            report.append(f"  - {cat}: {count}")
        
        # Top traits
        trait_counts = gwas_overlaps['disease_trait'].value_counts().head(15)
        report.append("\nTop 15 disease/traits in GWAS overlaps:")
        for trait, count in trait_counts.items():
            report.append(f"  - {trait}: {count}")
        
        # Distance distribution
        report.append("\nDistance to GWAS lead variant:")
        report.append(f"  - Exact match (0bp): {len(gwas_overlaps[gwas_overlaps['distance'] == 0]):,}")
        report.append(f"  - Within 100bp: {len(gwas_overlaps[gwas_overlaps['distance'] <= 100]):,}")
        report.append(f"  - Within 500bp: {len(gwas_overlaps[gwas_overlaps['distance'] <= 500]):,}")
        report.append(f"  - Within 1000bp: {len(gwas_overlaps[gwas_overlaps['distance'] <= 1000]):,}")
        
        # Significance
        significant = gwas_overlaps[gwas_overlaps['pvalue'] < 5e-8]
        report.append(f"\nGenome-wide significant associations (p<5e-8): {len(significant):,}")
        
        # High priority AP1 variants in GWAS
        high_priority = gwas_overlaps[
            gwas_overlaps['variant_id_str'].isin(
                ap1_df[ap1_df['quadrant'].str.contains('HIGH PRIORITY', na=False)]['variant_id_str']
            )
        ]
        report.append(f"\nHigh-priority AP1 variants near GWAS loci: {len(high_priority['variant_id_str'].unique()):,}")
    
    report.append("")
    
    # Key findings
    report.append("## KEY FINDINGS")
    
    # Calculate enrichment (if we had expected values)
    total_genome_size = 3e9  # approximate
    ap1_variant_count = len(ap1_df)
    clinvar_overlap_rate = len(clinvar_overlaps) / ap1_variant_count * 100 if ap1_variant_count > 0 else 0
    gwas_overlap_rate = len(gwas_unique) / ap1_variant_count * 100 if ap1_variant_count > 0 else 0
    
    report.append(f"ClinVar overlap rate: {clinvar_overlap_rate:.2f}%")
    report.append(f"GWAS overlap rate: {gwas_overlap_rate:.2f}%")
    
    if len(gwas_overlaps) > 0:
        # Highlight interesting hits
        report.append("\n### Notable Disease Associations:")
        
        # Get high-confidence hits
        high_conf = gwas_overlaps[
            (gwas_overlaps['pvalue'] < 1e-10) & 
            (gwas_overlaps['distance'] <= 100) &
            (gwas_overlaps['ap1_max_raw_score'] > 1000)
        ].head(10)
        
        for _, row in high_conf.iterrows():
            report.append(f"\n  {row['variant_id_str']}")
            report.append(f"    - Trait: {row['disease_trait']}")
            report.append(f"    - GWAS SNP: {row['snp_id']} (p={row['pvalue']:.2e})")
            report.append(f"    - Distance: {row['distance']}bp")
            report.append(f"    - AP1 raw score: {row['ap1_max_raw_score']:.1f}")
            report.append(f"    - Best TF: {row['ap1_best_tf']}")
    
    report.append("")
    report.append("=" * 80)
    
    return '\n'.join(report)


def create_visualizations(
    ap1_df: pd.DataFrame,
    clinvar_overlaps: pd.DataFrame,
    gwas_overlaps: pd.DataFrame,
    output_dir: str
):
    """
    Create visualization plots for disease overlap analysis.
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib
        matplotlib.use('Agg')
    except ImportError:
        logger.warning("matplotlib not available, skipping visualizations")
        return
    
    fig_dir = os.path.join(output_dir, 'figures')
    os.makedirs(fig_dir, exist_ok=True)
    
    # Figure 1: Disease category pie chart
    if len(gwas_overlaps) > 0:
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        
        # GWAS trait categories
        trait_cats, _ = categorize_diseases(gwas_overlaps['disease_trait'].dropna())
        
        # Sort and plot
        sorted_cats = sorted(trait_cats.items(), key=lambda x: -x[1])
        labels = [c[0] for c in sorted_cats]
        sizes = [c[1] for c in sorted_cats]
        
        colors = plt.cm.Set3(np.linspace(0, 1, len(labels)))
        
        wedges, texts, autotexts = axes[0].pie(
            sizes, labels=labels, autopct='%1.1f%%',
            colors=colors, pctdistance=0.85
        )
        axes[0].set_title('GWAS Trait Categories', fontsize=12, fontweight='bold')
        
        # Top traits bar chart
        trait_counts = gwas_overlaps['disease_trait'].value_counts().head(15)
        y_pos = np.arange(len(trait_counts))
        
        axes[1].barh(y_pos, trait_counts.values, color='steelblue', alpha=0.7)
        axes[1].set_yticks(y_pos)
        axes[1].set_yticklabels([t[:40] + '...' if len(t) > 40 else t for t in trait_counts.index])
        axes[1].invert_yaxis()
        axes[1].set_xlabel('Number of Associations')
        axes[1].set_title('Top 15 GWAS Traits', fontsize=12, fontweight='bold')
        
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, 'gwas_trait_distribution.png'), dpi=150)
        plt.close()
        
        logger.info("Created GWAS trait distribution figure")
    
    # Figure 2: Distance distribution for GWAS overlaps
    if len(gwas_overlaps) > 0:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        
        # Distance histogram
        distances = gwas_overlaps['distance'].values
        axes[0].hist(distances, bins=50, color='steelblue', alpha=0.7, edgecolor='black')
        axes[0].set_xlabel('Distance to GWAS Lead Variant (bp)')
        axes[0].set_ylabel('Count')
        axes[0].set_title('Distance Distribution')
        axes[0].axvline(x=0, color='red', linestyle='--', label='Exact match')
        axes[0].legend()
        
        # P-value distribution
        pvals = gwas_overlaps['pvalue'].dropna()
        log_pvals = -np.log10(pvals[pvals > 0])
        axes[1].hist(log_pvals, bins=50, color='coral', alpha=0.7, edgecolor='black')
        axes[1].axvline(x=-np.log10(5e-8), color='red', linestyle='--', 
                       label=f'Genome-wide significance')
        axes[1].set_xlabel('-log10(p-value)')
        axes[1].set_ylabel('Count')
        axes[1].set_title('GWAS P-value Distribution')
        axes[1].legend()
        
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, 'gwas_distance_pvalue.png'), dpi=150)
        plt.close()
        
        logger.info("Created GWAS distance/p-value figure")
    
    # Figure 3: AP1 score vs GWAS significance
    if len(gwas_overlaps) > 0:
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Filter for plottable data
        plot_df = gwas_overlaps[
            (gwas_overlaps['pvalue'] > 0) & 
            (gwas_overlaps['ap1_max_raw_score'] > 0)
        ].copy()
        
        if len(plot_df) > 0:
            x = np.log10(plot_df['ap1_max_raw_score'])
            y = -np.log10(plot_df['pvalue'])
            
            # Color by distance
            colors = plot_df['distance'].values
            
            scatter = ax.scatter(x, y, c=colors, cmap='viridis_r', 
                                alpha=0.6, s=30, edgecolors='none')
            
            plt.colorbar(scatter, label='Distance to GWAS SNP (bp)')
            
            ax.axhline(y=-np.log10(5e-8), color='red', linestyle='--', 
                      alpha=0.7, label='Genome-wide significance')
            
            ax.set_xlabel('log₁₀(AP1 Raw Score)', fontsize=12)
            ax.set_ylabel('-log₁₀(GWAS p-value)', fontsize=12)
            ax.set_title('AP1 Activation Score vs GWAS Significance', 
                        fontsize=14, fontweight='bold')
            ax.legend()
            
            plt.tight_layout()
            plt.savefig(os.path.join(fig_dir, 'ap1_vs_gwas_significance.png'), dpi=150)
            plt.close()
            
            logger.info("Created AP1 vs GWAS significance figure")
    
    # Figure 4: Overlap summary
    fig, ax = plt.subplots(figsize=(8, 6))
    
    total = len(ap1_df)
    clinvar_count = len(clinvar_overlaps['variant_id_str'].unique()) if len(clinvar_overlaps) > 0 else 0
    gwas_count = len(gwas_overlaps['variant_id_str'].unique()) if len(gwas_overlaps) > 0 else 0
    
    # Calculate overlaps between ClinVar and GWAS if both exist
    if len(clinvar_overlaps) > 0 and len(gwas_overlaps) > 0:
        both = len(set(clinvar_overlaps['variant_id_str']) & set(gwas_overlaps['variant_id_str']))
    else:
        both = 0
    
    categories = ['Total AP1\nVariants', 'ClinVar\nOverlap', 'GWAS\nOverlap', 'Both\nDatabases']
    counts = [total, clinvar_count, gwas_count, both]
    colors = ['steelblue', 'coral', 'mediumseagreen', 'gold']
    
    bars = ax.bar(categories, counts, color=colors, alpha=0.7, edgecolor='black')
    
    # Add count labels
    for bar, count in zip(bars, counts):
        height = bar.get_height()
        ax.annotate(f'{count:,}',
                   xy=(bar.get_x() + bar.get_width()/2, height),
                   xytext=(0, 3),
                   textcoords="offset points",
                   ha='center', va='bottom', fontweight='bold')
    
    ax.set_ylabel('Number of Variants', fontsize=12)
    ax.set_title('Disease Database Overlap Summary', fontsize=14, fontweight='bold')
    
    # Add percentage labels
    if total > 0:
        ax.text(1, clinvar_count + total*0.02, f'({clinvar_count/total*100:.1f}%)', 
               ha='center', fontsize=10)
        ax.text(2, gwas_count + total*0.02, f'({gwas_count/total*100:.1f}%)', 
               ha='center', fontsize=10)
    
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, 'overlap_summary.png'), dpi=150)
    plt.close()
    
    logger.info("Created overlap summary figure")


def main():
    parser = argparse.ArgumentParser(
        description='Disease Variant Overlap Analysis for Dormant Site Activation'
    )
    parser.add_argument(
        '--ap1-landscape',
        required=True,
        help='Path to AP1_activation_landscape.tsv'
    )
    parser.add_argument(
        '--clinvar-vcf',
        required=True,
        help='Path to ClinVar VCF file (can be gzipped)'
    )
    parser.add_argument(
        '--gwas-catalog',
        required=True,
        help='Path to GWAS Catalog associations TSV file'
    )
    parser.add_argument(
        '--output-dir',
        required=True,
        help='Output directory for results'
    )
    parser.add_argument(
        '--gwas-window',
        type=int,
        default=1000,
        help='Window size for GWAS position matching (default: 1000bp)'
    )
    parser.add_argument(
        '--all-clinvar',
        action='store_true',
        help='Include all ClinVar variants, not just pathogenic'
    )
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Set up file logging
    log_file = os.path.join(args.output_dir, 'disease_overlap.log')
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logger.addHandler(file_handler)
    
    logger.info("Starting Disease Variant Overlap Analysis")
    logger.info(f"Arguments: {args}")
    
    # Load data
    ap1_df = load_ap1_variants(args.ap1_landscape)
    clinvar_df = parse_clinvar_vcf(args.clinvar_vcf)
    gwas_df = parse_gwas_catalog(args.gwas_catalog)
    
    # Perform intersections
    clinvar_overlaps = intersect_with_clinvar(
        ap1_df, clinvar_df, 
        pathogenic_only=not args.all_clinvar
    )
    
    gwas_overlaps, gwas_unique = intersect_with_gwas(
        ap1_df, gwas_df,
        window_size=args.gwas_window
    )
    
    # Save results
    logger.info("Saving results...")
    
    # ClinVar overlaps
    clinvar_out = os.path.join(args.output_dir, 'clinvar_overlaps.tsv')
    clinvar_overlaps.to_csv(clinvar_out, sep='\t', index=False)
    logger.info(f"Saved ClinVar overlaps to {clinvar_out}")
    
    # GWAS overlaps
    gwas_out = os.path.join(args.output_dir, 'gwas_overlaps.tsv')
    gwas_overlaps.to_csv(gwas_out, sep='\t', index=False)
    logger.info(f"Saved GWAS overlaps to {gwas_out}")
    
    # GWAS unique variants
    gwas_unique_out = os.path.join(args.output_dir, 'gwas_unique_variants.tsv')
    gwas_unique.to_csv(gwas_unique_out, sep='\t', index=False)
    logger.info(f"Saved GWAS unique variants to {gwas_unique_out}")
    
    # Generate summary report
    report = generate_summary_report(
        ap1_df, clinvar_overlaps, gwas_overlaps, gwas_unique, args.output_dir
    )
    
    report_path = os.path.join(args.output_dir, 'disease_overlap_report.txt')
    with open(report_path, 'w') as f:
        f.write(report)
    logger.info(f"Saved summary report to {report_path}")
    
    # Print report to console
    print("\n" + report)
    
    # Create visualizations
    create_visualizations(ap1_df, clinvar_overlaps, gwas_overlaps, args.output_dir)
    
    logger.info("Disease overlap analysis complete!")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
