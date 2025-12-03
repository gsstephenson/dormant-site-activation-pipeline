#!/usr/bin/env python3
"""
Preprocess ClinVar VCF for enrichment analysis.

Extracts pathogenic and likely pathogenic variants,
creates BED files for interval analysis.
"""

import argparse
import gzip
import pandas as pd
from pathlib import Path
import logging
from collections import Counter

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser(description='Preprocess ClinVar VCF')
    parser.add_argument('--input', required=True, help='Path to clinvar.vcf.gz')
    parser.add_argument('--output', required=True, help='Output directory')
    return parser.parse_args()


def parse_info_field(info_str: str) -> dict:
    """Parse VCF INFO field into dictionary."""
    info = {}
    for item in info_str.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info[key] = value
        else:
            info[item] = True
    return info


def extract_clinvar_variants(vcf_path: str):
    """Extract variants from ClinVar VCF with clinical significance."""
    
    logger.info(f"Parsing ClinVar VCF: {vcf_path}")
    
    variants = []
    sig_counts = Counter()
    
    opener = gzip.open if vcf_path.endswith('.gz') else open
    
    with opener(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue
            
            chrom, pos, var_id, ref, alt, qual, filt, info_str = fields[:8]
            
            # Parse INFO field
            info = parse_info_field(info_str)
            
            # Get clinical significance
            clnsig = info.get('CLNSIG', 'not_provided')
            sig_counts[clnsig] += 1
            
            # Get disease name
            clndn = info.get('CLNDN', 'not_provided')
            
            # Get gene info
            geneinfo = info.get('GENEINFO', '')
            gene = geneinfo.split(':')[0] if geneinfo else ''
            
            # Get molecular consequence
            mc = info.get('MC', '')
            
            # Normalize chromosome
            if not chrom.startswith('chr'):
                chrom = f'chr{chrom}'
            
            variants.append({
                'chr': chrom,
                'pos': int(pos),
                'ref': ref,
                'alt': alt,
                'clnsig': clnsig,
                'clndn': clndn,
                'gene': gene,
                'mc': mc,
                'var_id': var_id,
            })
    
    logger.info(f"Parsed {len(variants):,} variants")
    logger.info("Clinical significance distribution (top 10):")
    for sig, count in sig_counts.most_common(10):
        logger.info(f"  {sig}: {count:,}")
    
    return pd.DataFrame(variants)


def categorize_significance(df: pd.DataFrame) -> pd.DataFrame:
    """Categorize clinical significance into major groups."""
    
    def get_category(clnsig):
        clnsig_lower = str(clnsig).lower()
        
        if 'pathogenic' in clnsig_lower and 'likely' not in clnsig_lower:
            if 'conflicting' in clnsig_lower or 'uncertain' in clnsig_lower:
                return 'conflicting'
            return 'pathogenic'
        elif 'likely_pathogenic' in clnsig_lower:
            return 'likely_pathogenic'
        elif 'benign' in clnsig_lower and 'likely' not in clnsig_lower:
            if 'conflicting' in clnsig_lower:
                return 'conflicting'
            return 'benign'
        elif 'likely_benign' in clnsig_lower:
            return 'likely_benign'
        elif 'uncertain' in clnsig_lower or 'vus' in clnsig_lower:
            return 'vus'
        elif 'conflicting' in clnsig_lower:
            return 'conflicting'
        else:
            return 'other'
    
    df['sig_category'] = df['clnsig'].apply(get_category)
    
    # Log category distribution
    cat_counts = df['sig_category'].value_counts()
    logger.info("Significance category distribution:")
    for cat, count in cat_counts.items():
        logger.info(f"  {cat}: {count:,}")
    
    return df


def create_bed_file(df: pd.DataFrame, output_path: str):
    """Create BED file from variants."""
    
    bed_df = pd.DataFrame({
        'chrom': df['chr'],
        'start': df['pos'] - 1,  # BED is 0-based
        'end': df['pos'],
        'name': df['var_id'],
        'score': 0,
        'strand': '.',
        'ref': df['ref'],
        'alt': df['alt'],
        'clnsig': df['clnsig'],
        'gene': df['gene'],
    })
    
    # Sort
    bed_df = bed_df.sort_values(['chrom', 'start'])
    
    bed_df.to_csv(output_path, sep='\t', index=False, header=False)
    logger.info(f"Wrote BED file to {output_path} ({len(bed_df):,} entries)")


def main():
    args = parse_args()
    
    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Extract variants
    df = extract_clinvar_variants(args.input)
    df = categorize_significance(df)
    
    # Filter to SNVs only (matching our forbidden variants)
    snvs = df[(df['ref'].str.len() == 1) & (df['alt'].str.len() == 1)].copy()
    logger.info(f"SNVs only: {len(snvs):,}")
    
    # Create output files
    
    # 1. All ClinVar variants
    create_bed_file(df, output_dir / 'clinvar_all.bed')
    
    # 2. Pathogenic only
    pathogenic = df[df['sig_category'] == 'pathogenic']
    create_bed_file(pathogenic, output_dir / 'clinvar_pathogenic.bed')
    
    # 3. Pathogenic + Likely pathogenic
    path_likely = df[df['sig_category'].isin(['pathogenic', 'likely_pathogenic'])]
    create_bed_file(path_likely, output_dir / 'clinvar_pathogenic_likely.bed')
    
    # 4. SNVs only (for exact matching with forbidden variants)
    create_bed_file(snvs, output_dir / 'clinvar_snvs.bed')
    
    # 5. Pathogenic SNVs
    path_snvs = snvs[snvs['sig_category'].isin(['pathogenic', 'likely_pathogenic'])]
    create_bed_file(path_snvs, output_dir / 'clinvar_pathogenic_snvs.bed')
    
    # 6. TSV with full information for downstream analysis
    df.to_csv(output_dir / 'clinvar_processed.tsv', sep='\t', index=False)
    
    # 7. Summary statistics
    summary = {
        'total_variants': len(df),
        'snvs': len(snvs),
        'pathogenic': len(df[df['sig_category'] == 'pathogenic']),
        'likely_pathogenic': len(df[df['sig_category'] == 'likely_pathogenic']),
        'vus': len(df[df['sig_category'] == 'vus']),
        'benign': len(df[df['sig_category'] == 'benign']),
        'likely_benign': len(df[df['sig_category'] == 'likely_benign']),
        'unique_genes': df['gene'].nunique(),
        'unique_diseases': df['clndn'].nunique(),
    }
    
    summary_df = pd.DataFrame([summary])
    summary_df.to_csv(output_dir / 'clinvar_summary.tsv', sep='\t', index=False)
    
    logger.info("ClinVar preprocessing complete!")
    logger.info(f"  Total variants: {summary['total_variants']:,}")
    logger.info(f"  Pathogenic: {summary['pathogenic']:,}")
    logger.info(f"  Likely pathogenic: {summary['likely_pathogenic']:,}")


if __name__ == '__main__':
    main()
