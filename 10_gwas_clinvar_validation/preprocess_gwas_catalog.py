#!/usr/bin/env python3
"""
Preprocess GWAS Catalog for enrichment analysis.

Filters to genome-wide significant associations (p < 5e-8),
extracts genomic coordinates, and creates BED files for interval analysis.
"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path
import logging
from collections import Counter

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser(description='Preprocess GWAS Catalog')
    parser.add_argument('--input', required=True, help='Path to gwas_catalog_associations.tsv')
    parser.add_argument('--output', required=True, help='Output directory')
    parser.add_argument('--pvalue-threshold', type=float, default=5e-8,
                        help='P-value threshold for genome-wide significance')
    return parser.parse_args()


def load_gwas_catalog(filepath: str) -> pd.DataFrame:
    """Load GWAS catalog with proper column handling."""
    logger.info(f"Loading GWAS catalog from {filepath}")
    
    # Read with proper column names
    df = pd.read_csv(filepath, sep='\t', low_memory=False)
    
    logger.info(f"Loaded {len(df):,} associations")
    logger.info(f"Columns: {list(df.columns)}")
    
    return df


def filter_significant(df: pd.DataFrame, pvalue_threshold: float) -> pd.DataFrame:
    """Filter to genome-wide significant associations."""
    
    # P-VALUE column
    df['P-VALUE'] = pd.to_numeric(df['P-VALUE'], errors='coerce')
    
    # Filter by significance
    significant = df[df['P-VALUE'] < pvalue_threshold].copy()
    
    logger.info(f"Filtered to {len(significant):,} genome-wide significant associations (p < {pvalue_threshold})")
    
    return significant


def extract_coordinates(df: pd.DataFrame) -> pd.DataFrame:
    """Extract clean chr:pos coordinates."""
    
    # Clean chromosome column
    df['chr'] = df['CHR_ID'].astype(str).str.strip()
    df['chr'] = 'chr' + df['chr'].str.replace('chr', '')
    
    # Clean position column
    df['pos'] = pd.to_numeric(df['CHR_POS'], errors='coerce')
    
    # Filter to valid coordinates
    valid = df[
        (df['chr'].str.match(r'^chr([1-9]|1[0-9]|2[0-2]|X|Y)$')) &
        (df['pos'].notna()) &
        (df['pos'] > 0)
    ].copy()
    
    logger.info(f"Extracted valid coordinates for {len(valid):,} associations")
    
    return valid


def categorize_traits(df: pd.DataFrame) -> pd.DataFrame:
    """Categorize traits into broad disease categories."""
    
    # Define category keywords
    categories = {
        'cancer': ['cancer', 'carcinoma', 'tumor', 'neoplasm', 'leukemia', 'lymphoma', 
                   'melanoma', 'myeloma', 'sarcoma', 'oncolog'],
        'immune': ['immune', 'autoimmune', 'inflammatory', 'arthritis', 'lupus', 
                   'crohn', 'colitis', 'asthma', 'allerg', 'celiac', 'psoriasis'],
        'cardiovascular': ['heart', 'cardiac', 'coronary', 'myocardial', 'atheroscler',
                          'hypertension', 'blood pressure', 'stroke', 'vascular'],
        'metabolic': ['diabetes', 'glucose', 'insulin', 'obesity', 'bmi', 'lipid',
                     'cholesterol', 'triglyceride', 'metabolic'],
        'neurological': ['alzheimer', 'parkinson', 'schizophren', 'bipolar', 'depress',
                        'autism', 'epilepsy', 'migraine', 'neurolog', 'cognitive'],
        'anthropometric': ['height', 'weight', 'bmi', 'body mass', 'waist', 'hip ratio'],
        'hematological': ['platelet', 'hemoglobin', 'blood cell', 'anemia', 'coagul',
                         'thrombos', 'erythrocyte', 'leukocyte'],
    }
    
    def get_category(trait):
        if pd.isna(trait):
            return 'other'
        trait_lower = str(trait).lower()
        for category, keywords in categories.items():
            if any(kw in trait_lower for kw in keywords):
                return category
        return 'other'
    
    df['trait_category'] = df['DISEASE/TRAIT'].apply(get_category)
    
    # Log category distribution
    cat_counts = df['trait_category'].value_counts()
    logger.info("Trait category distribution:")
    for cat, count in cat_counts.items():
        logger.info(f"  {cat}: {count:,}")
    
    return df


def create_bed_file(df: pd.DataFrame, output_path: str, window: int = 0):
    """Create BED file from coordinates with optional window extension."""
    
    bed_df = pd.DataFrame({
        'chrom': df['chr'],
        'start': (df['pos'] - 1 - window).clip(lower=0).astype(int),  # BED is 0-based
        'end': (df['pos'] + window).astype(int),
        'name': df['SNPS'].fillna('.'),
        'score': -np.log10(df['P-VALUE'].clip(lower=1e-300)),  # -log10(p) as score
        'strand': '.',
        'trait': df['DISEASE/TRAIT'].fillna('.').str.slice(0, 100),
        'category': df['trait_category']
    })
    
    # Sort by chromosome and position
    bed_df = bed_df.sort_values(['chrom', 'start'])
    
    # Write BED file
    bed_df.to_csv(output_path, sep='\t', index=False, header=False)
    logger.info(f"Wrote BED file to {output_path} ({len(bed_df):,} entries)")
    
    return bed_df


def create_unique_loci(df: pd.DataFrame, output_path: str, merge_window: int = 1000):
    """Create unique loci by merging nearby SNPs."""
    
    # Sort and deduplicate by position
    unique_pos = df.drop_duplicates(subset=['chr', 'pos']).sort_values(['chr', 'pos'])
    
    logger.info(f"Unique positions: {len(unique_pos):,} (from {len(df):,} associations)")
    
    # Create simple BED
    bed_df = pd.DataFrame({
        'chrom': unique_pos['chr'],
        'start': (unique_pos['pos'] - 1).astype(int),
        'end': unique_pos['pos'].astype(int),
        'name': unique_pos['SNPS'].fillna('.'),
    })
    
    bed_df.to_csv(output_path, sep='\t', index=False, header=False)
    logger.info(f"Wrote unique loci BED to {output_path}")
    
    return unique_pos


def write_efo_mapping(df: pd.DataFrame, output_path: str):
    """Write trait-to-category mapping file."""
    
    # Get unique traits with categories
    trait_mapping = df[['DISEASE/TRAIT', 'trait_category']].drop_duplicates()
    trait_mapping = trait_mapping.sort_values('trait_category')
    
    trait_mapping.to_csv(output_path, sep='\t', index=False)
    logger.info(f"Wrote EFO category mapping to {output_path}")


def main():
    args = parse_args()
    
    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load and process
    df = load_gwas_catalog(args.input)
    df = filter_significant(df, args.pvalue_threshold)
    df = extract_coordinates(df)
    df = categorize_traits(df)
    
    # Create output files
    
    # 1. Full BED with all significant associations
    create_bed_file(
        df, 
        output_dir / 'gwas_significant.bed'
    )
    
    # 2. Unique loci BED (deduplicated by position)
    unique_loci = create_unique_loci(
        df,
        output_dir / 'gwas_unique_loci.bed'
    )
    
    # 3. BED files with proximity windows
    for window in [50000, 100000, 500000]:
        create_bed_file(
            df,
            output_dir / f'gwas_significant_{window//1000}kb.bed',
            window=window
        )
    
    # 4. Category-specific BED files
    for category in df['trait_category'].unique():
        cat_df = df[df['trait_category'] == category]
        if len(cat_df) >= 100:  # Only create if substantial
            create_bed_file(
                cat_df,
                output_dir / f'gwas_{category}.bed'
            )
    
    # 5. EFO category mapping
    write_efo_mapping(df, output_dir / 'efo_categories.tsv')
    
    # 6. Summary statistics
    summary = {
        'total_associations': len(df),
        'unique_snps': df['SNPS'].nunique(),
        'unique_positions': len(unique_loci),
        'traits': df['DISEASE/TRAIT'].nunique(),
        'pvalue_threshold': args.pvalue_threshold,
    }
    
    # Add category counts
    for cat in df['trait_category'].value_counts().index:
        summary[f'n_{cat}'] = int(df[df['trait_category'] == cat].shape[0])
    
    summary_df = pd.DataFrame([summary])
    summary_df.to_csv(output_dir / 'gwas_summary.tsv', sep='\t', index=False)
    
    logger.info("GWAS preprocessing complete!")
    logger.info(f"  Total significant associations: {summary['total_associations']:,}")
    logger.info(f"  Unique positions: {summary['unique_positions']:,}")
    logger.info(f"  Unique traits: {summary['traits']:,}")


if __name__ == '__main__':
    main()
