#!/usr/bin/env python3
"""
Module 08: Compute AP1-Specific Scores from Forbidden Variant Predictions

This script post-processes the parquet partitions from score_forbidden_variants.py
to extract AP1-family-specific scores for downstream analysis and visualization.

The main scoring script saves all ~800M predictions to parquet, but the summary
only contains global max scores. This script computes AP1-specific metrics to
enable direct comparison with Module 05's observed variant AP1 scores.

Output: predictions_summary_ap1.tsv with AP1-family-specific columns

Author: George Stephenson
Date: December 2025
"""

import argparse
import logging
import sys
from pathlib import Path
import pandas as pd
import numpy as np
from tqdm import tqdm
import gc

# ============================================================================
# COMPREHENSIVE AP1-FAMILY CONFIGURATION
# ============================================================================

# Primary AP1 complex members (canonical FOS/JUN heterodimers)
AP1_CORE_TFS = [
    # Jun family - forms homodimers and heterodimers with Fos
    'JUN', 'c-Jun',           # Proto-oncogene c-Jun
    'JUNB', 'JunB',           # JunB
    'JUND', 'JunD',           # JunD
    
    # Fos family - requires Jun partner for DNA binding
    'FOS', 'c-Fos',           # Proto-oncogene c-Fos
    'FOSB', 'FosB',           # FosB
    'FOSL1', 'Fra-1', 'Fra1', # Fos-related antigen 1
    'FOSL2', 'Fra-2', 'Fra2', # Fos-related antigen 2
]

# ATF/CREB family - heterodimerize with Jun at AP1 sites
AP1_ATF_FAMILY = [
    'ATF1',                   # ATF1 - CREB family, can bind AP1
    'ATF2', 'CRE-BP1',        # ATF2 - binds AP1 and CRE sites
    'ATF3',                   # ATF3 - stress-induced, AP1 binding
    'ATF4', 'CREB2',          # ATF4 - integrated stress response
    'ATF5',                   # ATF5 - similar to ATF4
    'ATF6',                   # ATF6 - ER stress, AP1-like binding
    'ATF7',                   # ATF7 - ATF2 family member
    'CREB1', 'CREB',          # CREB - binds CRE, overlaps with AP1
    'CREM',                   # CREM - CREB family
]

# BATF family - AP1-like bZIP factors, critical for immune function
AP1_BATF_FAMILY = [
    'BATF',                   # BATF - partners with Jun, immune regulation
    'BATF2',                  # BATF2 - tumor suppressor
    'BATF3',                  # BATF3 - dendritic cell development
]

# Small Maf proteins - partner with AP1 at ARE/MARE elements
AP1_MAF_FAMILY = [
    'MAFK', 'MafK',           # Small Maf K
    'MAFF', 'MafF',           # Small Maf F  
    'MAFG', 'MafG',           # Small Maf G
    'MAF', 'c-Maf',           # Large Maf (can also bind AP1-like)
    'MAFB', 'MafB',           # MafB
]

# NRL/NRF family - bind AP1-like sites (ARE elements)
AP1_NRF_FAMILY = [
    'NFE2', 'NF-E2',          # NF-E2 (p45) - erythroid, binds MARE
    'NFE2L1', 'NRF1', 'Nrf1', # NRF1 - stress response
    'NFE2L2', 'NRF2', 'Nrf2', # NRF2 - antioxidant response, AP1-related
    'NFE2L3', 'NRF3', 'Nrf3', # NRF3
    'BACH1',                  # BACH1 - repressor at ARE/MARE
    'BACH2',                  # BACH2 - B cell development
]

# Combine all AP1-related TFs
AP1_ALL_TFS = (
    AP1_CORE_TFS + 
    AP1_ATF_FAMILY + 
    AP1_BATF_FAMILY + 
    AP1_MAF_FAMILY + 
    AP1_NRF_FAMILY
)

# Create case-insensitive lookup set
AP1_TFS_UPPER = set(tf.upper() for tf in AP1_ALL_TFS)

# Enhancer marks for validation
ENHANCER_HISTONE_MARKS = ['H3K27ac', 'H3K4me1', 'H3K4me3']

# Accessibility assays
ACCESSIBILITY_ASSAYS = ['ATAC', 'DNASE']


def setup_logger(name, log_file=None, level="INFO"):
    """Setup logging."""
    logger = logging.getLogger(name)
    logger.setLevel(getattr(logging, level))
    
    # Clear existing handlers
    logger.handlers = []
    
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    
    # Console handler
    ch = logging.StreamHandler()
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    
    # File handler
    if log_file:
        fh = logging.FileHandler(log_file)
        fh.setFormatter(formatter)
        logger.addHandler(fh)
    
    return logger


def is_ap1_tf(tf_name: str) -> bool:
    """Check if a TF name is in the AP1 family (case-insensitive)."""
    if pd.isna(tf_name):
        return False
    return tf_name.upper() in AP1_TFS_UPPER


def process_batch(parquet_path: Path, logger) -> pd.DataFrame:
    """
    Process a single parquet partition and extract AP1-specific scores.
    
    Returns a DataFrame with one row per variant containing:
    - AP1-family TF scores (raw and quantile)
    - Enhancer mark scores
    - Accessibility scores
    - Best TF and biosample info
    """
    df = pd.read_parquet(parquet_path)
    
    if len(df) == 0:
        return pd.DataFrame()
    
    variants = df['variant_id_str'].unique()
    
    # =========================================================================
    # AP1-FAMILY TF SCORES
    # =========================================================================
    
    # Filter to CHIP_TF output type and AP1-family TFs
    ap1_mask = (
        (df['output_type'] == 'CHIP_TF') & 
        (df['transcription_factor'].apply(is_ap1_tf))
    )
    ap1_data = df[ap1_mask]
    
    if len(ap1_data) > 0:
        # Aggregate scores per variant
        ap1_agg = ap1_data.groupby('variant_id_str').agg({
            'quantile_score': ['max', 'mean', 'std', 'count'],
            'raw_score': ['max', 'mean']
        }).reset_index()
        ap1_agg.columns = [
            'variant_id_str', 
            'ap1_quantile_max', 'ap1_quantile_mean', 'ap1_quantile_std', 'ap1_n_tracks',
            'ap1_raw_max', 'ap1_raw_mean'
        ]
        
        # Get best TF and biosample for each variant
        idx_max = ap1_data.groupby('variant_id_str')['raw_score'].idxmax()
        best_info = ap1_data.loc[idx_max, [
            'variant_id_str', 'transcription_factor', 'biosample_name'
        ]].copy()
        best_info.columns = ['variant_id_str', 'ap1_best_tf', 'ap1_best_biosample']
        
        ap1_scores = ap1_agg.merge(best_info, on='variant_id_str', how='left')
    else:
        ap1_scores = pd.DataFrame({
            'variant_id_str': variants,
            'ap1_quantile_max': np.nan,
            'ap1_quantile_mean': np.nan,
            'ap1_quantile_std': np.nan,
            'ap1_n_tracks': 0,
            'ap1_raw_max': np.nan,
            'ap1_raw_mean': np.nan,
            'ap1_best_tf': None,
            'ap1_best_biosample': None
        })
    
    # =========================================================================
    # ENHANCER HISTONE MARK SCORES
    # =========================================================================
    
    enhancer_mask = (
        (df['output_type'] == 'CHIP_HISTONE') &
        (df['histone_mark'].isin(ENHANCER_HISTONE_MARKS))
    )
    enhancer_data = df[enhancer_mask]
    
    if len(enhancer_data) > 0:
        enhancer_agg = enhancer_data.groupby('variant_id_str').agg({
            'quantile_score': 'max',
            'raw_score': 'max'
        }).reset_index()
        enhancer_agg.columns = ['variant_id_str', 'enhancer_quantile_max', 'enhancer_raw_max']
        
        # H3K27ac specifically
        h3k27ac_data = enhancer_data[enhancer_data['histone_mark'] == 'H3K27ac']
        if len(h3k27ac_data) > 0:
            h3k27ac_agg = h3k27ac_data.groupby('variant_id_str')['quantile_score'].max().reset_index()
            h3k27ac_agg.columns = ['variant_id_str', 'h3k27ac_quantile_max']
            enhancer_agg = enhancer_agg.merge(h3k27ac_agg, on='variant_id_str', how='left')
        else:
            enhancer_agg['h3k27ac_quantile_max'] = np.nan
    else:
        enhancer_agg = pd.DataFrame({
            'variant_id_str': variants,
            'enhancer_quantile_max': np.nan,
            'enhancer_raw_max': np.nan,
            'h3k27ac_quantile_max': np.nan
        })
    
    # =========================================================================
    # ACCESSIBILITY SCORES
    # =========================================================================
    
    access_mask = df['output_type'].isin(ACCESSIBILITY_ASSAYS)
    access_data = df[access_mask]
    
    if len(access_data) > 0:
        access_agg = access_data.groupby('variant_id_str').agg({
            'quantile_score': 'max',
            'raw_score': 'max'
        }).reset_index()
        access_agg.columns = ['variant_id_str', 'accessibility_quantile_max', 'accessibility_raw_max']
    else:
        access_agg = pd.DataFrame({
            'variant_id_str': variants,
            'accessibility_quantile_max': np.nan,
            'accessibility_raw_max': np.nan
        })
    
    # =========================================================================
    # GLOBAL SCORES (for comparison)
    # =========================================================================
    
    global_agg = df.groupby('variant_id_str').agg({
        'quantile_score': ['max', 'mean'],
        'raw_score': ['max', 'mean']
    }).reset_index()
    global_agg.columns = [
        'variant_id_str',
        'global_quantile_max', 'global_quantile_mean',
        'global_raw_max', 'global_raw_mean'
    ]
    
    # =========================================================================
    # MERGE ALL SCORES
    # =========================================================================
    
    result = ap1_scores.merge(enhancer_agg, on='variant_id_str', how='outer')
    result = result.merge(access_agg, on='variant_id_str', how='outer')
    result = result.merge(global_agg, on='variant_id_str', how='outer')
    
    return result


def main():
    parser = argparse.ArgumentParser(
        description='Extract AP1-specific scores from forbidden variant predictions'
    )
    parser.add_argument(
        '--partitions-dir', 
        required=True,
        help='Path to partitions directory containing batch parquet files'
    )
    parser.add_argument(
        '--original-summary',
        required=True,
        help='Path to original predictions_summary.tsv (for metadata)'
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Output path for predictions_summary_ap1.tsv'
    )
    parser.add_argument(
        '--log-file',
        help='Optional log file path'
    )
    args = parser.parse_args()
    
    # Setup
    partitions_dir = Path(args.partitions_dir)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    logger = setup_logger(__name__, args.log_file)
    
    logger.info("=" * 70)
    logger.info("Module 08: Extract AP1-Specific Scores from Forbidden Variants")
    logger.info("=" * 70)
    
    # List AP1 TFs being searched
    logger.info(f"\nAP1-family TFs included ({len(AP1_ALL_TFS)} total):")
    logger.info(f"  Core AP1 (FOS/JUN): {', '.join(AP1_CORE_TFS[:8])}...")
    logger.info(f"  ATF family: {', '.join([tf for tf in AP1_ATF_FAMILY if not tf.startswith('CRE')][:5])}...")
    logger.info(f"  BATF family: {', '.join(AP1_BATF_FAMILY)}")
    logger.info(f"  MAF family: {', '.join([tf for tf in AP1_MAF_FAMILY if len(tf) <= 5])}")
    logger.info(f"  NRF/ARE family: {', '.join([tf for tf in AP1_NRF_FAMILY if 'NRF' in tf or 'BACH' in tf])}")
    
    # Find partition files
    partition_files = sorted(partitions_dir.glob("batch_*.parquet"))
    logger.info(f"\nFound {len(partition_files)} partition files")
    
    if len(partition_files) == 0:
        logger.error(f"No partition files found in {partitions_dir}")
        sys.exit(1)
    
    # Process each partition
    logger.info("\nProcessing partitions...")
    all_results = []
    
    for pf in tqdm(partition_files, desc="Processing batches"):
        batch_result = process_batch(pf, logger)
        if len(batch_result) > 0:
            all_results.append(batch_result)
        gc.collect()
    
    # Combine all results
    logger.info("\nCombining results...")
    combined = pd.concat(all_results, ignore_index=True)
    logger.info(f"  Total variants with AP1 scores: {len(combined):,}")
    
    # Load original summary for metadata (AN, motif info, etc.)
    logger.info(f"\nMerging with original metadata from {args.original_summary}...")
    original = pd.read_csv(args.original_summary, sep='\t')
    
    # Rename variant_id column if needed
    if 'variant_id' in original.columns and 'variant_id_str' not in original.columns:
        original = original.rename(columns={'variant_id': 'variant_id_str'})
    
    # Merge metadata
    metadata_cols = ['variant_id_str', 'AN', 'motif_start', 'motif_end', 'strand', 'tier', 'pwm_score']
    metadata_cols = [c for c in metadata_cols if c in original.columns]
    
    final = combined.merge(original[metadata_cols], on='variant_id_str', how='left')
    
    # Sort by AP1 raw score (descending)
    final = final.sort_values('ap1_raw_max', ascending=False, na_position='last')
    
    # Save
    logger.info(f"\nSaving to {output_path}...")
    final.to_csv(output_path, sep='\t', index=False)
    
    # Summary statistics
    logger.info("\n" + "=" * 70)
    logger.info("SUMMARY")
    logger.info("=" * 70)
    
    n_with_ap1 = final['ap1_raw_max'].notna().sum()
    logger.info(f"  Total variants: {len(final):,}")
    logger.info(f"  Variants with AP1 TF predictions: {n_with_ap1:,} ({100*n_with_ap1/len(final):.1f}%)")
    
    if n_with_ap1 > 0:
        logger.info(f"\n  AP1 raw score statistics:")
        logger.info(f"    Mean: {final['ap1_raw_max'].mean():.2f}")
        logger.info(f"    Median: {final['ap1_raw_max'].median():.2f}")
        logger.info(f"    Max: {final['ap1_raw_max'].max():.2f}")
        logger.info(f"    Std: {final['ap1_raw_max'].std():.2f}")
        
        logger.info(f"\n  AP1 quantile score distribution:")
        for thresh in [0.9, 0.95, 0.99, 0.999]:
            count = (final['ap1_quantile_max'] > thresh).sum()
            pct = 100 * count / n_with_ap1
            logger.info(f"    > {thresh}: {count:,} ({pct:.1f}%)")
        
        # Top TFs
        tf_counts = final['ap1_best_tf'].value_counts().head(10)
        logger.info(f"\n  Top AP1-family TFs (by frequency as best responder):")
        for tf, count in tf_counts.items():
            logger.info(f"    {tf}: {count:,} ({100*count/n_with_ap1:.1f}%)")
        
        # Top 10 variants
        logger.info(f"\n  Top 10 forbidden variants by AP1 raw score:")
        for i, row in final.head(10).iterrows():
            logger.info(f"    {row['variant_id_str']}: ap1_raw_max={row['ap1_raw_max']:.1f}, "
                       f"tf={row['ap1_best_tf']}, biosample={row['ap1_best_biosample']}")
    
    logger.info(f"\nOutput saved to: {output_path}")
    logger.info("Done!")


if __name__ == '__main__':
    main()
