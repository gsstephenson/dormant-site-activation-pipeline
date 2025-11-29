#!/usr/bin/env python3
"""
Module 05: Compute Activation Landscape for AP1 Dormant Sites

This script creates a 2D "activation landscape" where:
- X-axis: Population accessibility (how easy to access via human variation)
- Y-axis: Functional impact (predicted change in AP1 binding and enhancer activity)

BIOLOGICAL RATIONALE FOR Y-AXIS DESIGN:
========================================
We're studying AP1 (FOS/JUN heterodimer) dormant site activation. Using max() 
across ALL AlphaGenome tracks would conflate unrelated signals (e.g., CTCF, 
splice sites). Instead, we use a biologically-focused approach:

PRIMARY METRIC - AP1-family TF binding:
  - JUND, JUN, JUNB (Jun family)
  - FOS, FOSL1, FOSL2 (Fos family)
  - ATF3, BATF (AP1-related bZIP factors)
  
  Rationale: These tracks directly measure predicted changes in AP1-family 
  binding, which is exactly what we hypothesize will change when a dormant 
  site becomes activated.

SECONDARY METRIC - Enhancer activation marks:
  - H3K27ac: Active enhancer mark
  - H3K4me1: Enhancer priming mark
  
  Rationale: AP1 sites function as enhancers. If a variant increases AP1 
  binding AND increases enhancer marks, this validates functional activation.

TERTIARY METRIC - Chromatin accessibility:
  - ATAC-seq, DNase-seq
  
  Rationale: AP1 is a pioneer factor that can open chromatin. Increased 
  accessibility supports the activation hypothesis.

Author: George Stephenson
Date: November 2025
"""

import pandas as pd
import numpy as np
from pathlib import Path
import argparse
import logging
from typing import Dict, List, Tuple, Optional

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# ============================================================================
# BIOLOGICAL CONFIGURATION
# ============================================================================

# AP1-family transcription factors (primary Y-axis)
AP1_FAMILY_TFS = [
    'JUND', 'JUN', 'JUNB',      # Jun family
    'FOS', 'FOSL1', 'FOSL2',    # Fos family  
    'ATF3', 'ATF2', 'ATF7',     # ATF family (heterodimerize with Jun)
    'BATF', 'BATF2',            # BATF family (AP1-like)
    'MAFK', 'MAFF', 'MAFG',     # Small Maf (partner with AP1)
]

# Enhancer-associated histone marks (secondary validation)
ENHANCER_HISTONE_MARKS = [
    'H3K27ac',   # Active enhancer mark
    'H3K4me1',   # Enhancer priming
    'H3K4me3',   # Promoter/active mark
]

# Accessibility tracks (tertiary)
ACCESSIBILITY_ASSAYS = ['ATAC', 'DNASE']


def load_alphagenome_predictions(parquet_path: Path) -> pd.DataFrame:
    """Load AlphaGenome predictions parquet file."""
    logger.info(f"Loading predictions from {parquet_path}")
    logger.info(f"  (This may take 1-2 minutes for large files...)")
    import sys
    sys.stdout.flush()
    df = pd.read_parquet(parquet_path)
    logger.info(f"Loaded {len(df):,} predictions for {df['variant_id_str'].nunique():,} variants")
    sys.stdout.flush()
    return df


def load_gnomad_metadata(gnomad_path: Path) -> pd.DataFrame:
    """Load gnomAD intersection data with path/step information."""
    logger.info(f"Loading gnomAD data from {gnomad_path}")
    df = pd.read_csv(gnomad_path, sep='\t', low_memory=False)
    
    # Handle VCF "." notation for missing values
    # bcftools outputs "." for missing INFO fields, which causes pandas to read as object type
    numeric_cols = ['AF', 'AC', 'AN', 'nhomalt']
    for col in numeric_cols:
        if col in df.columns:
            df[col] = df[col].replace('.', pd.NA)
            df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0.0)
    
    logger.info(f"Loaded {len(df):,} variant-path combinations")
    return df


def compute_ap1_impact_score(predictions: pd.DataFrame) -> pd.DataFrame:
    """
    Compute AP1-specific functional impact scores for each variant.
    
    Uses VECTORIZED groupby operations for ~100x speedup over loop-based approach.
    
    Returns DataFrame with columns:
    - variant_id_str: variant identifier
    - ap1_max_score: max quantile score across AP1-family TFs
    - ap1_mean_score: mean quantile score across AP1-family TFs
    - ap1_best_tf: which AP1-family TF showed highest signal
    - ap1_best_biosample: which cell type showed highest signal
    - enhancer_max_score: max score for H3K27ac/H3K4me1
    - accessibility_max_score: max score for ATAC/DNase
    """
    import sys
    from datetime import datetime
    
    start_time = datetime.now()
    n_variants = predictions['variant_id_str'].nunique()
    logger.info(f"Computing AP1 impact scores for {n_variants:,} variants using vectorized operations...")
    sys.stdout.flush()
    
    # =========================================================================
    # STEP 1: Get metadata for each variant (gnomAD AF, path_id, etc.)
    # =========================================================================
    logger.info("  [1/5] Extracting variant metadata...")
    sys.stdout.flush()
    
    variant_metadata = predictions.groupby('variant_id_str').agg({
        'gnomad_AF': 'first',
        'path_id': 'first',
    }).reset_index()
    
    # =========================================================================
    # STEP 2: Compute AP1-family TF scores (PRIMARY Y-axis)
    # =========================================================================
    logger.info("  [2/5] Computing AP1-family TF binding scores...")
    sys.stdout.flush()
    
    # Filter to AP1-family TFs only
    ap1_mask = (
        (predictions['output_type'] == 'CHIP_TF') & 
        (predictions['transcription_factor'].isin(AP1_FAMILY_TFS))
    )
    ap1_data = predictions[ap1_mask].copy()
    logger.info(f"        Found {len(ap1_data):,} AP1-family TF predictions")
    
    if len(ap1_data) > 0:
        # Get max and mean scores per variant
        ap1_agg = ap1_data.groupby('variant_id_str').agg({
            'quantile_score': ['max', 'mean', 'count']
        }).reset_index()
        ap1_agg.columns = ['variant_id_str', 'ap1_max_score', 'ap1_mean_score', 'ap1_n_tracks']
        
        # Get the best TF and biosample for each variant (where score is max)
        idx_max = ap1_data.groupby('variant_id_str')['quantile_score'].idxmax()
        best_tf_df = ap1_data.loc[idx_max, ['variant_id_str', 'transcription_factor', 'biosample_name']].copy()
        best_tf_df.columns = ['variant_id_str', 'ap1_best_tf', 'ap1_best_biosample']
        
        # Merge
        ap1_scores = ap1_agg.merge(best_tf_df, on='variant_id_str', how='left')
    else:
        ap1_scores = pd.DataFrame({
            'variant_id_str': predictions['variant_id_str'].unique(),
            'ap1_max_score': np.nan,
            'ap1_mean_score': np.nan,
            'ap1_n_tracks': 0,
            'ap1_best_tf': None,
            'ap1_best_biosample': None
        })
    
    # =========================================================================
    # STEP 3: Compute enhancer histone mark scores (SECONDARY validation)
    # =========================================================================
    logger.info("  [3/5] Computing enhancer mark scores (H3K27ac, H3K4me1)...")
    sys.stdout.flush()
    
    enhancer_mask = (
        (predictions['output_type'] == 'CHIP_HISTONE') &
        (predictions['histone_mark'].isin(ENHANCER_HISTONE_MARKS))
    )
    enhancer_data = predictions[enhancer_mask].copy()
    logger.info(f"        Found {len(enhancer_data):,} enhancer mark predictions")
    
    if len(enhancer_data) > 0:
        enhancer_agg = enhancer_data.groupby('variant_id_str')['quantile_score'].max().reset_index()
        enhancer_agg.columns = ['variant_id_str', 'enhancer_max_score']
        
        # Also get H3K27ac specifically
        h3k27ac_data = enhancer_data[enhancer_data['histone_mark'] == 'H3K27ac']
        if len(h3k27ac_data) > 0:
            h3k27ac_agg = h3k27ac_data.groupby('variant_id_str')['quantile_score'].max().reset_index()
            h3k27ac_agg.columns = ['variant_id_str', 'h3k27ac_max_score']
            enhancer_agg = enhancer_agg.merge(h3k27ac_agg, on='variant_id_str', how='left')
        else:
            enhancer_agg['h3k27ac_max_score'] = np.nan
    else:
        enhancer_agg = pd.DataFrame({
            'variant_id_str': predictions['variant_id_str'].unique(),
            'enhancer_max_score': np.nan,
            'h3k27ac_max_score': np.nan
        })
    
    # =========================================================================
    # STEP 4: Compute accessibility scores (TERTIARY)
    # =========================================================================
    logger.info("  [4/5] Computing chromatin accessibility scores (ATAC, DNase)...")
    sys.stdout.flush()
    
    access_mask = predictions['output_type'].isin(ACCESSIBILITY_ASSAYS)
    access_data = predictions[access_mask].copy()
    logger.info(f"        Found {len(access_data):,} accessibility predictions")
    
    if len(access_data) > 0:
        access_agg = access_data.groupby('variant_id_str')['quantile_score'].max().reset_index()
        access_agg.columns = ['variant_id_str', 'accessibility_max_score']
    else:
        access_agg = pd.DataFrame({
            'variant_id_str': predictions['variant_id_str'].unique(),
            'accessibility_max_score': np.nan
        })
    
    # =========================================================================
    # STEP 5: Compute global max (for comparison only)
    # =========================================================================
    logger.info("  [5/5] Computing global max scores (for comparison)...")
    sys.stdout.flush()
    
    global_agg = predictions.groupby('variant_id_str')['quantile_score'].max().reset_index()
    global_agg.columns = ['variant_id_str', 'global_max_score']
    
    # =========================================================================
    # MERGE ALL RESULTS
    # =========================================================================
    logger.info("  Merging all scores...")
    sys.stdout.flush()
    
    result = variant_metadata.merge(ap1_scores, on='variant_id_str', how='left')
    result = result.merge(enhancer_agg, on='variant_id_str', how='left')
    result = result.merge(access_agg, on='variant_id_str', how='left')
    result = result.merge(global_agg, on='variant_id_str', how='left')
    
    elapsed = (datetime.now() - start_time).total_seconds()
    logger.info(f"  DONE! Computed scores for {len(result):,} variants in {elapsed:.1f} seconds")
    sys.stdout.flush()
    
    return result


def compute_x_axis_accessibility(
    impact_scores: pd.DataFrame,
    gnomad_data: pd.DataFrame
) -> pd.DataFrame:
    """
    Compute X-axis: Population accessibility score.
    
    X = -log10(AF + 1e-12) × Hamming_distance
    
    Interpretation:
    - Low X = highly accessible (common variant, few mutations)
    - High X = hard to access (rare variant, many mutations)
    
    We need to join with gnomAD data to get Hamming distance (path length).
    
    IMPORTANT: 
    - step_num = position of this variant within a path (not what we want)
    - total_steps = total mutations in the path = Hamming distance (what we want)
    """
    import sys
    logger.info("Computing X-axis (population accessibility)...")
    sys.stdout.flush()
    
    gnomad_data = gnomad_data.copy()
    
    # Create variant_id matching format used in AlphaGenome predictions
    # Use genomic_position, ref_base, alt_base (our path coordinates)
    # NOT pos, ref, alt (which come from VCF matching and have NaN for non-matches)
    # AlphaGenome format: chr1:867861:A>C
    gnomad_data['variant_id_str'] = gnomad_data.apply(
        lambda row: f"{row['chr']}:{row['genomic_position']}:{row['ref_base']}>{row['alt_base']}", axis=1
    )
    
    logger.info(f"  gnomAD variants: {gnomad_data['variant_id_str'].nunique():,}")
    logger.info(f"  AlphaGenome variants: {impact_scores['variant_id_str'].nunique():,}")
    sys.stdout.flush()
    
    # For each variant, get:
    # - Minimum total_steps (shortest path = easiest activation route)
    # - Maximum AF across all observations
    variant_summary = gnomad_data.groupby('variant_id_str').agg({
        'total_steps': 'min',    # Shortest path to activation (Hamming distance)
        'AF': 'max',             # Highest observed AF
    }).reset_index()
    
    variant_summary.columns = ['variant_id_str', 'hamming_distance', 'gnomad_AF_from_file']
    
    logger.info(f"  Hamming distance range: {variant_summary['hamming_distance'].min()} - {variant_summary['hamming_distance'].max()}")
    logger.info(f"  Hamming distance distribution:")
    logger.info(f"    {variant_summary['hamming_distance'].value_counts().sort_index().to_dict()}")
    sys.stdout.flush()
    
    # Merge with impact scores
    merged = impact_scores.merge(variant_summary, on='variant_id_str', how='left')
    
    # Use gnomad_AF from predictions if available, otherwise from gnomad file
    merged['AF_final'] = merged['gnomad_AF'].fillna(merged['gnomad_AF_from_file'])
    
    # Fill missing hamming distance with median (shouldn't happen if data is complete)
    median_hamming = merged['hamming_distance'].median()
    merged['hamming_distance'] = merged['hamming_distance'].fillna(median_hamming)
    
    # Compute X-axis: -log10(AF) × Hamming_distance
    # Handle AF=0 by using 1e-12 floor
    merged['AF_safe'] = merged['AF_final'].clip(lower=1e-12)
    merged['neg_log10_AF'] = -np.log10(merged['AF_safe'])
    
    # X-axis: accessibility score (higher = harder to access)
    merged['x_accessibility'] = merged['neg_log10_AF'] * merged['hamming_distance']
    
    # Also keep just neg_log10_AF for simpler plotting
    merged['x_neg_log10_AF'] = merged['neg_log10_AF']
    
    logger.info(f"Computed accessibility for {len(merged):,} variants")
    logger.info(f"  X-axis (accessibility) range: {merged['x_accessibility'].min():.2f} to {merged['x_accessibility'].max():.2f}")
    logger.info(f"  -log10(AF) range: {merged['neg_log10_AF'].min():.2f} to {merged['neg_log10_AF'].max():.2f}")
    sys.stdout.flush()
    
    return merged


def create_activation_landscape(merged_data: pd.DataFrame) -> pd.DataFrame:
    """
    Create final activation landscape DataFrame with all coordinates and annotations.
    """
    landscape = merged_data.copy()
    
    # Define Y-axis options
    landscape['y_ap1_impact'] = landscape['ap1_max_score']
    landscape['y_enhancer_impact'] = landscape['enhancer_max_score']
    landscape['y_global_impact'] = landscape['global_max_score']
    
    # Categorize variants by quadrant
    # Use median as threshold (data-driven)
    x_median = landscape['x_accessibility'].median()
    y_median = landscape['y_ap1_impact'].median()
    
    def classify_quadrant(row):
        x = row['x_accessibility']
        y = row['y_ap1_impact']
        if pd.isna(x) or pd.isna(y):
            return 'Insufficient data'
        elif x < x_median and y > y_median:
            return 'HIGH PRIORITY: Accessible + High Impact'
        elif x >= x_median and y > y_median:
            return 'High Impact, Hard to Access'
        elif x < x_median and y <= y_median:
            return 'Accessible, Low Impact'
        else:
            return 'Low Priority'
    
    landscape['quadrant'] = landscape.apply(classify_quadrant, axis=1)
    
    # Add interpretable columns
    landscape['is_common'] = landscape['AF_final'] >= 0.01  # ≥1% MAF
    landscape['is_single_step'] = landscape['hamming_distance'] == 1
    landscape['shows_ap1_gain'] = landscape['ap1_max_score'] > 0.9  # Top 10% quantile
    landscape['shows_enhancer_gain'] = landscape['enhancer_max_score'] > 0.9
    
    return landscape


def save_results(
    landscape: pd.DataFrame,
    output_dir: Path,
    motif_name: str = 'AP1'
):
    """Save activation landscape results."""
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Full landscape data
    out_path = output_dir / f'{motif_name}_activation_landscape.tsv'
    landscape.to_csv(out_path, sep='\t', index=False)
    logger.info(f"Saved landscape data to {out_path}")
    
    # High-priority candidates (accessible + high impact)
    high_priority = landscape[
        (landscape['quadrant'] == 'HIGH PRIORITY: Accessible + High Impact')
    ].sort_values('y_ap1_impact', ascending=False)
    
    hp_path = output_dir / f'{motif_name}_high_priority_candidates.tsv'
    high_priority.to_csv(hp_path, sep='\t', index=False)
    logger.info(f"Saved {len(high_priority):,} high-priority candidates to {hp_path}")
    
    # Summary statistics
    summary = {
        'total_variants': len(landscape),
        'variants_with_ap1_data': landscape['ap1_max_score'].notna().sum(),
        'variants_with_enhancer_data': landscape['enhancer_max_score'].notna().sum(),
        'high_priority_count': len(high_priority),
        'common_variants': landscape['is_common'].sum(),
        'single_step_variants': landscape['is_single_step'].sum(),
        'strong_ap1_gain': landscape['shows_ap1_gain'].sum(),
        'strong_enhancer_gain': landscape['shows_enhancer_gain'].sum(),
        'x_axis_median': landscape['x_accessibility'].median(),
        'y_axis_median': landscape['y_ap1_impact'].median(),
        'x_axis_range': f"{landscape['x_accessibility'].min():.2f} - {landscape['x_accessibility'].max():.2f}",
        'y_axis_range': f"{landscape['y_ap1_impact'].min():.3f} - {landscape['y_ap1_impact'].max():.3f}",
    }
    
    # Quadrant breakdown
    quadrant_counts = landscape['quadrant'].value_counts().to_dict()
    summary.update({f'quadrant_{k}': v for k, v in quadrant_counts.items()})
    
    summary_path = output_dir / f'{motif_name}_landscape_summary.txt'
    with open(summary_path, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("AP1 DORMANT SITE ACTIVATION LANDSCAPE SUMMARY\n")
        f.write("=" * 60 + "\n\n")
        
        f.write("BIOLOGICAL RATIONALE:\n")
        f.write("-" * 40 + "\n")
        f.write("Y-axis uses AP1-family TF ChIP-seq predictions (JUND, JUN,\n")
        f.write("FOS, FOSL1, FOSL2, ATF3, BATF, etc.) rather than global max\n")
        f.write("across all tracks. This provides direct biological relevance\n")
        f.write("to the dormant AP1 site activation hypothesis.\n\n")
        
        f.write("SUMMARY STATISTICS:\n")
        f.write("-" * 40 + "\n")
        for key, value in summary.items():
            f.write(f"{key}: {value}\n")
    
    logger.info(f"Saved summary to {summary_path}")
    
    return summary


def main():
    parser = argparse.ArgumentParser(
        description='Compute AP1 Dormant Site Activation Landscape'
    )
    parser.add_argument(
        '--predictions',
        type=Path,
        default=Path('results/alphagenome/AP1/predictions.parquet'),
        help='Path to AlphaGenome predictions parquet file'
    )
    parser.add_argument(
        '--gnomad',
        type=Path,
        default=Path('results/gnomad_intersection/AP1/all_observed_variants.tsv'),
        help='Path to gnomAD intersection data'
    )
    parser.add_argument(
        '--output',
        type=Path,
        default=Path('results/landscape/AP1'),
        help='Output directory'
    )
    parser.add_argument(
        '--motif',
        type=str,
        default='AP1',
        help='Motif name for output files'
    )
    
    args = parser.parse_args()
    
    # Load data
    predictions = load_alphagenome_predictions(args.predictions)
    gnomad_data = load_gnomad_metadata(args.gnomad)
    
    # Compute AP1-specific impact scores (Y-axis)
    impact_scores = compute_ap1_impact_score(predictions)
    
    # Compute population accessibility (X-axis)
    merged = compute_x_axis_accessibility(impact_scores, gnomad_data)
    
    # Create final landscape
    landscape = create_activation_landscape(merged)
    
    # Save results
    summary = save_results(landscape, args.output, args.motif)
    
    logger.info("=" * 60)
    logger.info("ACTIVATION LANDSCAPE COMPLETE")
    logger.info("=" * 60)
    logger.info(f"Total variants: {summary['total_variants']:,}")
    logger.info(f"High-priority candidates: {summary['high_priority_count']:,}")
    logger.info(f"Strong AP1 gain (>90th percentile): {summary['strong_ap1_gain']:,}")


if __name__ == '__main__':
    main()
