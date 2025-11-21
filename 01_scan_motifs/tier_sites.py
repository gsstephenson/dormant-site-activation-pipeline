#!/usr/bin/env python3
"""
Tier motif sites by PWM score
Part of Dormant Site Activation Pipeline - Module 01
"""

import argparse
import sys
from pathlib import Path
import pandas as pd
import numpy as np

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils.io import load_config
from utils.logging_utils import setup_logger, log_step, log_parameters
from utils.constants import MOTIF_TIERS


def calculate_score_percentile(scores):
    """
    Calculate percentile rank for each score efficiently.
    
    For each score, computes the fraction of all scores that are <= this score.
    This is the proper percentile definition for tier assignment.
    
    Parameters
    ----------
    scores : array-like
        PWM scores
        
    Returns
    -------
    array
        Percentile ranks (0-1), where 1.0 = highest score
    """
    # Sort scores and get sorting indices
    sorted_indices = np.argsort(scores)
    
    # Create array to hold percentiles
    percentiles = np.zeros(len(scores))
    
    # Assign percentile based on position in sorted array
    # Handle ties by giving them the same percentile (average of their ranks)
    unique_scores, inverse, counts = np.unique(scores, return_inverse=True, return_counts=True)
    
    # For each unique score, calculate what fraction of scores are <= it
    cumsum = np.cumsum(counts)
    percentile_map = cumsum / len(scores)
    
    # Map back to original array
    percentiles = percentile_map[inverse]
    
    return percentiles


def tier_motif_sites(
    motif_file: str,
    tier_cutoffs: dict,
    output_dir: str,
    logger=None
):
    """
    Classify motif sites into tiers based on PWM score.
    
    Parameters
    ----------
    motif_file : str
        Path to motif hits TSV file
    tier_cutoffs : dict
        Dictionary of tier cutoffs (e.g., {'tier0': 0.95, 'tier1': 0.85})
    output_dir : str
        Output directory
    logger : logging.Logger
        Logger instance
        
    Returns
    -------
    pd.DataFrame
        Motif sites with tier assignments
    """
    logger.info("Loading motif hits...")
    df = pd.read_csv(motif_file, sep='\t', low_memory=False)
    
    logger.info(f"  Total sites: {len(df):,}")
    
    # Calculate score statistics
    scores = df['score'].values
    score_min = scores.min()
    score_max = scores.max()
    score_mean = scores.mean()
    score_std = scores.std()
    
    logger.info(f"  Score range: [{score_min:.3f}, {score_max:.3f}]")
    logger.info(f"  Score mean ± std: {score_mean:.3f} ± {score_std:.3f}")
    
    # Calculate percentile for each score
    logger.info("Calculating score percentiles...")
    percentiles = calculate_score_percentile(scores)
    df['score_percentile'] = percentiles
    
    # Assign tiers
    logger.info("Assigning tiers...")
    
    def assign_tier(percentile):
        if percentile >= tier_cutoffs.get('tier0', 0.95):
            return 0
        elif percentile >= tier_cutoffs.get('tier1', 0.85):
            return 1
        elif percentile >= tier_cutoffs.get('tier2', 0.70):
            return 2
        else:
            return 3
    
    df['tier'] = df['score_percentile'].apply(assign_tier)
    
    # Log tier statistics
    logger.info("\nTier distribution:")
    for tier in sorted(df['tier'].unique()):
        count = (df['tier'] == tier).sum()
        percent = 100 * count / len(df)
        tier_name = MOTIF_TIERS[tier]
        logger.info(f"  Tier {tier} ({tier_name}): {count:,} sites ({percent:.1f}%)")
    
    # Save tiered results
    output_file = Path(output_dir) / "motif_hits_tiered.tsv"
    df.to_csv(output_file, sep='\t', index=False)
    logger.info(f"\n✓ Saved tiered results: {output_file}")
    
    # Save separate BED files for each tier
    logger.info("\nSaving tier-specific BED files...")
    for tier in sorted(df['tier'].unique()):
        tier_df = df[df['tier'] == tier]
        tier_name = MOTIF_TIERS[tier]
        
        bed_file = Path(output_dir) / f"motif_hits_tier{tier}_{tier_name}.bed"
        
        bed_data = pd.DataFrame({
            'chr': tier_df['chr'],
            'start': tier_df['start'],
            'end': tier_df['stop'],
            'name': tier_df['motif_id'] + '_tier' + tier_df['tier'].astype(str) + '_' + tier_df.index.astype(str),
            'score': tier_df['score'],
            'strand': tier_df['strand'],
            'sequence': tier_df['matched_sequence'],
            'tier': tier_df['tier']
        })
        
        bed_data.to_csv(bed_file, sep='\t', header=False, index=False)
        logger.info(f"  ✓ Tier {tier}: {bed_file}")
    
    # Create summary statistics
    summary = {
        'total_sites': len(df),
        'score_min': float(score_min),
        'score_max': float(score_max),
        'score_mean': float(score_mean),
        'score_std': float(score_std),
        'tier_counts': {
            int(tier): int((df['tier'] == tier).sum())
            for tier in sorted(df['tier'].unique())
        },
        'tier_cutoffs': tier_cutoffs
    }
    
    import json
    summary_file = Path(output_dir) / "tier_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    logger.info(f"\n✓ Saved summary: {summary_file}")
    
    return df


def generate_tier_report(df, output_file: str, logger=None):
    """
    Generate a human-readable tier report.
    
    Parameters
    ----------
    df : pd.DataFrame
        Tiered motif sites
    output_file : str
        Output report file
    logger : logging.Logger
        Logger instance
    """
    with open(output_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("Motif Site Tiering Report\n")
        f.write("=" * 80 + "\n\n")
        
        f.write(f"Total motif sites: {len(df):,}\n\n")
        
        f.write("Tier Distribution:\n")
        f.write("-" * 80 + "\n")
        
        for tier in sorted(df['tier'].unique()):
            tier_df = df[df['tier'] == tier]
            count = len(tier_df)
            percent = 100 * count / len(df)
            tier_name = MOTIF_TIERS[tier]
            
            score_min = tier_df['score'].min()
            score_max = tier_df['score'].max()
            score_mean = tier_df['score'].mean()
            
            f.write(f"\nTier {tier} - {tier_name.upper()}\n")
            f.write(f"  Count: {count:,} sites ({percent:.1f}%)\n")
            f.write(f"  Score range: [{score_min:.3f}, {score_max:.3f}]\n")
            f.write(f"  Score mean: {score_mean:.3f}\n")
            
            # Show a few example sequences
            examples = tier_df.head(3)['matched_sequence'].tolist()
            f.write(f"  Example sequences:\n")
            for seq in examples:
                f.write(f"    {seq}\n")
        
        f.write("\n" + "=" * 80 + "\n")
    
    logger.info(f"✓ Generated report: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Tier motif sites by PWM score strength"
    )
    parser.add_argument(
        '--config',
        type=str,
        required=True,
        help='Path to pipeline configuration file'
    )
    parser.add_argument(
        '--input-file',
        type=str,
        help='Path to motif hits TSV (default: from Module 01 output)'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        help='Output directory (default: same as input)'
    )
    parser.add_argument(
        '--log-file',
        type=str,
        help='Log file path'
    )
    
    args = parser.parse_args()
    
    # Load config
    config = load_config(args.config)
    
    # Extract parameters
    tf_name = config['tf_name']
    tier_cutoffs = config['motif_scan']['tier_cutoffs']
    
    # Determine input file
    if args.input_file:
        input_file = args.input_file
        output_dir = args.output_dir or str(Path(input_file).parent)
    else:
        input_file = f"results/motif_scan/{tf_name}/motif_hits_full.tsv"
        output_dir = args.output_dir or f"results/motif_scan/{tf_name}"
    
    # Setup logger
    log_file = args.log_file or f"logs/01_tier_sites_{tf_name}.log"
    logger = setup_logger(
        name=__name__,
        log_file=log_file,
        level='INFO'
    )
    
    log_step(logger, "Module 01: Motif Site Tiering", start=True)
    
    # Log parameters
    params = {
        'TF': tf_name,
        'Input file': input_file,
        'Output directory': output_dir,
        'Tier cutoffs': tier_cutoffs
    }
    log_parameters(logger, params)
    
    # Check input file exists
    if not Path(input_file).exists():
        logger.error(f"Input file not found: {input_file}")
        logger.error("Run motif scanning first:")
        logger.error("  python 01_scan_motifs/scan_genome_fimo.py --config pipeline_config.yaml")
        sys.exit(1)
    
    # Tier motif sites
    try:
        df = tier_motif_sites(
            motif_file=input_file,
            tier_cutoffs=tier_cutoffs,
            output_dir=output_dir,
            logger=logger
        )
    except Exception as e:
        logger.error(f"Tiering failed: {e}")
        import traceback
        logger.error(traceback.format_exc())
        sys.exit(1)
    
    # Generate report
    report_file = Path(output_dir) / "tier_report.txt"
    generate_tier_report(df, str(report_file), logger)
    
    log_step(logger, "Module 01: Motif Site Tiering", start=False)
    logger.info("\nNext step: Module 02 - Generate mutation paths")
    logger.info("  python 02_generate_mutation_paths/enumerate_paths.py --config pipeline_config.yaml")


if __name__ == '__main__':
    main()
