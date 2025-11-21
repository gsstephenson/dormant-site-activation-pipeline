#!/usr/bin/env python3
"""
Compute Hamming distance between motif sites and consensus
Part of Dormant Site Activation Pipeline - Module 02
"""

import argparse
import sys
from pathlib import Path
import pandas as pd
import numpy as np

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils.io import load_config
from utils.logging_utils import setup_logger, log_step, ProgressLogger


def hamming_distance(seq1: str, seq2: str) -> int:
    """
    Calculate Hamming distance between two sequences.
    
    Parameters
    ----------
    seq1, seq2 : str
        DNA sequences (must be same length)
        
    Returns
    -------
    int
        Number of mismatches
    """
    if len(seq1) != len(seq2):
        raise ValueError(f"Sequences must be same length: {len(seq1)} vs {len(seq2)}")
    
    return sum(b1 != b2 for b1, b2 in zip(seq1, seq2))


def identify_mismatches(seq: str, consensus: str):
    """
    Identify positions where sequence differs from consensus.
    
    Parameters
    ----------
    seq : str
        Motif sequence
    consensus : str
        Consensus sequence
        
    Returns
    -------
    list
        List of (position, ref_base, consensus_base) tuples
    """
    mismatches = []
    for i, (s, c) in enumerate(zip(seq, consensus)):
        if s != c:
            mismatches.append((i, s, c))
    
    return mismatches


def compute_distances(
    motif_file: str,
    consensus: str,
    output_file: str,
    max_distance: int = None,
    logger=None
):
    """
    Compute Hamming distance for all motif sites.
    
    Parameters
    ----------
    motif_file : str
        Path to motif hits file
    consensus : str
        Consensus sequence
    output_file : str
        Output file path
    max_distance : int
        Maximum distance to keep (filter sites with distance > max_distance)
    logger : logging.Logger
        Logger instance
        
    Returns
    -------
    pd.DataFrame
        Motif sites with distance information
    """
    logger.info("Loading motif sites...")
    df = pd.read_csv(motif_file, sep='\t', low_memory=False)
    
    logger.info(f"  Total sites: {len(df):,}")
    logger.info(f"  Consensus length: {len(consensus)} bp")
    
    # Calculate Hamming distance
    logger.info("Computing Hamming distances...")
    
    distances = []
    mismatches_list = []
    
    progress = ProgressLogger(total=len(df), logger=logger)
    
    for idx, row in df.iterrows():
        seq = row['matched_sequence'].upper()
        
        # Handle sequences that don't match consensus length
        if len(seq) != len(consensus):
            logger.warning(f"Site {idx}: sequence length {len(seq)} != consensus {len(consensus)}")
            distances.append(np.nan)
            mismatches_list.append([])
            continue
        
        dist = hamming_distance(seq, consensus)
        mm = identify_mismatches(seq, consensus)
        
        distances.append(dist)
        mismatches_list.append(mm)
        
        progress.update()
    
    # Progress complete - final log is already done by update()
    
    # Add to dataframe
    df['hamming_distance'] = distances
    df['num_mismatches'] = [len(mm) for mm in mismatches_list]
    df['mismatch_positions'] = [
        ','.join([str(pos) for pos, _, _ in mm]) if mm else ''
        for mm in mismatches_list
    ]
    df['ref_bases'] = [
        ','.join([ref for _, ref, _ in mm]) if mm else ''
        for mm in mismatches_list
    ]
    df['consensus_bases'] = [
        ','.join([cons for _, _, cons in mm]) if mm else ''
        for mm in mismatches_list
    ]
    
    # Log distance distribution
    logger.info("\nHamming distance distribution:")
    dist_counts = df['hamming_distance'].value_counts().sort_index()
    for dist, count in dist_counts.items():
        if pd.isna(dist):
            continue
        percent = 100 * count / len(df)
        logger.info(f"  Distance {int(dist)}: {count:,} sites ({percent:.1f}%)")
    
    # Filter by max distance if specified
    if max_distance is not None:
        before = len(df)
        df = df[df['hamming_distance'] <= max_distance]
        after = len(df)
        logger.info(f"\nFiltered to distance ≤ {max_distance}: {after:,} / {before:,} sites")
    
    # Save results
    df.to_csv(output_file, sep='\t', index=False)
    logger.info(f"\n✓ Saved distances: {output_file}")
    
    return df


def main():
    parser = argparse.ArgumentParser(
        description="Compute Hamming distance to consensus for all motif sites"
    )
    parser.add_argument(
        '--config',
        type=str,
        required=True,
        help='Path to pipeline configuration file'
    )
    parser.add_argument(
        '--motif-file',
        type=str,
        help='Path to motif hits file (default: from Module 01)'
    )
    parser.add_argument(
        '--consensus-file',
        type=str,
        help='Path to consensus file (default: from consensus_from_pwm.py)'
    )
    parser.add_argument(
        '--output-file',
        type=str,
        help='Output file (default: results/mutation_paths/{TF}/distances.tsv)'
    )
    parser.add_argument(
        '--max-distance',
        type=int,
        help='Maximum Hamming distance to keep (default: from config)'
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
    max_distance_default = config['mutation_paths'].get('max_hamming_distance', 3)
    
    # Determine input files
    if args.motif_file:
        motif_file = args.motif_file
    else:
        motif_file = f"results/motif_scan/{tf_name}/motif_hits_tiered.tsv"
    
    if args.consensus_file:
        consensus_file = args.consensus_file
    else:
        consensus_file = f"data/motifs/consensus/{tf_name}_consensus.txt"
    
    # Determine output
    if args.output_file:
        output_file = args.output_file
    else:
        output_dir = Path(f"results/mutation_paths/{tf_name}")
        output_dir.mkdir(parents=True, exist_ok=True)
        output_file = output_dir / "distances.tsv"
    
    max_distance = args.max_distance or max_distance_default
    
    # Setup logger
    log_file = args.log_file or f"logs/02_compute_distance_{tf_name}.log"
    logger = setup_logger(
        name=__name__,
        log_file=log_file,
        level='INFO'
    )
    
    log_step(logger, "Module 02: Hamming Distance Computation", start=True)
    
    # Check input files exist
    if not Path(motif_file).exists():
        logger.error(f"Motif file not found: {motif_file}")
        logger.error("Run Module 01 first")
        sys.exit(1)
    
    if not Path(consensus_file).exists():
        logger.error(f"Consensus file not found: {consensus_file}")
        logger.error("Run consensus_from_pwm.py first")
        sys.exit(1)
    
    # Load consensus
    logger.info(f"Loading consensus from: {consensus_file}")
    with open(consensus_file, 'r') as f:
        lines = f.readlines()
        consensus = [line.strip() for line in lines if not line.startswith('#')][-1]
    
    logger.info(f"  Consensus: {consensus}")
    
    # Compute distances
    try:
        df = compute_distances(
            motif_file=motif_file,
            consensus=consensus,
            output_file=str(output_file),
            max_distance=max_distance,
            logger=logger
        )
    except Exception as e:
        logger.error(f"Distance computation failed: {e}")
        import traceback
        logger.error(traceback.format_exc())
        sys.exit(1)
    
    log_step(logger, "Module 02: Hamming Distance Computation", start=False)


if __name__ == '__main__':
    main()
