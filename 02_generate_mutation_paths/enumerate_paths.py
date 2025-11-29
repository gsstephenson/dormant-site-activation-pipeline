#!/usr/bin/env python3
"""
Enumerate minimal mutation paths from motif sites to consensus
Part of Dormant Site Activation Pipeline - Module 02
"""

import argparse
import sys
from pathlib import Path
import pandas as pd
import numpy as np
from itertools import combinations

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils.io import load_config
from utils.logging_utils import setup_logger, log_step, log_parameters, ProgressLogger


def reverse_complement(base: str) -> str:
    """
    Return the reverse complement of a single nucleotide base.
    
    For minus-strand motifs, the ref/alt bases are in motif orientation
    (5'->3' on the minus strand). To match VCF format (always plus strand),
    we need to reverse complement.
    
    Parameters
    ----------
    base : str
        Single nucleotide (A, T, G, or C)
        
    Returns
    -------
    str
        Reverse complement base
    """
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
    return complement.get(base, base)


def generate_mutation_paths(
    site_seq: str,
    consensus: str,
    mismatch_positions: list,
    max_steps: int = 3
):
    """
    Generate all minimal mutation paths from site to consensus.
    
    Parameters
    ----------
    site_seq : str
        Current motif sequence
    consensus : str
        Consensus sequence
    mismatch_positions : list
        List of (position, ref_base, consensus_base) tuples
    max_steps : int
        Maximum number of mutation steps
        
    Returns
    -------
    list
        List of mutation paths, each path is a list of steps
    """
    if not mismatch_positions or len(mismatch_positions) > max_steps:
        return []
    
    # For sites within max_steps, enumerate all paths
    # (order matters for intermediate alleles)
    num_mismatches = len(mismatch_positions)
    
    paths = []
    
    # Generate all permutations of mismatch order
    from itertools import permutations
    
    for perm in permutations(mismatch_positions):
        path = []
        current_seq = site_seq
        
        for pos, ref, alt in perm:
            # Create mutation step
            step = {
                'position': pos,
                'ref': ref,
                'alt': alt,
                'seq_before': current_seq,
                'seq_after': current_seq[:pos] + alt + current_seq[pos+1:]
            }
            path.append(step)
            current_seq = step['seq_after']
        
        paths.append(path)
    
    return paths


def enumerate_all_paths(
    distances_file: str,
    consensus: str,
    output_file: str,
    max_distance: int = 3,
    tier_filter: list = None,
    logger=None
):
    """
    Enumerate mutation paths for all motif sites.
    
    Parameters
    ----------
    distances_file : str
        Path to distances file from compute_distance.py
    consensus : str
        Consensus sequence
    output_file : str
        Output file path
    max_distance : int
        Maximum Hamming distance (number of mutations)
    tier_filter : list
        List of tiers to include (e.g., [1, 2] for Tier 1 & 2)
    logger : logging.Logger
        Logger instance
        
    Returns
    -------
    pd.DataFrame
        All mutation paths
    """
    logger.info("Loading distance data...")
    df = pd.read_csv(distances_file, sep='\t', low_memory=False)
    
    logger.info(f"  Total sites: {len(df):,}")
    
    # Filter by tier if specified
    if tier_filter is not None:
        before = len(df)
        df = df[df['tier'].isin(tier_filter)]
        after = len(df)
        logger.info(f"  Filtered to tiers {tier_filter}: {after:,} / {before:,} sites")
    
    # Filter by distance
    df = df[df['hamming_distance'] <= max_distance]
    logger.info(f"  Sites with distance ≤ {max_distance}: {len(df):,}")
    
    # Generate paths
    logger.info(f"\nGenerating mutation paths...")
    
    all_paths = []
    
    progress = ProgressLogger(total=len(df), logger=logger)
    
    for idx, row in df.iterrows():
        # Parse mismatches
        if pd.isna(row['mismatch_positions']) or row['mismatch_positions'] == '':
            # Perfect match (distance = 0)
            continue
        
        positions = [int(x) for x in row['mismatch_positions'].split(',')]
        ref_bases = row['ref_bases'].split(',')
        cons_bases = row['consensus_bases'].split(',')
        
        mismatches = list(zip(positions, ref_bases, cons_bases))
        
        # Generate paths
        paths = generate_mutation_paths(
            site_seq=row['matched_sequence'].upper(),
            consensus=consensus,
            mismatch_positions=mismatches,
            max_steps=max_distance
        )
        
        # Store each path
        for path_num, path in enumerate(paths):
            for step_num, step in enumerate(path):
                # For minus-strand motifs, convert ref/alt to genomic (plus strand) orientation
                # VCF format always uses plus strand reference, so we need to reverse complement
                # the bases when the motif is on the minus strand
                if row['strand'] == '-':
                    ref_genomic = reverse_complement(step['ref'])
                    alt_genomic = reverse_complement(step['alt'])
                else:
                    ref_genomic = step['ref']
                    alt_genomic = step['alt']
                
                path_record = {
                    'site_id': idx,
                    'chr': row['chr'],
                    'motif_start': row['start'],
                    'motif_end': row['stop'],
                    'strand': row['strand'],
                    'tier': row['tier'],
                    'pwm_score': row['score'],
                    'hamming_distance': row['hamming_distance'],
                    'path_id': f"{idx}_path{path_num}",
                    'path_num': path_num,
                    'step_num': step_num,
                    'total_steps': len(path),
                    'position_in_motif': step['position'],
                    'genomic_position': row['start'] + step['position'] if row['strand'] == '+' else row['stop'] - step['position'] - 1,
                    'ref_base': ref_genomic,        # Genomic orientation (plus strand)
                    'alt_base': alt_genomic,        # Genomic orientation (plus strand)
                    'ref_base_motif': step['ref'],  # Original motif orientation (for reference)
                    'alt_base_motif': step['alt'],  # Original motif orientation (for reference)
                    'seq_before': step['seq_before'],
                    'seq_after': step['seq_after']
                }
                all_paths.append(path_record)
        
        progress.update()
    
    # Progress complete - final log is already done by update()
    
    # Convert to dataframe
    paths_df = pd.DataFrame(all_paths)
    
    logger.info(f"\nGenerated {len(paths_df):,} mutation steps across {paths_df['path_id'].nunique():,} paths")
    logger.info(f"  From {paths_df['site_id'].nunique():,} unique motif sites")
    
    # Log path statistics
    logger.info("\nPath length distribution:")
    length_counts = paths_df.groupby('path_id')['total_steps'].first().value_counts().sort_index()
    for length, count in length_counts.items():
        percent = 100 * count / len(paths_df['path_id'].unique())
        logger.info(f"  {int(length)} steps: {count:,} paths ({percent:.1f}%)")
    
    # Save results
    paths_df.to_csv(output_file, sep='\t', index=False)
    logger.info(f"\n✓ Saved mutation paths: {output_file}")
    
    return paths_df


def main():
    parser = argparse.ArgumentParser(
        description="Enumerate mutation paths from motif sites to consensus"
    )
    parser.add_argument(
        '--config',
        type=str,
        required=True,
        help='Path to pipeline configuration file'
    )
    parser.add_argument(
        '--distances-file',
        type=str,
        help='Path to distances file (default: from compute_distance.py)'
    )
    parser.add_argument(
        '--consensus-file',
        type=str,
        help='Path to consensus file'
    )
    parser.add_argument(
        '--output-file',
        type=str,
        help='Output file'
    )
    parser.add_argument(
        '--max-distance',
        type=int,
        help='Maximum Hamming distance (default: from config)'
    )
    parser.add_argument(
        '--tiers',
        type=str,
        help='Comma-separated list of tiers to include (e.g., "1,2")'
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
    tier_filter_default = config['mutation_paths'].get('tier_filter', [1, 2])
    
    # Determine inputs
    if args.distances_file:
        distances_file = args.distances_file
    else:
        distances_file = f"results/mutation_paths/{tf_name}/distances.tsv"
    
    if args.consensus_file:
        consensus_file = args.consensus_file
    else:
        consensus_file = f"data/motifs/consensus/{tf_name}_consensus.txt"
    
    if args.output_file:
        output_file = args.output_file
    else:
        output_file = f"results/mutation_paths/{tf_name}/paths.tsv"
    
    max_distance = args.max_distance or max_distance_default
    
    if args.tiers:
        tier_filter = [int(x) for x in args.tiers.split(',')]
    else:
        tier_filter = tier_filter_default
    
    # Setup logger
    log_file = args.log_file or f"logs/02_enumerate_paths_{tf_name}.log"
    logger = setup_logger(
        name=__name__,
        log_file=log_file,
        level='INFO'
    )
    
    log_step(logger, "Module 02: Mutation Path Enumeration", start=True)
    
    # Log parameters
    params = {
        'TF': tf_name,
        'Distances file': distances_file,
        'Max distance': max_distance,
        'Tier filter': tier_filter,
        'Output file': output_file
    }
    log_parameters(logger, params)
    
    # Check inputs
    if not Path(distances_file).exists():
        logger.error(f"Distances file not found: {distances_file}")
        logger.error("Run compute_distance.py first")
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
    
    # Enumerate paths
    try:
        paths_df = enumerate_all_paths(
            distances_file=distances_file,
            consensus=consensus,
            output_file=output_file,
            max_distance=max_distance,
            tier_filter=tier_filter,
            logger=logger
        )
    except Exception as e:
        logger.error(f"Path enumeration failed: {e}")
        import traceback
        logger.error(traceback.format_exc())
        sys.exit(1)
    
    log_step(logger, "Module 02: Mutation Path Enumeration", start=False)
    logger.info("\nNext step: Module 03 - Intersect with gnomAD")
    logger.info("  python 03_intersect_gnomad/query_gnomad_vcfs.py --config pipeline_config.yaml")


if __name__ == '__main__':
    main()
