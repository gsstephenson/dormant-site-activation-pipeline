#!/usr/bin/env python3
"""
Infer consensus sequence from PWM
Part of Dormant Site Activation Pipeline - Module 02
"""

import argparse
import sys
from pathlib import Path
import numpy as np

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils.io import load_config
from utils.logging_utils import setup_logger, log_step


def parse_meme_pwm(meme_file: str, logger=None):
    """
    Parse PWM from MEME format file.
    
    Parameters
    ----------
    meme_file : str
        Path to MEME format motif file
    logger : logging.Logger
        Logger instance
        
    Returns
    -------
    tuple
        (motif_id, pwm_matrix) where pwm_matrix is numpy array [positions x 4]
    """
    with open(meme_file, 'r') as f:
        lines = f.readlines()
    
    # Find MOTIF line
    motif_id = None
    pwm_lines = []
    in_matrix = False
    
    for i, line in enumerate(lines):
        if line.startswith('MOTIF'):
            motif_id = line.split()[1]
        
        if line.startswith('letter-probability matrix'):
            in_matrix = True
            continue
        
        if in_matrix:
            # Read PWM rows until blank line
            stripped = line.strip()
            if not stripped:
                break
            
            try:
                values = [float(x) for x in stripped.split()]
                if len(values) == 4:  # A C G T
                    pwm_lines.append(values)
            except ValueError:
                continue
    
    pwm = np.array(pwm_lines)
    
    if logger:
        logger.info(f"Parsed motif: {motif_id}")
        logger.info(f"  PWM shape: {pwm.shape[0]} positions x 4 bases")
    
    return motif_id, pwm


def pwm_to_consensus(pwm: np.ndarray, logger=None):
    """
    Convert PWM to consensus sequence.
    
    Takes the base with maximum probability at each position.
    
    Parameters
    ----------
    pwm : np.ndarray
        PWM matrix [positions x 4], columns in order [A, C, G, T]
    logger : logging.Logger
        Logger instance
        
    Returns
    -------
    str
        Consensus sequence
    """
    bases = ['A', 'C', 'G', 'T']
    
    consensus = []
    for pos in range(pwm.shape[0]):
        max_idx = np.argmax(pwm[pos, :])
        max_prob = pwm[pos, max_idx]
        consensus.append(bases[max_idx])
        
        if logger and max_prob < 0.5:
            logger.warning(f"  Position {pos+1}: max prob = {max_prob:.3f} (ambiguous)")
    
    consensus_seq = ''.join(consensus)
    
    if logger:
        logger.info(f"Consensus sequence: {consensus_seq}")
        
        # Calculate information content
        ic_per_pos = []
        for pos in range(pwm.shape[0]):
            # IC = 2 + sum(p_i * log2(p_i))
            probs = pwm[pos, :]
            # Avoid log(0)
            probs_safe = probs[probs > 0]
            ic = 2 + np.sum(probs_safe * np.log2(probs_safe))
            ic_per_pos.append(ic)
        
        total_ic = sum(ic_per_pos)
        logger.info(f"  Total information content: {total_ic:.2f} bits")
        logger.info(f"  Average IC per position: {total_ic/len(ic_per_pos):.2f} bits")
    
    return consensus_seq


def calculate_max_pwm_score(pwm: np.ndarray):
    """
    Calculate maximum possible PWM score (consensus score).
    
    Parameters
    ----------
    pwm : np.ndarray
        PWM matrix
        
    Returns
    -------
    float
        Maximum PWM score
    """
    # Max score = sum of max probability at each position
    # PWM scores in FIMO are typically log-odds scores
    # For simplicity, we use sum of max probabilities
    max_score = np.sum(np.max(pwm, axis=1))
    return max_score


def main():
    parser = argparse.ArgumentParser(
        description="Infer consensus sequence from PWM"
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
        help='Path to MEME format motif file (default: from config)'
    )
    parser.add_argument(
        '--output-file',
        type=str,
        help='Output consensus file (default: data/motifs/consensus/{TF}_consensus.txt)'
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
    motif_id = config['motif_id']
    
    # Determine motif file
    if args.motif_file:
        motif_file = args.motif_file
    else:
        motif_file = f"data/motifs/JASPAR/{motif_id}.meme"
    
    # Determine output file
    if args.output_file:
        output_file = args.output_file
    else:
        output_dir = Path("data/motifs/consensus")
        output_dir.mkdir(parents=True, exist_ok=True)
        output_file = output_dir / f"{tf_name}_consensus.txt"
    
    # Setup logger
    log_file = args.log_file or f"logs/02_consensus_{tf_name}.log"
    logger = setup_logger(
        name=__name__,
        log_file=log_file,
        level='INFO'
    )
    
    log_step(logger, "Module 02: Consensus Sequence Inference", start=True)
    
    # Check input file exists
    if not Path(motif_file).exists():
        logger.error(f"Motif file not found: {motif_file}")
        logger.error("Run Module 00 first to download motifs")
        sys.exit(1)
    
    # Parse PWM
    try:
        motif_id_parsed, pwm = parse_meme_pwm(motif_file, logger)
    except Exception as e:
        logger.error(f"Failed to parse PWM: {e}")
        import traceback
        logger.error(traceback.format_exc())
        sys.exit(1)
    
    # Generate consensus
    consensus = pwm_to_consensus(pwm, logger)
    
    # Calculate max PWM score
    max_score = calculate_max_pwm_score(pwm)
    logger.info(f"Maximum PWM score: {max_score:.3f}")
    
    # Save consensus
    with open(output_file, 'w') as f:
        f.write(f"# Consensus sequence for {tf_name} ({motif_id})\n")
        f.write(f"# Motif length: {len(consensus)} bp\n")
        f.write(f"# Maximum PWM score: {max_score:.3f}\n")
        f.write(f"{consensus}\n")
    
    logger.info(f"\n✓ Saved consensus: {output_file}")
    
    log_step(logger, "Module 02: Consensus Sequence Inference", start=False)


if __name__ == '__main__':
    main()
