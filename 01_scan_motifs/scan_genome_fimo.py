#!/usr/bin/env python3
"""
Genome-wide motif scanning using FIMO
Part of Dormant Site Activation Pipeline - Module 01
"""

import argparse
import subprocess
import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils.io import load_config, ensure_dir
from utils.logging_utils import setup_logger, log_step, log_parameters
from utils.constants import STANDARD_CHROMOSOMES


def run_fimo(
    motif_file: str,
    fasta_file: str,
    output_dir: str,
    pvalue: float = 1e-3,
    background: str = "uniform",
    logger=None
):
    """
    Run FIMO to scan genome for motif matches.
    
    Parameters
    ----------
    motif_file : str
        Path to motif file in MEME format
    fasta_file : str
        Path to reference genome FASTA
    output_dir : str
        Output directory for FIMO results
    pvalue : float
        P-value threshold for motif matches
    background : str
        Background model ("uniform" or path to background file)
    logger : logging.Logger
        Logger instance
        
    Returns
    -------
    str
        Path to FIMO output file (tsv.gz)
    """
    ensure_dir(output_dir)
    
    cmd = [
        "fimo",
        "--thresh", str(pvalue),
        "--verbosity", "2",
        "--text",  # Output to stdout
        motif_file,
        fasta_file
    ]
    
    # Add background if not uniform
    if background != "uniform":
        cmd.extend(["--bgfile", background])
    
    logger.info("Running FIMO...")
    logger.info(f"Command: {' '.join(cmd)}")
    
    output_file = Path(output_dir) / "fimo_output.tsv"
    
    try:
        with open(output_file, 'w') as f:
            result = subprocess.run(
                cmd,
                stdout=f,
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )
        
        logger.info(f"✓ FIMO completed successfully")
        logger.info(f"  Output: {output_file}")
        
        # Count hits
        with open(output_file, 'r') as f:
            n_lines = sum(1 for line in f if not line.startswith('#'))
        n_hits = n_lines - 1  # Subtract header
        
        logger.info(f"  Total motif hits: {n_hits:,}")
        
        return str(output_file)
        
    except subprocess.CalledProcessError as e:
        logger.error(f"FIMO failed: {e}")
        logger.error(f"stderr: {e.stderr}")
        raise
    except FileNotFoundError:
        logger.error("FIMO not found. Please install MEME Suite:")
        logger.error("  conda install -c bioconda meme")
        logger.error("  or visit: https://meme-suite.org/meme/doc/download.html")
        raise


def parse_fimo_output(fimo_file: str, logger=None):
    """
    Parse FIMO output TSV file.
    
    Parameters
    ----------
    fimo_file : str
        Path to FIMO output file
    logger : logging.Logger
        Logger instance
        
    Returns
    -------
    pd.DataFrame
        Parsed FIMO results
    """
    import pandas as pd
    
    logger.info("Parsing FIMO output...")
    
    # FIMO --text output has NO header row, just data
    # Columns: motif_id, sequence_name, start, stop, strand, score, p-value, q-value, matched_sequence
    df = pd.read_csv(
        fimo_file, 
        sep='\t', 
        comment='#',
        header=None,
        names=['motif_id', 'sequence_name', 'start', 'stop', 'strand', 'score', 'pvalue', 'qvalue', 'matched_sequence'],
        low_memory=False
    )
    
    logger.info(f"  Parsed {len(df):,} motif hits")
    logger.info(f"  Score range: [{df['score'].min():.2f}, {df['score'].max():.2f}]")
    logger.info(f"  Score mean: {df['score'].mean():.2f}")
    
    # Convert to 0-based coordinates (FIMO uses 1-based)
    df['start'] = df['start'] - 1
    
    # Add chromosome column (rename sequence_name)
    df['chr'] = df['sequence_name']
    
    return df


def convert_to_bed(df, output_file: str, logger=None):
    """
    Convert FIMO results to BED format.
    
    Parameters
    ----------
    df : pd.DataFrame
        FIMO results DataFrame
    output_file : str
        Output BED file path
    logger : logging.Logger
        Logger instance
    """
    import pandas as pd
    
    logger.info("Converting to BED format...")
    
    # Create BED format DataFrame
    bed_df = pd.DataFrame({
        'chr': df['chr'],
        'start': df['start'],
        'end': df['stop'],
        'name': df['motif_id'] + '_' + df.index.astype(str),
        'score': df['score'],
        'strand': df['strand'],
        'sequence': df['matched_sequence'],
        'pvalue': df['pvalue']
    })
    
    # Sort by chromosome and position
    bed_df = bed_df.sort_values(['chr', 'start'])
    
    bed_df.to_csv(output_file, sep='\t', header=False, index=False)
    
    logger.info(f"✓ Saved BED file: {output_file}")
    logger.info(f"  {len(bed_df):,} motif sites")


def main():
    parser = argparse.ArgumentParser(
        description="Scan genome for transcription factor motif matches using FIMO"
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
        help='Path to motif file (overrides config)'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        help='Output directory (default: results/motif_scan/{TF})'
    )
    parser.add_argument(
        '--pvalue',
        type=float,
        help='P-value threshold (overrides config)'
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
    reference_genome = config['reference_genome']
    pvalue = args.pvalue or config['motif_scan']['pvalue_threshold']
    background = config['motif_scan'].get('background', 'uniform')
    
    # Determine motif file
    if args.motif_file:
        motif_file = args.motif_file
    else:
        motif_file = f"data/motifs/JASPAR/{motif_id}.meme"
    
    # Determine output directory
    if args.output_dir:
        output_dir = args.output_dir
    else:
        output_dir = f"results/motif_scan/{tf_name}"
    
    # Setup logger
    log_file = args.log_file or f"logs/01_scan_motifs_{tf_name}.log"
    logger = setup_logger(
        name=__name__,
        log_file=log_file,
        level='INFO'
    )
    
    log_step(logger, "Module 01: Genome-wide Motif Scanning", start=True)
    
    # Log parameters
    params = {
        'TF': tf_name,
        'Motif ID': motif_id,
        'Motif file': motif_file,
        'Reference genome': reference_genome,
        'P-value threshold': pvalue,
        'Background': background,
        'Output directory': output_dir
    }
    log_parameters(logger, params)
    
    # Check input files exist
    if not Path(motif_file).exists():
        logger.error(f"Motif file not found: {motif_file}")
        logger.error("Run Module 00 to download motifs:")
        logger.error("  python 00_fetch_data/download_motifs.py")
        sys.exit(1)
    
    if not Path(reference_genome).exists():
        logger.error(f"Reference genome not found: {reference_genome}")
        logger.error("Run Module 00 to download reference:")
        logger.error("  bash 00_fetch_data/download_reference.sh")
        sys.exit(1)
    
    # Run FIMO
    try:
        fimo_output = run_fimo(
            motif_file=motif_file,
            fasta_file=reference_genome,
            output_dir=output_dir,
            pvalue=pvalue,
            background=background,
            logger=logger
        )
    except Exception as e:
        logger.error(f"FIMO scanning failed: {e}")
        sys.exit(1)
    
    # Parse FIMO output
    df = parse_fimo_output(fimo_output, logger)
    
    # Convert to BED format
    bed_file = Path(output_dir) / "motif_hits.bed"
    convert_to_bed(df, str(bed_file), logger)
    
    # Save full results as TSV
    tsv_file = Path(output_dir) / "motif_hits_full.tsv"
    df.to_csv(tsv_file, sep='\t', index=False)
    logger.info(f"✓ Saved full results: {tsv_file}")
    
    log_step(logger, "Module 01: Genome-wide Motif Scanning", start=False)
    logger.info("\nNext step: Module 02 - Tier motif sites and generate mutation paths")
    logger.info("  python 01_scan_motifs/tier_sites.py --config pipeline_config.yaml")


if __name__ == '__main__':
    main()
