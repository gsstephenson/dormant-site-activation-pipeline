#!/usr/bin/env python3
"""
Download motif matrices from JASPAR database
Part of Dormant Site Activation Pipeline - Module 00
"""

import argparse
import requests
import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils.io import ensure_dir, save_config
from utils.logging_utils import setup_logger
from utils.constants import JASPAR_BASE_URL, JASPAR_MATRIX_URL


def download_jaspar_motif(motif_id: str, output_dir: str, logger):
    """
    Download a motif from JASPAR database.
    
    Parameters
    ----------
    motif_id : str
        JASPAR motif ID (e.g., MA0099.3 for AP1)
    output_dir : str
        Output directory for motif files
    logger : logging.Logger
        Logger instance
    """
    logger.info(f"Downloading motif {motif_id} from JASPAR...")
    
    # Fetch motif data
    url = f"{JASPAR_MATRIX_URL}/{motif_id}/"
    response = requests.get(url)
    
    if response.status_code != 200:
        logger.error(f"Failed to download motif {motif_id}: HTTP {response.status_code}")
        return None
    
    motif_data = response.json()
    
    # Extract information
    tf_name = motif_data.get('name', 'Unknown')
    tf_class = motif_data.get('class', 'Unknown')
    species = motif_data.get('species', [{}])[0].get('name', 'Unknown')
    
    logger.info(f"  TF Name: {tf_name}")
    logger.info(f"  Class: {tf_class}")
    logger.info(f"  Species: {species}")
    
    # PFM is included in the main motif data
    pfm_data = motif_data.get('pfm')
    
    if not pfm_data:
        logger.error("PFM data not found in motif response")
        return None
    
    # Save PFM in MEME format
    output_dir_path = Path(output_dir)
    ensure_dir(output_dir_path)
    
    pfm_file = output_dir_path / f"{motif_id}.meme"
    
    with open(pfm_file, 'w') as f:
        f.write("MEME version 4\n\n")
        f.write(f"MOTIF {motif_id} {tf_name}\n")
        f.write("letter-probability matrix: ")
        f.write(f"alength= 4 w= {len(pfm_data['A'])} nsites= 20\n")
        
        # Convert counts to frequencies (normalize)
        for i in range(len(pfm_data['A'])):
            total = pfm_data['A'][i] + pfm_data['C'][i] + pfm_data['G'][i] + pfm_data['T'][i]
            if total > 0:
                a_freq = pfm_data['A'][i] / total
                c_freq = pfm_data['C'][i] / total
                g_freq = pfm_data['G'][i] / total
                t_freq = pfm_data['T'][i] / total
            else:
                a_freq = c_freq = g_freq = t_freq = 0.25
            
            f.write(f" {a_freq:.6f} {c_freq:.6f} {g_freq:.6f} {t_freq:.6f}\n")
    
    logger.info(f"✓ Saved PFM: {pfm_file}")
    
    # Save raw JSON
    json_file = output_dir_path / f"{motif_id}.json"
    with open(json_file, 'w') as f:
        import json
        json.dump(motif_data, f, indent=2)
    
    logger.info(f"✓ Saved metadata: {json_file}")
    
    # Save metadata for pipeline
    metadata = {
        'motif_id': motif_id,
        'tf_name': tf_name,
        'class': tf_class,
        'species': species,
        'pfm_file': str(pfm_file),
        'length': len(pfm_data['A'])
    }
    
    return metadata


def main():
    parser = argparse.ArgumentParser(
        description="Download transcription factor motifs from JASPAR"
    )
    parser.add_argument(
        '--motif-id',
        type=str,
        default='MA0099.3',
        help='JASPAR motif ID (default: MA0099.3 for AP1)'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default='../data/motifs/JASPAR',
        help='Output directory for motif files'
    )
    parser.add_argument(
        '--log-file',
        type=str,
        default='../logs/download_motifs.log',
        help='Log file path'
    )
    
    args = parser.parse_args()
    
    # Setup logger
    logger = setup_logger(
        name=__name__,
        log_file=args.log_file,
        level='INFO'
    )
    
    logger.info("=" * 80)
    logger.info("JASPAR Motif Download")
    logger.info("=" * 80)
    
    # Download motif
    metadata = download_jaspar_motif(args.motif_id, args.output_dir, logger)
    
    if metadata:
        # Save metadata
        metadata_file = Path(args.output_dir) / 'motif_metadata.yaml'
        save_config(metadata, str(metadata_file))
        logger.info(f"✓ Saved pipeline metadata: {metadata_file}")
        
        logger.info("=" * 80)
        logger.info("Download complete!")
        logger.info("=" * 80)
    else:
        logger.error("Failed to download motif")
        sys.exit(1)


if __name__ == '__main__':
    main()
