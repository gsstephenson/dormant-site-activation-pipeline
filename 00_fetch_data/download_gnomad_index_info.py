#!/usr/bin/env python3
"""
Verify gnomAD installation and create symlink
Part of Dormant Site Activation Pipeline - Module 00
"""

import argparse
import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils.logging_utils import setup_logger
from utils.io import load_config


def check_file_exists(file_path: Path, logger) -> bool:
    """Check if file exists and log its size."""
    if file_path.exists():
        size_gb = file_path.stat().st_size / (1024**3)
        logger.info(f"✓ Found: {file_path} ({size_gb:.2f} GB)")
        return True
    else:
        logger.warning(f"✗ Missing: {file_path}")
        return False


def verify_gnomad_installation(config: dict, logger) -> bool:
    """
    Verify that gnomAD files are present.
    
    Parameters
    ----------
    config : dict
        Pipeline configuration
    logger : logging.Logger
        Logger instance
        
    Returns
    -------
    bool
        True if all required files are present
    """
    logger.info("=" * 80)
    logger.info("Verifying gnomAD v4.1 Installation")
    logger.info("=" * 80)
    
    gnomad_config = config.get('gnomad', {})
    
    # Check required files
    required_files = {
        'VCF': gnomad_config.get('vcf'),
        'VCF Index': gnomad_config.get('vcf_tbi'),
    }
    
    optional_files = {
        'Coverage': gnomad_config.get('coverage'),
        'Constraint': gnomad_config.get('constraint'),
    }
    
    logger.info("\nRequired Files:")
    logger.info("-" * 80)
    
    all_required_present = True
    for name, path in required_files.items():
        if path:
            file_path = Path(path)
            if not check_file_exists(file_path, logger):
                all_required_present = False
        else:
            logger.warning(f"✗ {name}: Path not configured")
            all_required_present = False
    
    logger.info("\nOptional Files:")
    logger.info("-" * 80)
    
    for name, path in optional_files.items():
        if path:
            path_obj = Path(path)
            if path_obj.exists():
                if path_obj.is_dir():
                    logger.info(f"✓ Found directory: {path}")
                else:
                    check_file_exists(path_obj, logger)
            else:
                logger.warning(f"✗ Not found: {path}")
        else:
            logger.info(f"○ {name}: Not configured (optional)")
    
    return all_required_present


def create_symlink(target: str, link_path: str, logger):
    """
    Create symbolic link to gnomAD data.
    
    Parameters
    ----------
    target : str
        Target directory (gnomAD root)
    link_path : str
        Path where symlink should be created
    logger : logging.Logger
        Logger instance
    """
    target_path = Path(target)
    link = Path(link_path)
    
    if not target_path.exists():
        logger.error(f"Target directory does not exist: {target}")
        return False
    
    # Create parent directory
    link.parent.mkdir(parents=True, exist_ok=True)
    
    # Remove existing symlink if present
    if link.is_symlink():
        logger.info(f"Removing existing symlink: {link}")
        link.unlink()
    elif link.exists():
        logger.error(f"Path exists but is not a symlink: {link}")
        return False
    
    # Create symlink
    link.symlink_to(target_path)
    logger.info(f"✓ Created symlink: {link} -> {target}")
    
    return True


def print_download_instructions(logger):
    """Print instructions for downloading gnomAD data."""
    logger.info("\n" + "=" * 80)
    logger.info("gnomAD Download Instructions")
    logger.info("=" * 80)
    logger.info("""
To download gnomAD v4.1 data:

1. Create the directory structure:
   sudo mkdir -p /mnt/data_1/gnomAD_data/raw/gnomad_v4.1
   sudo chown $USER:$USER /mnt/data_1/gnomAD_data

2. Download required files:
   cd /mnt/data_1/gnomAD_data/raw/gnomad_v4.1
   
   mkdir -p genomes_vcf
   gsutil cp gs://gcp-public-data--gnomad/release/4.1/genomes_vcf/genomes.vcf.bgz genomes_vcf/
   gsutil cp gs://gcp-public-data--gnomad/release/4.1/genomes_vcf/genomes.vcf.bgz.tbi genomes_vcf/

3. Download optional files (recommended):
   mkdir -p coverage
   gsutil -m cp -r gs://gcp-public-data--gnomad/release/4.1/coverage/genomes/* coverage/
   
   mkdir -p constraint
   gsutil cp gs://gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint.json.bgz constraint/

Note: Total download size is approximately 750 GB
    """)


def main():
    parser = argparse.ArgumentParser(
        description="Verify gnomAD installation and setup symlinks"
    )
    parser.add_argument(
        '--config',
        type=str,
        default='../pipeline_config.yaml',
        help='Path to pipeline configuration file'
    )
    parser.add_argument(
        '--create-symlink',
        action='store_true',
        help='Create symlink from data/gnomad to gnomAD root'
    )
    parser.add_argument(
        '--show-download-instructions',
        action='store_true',
        help='Show instructions for downloading gnomAD data'
    )
    parser.add_argument(
        '--log-file',
        type=str,
        default='../logs/gnomad_setup.log',
        help='Log file path'
    )
    
    args = parser.parse_args()
    
    # Setup logger
    logger = setup_logger(
        name=__name__,
        log_file=args.log_file,
        level='INFO'
    )
    
    # Show download instructions if requested
    if args.show_download_instructions:
        print_download_instructions(logger)
        return
    
    # Load config
    try:
        config = load_config(args.config)
    except Exception as e:
        logger.error(f"Failed to load config: {e}")
        sys.exit(1)
    
    # Verify installation
    all_present = verify_gnomad_installation(config, logger)
    
    # Create symlink if requested
    if args.create_symlink:
        logger.info("\n" + "=" * 80)
        logger.info("Creating Symlink")
        logger.info("=" * 80)
        
        gnomad_root = config['gnomad']['root']
        symlink_path = Path(__file__).parent.parent / 'data' / 'gnomad'
        
        create_symlink(gnomad_root, str(symlink_path), logger)
    
    # Summary
    logger.info("\n" + "=" * 80)
    logger.info("Summary")
    logger.info("=" * 80)
    
    if all_present:
        logger.info("✓ All required gnomAD files are present")
        logger.info("  You can proceed with the pipeline")
    else:
        logger.warning("✗ Some required files are missing")
        logger.info("  Run with --show-download-instructions for help")
        sys.exit(1)


if __name__ == '__main__':
    main()
