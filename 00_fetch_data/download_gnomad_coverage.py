#!/usr/bin/env python3
"""
Download and setup gnomAD v4.1 allele number (coverage) data with symlinks
Part of Dormant Site Activation Pipeline - Module 00

From gnomAD v4.1 release notes (April 2024):
  "We have calculated AN information for loci that do not yet have variation 
   observed within gnomAD. This means that even for loci where zero gnomAD 
   samples have a non-reference genotype call, we report the exact AN based 
   on the total number of defined sample genotype calls at that site."

This file provides AN (allele number = 2 × sample count) for ALL callable sites,
which is essential for distinguishing:
  - "Truly constrained" (high AN, variant absent → strong purifying selection)
  - "Uncertain" (low AN or missing → variant might exist but undetected)

File: gnomad.genomes.v4.1.allele_number_all_sites.tsv.bgz (~12 GB)
Location: gs://gcp-public-data--gnomad/release/4.1/tsv/genomes/
"""

import argparse
import subprocess
import sys
from pathlib import Path
import json

# Storage locations (following same pattern as download_gnomad.py)
STORAGE_ROOT = Path("/mnt/data_1/gnomAD_data/raw/gnomad_v4.1")  # Large HDD storage
WORK_ROOT = Path("/mnt/work_1/gest9386/CU_Boulder/rotations/LAYER/dormant_site_activation_pipeline/data/gnomad_coverage")  # SSD symlinks

# Google Cloud Storage path for gnomAD v4.1 allele number data
GCS_BASE = "gs://gcp-public-data--gnomad/release/4.1/tsv/genomes"
AN_FILE = "gnomad.genomes.v4.1.allele_number_all_sites.tsv.bgz"

# Expected file size (~12 GB)
EXPECTED_SIZE_GB = 12


def check_gsutil():
    """Check if gsutil is installed."""
    try:
        result = subprocess.run(['gsutil', 'version'], capture_output=True, text=True)
        print(f"✓ gsutil found: {result.stdout.split()[2]}")
        return True
    except FileNotFoundError:
        print("✗ gsutil not found")
        print("\nTo install gsutil:")
        print("  1. Install Google Cloud SDK: https://cloud.google.com/sdk/docs/install")
        print("  2. Or use conda: conda install -c conda-forge google-cloud-sdk")
        return False


def check_file_exists():
    """Check if the allele number file exists on GCS and get its size."""
    print("\n" + "="*80)
    print("Checking gnomAD v4.1 allele number file on Google Cloud Storage...")
    print("="*80)
    
    gcs_path = f"{GCS_BASE}/{AN_FILE}"
    cmd = ['gsutil', 'ls', '-l', gcs_path]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        lines = result.stdout.strip().split('\n')
        if lines:
            parts = lines[0].split()
            if len(parts) >= 2:
                size_bytes = int(parts[0])
                size_gb = size_bytes / (1024**3)
                print(f"\n✓ File found: {AN_FILE}")
                print(f"  Size: {size_gb:.2f} GB")
                print(f"  Source: {gcs_path}")
                return True, size_gb
    except subprocess.CalledProcessError as e:
        print(f"✗ File not found: {e.stderr}")
        return False, 0
    
    return False, 0


def download_an_file(storage_dir, dry_run=False):
    """Download the allele number file."""
    gcs_path = f"{GCS_BASE}/{AN_FILE}"
    local_path = storage_dir / AN_FILE
    
    print(f"\n{'[DRY RUN] ' if dry_run else ''}Downloading allele number file:")
    print(f"  Source: {gcs_path}")
    print(f"  Destination: {local_path}")
    print(f"  Expected size: ~{EXPECTED_SIZE_GB} GB")
    
    # Check if already downloaded
    if local_path.exists():
        existing_size = local_path.stat().st_size / (1024**3)
        print(f"\n✓ File already exists ({existing_size:.2f} GB)")
        return True
    
    if dry_run:
        print(f"\n→ Would download ~{EXPECTED_SIZE_GB} GB file")
        return True
    
    # Download
    print(f"\nDownloading... (this may take 10-30 minutes)")
    cmd = ['gsutil', '-m', 'cp', gcs_path, str(local_path)]
    
    try:
        subprocess.run(cmd, check=True)
        print(f"✓ Download complete")
        return True
    except subprocess.CalledProcessError as e:
        print(f"✗ Download failed: {e}")
        return False


def create_symlink(storage_dir, work_dir, dry_run=False):
    """Create symlink from work directory to storage."""
    print("\n" + "="*80)
    print(f"{'[DRY RUN] ' if dry_run else ''}Creating symlink...")
    print("="*80)
    
    source = storage_dir / AN_FILE
    link = work_dir / AN_FILE
    
    print(f"\nStorage location: {source}")
    print(f"Work location: {link}")
    
    if not source.exists() and not dry_run:
        print(f"✗ Source file not found: {source}")
        return False
    
    if dry_run:
        print(f"→ Would create symlink: {link} -> {source}")
        return True
    
    # Create work directory
    work_dir.mkdir(parents=True, exist_ok=True)
    
    # Create symlink
    try:
        if link.exists() or link.is_symlink():
            link.unlink()
        link.symlink_to(source)
        print(f"✓ Symlink created: {link.name}")
        return True
    except Exception as e:
        print(f"✗ Failed to create symlink: {e}")
        return False


def create_tabix_index(storage_dir, dry_run=False):
    """Create tabix index for the file (enables fast random access)."""
    file_path = storage_dir / AN_FILE
    index_path = storage_dir / f"{AN_FILE}.tbi"
    
    print("\n" + "="*80)
    print(f"{'[DRY RUN] ' if dry_run else ''}Creating tabix index...")
    print("="*80)
    
    if index_path.exists():
        print(f"✓ Index already exists: {index_path.name}")
        return True
    
    if dry_run:
        print(f"→ Would create index with: tabix -s1 -b2 -e2 {file_path}")
        return True
    
    # Check if tabix is available
    try:
        subprocess.run(['tabix', '--version'], capture_output=True, check=True)
    except (FileNotFoundError, subprocess.CalledProcessError):
        print("⚠ tabix not found - skipping index creation")
        print("  Install with: conda install -c bioconda htslib")
        print("  Index enables fast random access by genomic position")
        return False
    
    print(f"Creating index (this may take several minutes)...")
    cmd = ['tabix', '-s1', '-b2', '-e2', str(file_path)]
    
    try:
        subprocess.run(cmd, check=True)
        print(f"✓ Index created: {index_path.name}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"✗ Index creation failed: {e}")
        return False


def save_manifest(storage_dir, work_dir):
    """Save manifest with usage documentation."""
    file_path = storage_dir / AN_FILE
    
    manifest = {
        'storage_location': str(storage_dir),
        'work_location': str(work_dir),
        'gnomad_version': '4.1',
        'data_type': 'allele_number_all_sites',
        'file': AN_FILE,
        'description': (
            'Allele number (AN) for ALL callable sites in gnomAD v4.1 genomes. '
            'From v4.1 release: "We have calculated AN information for loci that '
            'do not yet have variation observed within gnomAD."'
        ),
        'columns': {
            'locus': 'Genomic position (chr:pos format)',
            'AN': 'Allele number = 2 × number of samples with genotype call'
        },
        'interpretation': {
            'high_AN': 'AN >= 100,000 (~50K samples) = high confidence coverage',
            'medium_AN': 'AN >= 50,000 (~25K samples) = medium confidence',
            'low_AN': 'AN < 50,000 = low confidence, interpret with caution',
            'missing': 'Position not in file = no genotype calls (uncallable region)'
        },
        'usage': {
            'constraint_analysis': (
                'If a position has high AN AND no observed variant in gnomAD VCF, '
                'this is strong evidence of purifying selection (truly constrained). '
                'If a position has low AN or is missing, the absence of a variant '
                'could be due to poor coverage rather than selection.'
            ),
            'recommended_threshold': 'AN >= 100,000 for high-confidence constraint calls'
        }
    }
    
    if file_path.exists():
        manifest['size_gb'] = round(file_path.stat().st_size / (1024**3), 2)
    
    manifest_file = storage_dir / 'allele_number_manifest.json'
    with open(manifest_file, 'w') as f:
        json.dump(manifest, f, indent=2)
    
    print(f"\n✓ Saved manifest: {manifest_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Download gnomAD v4.1 allele number (coverage) data and create symlinks",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Check if file exists on GCS
  python download_gnomad_coverage.py --list-only
  
  # Preview what would be downloaded
  python download_gnomad_coverage.py --dry-run
  
  # Download the file (~12 GB)
  python download_gnomad_coverage.py
  
  # Create symlink only (if file already downloaded)
  python download_gnomad_coverage.py --skip-download
"""
    )
    parser.add_argument(
        '--list-only',
        action='store_true',
        help='Only check if file exists on GCS, do not download'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Show what would be downloaded without downloading'
    )
    parser.add_argument(
        '--skip-download',
        action='store_true',
        help='Skip download, only create symlinks (if file already exists)'
    )
    parser.add_argument(
        '--skip-index',
        action='store_true',
        help='Skip tabix index creation'
    )
    parser.add_argument(
        '--storage-dir',
        type=str,
        default=str(STORAGE_ROOT),
        help=f'Storage directory for downloaded file (default: {STORAGE_ROOT})'
    )
    parser.add_argument(
        '--work-dir',
        type=str,
        default=str(WORK_ROOT),
        help=f'Work directory for symlinks (default: {WORK_ROOT})'
    )
    
    args = parser.parse_args()
    
    print("="*80)
    print("gnomAD v4.1 Allele Number (Coverage) Data Download")
    print("="*80)
    print("\nThis data enables HIGH-CONFIDENCE constraint analysis.")
    print("\nFrom gnomAD v4.1 release notes:")
    print('  "We have calculated AN information for loci that do not yet have')
    print('   variation observed within gnomAD... even for loci where zero samples')
    print('   have a non-reference genotype call."')
    print(f"\nFile: {AN_FILE}")
    print(f"Size: ~{EXPECTED_SIZE_GB} GB")
    
    # Check gsutil
    if not check_gsutil():
        sys.exit(1)
    
    # Check file exists on GCS
    exists, size_gb = check_file_exists()
    if not exists:
        print("\n✗ Could not find gnomAD v4.1 allele number file on GCS")
        sys.exit(1)
    
    if args.list_only:
        print("\n[LIST ONLY] File verified on GCS, exiting")
        sys.exit(0)
    
    # Setup directories
    storage_dir = Path(args.storage_dir) / "coverage"
    work_dir = Path(args.work_dir)
    
    if not args.dry_run and not args.skip_download:
        storage_dir.mkdir(parents=True, exist_ok=True)
        print(f"\n✓ Storage directory: {storage_dir}")
    
    # Confirm download
    if not args.skip_download and not args.dry_run:
        print("\n" + "="*80)
        print("DOWNLOAD CONFIRMATION")
        print("="*80)
        print(f"File: {AN_FILE}")
        print(f"Size: {size_gb:.2f} GB")
        print(f"Storage: {storage_dir}")
        print("\nThis download may take 10-30 minutes.")
        
        response = input("\nProceed with download? (yes/no): ")
        if response.lower() not in ['yes', 'y']:
            print("Download cancelled")
            sys.exit(0)
    
    # Download file
    if not args.skip_download:
        success = download_an_file(storage_dir, dry_run=args.dry_run)
        if not success and not args.dry_run:
            print("\n✗ Download failed")
            sys.exit(1)
    
    # Create tabix index
    if not args.skip_index and not args.dry_run:
        create_tabix_index(storage_dir, dry_run=args.dry_run)
    
    # Create symlink
    create_symlink(storage_dir, work_dir, dry_run=args.dry_run)
    
    # Save manifest
    if not args.dry_run:
        save_manifest(storage_dir, work_dir)
    
    print("\n" + "="*80)
    print("SETUP COMPLETE")
    print("="*80)
    print(f"\nStorage: {storage_dir / AN_FILE}")
    print(f"Symlink: {work_dir / AN_FILE}")
    
    if args.dry_run:
        print("\n[DRY RUN] No files were downloaded or modified")
    
    print("\n" + "-"*80)
    print("NEXT STEPS")
    print("-"*80)
    print("\n1. Update pipeline_config.yaml:")
    print("   gnomad:")
    print(f"     allele_number_file: \"{work_dir / AN_FILE}\"")
    print("\n2. Update module 03 (gnomAD intersection) to use AN data:")
    print("   - Query AN for each mutation position")
    print("   - Flag positions with high AN (≥100K) as 'high_confidence'")
    print("   - Flag positions with low AN (<50K) as 'low_confidence'")
    print("\n3. Re-run pipeline for HIGH-CONFIDENCE constraint analysis:")
    print("   - H=1 variants absent with high AN = truly under purifying selection")
    print("   - H=1 variants absent with low AN = uncertain (may be coverage artifact)")


if __name__ == '__main__':
    main()
