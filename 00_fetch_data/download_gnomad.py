#!/usr/bin/env python3
"""
Download and setup gnomAD v4.1 data with symlinks
Part of Dormant Site Activation Pipeline - Module 00
"""

import argparse
import subprocess
import sys
from pathlib import Path
import json

# Storage locations
STORAGE_ROOT = Path("/mnt/data_1/gnomAD_data/raw/gnomad_v4.1")  # Large HDD storage
WORK_ROOT = Path("/mnt/work_1/gest9386/CU_Boulder/rotations/LAYER/dormant_site_activation_pipeline/data/gnomad")  # SSD symlinks

# Google Cloud Storage paths
GCS_BASE = "gs://gcp-public-data--gnomad/release/4.1/vcf/genomes"

# Chromosomes to download
CHROMOSOMES = [str(i) for i in range(1, 23)] + ['X', 'Y']


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


def list_available_files():
    """List available gnomAD files on GCS."""
    print("\n" + "="*80)
    print("Checking available gnomAD v4.1 files on Google Cloud Storage...")
    print("="*80)
    
    cmd = ['gsutil', 'ls', f'{GCS_BASE}/']
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        files = [line.strip() for line in result.stdout.strip().split('\n') if line.strip()]
        
        # Filter for chr files
        chr_files = [f for f in files if 'chr' in f and 'vcf.bgz' in f]
        
        print(f"\nFound {len(chr_files)} chromosome VCF files")
        
        # Show sample
        print("\nSample files:")
        for f in chr_files[:5]:
            print(f"  {Path(f).name}")
        
        if len(chr_files) > 5:
            print(f"  ... and {len(chr_files) - 5} more")
        
        return chr_files
        
    except subprocess.CalledProcessError as e:
        print(f"✗ Failed to list files: {e}")
        print(f"stderr: {e.stderr}")
        return []


def get_file_size(gcs_path):
    """Get size of file on GCS."""
    cmd = ['gsutil', 'du', '-s', gcs_path]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        size_bytes = int(result.stdout.split()[0])
        size_gb = size_bytes / (1024**3)
        return size_gb
    except:
        return None


def estimate_total_size(chromosomes):
    """Estimate total download size."""
    print("\n" + "="*80)
    print("Estimating download sizes...")
    print("="*80)
    
    # Check size of chr1 as reference
    chr1_file = f"{GCS_BASE}/gnomad.genomes.v4.1.sites.chr1.vcf.bgz"
    chr1_size = get_file_size(chr1_file)
    
    if chr1_size:
        print(f"\nChromosome 1 VCF size: {chr1_size:.2f} GB")
        
        # Estimate total (rough approximation)
        estimated_total = chr1_size * len(chromosomes) * 0.8  # chr1 is usually larger
        print(f"Estimated total download: {estimated_total:.2f} GB")
        print(f"  (for {len(chromosomes)} chromosomes)")
        
        return estimated_total
    else:
        print("\n⚠ Could not estimate size")
        return None


def download_chromosome(chrom, storage_dir, dry_run=False):
    """Download VCF and index for one chromosome."""
    vcf_file = f"gnomad.genomes.v4.1.sites.chr{chrom}.vcf.bgz"
    tbi_file = f"{vcf_file}.tbi"
    
    vcf_gcs = f"{GCS_BASE}/{vcf_file}"
    tbi_gcs = f"{GCS_BASE}/{tbi_file}"
    
    vcf_local = storage_dir / vcf_file
    tbi_local = storage_dir / tbi_file
    
    print(f"\n{'[DRY RUN] ' if dry_run else ''}Chromosome {chrom}:")
    print(f"  VCF: {vcf_file}")
    
    # Check if already downloaded
    if vcf_local.exists() and tbi_local.exists():
        print(f"  ✓ Already exists, skipping")
        return True
    
    if dry_run:
        print(f"  → Would download from: {vcf_gcs}")
        print(f"  → Would save to: {vcf_local}")
        return True
    
    # Download VCF
    print(f"  Downloading VCF...")
    cmd_vcf = ['gsutil', '-m', 'cp', vcf_gcs, str(vcf_local)]
    
    try:
        subprocess.run(cmd_vcf, check=True)
        print(f"  ✓ VCF downloaded")
    except subprocess.CalledProcessError as e:
        print(f"  ✗ VCF download failed: {e}")
        return False
    
    # Download index
    print(f"  Downloading index...")
    cmd_tbi = ['gsutil', '-m', 'cp', tbi_gcs, str(tbi_local)]
    
    try:
        subprocess.run(cmd_tbi, check=True)
        print(f"  ✓ Index downloaded")
    except subprocess.CalledProcessError as e:
        print(f"  ✗ Index download failed: {e}")
        return False
    
    return True


def create_symlinks(storage_dir, work_dir, chromosomes, dry_run=False):
    """Create symlinks from work directory to storage."""
    print("\n" + "="*80)
    print(f"{'[DRY RUN] ' if dry_run else ''}Creating symlinks...")
    print("="*80)
    
    print(f"\nStorage location: {storage_dir}")
    print(f"Work location: {work_dir}")
    
    if not dry_run:
        work_dir.mkdir(parents=True, exist_ok=True)
    
    success_count = 0
    
    for chrom in chromosomes:
        vcf_file = f"gnomad.genomes.v4.1.sites.chr{chrom}.vcf.bgz"
        tbi_file = f"{vcf_file}.tbi"
        
        vcf_source = storage_dir / vcf_file
        tbi_source = storage_dir / tbi_file
        
        vcf_link = work_dir / vcf_file
        tbi_link = work_dir / tbi_file
        
        if not vcf_source.exists():
            print(f"  ⚠ chr{chrom}: Source file not found, skipping")
            continue
        
        if dry_run:
            print(f"  chr{chrom}: Would create symlink")
            print(f"    {vcf_link} -> {vcf_source}")
            success_count += 1
            continue
        
        # Create VCF symlink
        try:
            if vcf_link.exists() or vcf_link.is_symlink():
                vcf_link.unlink()
            vcf_link.symlink_to(vcf_source)
            
            # Create index symlink
            if tbi_link.exists() or tbi_link.is_symlink():
                tbi_link.unlink()
            tbi_link.symlink_to(tbi_source)
            
            print(f"  ✓ chr{chrom}: Symlinks created")
            success_count += 1
            
        except Exception as e:
            print(f"  ✗ chr{chrom}: Failed to create symlink: {e}")
    
    print(f"\nCreated symlinks for {success_count}/{len(chromosomes)} chromosomes")
    
    return success_count


def save_manifest(storage_dir, work_dir, chromosomes):
    """Save manifest of downloaded files."""
    manifest = {
        'storage_location': str(storage_dir),
        'work_location': str(work_dir),
        'gnomad_version': '4.1',
        'chromosomes': chromosomes,
        'files': []
    }
    
    for chrom in chromosomes:
        vcf_file = f"gnomad.genomes.v4.1.sites.chr{chrom}.vcf.bgz"
        vcf_path = storage_dir / vcf_file
        
        if vcf_path.exists():
            size_gb = vcf_path.stat().st_size / (1024**3)
            manifest['files'].append({
                'chromosome': chrom,
                'file': vcf_file,
                'size_gb': round(size_gb, 2),
                'path': str(vcf_path)
            })
    
    manifest_file = storage_dir / 'download_manifest.json'
    with open(manifest_file, 'w') as f:
        json.dump(manifest, f, indent=2)
    
    print(f"\n✓ Saved manifest: {manifest_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Download gnomAD v4.1 data and create symlinks"
    )
    parser.add_argument(
        '--chromosomes',
        type=str,
        default='21,22',
        help='Comma-separated list of chromosomes (default: 21,22 for testing)'
    )
    parser.add_argument(
        '--all-chromosomes',
        action='store_true',
        help='Download all chromosomes (1-22,X,Y)'
    )
    parser.add_argument(
        '--list-only',
        action='store_true',
        help='Only list available files, do not download'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Show what would be downloaded without downloading'
    )
    parser.add_argument(
        '--skip-download',
        action='store_true',
        help='Skip download, only create symlinks'
    )
    parser.add_argument(
        '--storage-dir',
        type=str,
        default=str(STORAGE_ROOT),
        help=f'Storage directory (default: {STORAGE_ROOT})'
    )
    parser.add_argument(
        '--work-dir',
        type=str,
        default=str(WORK_ROOT),
        help=f'Work directory for symlinks (default: {WORK_ROOT})'
    )
    
    args = parser.parse_args()
    
    print("="*80)
    print("gnomAD v4.1 Download and Setup")
    print("="*80)
    
    # Check gsutil
    if not check_gsutil():
        sys.exit(1)
    
    # List files
    available_files = list_available_files()
    if not available_files:
        print("\n✗ Could not access gnomAD files on Google Cloud Storage")
        sys.exit(1)
    
    if args.list_only:
        sys.exit(0)
    
    # Determine chromosomes
    if args.all_chromosomes:
        chromosomes = CHROMOSOMES
    else:
        chromosomes = args.chromosomes.split(',')
    
    print(f"\nChromosomes to process: {', '.join(chromosomes)}")
    
    # Estimate size
    estimated_size = estimate_total_size(chromosomes)
    
    # Setup directories
    storage_dir = Path(args.storage_dir) / "vcf"
    work_dir = Path(args.work_dir)
    
    if not args.dry_run and not args.skip_download:
        storage_dir.mkdir(parents=True, exist_ok=True)
        print(f"\n✓ Created storage directory: {storage_dir}")
    
    # Confirm download
    if not args.skip_download and not args.dry_run:
        print("\n" + "="*80)
        print("DOWNLOAD CONFIRMATION")
        print("="*80)
        print(f"Chromosomes: {', '.join(chromosomes)}")
        print(f"Storage location: {storage_dir}")
        if estimated_size:
            print(f"Estimated size: {estimated_size:.2f} GB")
        print("\nThis will download large files and may take hours.")
        
        response = input("\nProceed with download? (yes/no): ")
        if response.lower() not in ['yes', 'y']:
            print("Download cancelled")
            sys.exit(0)
    
    # Download files
    if not args.skip_download:
        print("\n" + "="*80)
        print(f"{'[DRY RUN] ' if args.dry_run else ''}DOWNLOADING FILES")
        print("="*80)
        
        for i, chrom in enumerate(chromosomes, 1):
            print(f"\n[{i}/{len(chromosomes)}] ", end='')
            success = download_chromosome(chrom, storage_dir, dry_run=args.dry_run)
            if not success and not args.dry_run:
                print(f"\n⚠ Failed to download chr{chrom}, continuing with others...")
    
    # Create symlinks
    success_count = create_symlinks(storage_dir, work_dir, chromosomes, dry_run=args.dry_run)
    
    # Save manifest
    if not args.dry_run:
        save_manifest(storage_dir, work_dir, chromosomes)
    
    print("\n" + "="*80)
    print("SETUP COMPLETE")
    print("="*80)
    print(f"\nStorage: {storage_dir}")
    print(f"Symlinks: {work_dir}")
    print(f"Chromosomes: {success_count} configured")
    
    if args.dry_run:
        print("\n[DRY RUN] No files were downloaded or modified")
    
    print("\nNext: Update pipeline_config.yaml with:")
    print(f"  gnomad:")
    print(f"    vcf_dir: \"{work_dir}\"")
    print(f"    vcf_pattern: \"{work_dir}/gnomad.genomes.v4.1.sites.chr{{chr}}.vcf.bgz\"")


if __name__ == '__main__':
    main()
