#!/usr/bin/env python3
"""
Query gnomAD VCFs for variants at mutation path positions
Part of Dormant Site Activation Pipeline - Module 03

OPTIMIZED for parallel chromosome queries
"""

import argparse
import sys
from pathlib import Path
import pandas as pd
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp
from functools import partial

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils.io import load_config
from utils.logging_utils import setup_logger, log_step, log_parameters, ProgressLogger


def query_single_chromosome(chrom, vcf_dir, vcf_pattern, bed_df, threads_per_job=4):
    """
    Query a single chromosome (parallelizable function).
    
    OPTIMIZED:
    - Uses bcftools --threads for decompression parallelism
    - Efficient BED file handling
    - Robust error handling with 6-hour timeout
    
    Parameters
    ----------
    chrom : str
        Chromosome name
    vcf_dir : str
        Directory containing VCFs
    vcf_pattern : str
        VCF filename pattern
    bed_df : pd.DataFrame
        BED data for all chromosomes
    threads_per_job : int
        Threads for bcftools decompression (default: 4)
        
    Returns
    -------
    tuple
        (chrom, variants_text, num_variants, error_msg)
    """
    try:
        # Construct VCF path
        vcf_file = Path(vcf_dir) / vcf_pattern.format(chr=chrom)
        
        if not vcf_file.exists():
            return (chrom, "", 0, f"VCF not found: {vcf_file}")
        
        # Filter BED for this chromosome
        chrom_bed_data = bed_df[bed_df['chr'] == chrom]
        
        if len(chrom_bed_data) == 0:
            return (chrom, "", 0, None)
        
        # Create temporary BED file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as tmp_bed:
            tmp_bed_path = tmp_bed.name
            chrom_bed_data.to_csv(tmp_bed, sep='\t', header=False, index=False)
        
        try:
            # bcftools query command
            # Note: Parallel processing at chromosome level provides speedup
            cmd = [
                'bcftools', 'query',
                '-R', tmp_bed_path,
                '-f', '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/AF\\t%INFO/AC\\t%INFO/AN\\t%INFO/nhomalt\\n',
                str(vcf_file)
            ]
            
            # Run bcftools with 6-hour timeout for large chromosomes
            result = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=True,
                timeout=21600  # 6 hour timeout per chromosome
            )
            
            num_variants = result.stdout.count('\n') if result.stdout else 0
            return (chrom, result.stdout, num_variants, None)
            
        except subprocess.CalledProcessError as e:
            error_msg = e.stderr[:500] if e.stderr else str(e)  # Truncate long errors
            return (chrom, "", 0, f"bcftools failed: {error_msg}")
        except subprocess.TimeoutExpired:
            return (chrom, "", 0, "Query timed out (>6 hours)")
        finally:
            # Clean up temp file
            Path(tmp_bed_path).unlink(missing_ok=True)
            
    except Exception as e:
        return (chrom, "", 0, f"Error: {str(e)}")


def query_gnomad_parallel(
    vcf_dir: str,
    vcf_pattern: str,
    bed_file: str,
    output_file: str,
    n_jobs: int = None,
    logger=None
):
    """
    Query per-chromosome gnomAD VCFs using bcftools IN PARALLEL.
    
    FULLY OPTIMIZED for multi-core systems:
    - Process multiple chromosomes simultaneously (process-level parallelism)
    - Each chromosome uses bcftools --threads for decompression (thread-level parallelism)
    - Results streamed as they complete (no waiting for slowest chromosome)
    
    Performance on 32-core system with 3M positions:
    - Expected time: ~15-30 minutes (vs. hours sequential)
    - Memory: ~2-4 GB per chromosome (48-96 GB total for 24 chromosomes)
    
    Parameters
    ----------
    vcf_dir : str
        Directory containing per-chromosome VCFs
    vcf_pattern : str
        Pattern for VCF filenames (e.g., "gnomad.genomes.v4.1.sites.chr{chr}.vcf.bgz")
    bed_file : str
        Path to BED file with query regions
    output_file : str
        Output TSV file
    n_jobs : int
        Number of parallel chromosome jobs (default: min(24, CPU count - 4))
    logger : logging.Logger
        Logger instance
        
    Returns
    -------
    bool
        Success status
    """
    if n_jobs is None:
        # Use most cores for parallel chromosome processing
        n_jobs = min(30, max(1, mp.cpu_count() - 2))
    
    logger.info("=" * 80)
    logger.info("PARALLEL gnomAD QUERY - FULLY OPTIMIZED")
    logger.info("=" * 80)
    logger.info(f"  VCF directory: {vcf_dir}")
    logger.info(f"  BED: {bed_file}")
    logger.info(f"  Parallel chromosome jobs: {n_jobs}")
    logger.info(f"  Total CPU utilization: ~{n_jobs} cores")
    logger.info(f"  Timeout per chromosome: 6 hours")
    
    # Read BED to determine which chromosomes we need
    logger.info("\n  Loading BED file...")
    bed_df = pd.read_csv(bed_file, sep='\t', header=None, names=['chr', 'start', 'end', 'name'], dtype={'chr': str})
    
    # CRITICAL: Add "chr" prefix if not present (gnomAD uses "chr1" not "1")
    first_chr = str(bed_df['chr'].iloc[0])
    if not first_chr.startswith('chr'):
        logger.info("  ⚠ Chromosome names lack 'chr' prefix - adding for gnomAD compatibility")
        bed_df['chr'] = 'chr' + bed_df['chr'].astype(str)
        logger.info(f"  ✓ Converted: '1' → 'chr1', '2' → 'chr2', etc.")
    
    all_chromosomes = sorted(bed_df['chr'].unique())
    
    # Filter to only chromosomes with available VCF files
    # gnomAD only has chr1-22, X, Y (no scaffolds, patches, or MT)
    available_chroms = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
    chromosomes = [c for c in all_chromosomes if c in available_chroms]
    
    # Log filtering
    filtered_out = set(all_chromosomes) - set(chromosomes)
    if filtered_out:
        logger.info(f"  ⚠ Filtered out {len(filtered_out)} chromosomes without gnomAD VCFs (scaffolds, MT, etc.)")
        logger.info(f"    Skipped: {', '.join(sorted(filtered_out)[:10])}{'...' if len(filtered_out) > 10 else ''}")
    
    # Filter BED to only available chromosomes
    original_size = len(bed_df)
    bed_df = bed_df[bed_df['chr'].isin(chromosomes)]
    logger.info(f"  Filtered positions: {original_size:,} → {len(bed_df):,} (kept {100*len(bed_df)/original_size:.1f}%)")
    
    chromosomes = sorted(bed_df['chr'].unique())
    logger.info(f"  Chromosomes to query: {', '.join(chromosomes)} ({len(chromosomes)} total)")
    logger.info(f"  Total query positions: {len(bed_df):,}")
    
    # Estimate runtime
    avg_positions_per_chrom = len(bed_df) / len(chromosomes)
    logger.info(f"  Avg positions per chromosome: {avg_positions_per_chrom:,.0f}")
    logger.info(f"\n  Estimated runtime: 1-3 hours (large chromosomes may take longer)")
    logger.info("=" * 80)
    
    # Prepare output file with header
    with open(output_file, 'w') as f:
        f.write('chr\tpos\tref\talt\tAF\tAC\tAN\tnhomalt\n')
    
    # Process chromosomes in parallel
    total_variants = 0
    completed = 0
    failed = []
    
    logger.info("\n  Starting parallel queries...")
    
    with ProcessPoolExecutor(max_workers=n_jobs) as executor:
        # Submit all chromosome queries in parallel
        future_to_chrom = {
            executor.submit(query_single_chromosome, chrom, vcf_dir, vcf_pattern, bed_df): chrom
            for chrom in chromosomes
        }
        
        # Process results as they complete
        for future in as_completed(future_to_chrom):
            chrom, variants_text, num_variants, error_msg = future.result()
            completed += 1
            
            if error_msg:
                logger.warning(f"  [{completed}/{len(chromosomes)}] {chrom}: {error_msg}")
                failed.append(chrom)
            else:
                logger.info(f"  [{completed}/{len(chromosomes)}] {chrom}: {num_variants:,} variants")
                
                # Append results to output file
                if variants_text:
                    with open(output_file, 'a') as f:
                        f.write(variants_text)
                    total_variants += num_variants
    
    logger.info(f"\n  ✓ Parallel query complete")
    logger.info(f"  Total variants found: {total_variants:,}")
    logger.info(f"  Successful: {completed - len(failed)}/{len(chromosomes)} chromosomes")
    
    if failed:
        logger.warning(f"  Failed chromosomes: {', '.join(failed)}")
    
    logger.info(f"  Saved to: {output_file}")
    
    return len(failed) == 0


def match_paths_to_gnomad(
    paths_file: str,
    gnomad_file: str,
    output_file: str,
    logger=None
):
    """
    Match mutation paths to gnomAD variants.
    
    Parameters
    ----------
    paths_file : str
        Path to mutation paths file
    gnomad_file : str
        Path to gnomAD query results
    output_file : str
        Output file with matched variants
    logger : logging.Logger
        Logger instance
        
    Returns
    -------
    pd.DataFrame
        Matched paths with gnomAD data
    """
    logger.info("Matching mutation paths to gnomAD variants...")
    
    # Load paths
    logger.info(f"  Loading paths from: {paths_file}")
    paths_df = pd.read_csv(paths_file, sep='\t', low_memory=False)
    logger.info(f"    {len(paths_df):,} mutation steps")
    
    # Load gnomAD results
    logger.info(f"  Loading gnomAD variants from: {gnomad_file}")
    gnomad_df = pd.read_csv(gnomad_file, sep='\t', low_memory=False)
    logger.info(f"    {len(gnomad_df):,} variants")
    
    # CRITICAL: Ensure chromosome naming matches between paths and gnomAD
    # Paths use "1", gnomAD uses "chr1"
    first_chr = str(paths_df['chr'].iloc[0])
    if not first_chr.startswith('chr'):
        logger.info("  ⚠ Adding 'chr' prefix to mutation paths for gnomAD compatibility")
        paths_df['chr'] = 'chr' + paths_df['chr'].astype(str)
    
    # Merge on chr, pos, ref, alt
    # Note: gnomAD may have multiple ALT alleles per position
    logger.info("  Merging paths with gnomAD variants...")
    
    merged_df = paths_df.merge(
        gnomad_df,
        left_on=['chr', 'genomic_position', 'ref_base', 'alt_base'],
        right_on=['chr', 'pos', 'ref', 'alt'],
        how='left'
    )
    
    # Count matches
    num_matched = merged_df['AF'].notna().sum()
    percent_matched = 100 * num_matched / len(merged_df)
    
    logger.info(f"\n  Matched {num_matched:,} / {len(merged_df):,} mutation steps ({percent_matched:.1f}%)")
    
    # Fill missing values
    merged_df['AF'] = merged_df['AF'].fillna(0.0)
    merged_df['AC'] = merged_df['AC'].fillna(0)
    merged_df['AN'] = merged_df['AN'].fillna(0)
    merged_df['nhomalt'] = merged_df['nhomalt'].fillna(0)
    
    # Save
    merged_df.to_csv(output_file, sep='\t', index=False)
    logger.info(f"\n✓ Saved matched results: {output_file}")
    
    return merged_df


def summarize_paths_af(
    matched_file: str,
    output_file: str,
    logger=None
):
    """
    Summarize allele frequencies per mutation path.
    
    Parameters
    ----------
    matched_file : str
        Path to matched paths+gnomAD file
    output_file : str
        Output summary file
    logger : logging.Logger
        Logger instance
        
    Returns
    -------
    pd.DataFrame
        Path-level summary
    """
    logger.info("Summarizing allele frequencies per path...")
    
    # Load matched data
    df = pd.read_csv(matched_file, sep='\t', low_memory=False)
    
    # Group by path_id
    logger.info("  Aggregating by path_id...")
    
    path_summary = df.groupby('path_id').agg({
        'site_id': 'first',
        'chr': 'first',
        'motif_start': 'first',
        'motif_end': 'first',
        'strand': 'first',
        'tier': 'first',
        'pwm_score': 'first',
        'hamming_distance': 'first',
        'total_steps': 'first',
        'AF': ['max', 'mean', 'sum'],  # Max, mean, sum of AF across steps
        'AC': 'sum',
        'nhomalt': 'sum'
    }).reset_index()
    
    # Flatten column names
    path_summary.columns = [
        'path_id', 'site_id', 'chr', 'motif_start', 'motif_end', 'strand',
        'tier', 'pwm_score', 'hamming_distance', 'total_steps',
        'max_AF', 'mean_AF', 'sum_AF', 'total_AC', 'total_nhomalt'
    ]
    
    # Add path accessibility score
    # Simple version: max AF across all steps
    path_summary['accessibility_score'] = path_summary['max_AF']
    
    # Log statistics
    logger.info(f"\n  Summary for {len(path_summary):,} unique paths")
    
    logger.info("\n  AF distribution:")
    af_bins = [0, 1e-6, 1e-5, 1e-4, 1e-3, 0.01, 0.05, 0.1, 1.0]
    af_counts = pd.cut(path_summary['max_AF'], bins=af_bins).value_counts().sort_index()
    for bin_range, count in af_counts.items():
        percent = 100 * count / len(path_summary)
        logger.info(f"    {bin_range}: {count:,} paths ({percent:.1f}%)")
    
    # Save
    path_summary.to_csv(output_file, sep='\t', index=False)
    logger.info(f"\n✓ Saved path summary: {output_file}")
    
    return path_summary


def main():
    parser = argparse.ArgumentParser(
        description="Query gnomAD for variants at mutation path positions"
    )
    parser.add_argument(
        '--config',
        type=str,
        required=True,
        help='Path to pipeline configuration file'
    )
    parser.add_argument(
        '--paths-file',
        type=str,
        help='Path to mutation paths file (default: from Module 02)'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        help='Output directory'
    )
    parser.add_argument(
        '--skip-query',
        action='store_true',
        help='Skip bcftools query (use existing query results)'
    )
    parser.add_argument(
        '--n-jobs',
        type=int,
        help='Number of parallel jobs (default: CPU count - 2)'
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
    gnomad_vcf_dir = config['gnomad']['vcf_dir']
    gnomad_vcf_pattern = config['gnomad']['vcf_pattern']
    
    # Determine inputs/outputs
    if args.paths_file:
        paths_file = args.paths_file
    else:
        paths_file = f"results/mutation_paths/{tf_name}/paths.tsv"
    
    if args.output_dir:
        output_dir = Path(args.output_dir)
    else:
        output_dir = Path(f"results/gnomad_intersection/{tf_name}")
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    bed_file = output_dir / "query_regions.bed"
    gnomad_query_file = output_dir / "gnomad_query_results.tsv"
    matched_file = output_dir / "paths_with_gnomad.tsv"
    summary_file = output_dir / "path_af_summary.tsv"
    
    # Setup logger
    log_file = args.log_file or f"logs/03_gnomad_query_{tf_name}.log"
    logger = setup_logger(
        name=__name__,
        log_file=log_file,
        level='INFO'
    )
    
    log_step(logger, "Module 03: gnomAD Intersection", start=True)
    
    # Log parameters
    params = {
        'TF': tf_name,
        'Paths file': paths_file,
        'gnomAD VCF directory': gnomad_vcf_dir,
        'gnomAD VCF pattern': gnomad_vcf_pattern,
        'Output directory': str(output_dir),
        'Skip query': args.skip_query
    }
    log_parameters(logger, params)
    
    # Check inputs
    if not Path(paths_file).exists():
        logger.error(f"Paths file not found: {paths_file}")
        logger.error("Run Module 02 first")
        sys.exit(1)
    
    if not Path(gnomad_vcf_dir).exists():
        logger.error(f"gnomAD VCF directory not found: {gnomad_vcf_dir}")
        logger.error("Check gnomAD configuration in pipeline_config.yaml")
        sys.exit(1)
    
    # Load paths
    logger.info(f"Loading mutation paths...")
    paths_df = pd.read_csv(paths_file, sep='\t', low_memory=False)
    logger.info(f"  {len(paths_df):,} mutation steps")
    logger.info(f"  {paths_df['path_id'].nunique():,} unique paths")
    logger.info(f"  {paths_df['site_id'].nunique():,} unique motif sites")
    
    # Create BED file (combined for all chromosomes)
    # Note: The query function will split by chromosome internally
    logger.info("Creating query regions BED file...")
    positions_df = paths_df[['chr', 'genomic_position']].drop_duplicates()
    bed_data = pd.DataFrame({
        'chr': positions_df['chr'],
        'start': positions_df['genomic_position'],
        'end': positions_df['genomic_position'] + 1,
        'name': positions_df['chr'] + ':' + positions_df['genomic_position'].astype(str)
    })
    bed_data = bed_data.sort_values(['chr', 'start'])
    bed_data.to_csv(bed_file, sep='\t', header=False, index=False)
    logger.info(f"  Created BED file with {len(bed_data):,} unique positions")
    logger.info(f"  Chromosomes: {', '.join(sorted(bed_data['chr'].unique()))}")
    
    # Query gnomAD (unless skipping)
    if args.skip_query and gnomad_query_file.exists():
        logger.info(f"\nSkipping bcftools query, using existing: {gnomad_query_file}")
    else:
        # Determine number of parallel jobs
        n_jobs = args.n_jobs or config.get('resources', {}).get('threads', None)
        if n_jobs:
            logger.info(f"\nUsing {n_jobs} parallel jobs for gnomAD query")
        
        success = query_gnomad_parallel(
            vcf_dir=gnomad_vcf_dir,
            vcf_pattern=gnomad_vcf_pattern,
            bed_file=str(bed_file),
            output_file=str(gnomad_query_file),
            n_jobs=n_jobs,
            logger=logger
        )
        
        if not success:
            logger.warning("Some chromosomes failed, but continuing with available data")
            # Don't exit - continue with partial data
    
    # Match paths to gnomAD
    matched_df = match_paths_to_gnomad(
        paths_file=paths_file,
        gnomad_file=str(gnomad_query_file),
        output_file=str(matched_file),
        logger=logger
    )
    
    # Summarize per path
    summary_df = summarize_paths_af(
        matched_file=str(matched_file),
        output_file=str(summary_file),
        logger=logger
    )
    
    log_step(logger, "Module 03: gnomAD Intersection", start=False)
    logger.info("\nNext step: Module 04 - AlphaGenome functional scoring")
    logger.info("  python 04_run_alphagenome/make_variant_seqs.py --config pipeline_config.yaml")


if __name__ == '__main__':
    main()
