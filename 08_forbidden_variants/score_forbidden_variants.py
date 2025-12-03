#!/usr/bin/env python3
"""
Module 08: Score Forbidden Variants with AlphaGenome (Streaming Version)

Score ~30,000 "forbidden variants" using a streaming/chunked approach to avoid
memory exhaustion. Each batch is written to disk immediately.

Memory-safe design:
- Process variants in batches of ~500
- Write each batch to parquet immediately
- Never hold more than ~30 GB in RAM
- Final output is partitioned parquet directory
"""

import argparse
import logging
import os
import sys
from pathlib import Path
import pandas as pd
import numpy as np
from tqdm import tqdm
import time
import gc

# Add utils to path
sys.path.insert(0, str(Path(__file__).parent.parent / "utils"))


def setup_logger(name, log_file=None, level="INFO"):
    """Setup logging."""
    logger = logging.getLogger(name)
    logger.setLevel(getattr(logging, level))
    
    # Clear existing handlers
    logger.handlers = []
    
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    
    # Console handler
    ch = logging.StreamHandler()
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    
    # File handler
    if log_file:
        fh = logging.FileHandler(log_file)
        fh.setFormatter(formatter)
        logger.addHandler(fh)
    
    return logger


def score_batch_and_save(batch_df, batch_idx, client, scorers, organism, output_dir, logger):
    """
    Score a batch of variants and save immediately to parquet.
    
    Returns summary statistics for the batch (not full predictions).
    """
    from alphagenome.data import genome
    from alphagenome.models import dna_client, variant_scorers
    
    batch_results = []
    batch_summaries = []
    failed = []
    
    for _, row in batch_df.iterrows():
        try:
            chrom = str(row['chr'])
            if not chrom.startswith('chr'):
                chrom = f"chr{chrom}"
            
            pos = int(row['pos'])
            ref = str(row['ref'])
            alt = str(row['alt'])
            
            variant_name = f"{chrom}:{pos}:{ref}>{alt}"
            
            variant = genome.Variant(
                chromosome=chrom,
                position=pos,
                reference_bases=ref,
                alternate_bases=alt,
                name=variant_name
            )
            
            interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
            
            variant_scores = client.score_variant(
                interval=interval,
                variant=variant,
                variant_scorers=scorers,
                organism=organism
            )
            
            batch_results.extend(variant_scores)
            
        except Exception as e:
            failed.append({
                'variant': f"{row['chr']}:{row['pos']}:{row['ref']}>{row['alt']}",
                'error': str(e)
            })
            continue
    
    if not batch_results:
        return [], failed
    
    # Convert batch to DataFrame
    df = variant_scorers.tidy_scores(batch_results)
    
    # Convert Variant objects to strings BEFORE saving to parquet
    # The variant_id column contains Variant objects which pyarrow can't serialize
    df['variant_id_str'] = df['variant_id'].astype(str)
    df['scored_interval'] = df['scored_interval'].astype(str)
    
    # Drop object columns that can't be serialized to parquet
    # (variant_id is a Variant object, scorer may be a scorer object)
    cols_to_drop = [c for c in ['variant_id', 'scorer'] if c in df.columns]
    if cols_to_drop:
        df = df.drop(columns=cols_to_drop)
    
    # Save batch to parquet
    batch_path = output_dir / "partitions" / f"batch_{batch_idx:04d}.parquet"
    batch_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(batch_path, index=False)
    
    # Compute summary for this batch (keeps only essential info)
    summary = df.groupby('variant_id_str').agg({
        'quantile_score': ['mean', 'max', 'std', 'count'],
        'raw_score': ['mean', 'max']
    }).reset_index()
    summary.columns = ['variant_id', 'quantile_mean', 'quantile_max', 'quantile_std', 
                       'n_tracks', 'raw_mean', 'raw_max']
    
    # Extract AP1-specific scores if available
    ap1_tfs = ['JUN', 'JUND', 'FOS', 'FOSL1', 'FOSL2', 'ATF3', 'ATF2', 'BATF', 'MAFK']
    if 'experiment_target' in df.columns:
        ap1_pattern = '|'.join(ap1_tfs)
        ap1_mask = df['experiment_target'].str.upper().str.contains(ap1_pattern, na=False)
        ap1_df = df[ap1_mask]
        
        if len(ap1_df) > 0:
            ap1_summary = ap1_df.groupby('variant_id_str').agg({
                'quantile_score': ['mean', 'max'],
                'raw_score': ['mean', 'max']
            }).reset_index()
            ap1_summary.columns = ['variant_id', 'ap1_quantile_mean', 'ap1_quantile_max', 
                                   'ap1_raw_mean', 'ap1_raw_max']
            summary = summary.merge(ap1_summary, on='variant_id', how='left')
    
    # Free memory
    del df, batch_results
    gc.collect()
    
    return summary.to_dict('records'), failed


def main():
    parser = argparse.ArgumentParser(description='Score forbidden variants (streaming)')
    parser.add_argument('--input', required=True, help='Path to forbidden_variants.tsv')
    parser.add_argument('--output', required=True, help='Output directory')
    parser.add_argument('--api-key', help='AlphaGenome API key (or set ALPHA_GENOME_KEY)')
    parser.add_argument('--batch-size', type=int, default=500, 
                        help='Variants per batch (default: 500, ~30GB RAM)')
    parser.add_argument('--limit', type=int, help='Limit variants (for testing)')
    parser.add_argument('--resume-from', type=int, default=0,
                        help='Resume from batch N (for crash recovery)')
    args = parser.parse_args()
    
    # Setup
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    log_dir = Path(args.output).parent.parent / "logs"
    log_dir.mkdir(exist_ok=True)
    
    logger = setup_logger(
        __name__,
        log_file=str(log_dir / "08_forbidden_variants_scoring.log")
    )
    
    logger.info("=" * 70)
    logger.info("Module 08: Score Forbidden Variants (STREAMING MODE)")
    logger.info("=" * 70)
    logger.info("Memory-safe: each batch saved immediately, ~30GB RAM max")
    
    # Get API key
    api_key = args.api_key or os.environ.get('ALPHA_GENOME_KEY')
    if not api_key:
        raise ValueError("Must provide --api-key or set ALPHA_GENOME_KEY")
    
    # Load forbidden variants
    logger.info(f"\n[1/5] Loading forbidden variants from {args.input}...")
    variants_df = pd.read_csv(args.input, sep='\t')
    logger.info(f"  Loaded {len(variants_df):,} forbidden variants")
    
    if args.limit:
        variants_df = variants_df.head(args.limit)
        logger.info(f"  Limited to {args.limit} variants for testing")
    
    # Calculate batches
    n_variants = len(variants_df)
    n_batches = (n_variants + args.batch_size - 1) // args.batch_size
    
    logger.info(f"\n[2/5] Batch configuration:")
    logger.info(f"  Total variants: {n_variants:,}")
    logger.info(f"  Batch size: {args.batch_size}")
    logger.info(f"  Total batches: {n_batches}")
    logger.info(f"  Est. RAM per batch: ~{args.batch_size * 28000 * 270 / 1e9:.1f} GB")
    
    if args.resume_from > 0:
        logger.info(f"  Resuming from batch {args.resume_from}")
    
    # Initialize AlphaGenome
    logger.info("\n[3/5] Initializing AlphaGenome client...")
    
    try:
        from alphagenome.data import genome
        from alphagenome.models import dna_client, variant_scorers
    except ImportError as e:
        logger.error("Failed to import AlphaGenome. Activate alphagenome-env.")
        raise
    
    client = dna_client.create(api_key=api_key)
    scorers = list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.values())
    organism = dna_client.Organism.HOMO_SAPIENS
    
    logger.info(f"  Using {len(scorers)} scorers")
    logger.info(f"  Context window: {dna_client.SEQUENCE_LENGTH_1MB:,} bp")
    
    # Score variants in batches
    logger.info(f"\n[4/5] Scoring variants in {n_batches} batches...")
    start_time = time.time()
    
    all_summaries = []
    all_failed = []
    
    for batch_idx in range(args.resume_from, n_batches):
        batch_start = batch_idx * args.batch_size
        batch_end = min((batch_idx + 1) * args.batch_size, n_variants)
        batch_df = variants_df.iloc[batch_start:batch_end]
        
        logger.info(f"\n  Batch {batch_idx + 1}/{n_batches}: variants {batch_start}-{batch_end}")
        
        batch_time = time.time()
        summaries, failed = score_batch_and_save(
            batch_df, batch_idx, client, scorers, organism, output_dir, logger
        )
        batch_elapsed = time.time() - batch_time
        
        all_summaries.extend(summaries)
        all_failed.extend(failed)
        
        # Progress update
        variants_done = batch_end
        elapsed = time.time() - start_time
        rate = variants_done / elapsed if elapsed > 0 else 0
        eta = (n_variants - variants_done) / rate if rate > 0 else 0
        
        logger.info(f"    Scored: {len(batch_df) - len(failed)}, Failed: {len(failed)}")
        logger.info(f"    Batch time: {batch_elapsed:.1f}s, Rate: {rate:.1f} var/s, ETA: {eta/3600:.1f}h")
        
        # Save running summary every 10 batches
        if (batch_idx + 1) % 10 == 0:
            checkpoint_df = pd.DataFrame(all_summaries)
            checkpoint_path = output_dir / "summary_checkpoint.tsv"
            checkpoint_df.to_csv(checkpoint_path, sep='\t', index=False)
            logger.info(f"    Checkpoint saved: {len(all_summaries)} variants")
    
    total_time = time.time() - start_time
    
    # Final outputs
    logger.info(f"\n[5/5] Finalizing outputs...")
    
    # Save final summary
    summary_df = pd.DataFrame(all_summaries)
    
    # Merge with original metadata
    variants_df['variant_id'] = variants_df.apply(
        lambda r: f"chr{r['chr']}:{r['pos']}:{r['ref']}>{r['alt']}" 
        if not str(r['chr']).startswith('chr') 
        else f"{r['chr']}:{r['pos']}:{r['ref']}>{r['alt']}", axis=1
    )
    summary_df = summary_df.merge(
        variants_df[['variant_id', 'AN', 'motif_start', 'motif_end', 'strand', 'tier', 'pwm_score']],
        left_on='variant_id', right_on='variant_id', how='left'
    )
    
    # Sort by predicted impact
    summary_df = summary_df.sort_values('raw_max', ascending=False)
    
    summary_path = output_dir / "predictions_summary.tsv"
    summary_df.to_csv(summary_path, sep='\t', index=False)
    logger.info(f"  Saved summary: {summary_path}")
    
    # Top candidates
    top_n = min(1000, len(summary_df))
    top_path = output_dir / "top_candidates.tsv"
    summary_df.head(top_n).to_csv(top_path, sep='\t', index=False)
    logger.info(f"  Saved top {top_n}: {top_path}")
    
    # Save failures
    if all_failed:
        failed_path = output_dir / "failed_variants.tsv"
        pd.DataFrame(all_failed).to_csv(failed_path, sep='\t', index=False)
        logger.info(f"  Saved {len(all_failed)} failures: {failed_path}")
    
    # List partition files
    partition_dir = output_dir / "partitions"
    partition_files = list(partition_dir.glob("*.parquet"))
    total_size = sum(f.stat().st_size for f in partition_files) / 1e9
    
    logger.info("\n" + "=" * 70)
    logger.info("COMPLETE")
    logger.info("=" * 70)
    logger.info(f"  Total time: {total_time/3600:.1f} hours")
    logger.info(f"  Variants scored: {len(summary_df):,}")
    logger.info(f"  Failed: {len(all_failed)}")
    logger.info(f"  Partition files: {len(partition_files)}")
    logger.info(f"  Total disk: {total_size:.2f} GB")
    logger.info(f"\nOutputs:")
    logger.info(f"  - {summary_path}")
    logger.info(f"  - {top_path}")
    logger.info(f"  - {partition_dir}/ (full predictions)")
    
    # Quick stats
    if len(summary_df) > 0:
        logger.info(f"\nTop 10 forbidden variants by predicted impact:")
        for i, row in summary_df.head(10).iterrows():
            ap1_max = row.get('ap1_raw_max', 'N/A')
            logger.info(f"  {row['variant_id']}: raw_max={row['raw_max']:.1f}, ap1_max={ap1_max}")


if __name__ == '__main__':
    main()
