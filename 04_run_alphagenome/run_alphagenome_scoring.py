#!/usr/bin/env python3
"""
Module 04: AlphaGenome Variant Scoring

Score functional impact of activating mutations using AlphaGenome.
Uses 1MB context windows following official AlphaGenome methodology.

Based on working implementation from alphagenome-qtl-validation repository.

Input:
- results/gnomad_intersection/{TF}/paths_with_gnomad.tsv
- AlphaGenome API key

Output:
- results/alphagenome/{TF}/predictions.parquet - Raw scores per variant-track
- results/alphagenome/{TF}/predictions_summary.tsv - Mean scores per variant
"""

import argparse
import logging
import sys
import os
from pathlib import Path
import pandas as pd
import yaml
from tqdm import tqdm

# Add utils to path
sys.path.insert(0, str(Path(__file__).parent.parent / "utils"))
from logging_utils import setup_logger


def load_config(config_path: Path) -> dict:
    """Load pipeline configuration."""
    with open(config_path) as f:
        return yaml.safe_load(f)


def score_variants_alphagenome(variants_df, api_key, output_dir, limit=None):
    """
    Score variants using AlphaGenome with 1MB context windows.
    
    Follows methodology from alphagenome-qtl-validation/scripts/02_predict.py
    
    Args:
        variants_df: DataFrame with columns [chrom, mutation_pos or start, ref, alt, ...]
        api_key: AlphaGenome API key
        output_dir: Path to output directory (for checkpoint saving)
        limit: Optional limit on number of variants to process (for testing)
        
    Returns:
        DataFrame with prediction results (variant-track level)
    """
    logger = logging.getLogger(__name__)
    
    # Import AlphaGenome modules
    try:
        from alphagenome.data import genome
        from alphagenome.models import dna_client, variant_scorers
    except ImportError as e:
        logger.error("Failed to import AlphaGenome modules. Ensure alphagenome package is installed.")
        logger.error("Activate the conda environment: conda activate alphagenome-env")
        raise
    
    # Initialize client
    logger.info("Initializing AlphaGenome client with API key...")
    client = dna_client.create(api_key=api_key)
    
    # Define scorers - USE ALL SCORERS to capture multi-modal features
    # This captures: expression (RNA_SEQ), chromatin accessibility (ATAC, DNASE),
    # TF binding (CHIP_TF), histone modifications (CHIP_HISTONE), 3D contacts (CONTACT_MAPS),
    # transcript features (CAGE, PROCAP, SPLICE_*), and more
    # Total: 19 scorers × ~5K tracks each = ~89K tracks per variant
    scorers = list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.values())
    organism = dna_client.Organism.HOMO_SAPIENS
    
    logger.info(f"Scoring {len(variants_df)} variants with 1MB context windows")
    logger.info(f"Context window size: {dna_client.SEQUENCE_LENGTH_1MB:,} bp")
    logger.info(f"Using {len(scorers)} scorers to capture multi-modal features:")
    logger.info(f"  Output types: RNA_SEQ, ATAC, DNASE, CHIP_TF, CHIP_HISTONE, CONTACT_MAPS, CAGE, PROCAP, SPLICE_*")
    logger.info(f"  Expected ~89K tracks per variant (vs 1.1K with single scorer)")
    
    if limit:
        variants_df = variants_df.head(limit)
        logger.info(f"Limited to {limit} variants for testing")
    
    results = []
    failed_variants = []
    variant_metadata = {}  # Store metadata keyed by variant_name
    
    # Score each variant
    for idx, row in tqdm(variants_df.iterrows(), total=len(variants_df), desc="Scoring variants"):
        chrom = None
        pos = None
        try:
            # Get genomic coordinates (using actual column names from Module 03 output)
            chrom = str(row['chr'])
            if not chrom.startswith('chr'):
                chrom = f"chr{chrom}"
            
            # Use the mutation path columns (ref_base, alt_base, genomic_position)
            pos = int(row['genomic_position'])
            ref = str(row['ref_base'])
            alt = str(row['alt_base'])
            
            # Create variant ID matching AlphaGenome's Variant.__str__ format: chr:pos:ref>alt
            variant_name = f"{chrom}:{pos}:{ref}>{alt}"
            
            # Create Variant object
            variant = genome.Variant(
                chromosome=chrom,
                position=pos,
                reference_bases=ref,
                alternate_bases=alt,
                name=variant_name
            )
            
            # Get 1MB interval around variant (same as qtl-validation implementation)
            interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
            
            # Score variant (returns list of VariantScores)
            variant_scores = client.score_variant(
                interval=interval,
                variant=variant,
                variant_scorers=scorers,
                organism=organism
            )
            
            # Extend results with all scores from this variant
            results.extend(variant_scores)
            
            # Store metadata keyed by variant_id (which matches variant.name)
            variant_metadata[variant_name] = {
                'gnomad_AF': row.get('AF'),
                'gnomad_AC': row.get('AC'),
                'gnomad_AN': row.get('AN'),
                'path_id': row.get('path_id'),
                'step_num': row.get('step_num'),
                'site_id': row.get('site_id'),
                'tier': row.get('tier')
            }
            
        except Exception as e:
            logger.error(f"Error scoring variant at {chrom}:{pos}: {e}")
            failed_variants.append({
                'chrom': chrom,
                'pos': pos,
                'ref': ref,
                'alt': alt,
                'error': str(e)
            })
            continue
    
    if not results:
        logger.error("No variants were successfully scored!")
        return None
    
    # Convert to tidy dataframe (one row per variant-track combination)
    logger.info("Converting results to dataframe...")
    df = variant_scorers.tidy_scores(results)
    
    logger.info(f"Tidy scores shape: {df.shape}")
    logger.info(f"Tidy scores columns: {df.columns.tolist()}")
    
    # Check variant_id types before conversion
    logger.info(f"Sample variant_id object: {repr(df['variant_id'].iloc[0])}")
    # NOTE: Cannot call .nunique() on Variant objects (unhashable) - convert first!
    
    # Convert object columns to strings (following qtl-validation pattern)
    # The Variant object's __str__ method returns the proper format: chr:pos:ref>alt
    df['variant_id_str'] = df['variant_id'].astype(str)
    
    # Also convert scored_interval to string
    df['scored_interval'] = df['scored_interval'].astype(str)
    
    # Check conversion worked
    logger.info(f"Sample variant_id_str after conversion: {df['variant_id_str'].iloc[0]}")
    logger.info(f"Unique variants after conversion: {df['variant_id_str'].nunique()}")
    
    # Create a metadata DataFrame from the dict
    # The keys in variant_metadata are the variant_name strings we created (chr_pos_ref_alt format)
    # The variant_id_str from AlphaGenome will be in format "chr:pos:ref>alt"
    # We need to convert our metadata keys to match
    metadata_df = pd.DataFrame([
        {'variant_id_str': vid, **meta} 
        for vid, meta in variant_metadata.items()
    ])
    
    logger.info(f"Metadata entries: {len(metadata_df)}")
    logger.info(f"Sample metadata variant_id: {metadata_df['variant_id_str'].iloc[0] if len(metadata_df) > 0 else 'N/A'}")
    
    # SKIP CHECKPOINT: System doesn't have enough RAM to save 199M rows
    # Instead, proceed directly to merge which is more memory-efficient
    # If merge fails, user will need to re-run API calls (~2.5 hours)
    logger.info(f"Skipping checkpoint save (dataset too large: {len(df):,} rows)")
    logger.info(f"Proceeding directly to memory-efficient merge...")
    
    # Merge metadata into scores using memory-efficient approach
    logger.info("Merging metadata (memory-efficient left join)...")
    df = df.merge(metadata_df, on='variant_id_str', how='left', copy=False)
    
    # Clean up metadata_df to free memory before continuing
    del metadata_df
    import gc
    gc.collect()
    logger.info("Metadata merged and temporary objects cleaned up")
    
    logger.info(f"Successfully scored {len(df)} variant-track pairs")
    logger.info(f"Unique variants: {df['variant_id_str'].nunique()}")
    logger.info(f"Failed variants: {len(failed_variants)}")
    logger.info(f"Rows with metadata: {df['gnomad_AF'].notna().sum()}/{len(df)}")
    
    if failed_variants:
        logger.warning(f"First 5 failures: {failed_variants[:5]}")
    
    return df


def compute_summary_scores(predictions_df):
    """
    Compute summary statistics per variant and per variant-output_type.
    
    With multi-modal AlphaGenome output, we get predictions for:
    - RNA_SEQ (gene expression)
    - ATAC, DNASE (chromatin accessibility)
    - CHIP_TF (TF binding)
    - CHIP_HISTONE (histone modifications)
    - CONTACT_MAPS (3D genome contacts)
    - CAGE, PROCAP (transcription start sites)
    - SPLICE_* (splicing)
    
    Args:
        predictions_df: Tidy scores from AlphaGenome (variant-track level)
        
    Returns:
        tuple: (per_variant_summary, per_variant_output_type_summary)
    """
    logger = logging.getLogger(__name__)
    
    logger.info("Computing summary scores per variant...")
    
    # Overall summary: mean across ALL tracks
    summary = predictions_df.groupby('variant_id_str').agg({
        'quantile_score': 'mean',  # Mean across all tracks
        'gnomad_AF': 'first',
        'gnomad_AC': 'first',
        'gnomad_AN': 'first',
        'path_id': 'first',
        'step_num': 'first'
    }).reset_index()
    
    # Also get count of tracks contributing
    track_counts = predictions_df.groupby('variant_id_str').size().rename('n_tracks')
    summary = summary.merge(track_counts, on='variant_id_str')
    
    logger.info(f"Summary: {len(summary)} variants with mean scores across all tracks")
    
    # Feature-specific summary: mean per output_type (enables feature-specific analysis)
    if 'output_type' in predictions_df.columns:
        logger.info("Computing per-output_type summaries...")
        feature_summary = predictions_df.groupby(['variant_id_str', 'output_type']).agg({
            'quantile_score': 'mean',
            'gnomad_AF': 'first',
            'gnomad_AC': 'first',
            'gnomad_AN': 'first'
        }).reset_index()
        
        # Count tracks per feature
        feature_track_counts = predictions_df.groupby(['variant_id_str', 'output_type']).size().rename('n_tracks')
        feature_summary = feature_summary.merge(
            feature_track_counts.reset_index(),
            on=['variant_id_str', 'output_type']
        )
        
        logger.info(f"Feature-specific summary: {len(feature_summary)} variant-output_type pairs")
        logger.info(f"Output types: {sorted(feature_summary['output_type'].unique())}")
        
        return summary, feature_summary
    else:
        logger.warning("No output_type column found - skipping feature-specific summary")
        return summary, None


def main():
    parser = argparse.ArgumentParser(
        description="Score variants with AlphaGenome using 1MB context windows"
    )
    parser.add_argument('--config', type=Path, required=True,
                       help='Pipeline config YAML')
    parser.add_argument('--api-key', type=str,
                       help='AlphaGenome API key (or set ALPHA_GENOME_KEY env var)')
    parser.add_argument('--limit', type=int,
                       help='Limit to first N variants (for testing)')
    parser.add_argument('--clean', action='store_true',
                       help='Remove checkpoint files before starting (forces re-scoring)')
    args = parser.parse_args()
    
    # Get API key
    api_key = args.api_key or os.environ.get('ALPHA_GENOME_KEY')
    if not api_key:
        raise ValueError("Must provide --api-key or set ALPHA_GENOME_KEY environment variable")
    
    # Load config
    config = load_config(args.config)
    tf_name = config['tf_name']
    
    # Setup paths
    base_dir = Path(__file__).parent.parent
    gnomad_dir = base_dir / "results" / "gnomad_intersection" / tf_name
    output_dir = base_dir / "results" / "alphagenome" / tf_name
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Setup logging
    log_dir = base_dir / "logs"
    log_dir.mkdir(exist_ok=True)
    logger = setup_logger(
        __name__,
        log_file=str(log_dir / f"04_alphagenome_scoring_{tf_name}.log"),
        level="INFO"
    )
    
    logging.info("="*80)
    logging.info("Starting: Module 04: AlphaGenome Variant Scoring")
    logging.info("="*80)
    logging.info(f"Parameters:")
    logging.info(f"  TF: {tf_name}")
    logging.info(f"  Input: {gnomad_dir / 'paths_with_gnomad.tsv'}")
    logging.info(f"  Output: {output_dir}")
    if args.limit:
        logging.info(f"  Limit: {args.limit} variants (testing mode)")
    
    # Load variants - use unique variants file if available, otherwise full paths
    unique_variants_path = gnomad_dir / "unique_variants.tsv"
    paths_variants_path = gnomad_dir / "paths_with_gnomad.tsv"
    
    if unique_variants_path.exists():
        variants_path = unique_variants_path
        logging.info(f"Loading unique variants from {variants_path}...")
    else:
        variants_path = paths_variants_path
        logging.info(f"Loading all mutation steps from {variants_path}...")
        logging.warning("Consider running prepare_unique_variants.py to deduplicate variants")
    
    if not variants_path.exists():
        raise FileNotFoundError(f"Input not found: {variants_path}")
    
    variants_df = pd.read_csv(variants_path, sep='\t', low_memory=False)
    logging.info(f"Loaded {len(variants_df)} variants")
    
    # Check if we have a raw predictions checkpoint (from previous crash during merge)
    raw_checkpoint = output_dir / 'predictions_raw.parquet'
    
    # Clean checkpoint if requested
    if args.clean and raw_checkpoint.exists():
        checkpoint_size = raw_checkpoint.stat().st_size / 1e9
        logging.info(f"--clean flag: Removing checkpoint file: {raw_checkpoint} ({checkpoint_size:.2f} GB)")
        raw_checkpoint.unlink()
        logging.info("Checkpoint removed - will re-run AlphaGenome scoring")
    
    if raw_checkpoint.exists():
        checkpoint_size = raw_checkpoint.stat().st_size / 1e9
        logging.warning(f"Found raw predictions checkpoint: {raw_checkpoint} ({checkpoint_size:.2f} GB)")
        logging.warning("This suggests a previous run crashed during metadata merge")
        logging.info("Loading checkpoint instead of re-running AlphaGenome scoring...")
        logging.info("(Use --clean flag to force re-scoring)")
        predictions_df = pd.read_parquet(raw_checkpoint)
        logging.info(f"Loaded {len(predictions_df):,} variant-track pairs from checkpoint")
    else:
        # Score variants normally
        predictions_df = score_variants_alphagenome(
            variants_df=variants_df,
            api_key=api_key,
            output_dir=output_dir,
            limit=args.limit
        )
    
    if predictions_df is None or len(predictions_df) == 0:
        logging.error("No predictions generated. Exiting.")
        return 1
    
    # Drop object columns (Variant and Interval objects) before saving
    # We've already converted them to variant_id_str and scored_interval strings
    columns_to_drop = ['variant_id', 'scored_interval']
    predictions_df = predictions_df.drop(columns=columns_to_drop, errors='ignore')
    
    # Save raw predictions
    predictions_path = output_dir / "predictions.parquet"
    logging.info(f"Saving final predictions to {predictions_path}...")
    predictions_df.to_parquet(predictions_path, index=False, compression='snappy')
    file_size_gb = predictions_path.stat().st_size / 1e9
    logging.info(f"Saved raw predictions: {predictions_path}")
    logging.info(f"  {len(predictions_df):,} variant-track pairs")
    logging.info(f"  File size: {file_size_gb:.2f} GB")
    logging.info(f"  Unique variants: {predictions_df['variant_id_str'].nunique():,}")
    
    # Compute and save summary scores (both overall and per-feature)
    summary_result = compute_summary_scores(predictions_df)
    
    # Handle both single and tuple returns for backwards compatibility
    if isinstance(summary_result, tuple):
        summary_df, feature_summary_df = summary_result
    else:
        summary_df = summary_result
        feature_summary_df = None
    
    # Save overall summary
    summary_path = output_dir / "predictions_summary.tsv"
    summary_df.to_csv(summary_path, sep='\t', index=False)
    logging.info(f"Saved overall summary scores: {summary_path}")
    logging.info(f"  {len(summary_df)} variants")
    
    # Save feature-specific summary
    if feature_summary_df is not None:
        feature_summary_path = output_dir / "predictions_summary_by_feature.tsv"
        feature_summary_df.to_csv(feature_summary_path, sep='\t', index=False)
        logging.info(f"Saved feature-specific summary: {feature_summary_path}")
        logging.info(f"  {len(feature_summary_df)} variant-output_type pairs")
        
        # Log feature distribution
        feature_counts = feature_summary_df.groupby('output_type').size().sort_values(ascending=False)
        logging.info(f"  Feature distribution:")
        for feature, count in feature_counts.head(10).items():
            logging.info(f"    {feature}: {count} variants")
    
    # Clean up checkpoint file on successful completion
    raw_checkpoint = output_dir / 'predictions_raw.parquet'
    if raw_checkpoint.exists():
        checkpoint_size_gb = raw_checkpoint.stat().st_size / 1e9
        logging.info(f"Removing checkpoint file (no longer needed): {raw_checkpoint} ({checkpoint_size_gb:.2f} GB)")
        raw_checkpoint.unlink()
        logging.info("Checkpoint removed successfully")
    else:
        logging.info("No checkpoint file to clean up")
    
    logging.info("="*80)
    logging.info("Module 04 complete!")
    logging.info("="*80)
    logging.info(f"Next step: Module 05 - Compute activation landscape")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
