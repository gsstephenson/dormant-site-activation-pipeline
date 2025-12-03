#!/usr/bin/env python3
"""
Module 08: Extract Forbidden Variants

Extract H=1 mutations that are absent from gnomAD despite high coverage.
These are the "forbidden variants" under strongest purifying selection.

Input:
- paths.tsv from Module 02 (all mutation paths)
- unique_variants.tsv from Module 03 (observed gnomAD variants)
- Coverage file with AN data

Output:
- forbidden_variants.tsv with columns for AlphaGenome scoring
"""

import argparse
import pandas as pd
import pysam
from pathlib import Path
from tqdm import tqdm
import multiprocessing as mp
from functools import partial


def query_coverage_batch(positions, coverage_file, min_an=100000):
    """Query coverage for a batch of positions."""
    try:
        tbx = pysam.TabixFile(coverage_file)
    except Exception as e:
        print(f"Error opening tabix file: {e}")
        return []
    
    results = []
    for chrom, pos, ref, alt in positions:
        an = 0
        try:
            chrom_query = chrom if chrom.startswith('chr') else f'chr{chrom}'
            for row in tbx.fetch(chrom_query, pos-1, pos):
                fields = row.split('\t')
                an = int(fields[2])
                break
        except:
            pass
        
        if an >= min_an:
            results.append((chrom, pos, ref, alt, an))
    
    tbx.close()
    return results


def main():
    parser = argparse.ArgumentParser(description='Extract forbidden variants')
    parser.add_argument('--paths', required=True, help='Path to paths.tsv from Module 02')
    parser.add_argument('--observed', required=True, help='Path to unique_variants.tsv from Module 03')
    parser.add_argument('--coverage', required=True, help='Path to gnomAD AN tabix file')
    parser.add_argument('--output', required=True, help='Output path for forbidden_variants.tsv')
    parser.add_argument('--min-an', type=int, default=100000, help='Minimum AN for high confidence')
    parser.add_argument('--workers', type=int, default=32, help='Number of parallel workers')
    args = parser.parse_args()
    
    print("=" * 60)
    print("Module 08: Extract Forbidden Variants")
    print("=" * 60)
    
    # Create output directory
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    
    # Load H=1 mutations from paths
    print("\n[1/4] Loading H=1 mutations from paths.tsv...")
    
    # Read with proper column parsing (low_memory=False silences dtype warning)
    paths_df = pd.read_csv(args.paths, sep='\t', low_memory=False)
    h1_paths = paths_df[paths_df['hamming_distance'] == 1].copy()
    
    print(f"  Total H=1 mutations: {len(h1_paths):,}")
    
    # Extract unique variants (chr, pos, ref, alt)
    h1_variants = h1_paths[['chr', 'genomic_position', 'ref_base', 'alt_base']].drop_duplicates()
    h1_variants.columns = ['chrom', 'pos', 'ref', 'alt']
    print(f"  Unique H=1 positions: {len(h1_variants):,}")
    
    # Load observed gnomAD variants
    print("\n[2/4] Loading observed gnomAD variants...")
    observed_df = pd.read_csv(args.observed, sep='\t')
    
    # Create set of observed variant keys
    observed_set = set()
    for _, row in observed_df.iterrows():
        chrom = str(row.get('chr', row.get('chrom', '')))
        pos = int(row.get('pos', row.get('genomic_position', 0)))
        ref = str(row.get('ref', row.get('ref_base', '')))
        alt = str(row.get('alt', row.get('alt_base', '')))
        observed_set.add((chrom, pos, ref, alt))
    
    print(f"  Observed variants: {len(observed_set):,}")
    
    # Filter out observed variants
    print("\n[3/4] Filtering to absent variants...")
    absent_variants = []
    for _, row in h1_variants.iterrows():
        key = (str(row['chrom']), int(row['pos']), row['ref'], row['alt'])
        if key not in observed_set:
            absent_variants.append(key)
    
    print(f"  Absent H=1 variants: {len(absent_variants):,}")
    
    # Query coverage to filter for high-confidence positions
    print(f"\n[4/4] Querying coverage for {len(absent_variants):,} positions...")
    print(f"  Using {args.workers} workers, min AN = {args.min_an:,}")
    
    # Split into batches
    batch_size = len(absent_variants) // args.workers + 1
    batches = [absent_variants[i:i+batch_size] for i in range(0, len(absent_variants), batch_size)]
    
    # Query in parallel
    query_func = partial(query_coverage_batch, coverage_file=args.coverage, min_an=args.min_an)
    
    high_conf_variants = []
    with mp.Pool(args.workers) as pool:
        for batch_results in tqdm(pool.imap(query_func, batches), total=len(batches), desc="Querying coverage"):
            high_conf_variants.extend(batch_results)
    
    print(f"  High-confidence forbidden variants: {len(high_conf_variants):,}")
    
    # Create output DataFrame
    forbidden_df = pd.DataFrame(high_conf_variants, columns=['chr', 'pos', 'ref', 'alt', 'AN'])
    
    # Add additional columns needed for AlphaGenome scoring
    forbidden_df['variant_id'] = forbidden_df.apply(
        lambda r: f"{r['chr']}:{r['pos']}:{r['ref']}>{r['alt']}", axis=1
    )
    
    # Merge with original paths data for motif context
    h1_paths_unique = h1_paths.drop_duplicates(subset=['chr', 'genomic_position', 'ref_base', 'alt_base'])
    h1_paths_unique = h1_paths_unique.rename(columns={
        'genomic_position': 'pos',
        'ref_base': 'ref', 
        'alt_base': 'alt'
    })
    
    forbidden_df = forbidden_df.merge(
        h1_paths_unique[['chr', 'pos', 'ref', 'alt', 'motif_start', 'motif_end', 'strand', 'tier', 'pwm_score', 'seq_before', 'seq_after']],
        on=['chr', 'pos', 'ref', 'alt'],
        how='left'
    )
    
    # Sort by AN (highest confidence first)
    forbidden_df = forbidden_df.sort_values('AN', ascending=False)
    
    # Save
    forbidden_df.to_csv(args.output, sep='\t', index=False)
    
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"  Total H=1 mutations:              {len(h1_variants):,}")
    print(f"  Observed in gnomAD:               {len(h1_variants) - len(absent_variants):,}")
    print(f"  Absent from gnomAD:               {len(absent_variants):,}")
    print(f"  High-confidence (AN >= {args.min_an:,}): {len(high_conf_variants):,}")
    print(f"\nOutput saved to: {args.output}")
    print(f"\nThese {len(high_conf_variants):,} forbidden variants are ready for AlphaGenome scoring.")


if __name__ == '__main__':
    main()
