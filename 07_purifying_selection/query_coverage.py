#!/usr/bin/env python3
"""
Module 07: Query Coverage and Compute Constraint Statistics

This script queries the tabix-indexed gnomAD coverage file for all dormant
site positions and computes purifying selection statistics.

Prerequisites:
  - Run the bash script first to create the reformatted + indexed coverage file
  - Or ensure gnomad_AN_tabix.tsv.bgz and .tbi exist

Usage:
  conda activate alphagenome-env
  python query_coverage.py

Author: George Stephenson
Date: December 2025
"""

import pysam
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict
import time
from scipy import stats
import sys
import multiprocessing as mp
from functools import partial
import os

# Configuration
COVERAGE_FILE = "/mnt/data_1/gnomAD_data/raw/gnomad_v4.1/coverage/gnomad_AN_tabix.tsv.bgz"
PATHS_FILE = "/mnt/work_1/gest9386/CU_Boulder/rotations/LAYER/dormant_site_activation_pipeline/results/mutation_paths/AP1/paths.tsv"
LANDSCAPE_FILE = "/mnt/work_1/gest9386/CU_Boulder/rotations/LAYER/dormant_site_activation_pipeline/results/landscape/AP1/AP1_activation_landscape.tsv"
OUTPUT_DIR = Path("/mnt/work_1/gest9386/CU_Boulder/rotations/LAYER/dormant_site_activation_pipeline/results/purifying_selection/AP1")

AN_HIGH = 100_000
AN_MED = 50_000
MAX_AN = 1_614_006

# Use all available CPUs
N_WORKERS = os.cpu_count() or 32


def query_batch(positions_batch, coverage_file):
    """Query a batch of positions using a thread-local tabix file."""
    tbx = pysam.TabixFile(coverage_file)
    results = []
    
    for chrom, pos in positions_batch:
        try:
            an = 0
            for row in tbx.fetch(chrom, pos-1, pos):
                parts = row.split('\t')
                an = int(parts[2])
                break
            results.append((chrom, pos, an))
        except:
            results.append((chrom, pos, 0))
    
    tbx.close()
    return results


def main():
    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Check that indexed file exists
    if not Path(COVERAGE_FILE).exists():
        print(f"ERROR: Coverage file not found: {COVERAGE_FILE}")
        print("Run the bash script first to create the reformatted file.")
        sys.exit(1)
    
    if not Path(COVERAGE_FILE + ".tbi").exists():
        print(f"ERROR: Tabix index not found: {COVERAGE_FILE}.tbi")
        print("Run the bash script first to create the index.")
        sys.exit(1)
    
    print("="*60)
    print("Querying Coverage for Dormant Site Positions")
    print("="*60)
    
    print("\nLoading mutation positions from paths.tsv...")
    start = time.time()
    
    positions_by_hamming = defaultdict(set)
    with open(PATHS_FILE, 'r') as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 14:
                try:
                    chrom = parts[1]
                    pos = int(parts[13])
                    hamming = int(parts[7])
                    total_steps = int(parts[11])
                    if hamming == total_steps:  # single-step path
                        positions_by_hamming[hamming].add((f"chr{chrom}", pos))
                except:
                    continue
    
    for h in sorted(positions_by_hamming.keys()):
        print(f"  H={h}: {len(positions_by_hamming[h]):,} positions")
    
    # Load observed variants
    print("\nLoading observed variants...")
    landscape_df = pd.read_csv(LANDSCAPE_FILE, sep='\t')
    observed_by_hamming = defaultdict(set)
    for _, row in landscape_df.iterrows():
        parts = row['variant_id_str'].split(':')
        if len(parts) >= 2:
            observed_by_hamming[int(row['hamming_distance'])].add((parts[0], int(parts[1])))
    
    for h in sorted(observed_by_hamming.keys()):
        print(f"  H={h}: {len(observed_by_hamming[h]):,} observed")
    
    # Query coverage using tabix - PARALLELIZED
    print(f"\nQuerying coverage file with tabix (using {N_WORKERS} workers)...")
    print(f"  File: {COVERAGE_FILE}")
    
    results = []
    for hamming in sorted(positions_by_hamming.keys()):
        positions = list(positions_by_hamming[hamming])
        print(f"\n  Querying H={hamming} ({len(positions):,} positions)...")
        
        # Split positions into batches for parallel processing
        batch_size = max(1, len(positions) // N_WORKERS)
        batches = [positions[i:i+batch_size] for i in range(0, len(positions), batch_size)]
        print(f"    Split into {len(batches)} batches of ~{batch_size} positions each")
        
        # Query in parallel
        query_func = partial(query_batch, coverage_file=COVERAGE_FILE)
        
        with mp.Pool(N_WORKERS) as pool:
            batch_results = pool.map(query_func, batches)
        
        # Flatten results
        all_results = [r for batch in batch_results for r in batch]
        an_values = [r[2] for r in all_results]
        
        # Compute statistics
        an_array = np.array(an_values)
        high_conf = np.sum(an_array >= AN_HIGH)
        med_conf = np.sum((an_array >= AN_MED) & (an_array < AN_HIGH))
        low_conf = np.sum((an_array > 0) & (an_array < AN_MED))
        missing = np.sum(an_array == 0)
        
        n_observed = len(observed_by_hamming.get(hamming, set()))
        
        results.append({
            'hamming': hamming,
            'total_possible': len(positions),
            'high_conf': high_conf,
            'med_conf': med_conf,
            'low_conf': low_conf,
            'missing': missing,
            'observed': n_observed,
            'mean_an': np.mean(an_array[an_array > 0]) if np.any(an_array > 0) else 0,
            'median_an': np.median(an_array[an_array > 0]) if np.any(an_array > 0) else 0,
        })
        
        print(f"    High-conf (AN>=100K): {high_conf:,}")
        print(f"    Medium-conf: {med_conf:,}")
        print(f"    Low-conf: {low_conf:,}")
        print(f"    Missing: {missing:,}")
        print(f"    Observed in gnomAD: {n_observed:,}")
    
    # Compute constraint metrics
    print("\n" + "="*60)
    print("CONSTRAINT ANALYSIS RESULTS")
    print("="*60)
    
    results_df = pd.DataFrame(results)
    
    # Use H=3 as baseline rate
    h3 = results_df[results_df['hamming'] == 3].iloc[0]
    baseline_rate = h3['observed'] / h3['total_possible']
    
    for i, row in results_df.iterrows():
        expected = row['total_possible'] * baseline_rate
        fold_depletion = expected / row['observed'] if row['observed'] > 0 else float('inf')
        results_df.loc[i, 'expected'] = expected
        results_df.loc[i, 'fold_depletion'] = fold_depletion
        
        # Binomial test
        if row['total_possible'] > 0:
            p = stats.binom.cdf(row['observed'], row['total_possible'], baseline_rate)
            results_df.loc[i, 'binom_p'] = p
    
    # Print summary
    print(f"\nBaseline rate (H=3): {baseline_rate:.6f}")
    print(f"  ({h3['observed']:,} / {h3['total_possible']:,})")
    print()
    
    for _, row in results_df.iterrows():
        print(f"Hamming = {int(row['hamming'])}")
        print(f"  Possible:     {int(row['total_possible']):>12,}")
        print(f"  High-conf:    {int(row['high_conf']):>12,} ({row['high_conf']/row['total_possible']*100:.1f}%)")
        print(f"  Observed:     {int(row['observed']):>12,}")
        print(f"  Expected:     {row['expected']:>12,.1f}")
        print(f"  Fold depletion: {row['fold_depletion']:>10.1f}x")
        print(f"  Binomial p:   {row['binom_p']:>12.2e}")
        print()
    
    # Save results
    results_df.to_csv(OUTPUT_DIR / 'constraint_by_hamming.tsv', sep='\t', index=False)
    print(f"\nSaved: {OUTPUT_DIR / 'constraint_by_hamming.tsv'}")
    
    # Generate summary report
    with open(OUTPUT_DIR / 'purifying_selection_summary.txt', 'w') as f:
        f.write("="*70 + "\n")
        f.write("PURIFYING SELECTION ANALYSIS: AP1 DORMANT SITE ACTIVATION\n")
        f.write("="*70 + "\n\n")
        
        h1 = results_df[results_df['hamming'] == 1].iloc[0]
        
        f.write("KEY FINDING:\n")
        f.write("-"*40 + "\n")
        f.write(f"Of {int(h1['total_possible']):,} possible single-nucleotide mutations\n")
        f.write(f"that would create functional AP1 binding sites:\n\n")
        f.write(f"  * {int(h1['high_conf']):,} are at high-confidence positions (AN >= 100K)\n")
        f.write(f"  * Only {int(h1['observed'])} variant(s) observed in 807K individuals\n")
        f.write(f"  * Expected under neutrality: {h1['expected']:.1f}\n")
        f.write(f"  * Fold depletion: {h1['fold_depletion']:.1f}x\n")
        f.write(f"  * Binomial p-value: {h1['binom_p']:.2e}\n\n")
        
        if h1['fold_depletion'] > 10:
            f.write("-> STRONG EVIDENCE of purifying selection against AP1 site creation\n")
        
        f.write("\n" + "="*70 + "\n")
    
    print(f"Saved: {OUTPUT_DIR / 'purifying_selection_summary.txt'}")
    
    elapsed = time.time() - start
    print(f"\nTotal time: {elapsed:.1f} seconds")


if __name__ == '__main__':
    main()
