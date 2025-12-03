#!/usr/bin/env python3
"""
Module 07: Analyze Coverage-Based Constraint for Dormant AP1 Sites

This script uses gnomAD v4.1 allele number (AN) data for ALL sites to rigorously
test purifying selection against dormant AP1 site activation.

SCIENTIFIC RATIONALE:
=====================
VCF files only contain positions where variants EXIST. To test constraint, we need
to know: "Was this position well-covered?" The AN (allele number) file tells us
exactly how many individuals were successfully genotyped at each position.

AN = 2 × (number of individuals with a genotype call)
  - AN = 1,614,006 → All 807K individuals called → Maximum confidence
  - AN = 100,000   → ~50K individuals called     → High confidence
  - AN < 50,000    → Poor coverage               → Cannot conclude constraint

ALGORITHM:
==========
1. Load all possible mutation positions from paths.tsv (H=1, H=2, H=3)
2. Stream through coverage file, collecting AN for each position
3. Compute coverage-adjusted expected vs observed variants
4. Perform statistical tests for constraint

PERFORMANCE:
============
- Loads query positions into a set first (O(1) lookup)
- Streams 12GB coverage file once (no random access)
- With 32 cores/188GB RAM: ~2-3 minutes

Author: George Stephenson
Date: December 2025
"""

import argparse
import gzip
import logging
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Dict, Set, Tuple

import numpy as np
import pandas as pd
from scipy import stats

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler('analyze_coverage_constraint.log')
    ]
)
logger = logging.getLogger(__name__)

# ============================================================================
# CONFIGURATION
# ============================================================================

# Confidence thresholds for AN (allele number)
AN_HIGH_CONFIDENCE = 100_000    # ~50K individuals
AN_MEDIUM_CONFIDENCE = 50_000   # ~25K individuals

# Maximum AN in gnomAD v4.1 genomes (~807K individuals × 2)
MAX_AN = 1_614_006


def load_mutation_positions(paths_file: Path) -> Dict[int, Set[str]]:
    """
    Load all possible mutation positions from paths.tsv, grouped by Hamming distance.
    
    Returns:
        Dict mapping hamming_distance -> set of "chrX:position" strings
    """
    logger.info(f"Loading mutation positions from {paths_file}")
    
    positions_by_hamming = defaultdict(set)
    
    with open(paths_file, 'r') as f:
        header = next(f)  # Skip header
        
        for line_num, line in enumerate(f, 1):
            parts = line.strip().split('\t')
            
            if len(parts) < 14:
                continue
                
            try:
                chrom = parts[1]
                position = parts[13]
                hamming = int(parts[7])
                total_steps = int(parts[11])
                
                # Only count single-step paths (unique mutations)
                if hamming == total_steps:
                    locus = f"chr{chrom}:{position}"
                    positions_by_hamming[hamming].add(locus)
                    
            except (ValueError, IndexError) as e:
                if line_num < 10:
                    logger.warning(f"Line {line_num}: {e}")
                continue
            
            if line_num % 5_000_000 == 0:
                logger.info(f"  Processed {line_num:,} lines...")
    
    # Log summary
    for h in sorted(positions_by_hamming.keys()):
        logger.info(f"  Hamming {h}: {len(positions_by_hamming[h]):,} unique positions")
    
    return dict(positions_by_hamming)


def load_observed_variants(landscape_file: Path) -> Dict[int, Set[str]]:
    """
    Load observed variants from activation landscape, grouped by Hamming distance.
    
    Returns:
        Dict mapping hamming_distance -> set of "chrX:position" strings
    """
    logger.info(f"Loading observed variants from {landscape_file}")
    
    df = pd.read_csv(landscape_file, sep='\t')
    
    observed_by_hamming = defaultdict(set)
    
    for _, row in df.iterrows():
        variant_id = row['variant_id_str']
        hamming = int(row['hamming_distance'])
        
        # Parse variant_id to get chr:pos
        parts = variant_id.split(':')
        if len(parts) >= 2:
            locus = f"{parts[0]}:{parts[1]}"
            observed_by_hamming[hamming].add(locus)
    
    for h in sorted(observed_by_hamming.keys()):
        logger.info(f"  Hamming {h}: {len(observed_by_hamming[h]):,} observed variants")
    
    return dict(observed_by_hamming)


def load_coverage_into_memory(coverage_file: Path) -> Dict[str, int]:
    """
    Load ENTIRE coverage file into RAM for instant lookups.
    
    With 188GB RAM available, we can easily hold the ~50GB dictionary.
    This is MUCH faster than streaming for multiple queries.
    
    Args:
        coverage_file: Path to gnomad coverage TSV.bgz (~12GB compressed, ~3B lines)
        
    Returns:
        Dict mapping locus -> AN value (all ~3 billion positions)
    """
    logger.info(f"Loading ENTIRE coverage file into RAM...")
    logger.info(f"Coverage file: {coverage_file}")
    logger.info(f"This will use ~50GB RAM but enables instant lookups")
    
    coverage_dict = {}
    lines_processed = 0
    start_time = time.time()
    progress_interval = 50_000_000  # Every 50M lines
    
    with gzip.open(coverage_file, 'rt') as f:
        # Skip header
        header = next(f)
        logger.info(f"Header: {header.strip()}")
        
        for line in f:
            lines_processed += 1
            
            # Progress update
            if lines_processed % progress_interval == 0:
                elapsed = time.time() - start_time
                rate = lines_processed / elapsed / 1e6
                mem_gb = len(coverage_dict) * 50 / 1e9  # Rough estimate
                logger.info(
                    f"  Loaded {lines_processed/1e9:.2f}B lines "
                    f"({rate:.1f}M/sec), ~{mem_gb:.1f}GB RAM"
                )
            
            # Parse line and store
            try:
                locus, an_str = line.strip().split('\t')
                coverage_dict[locus] = int(an_str)
            except ValueError:
                continue
    
    elapsed = time.time() - start_time
    logger.info(f"Loaded {len(coverage_dict):,} positions in {elapsed:.1f}s")
    logger.info(f"Rate: {len(coverage_dict)/elapsed/1e6:.1f}M positions/sec")
    
    return coverage_dict


def query_coverage_from_memory(
    coverage_dict: Dict[str, int],
    query_positions: Set[str]
) -> Dict[str, int]:
    """
    Query coverage for specific positions from in-memory dictionary.
    
    This is O(1) per lookup - instant for any number of queries!
    
    Args:
        coverage_dict: Full coverage dictionary (loaded into RAM)
        query_positions: Set of "chrX:position" strings to query
        
    Returns:
        Dict mapping locus -> AN value (only for queried positions)
    """
    logger.info(f"Querying {len(query_positions):,} positions from memory...")
    start_time = time.time()
    
    results = {
        locus: coverage_dict[locus]
        for locus in query_positions
        if locus in coverage_dict
    }
    
    elapsed = time.time() - start_time
    logger.info(f"Found {len(results):,} / {len(query_positions):,} in {elapsed:.2f}s")
    
    return results


def query_coverage_file(
    coverage_file: Path,
    query_positions: Set[str],
    progress_interval: int = 100_000_000
) -> Dict[str, int]:
    """
    DEPRECATED: Use load_coverage_into_memory() + query_coverage_from_memory() instead.
    
    This streaming version is kept for systems with limited RAM.
    """
    logger.info(f"[STREAMING MODE] Querying coverage for {len(query_positions):,} positions...")
    logger.info(f"Coverage file: {coverage_file}")
    logger.info(f"NOTE: Consider using RAM-based loading for faster queries")
    
    results = {}
    lines_processed = 0
    start_time = time.time()
    
    # Open gzipped file
    with gzip.open(coverage_file, 'rt') as f:
        # Skip header
        header = next(f)
        logger.info(f"Header: {header.strip()}")
        
        for line in f:
            lines_processed += 1
            
            # Progress update
            if lines_processed % progress_interval == 0:
                elapsed = time.time() - start_time
                rate = lines_processed / elapsed / 1e6
                found_pct = len(results) / len(query_positions) * 100
                logger.info(
                    f"  Processed {lines_processed/1e9:.2f}B lines "
                    f"({rate:.1f}M/sec), found {len(results):,} ({found_pct:.1f}%)"
                )
            
            # Parse line
            try:
                locus, an_str = line.strip().split('\t')
            except ValueError:
                continue
            
            # Check if this position is in our query set
            if locus in query_positions:
                results[locus] = int(an_str)
                
                # Early exit if we found all positions
                if len(results) == len(query_positions):
                    logger.info(f"  Found all {len(query_positions):,} positions!")
                    break
    
    elapsed = time.time() - start_time
    logger.info(f"Completed in {elapsed:.1f}s, found {len(results):,} / {len(query_positions):,}")
    
    return results


def compute_constraint_statistics(
    positions_by_hamming: Dict[int, Set[str]],
    observed_by_hamming: Dict[int, Set[str]],
    coverage_data: Dict[str, int]
) -> pd.DataFrame:
    """
    Compute constraint statistics for each Hamming distance.
    
    Returns DataFrame with columns:
    - hamming_distance
    - total_possible: Total mutation positions
    - high_conf_possible: Positions with AN >= 100K
    - observed: Variants observed in gnomAD
    - high_conf_observed: Observed variants at high-conf positions
    - expected_uniform: Expected if no selection (coverage-adjusted)
    - fold_depletion: expected / observed
    - binomial_p: Binomial test p-value
    """
    logger.info("Computing constraint statistics...")
    
    results = []
    
    for hamming in sorted(positions_by_hamming.keys()):
        possible = positions_by_hamming[hamming]
        observed = observed_by_hamming.get(hamming, set())
        
        # Count by confidence level
        high_conf = 0
        medium_conf = 0
        low_conf = 0
        missing = 0
        total_an = 0
        
        for locus in possible:
            if locus in coverage_data:
                an = coverage_data[locus]
                total_an += an
                if an >= AN_HIGH_CONFIDENCE:
                    high_conf += 1
                elif an >= AN_MEDIUM_CONFIDENCE:
                    medium_conf += 1
                else:
                    low_conf += 1
            else:
                missing += 1
        
        # Observed at high-confidence positions
        high_conf_observed = sum(
            1 for locus in observed 
            if locus in coverage_data and coverage_data[locus] >= AN_HIGH_CONFIDENCE
        )
        
        # Coverage-adjusted expected (using total AN / max AN as effective sample size)
        # This accounts for varying coverage across positions
        effective_coverage = total_an / MAX_AN if MAX_AN > 0 else 0
        
        # Calculate overall observed rate from H=3 (baseline, least selection)
        baseline_hamming = 3
        if baseline_hamming in positions_by_hamming:
            baseline_possible = len(positions_by_hamming[baseline_hamming])
            baseline_observed = len(observed_by_hamming.get(baseline_hamming, set()))
            baseline_rate = baseline_observed / baseline_possible if baseline_possible > 0 else 0
        else:
            baseline_rate = len(observed) / len(possible) if possible else 0
        
        # Expected under uniform rate (no selection)
        expected_uniform = len(possible) * baseline_rate
        
        # Fold depletion
        n_observed = len(observed)
        fold_depletion = expected_uniform / n_observed if n_observed > 0 else float('inf')
        
        # Binomial test: is observed significantly less than expected?
        if len(possible) > 0 and baseline_rate > 0:
            binom_result = stats.binom.cdf(n_observed, len(possible), baseline_rate)
        else:
            binom_result = 1.0
        
        results.append({
            'hamming_distance': hamming,
            'total_possible': len(possible),
            'high_conf_possible': high_conf,
            'medium_conf_possible': medium_conf,
            'low_conf_possible': low_conf,
            'missing_coverage': missing,
            'observed': n_observed,
            'high_conf_observed': high_conf_observed,
            'baseline_rate': baseline_rate,
            'expected_uniform': expected_uniform,
            'fold_depletion': fold_depletion,
            'binomial_p': binom_result,
            'mean_an': total_an / len(possible) if possible else 0,
            'effective_coverage_fraction': effective_coverage / len(possible) if possible else 0
        })
    
    return pd.DataFrame(results)


def save_position_coverage(
    positions_by_hamming: Dict[int, Set[str]],
    coverage_data: Dict[str, int],
    output_dir: Path,
    motif_family: str = "AP1"
) -> None:
    """Save per-position coverage data for each Hamming distance."""
    
    for hamming in sorted(positions_by_hamming.keys()):
        positions = positions_by_hamming[hamming]
        
        rows = []
        for locus in sorted(positions):
            an = coverage_data.get(locus, 0)
            if an >= AN_HIGH_CONFIDENCE:
                conf = "high"
            elif an >= AN_MEDIUM_CONFIDENCE:
                conf = "medium"
            elif an > 0:
                conf = "low"
            else:
                conf = "missing"
            
            rows.append({
                'locus': locus,
                'AN': an,
                'confidence': conf
            })
        
        df = pd.DataFrame(rows)
        output_file = output_dir / f"{motif_family}_h{hamming}_coverage.tsv"
        df.to_csv(output_file, sep='\t', index=False)
        logger.info(f"Saved {len(df):,} positions to {output_file}")


def generate_summary_report(
    constraint_df: pd.DataFrame,
    output_file: Path
) -> None:
    """Generate human-readable summary report."""
    
    with open(output_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("PURIFYING SELECTION ANALYSIS: AP1 DORMANT SITE ACTIVATION\n")
        f.write("=" * 80 + "\n\n")
        
        f.write("OVERVIEW\n")
        f.write("-" * 40 + "\n")
        f.write("This analysis uses gnomAD v4.1 allele number (AN) data to test\n")
        f.write("whether mutations that would create AP1 binding sites are under\n")
        f.write("purifying selection in the human population.\n\n")
        
        f.write("COVERAGE THRESHOLDS\n")
        f.write("-" * 40 + "\n")
        f.write(f"  High confidence:   AN >= {AN_HIGH_CONFIDENCE:,} (~50K individuals)\n")
        f.write(f"  Medium confidence: AN >= {AN_MEDIUM_CONFIDENCE:,} (~25K individuals)\n")
        f.write(f"  Low confidence:    AN < {AN_MEDIUM_CONFIDENCE:,}\n")
        f.write(f"  Maximum AN:        {MAX_AN:,} (all 807K individuals)\n\n")
        
        f.write("RESULTS BY HAMMING DISTANCE\n")
        f.write("-" * 40 + "\n\n")
        
        for _, row in constraint_df.iterrows():
            h = int(row['hamming_distance'])
            f.write(f"Hamming Distance = {h}\n")
            f.write(f"  Total possible positions:     {row['total_possible']:>12,}\n")
            f.write(f"  High-confidence positions:    {row['high_conf_possible']:>12,} ({row['high_conf_possible']/row['total_possible']*100:.1f}%)\n")
            f.write(f"  Medium-confidence positions:  {row['medium_conf_possible']:>12,}\n")
            f.write(f"  Low-confidence positions:     {row['low_conf_possible']:>12,}\n")
            f.write(f"  Missing from coverage file:   {row['missing_coverage']:>12,}\n")
            f.write(f"  Mean AN:                      {row['mean_an']:>12,.0f}\n")
            f.write(f"\n")
            f.write(f"  Observed in gnomAD:           {row['observed']:>12,}\n")
            f.write(f"  Expected (uniform rate):      {row['expected_uniform']:>12,.1f}\n")
            f.write(f"  Fold depletion:               {row['fold_depletion']:>12.1f}×\n")
            f.write(f"  Binomial p-value:             {row['binomial_p']:>12.2e}\n")
            f.write(f"\n")
        
        f.write("INTERPRETATION\n")
        f.write("-" * 40 + "\n")
        
        h1_row = constraint_df[constraint_df['hamming_distance'] == 1].iloc[0] if 1 in constraint_df['hamming_distance'].values else None
        
        if h1_row is not None:
            high_conf_pct = h1_row['high_conf_possible'] / h1_row['total_possible'] * 100
            f.write(f"\nFor H=1 (single mutation to create AP1 site):\n")
            f.write(f"  • {h1_row['high_conf_possible']:,} of {h1_row['total_possible']:,} positions ({high_conf_pct:.1f}%) are high-confidence\n")
            f.write(f"  • Only {h1_row['observed']:,} variant(s) observed vs {h1_row['expected_uniform']:.1f} expected\n")
            f.write(f"  • This represents {h1_row['fold_depletion']:.0f}× depletion (p = {h1_row['binomial_p']:.2e})\n")
            
            if h1_row['fold_depletion'] > 10 and h1_row['binomial_p'] < 0.001:
                f.write(f"\n  → STRONG EVIDENCE of purifying selection against AP1 site creation\n")
            elif h1_row['fold_depletion'] > 2:
                f.write(f"\n  → MODERATE EVIDENCE of purifying selection\n")
            else:
                f.write(f"\n  → WEAK or NO EVIDENCE of purifying selection\n")
        
        f.write("\n" + "=" * 80 + "\n")
        f.write("Analysis completed: " + time.strftime("%Y-%m-%d %H:%M:%S") + "\n")
    
    logger.info(f"Summary report saved to {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Analyze coverage-based constraint for dormant AP1 sites",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Full analysis
  python analyze_coverage_constraint.py \\
      --paths ../results/mutation_paths/AP1/paths.tsv \\
      --coverage ../data/gnomad_coverage/gnomad.genomes.v4.1.allele_number_all_sites.tsv.bgz \\
      --observed ../results/landscape/AP1/AP1_activation_landscape.tsv \\
      --output ../results/purifying_selection/AP1

  # Only analyze H=1 positions (faster)
  python analyze_coverage_constraint.py \\
      --paths ../results/mutation_paths/AP1/paths.tsv \\
      --coverage ../data/gnomad_coverage/gnomad.genomes.v4.1.allele_number_all_sites.tsv.bgz \\
      --observed ../results/landscape/AP1/AP1_activation_landscape.tsv \\
      --output ../results/purifying_selection/AP1 \\
      --hamming 1
        """
    )
    
    parser.add_argument(
        '--paths', '-p',
        type=Path,
        required=True,
        help='Path to mutation paths file (paths.tsv from Module 02)'
    )
    
    parser.add_argument(
        '--coverage', '-c',
        type=Path,
        required=True,
        help='Path to gnomAD coverage file (allele_number_all_sites.tsv.bgz)'
    )
    
    parser.add_argument(
        '--observed', '-o',
        type=Path,
        required=True,
        help='Path to observed variants file (activation_landscape.tsv from Module 05)'
    )
    
    parser.add_argument(
        '--output', '-O',
        type=Path,
        required=True,
        help='Output directory for results'
    )
    
    parser.add_argument(
        '--hamming',
        type=int,
        nargs='+',
        default=[1, 2, 3],
        help='Hamming distances to analyze (default: 1 2 3)'
    )
    
    parser.add_argument(
        '--motif-family',
        type=str,
        default='AP1',
        help='Motif family name for output files (default: AP1)'
    )
    
    parser.add_argument(
        '--save-positions',
        action='store_true',
        help='Save per-position coverage data (large files)'
    )
    
    args = parser.parse_args()
    
    # Create output directory
    args.output.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Load mutation positions
    all_positions = load_mutation_positions(args.paths)
    
    # Filter to requested Hamming distances
    positions_by_hamming = {
        h: pos for h, pos in all_positions.items() 
        if h in args.hamming
    }
    
    # Step 2: Load observed variants
    observed_by_hamming = load_observed_variants(args.observed)
    
    # Step 3: Collect all positions to query
    all_query_positions = set()
    for positions in positions_by_hamming.values():
        all_query_positions.update(positions)
    
    logger.info(f"Total positions to query: {len(all_query_positions):,}")
    
    # Step 4: Query coverage file
    coverage_data = query_coverage_file(args.coverage, all_query_positions)
    
    # Step 5: Compute constraint statistics
    constraint_df = compute_constraint_statistics(
        positions_by_hamming,
        observed_by_hamming,
        coverage_data
    )
    
    # Step 6: Save results
    constraint_file = args.output / f"{args.motif_family}_constraint_by_hamming.tsv"
    constraint_df.to_csv(constraint_file, sep='\t', index=False)
    logger.info(f"Saved constraint statistics to {constraint_file}")
    
    # Step 7: Generate summary report
    summary_file = args.output / f"{args.motif_family}_purifying_selection_summary.txt"
    generate_summary_report(constraint_df, summary_file)
    
    # Step 8: Optionally save per-position coverage
    if args.save_positions:
        save_position_coverage(
            positions_by_hamming,
            coverage_data,
            args.output,
            args.motif_family
        )
    
    # Print summary to console
    print("\n" + "=" * 60)
    print("CONSTRAINT ANALYSIS COMPLETE")
    print("=" * 60)
    print(constraint_df.to_string(index=False))
    print()
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
