#!/usr/bin/env python3
"""
Per-variant analysis: Identify interesting variants for dormant site activation.

Instead of analyzing paths, analyze each observed variant individually:
- Variant-level constraint (AF, AN)
- Variant-level impact (AlphaGenome predictions)
- No assumptions about path completeness

This approach asks: "Which individual variants that could activate motifs
show strong constraint and/or high predicted impact?"
"""

import pandas as pd
from pathlib import Path
import sys

def main():
    input_file = Path("results/gnomad_intersection/AP1/paths_with_gnomad.tsv")
    
    print(f"Loading {input_file}...")
    df = pd.read_csv(input_file, sep='\t', 
                     usecols=['path_id', 'step_num', 'chr', 'pos', 
                              'ref', 'alt', 'AF', 'AN', 'is_missing', 'total_steps'],
                     low_memory=False)
    
    # Convert types after loading
    df['chr'] = df['chr'].astype(str)
    df['pos'] = pd.to_numeric(df['pos'], errors='coerce').fillna(0).astype(int)
    df['AF'] = pd.to_numeric(df['AF'], errors='coerce').fillna(0.0)
    df['AN'] = pd.to_numeric(df['AN'], errors='coerce').fillna(0.0)
    df['is_missing'] = pd.to_numeric(df['is_missing'], errors='coerce').fillna(1).astype(int)
    
    print(f"Loaded {len(df):,} mutation steps")
    
    # Create variant_id from chr:pos:ref:alt
    df['variant_id'] = df['chr'] + ':' + df['pos'].astype(str) + ':' + df['ref'] + ':' + df['alt']
    
    # Filter to observed variants only
    observed = df[df['is_missing'] == 0].copy()
    print(f"\nObserved variants in gnomAD: {len(observed):,}")
    print(f"Unique observed variants: {observed['variant_id'].nunique():,}")
    
    # Get unique variants (some variants appear in multiple paths)
    unique_variants = observed.drop_duplicates(subset='variant_id').copy()
    print(f"\nAnalyzing {len(unique_variants):,} unique observed variants...")
    
    # Categorize by constraint level
    print(f"\n{'='*70}")
    print(f"CONSTRAINT STRATIFICATION")
    print(f"{'='*70}\n")
    
    # AF=0 with high confidence
    af0_high = unique_variants[(unique_variants['AF'] == 0) & (unique_variants['AN'] >= 50000)]
    af0_med = unique_variants[(unique_variants['AF'] == 0) & 
                              (unique_variants['AN'] >= 10000) & 
                              (unique_variants['AN'] < 50000)]
    af0_low = unique_variants[(unique_variants['AF'] == 0) & (unique_variants['AN'] < 10000)]
    
    print(f"AF=0 (Never observed):")
    print(f"  High confidence (AN≥50K):  {len(af0_high):>6,} variants")
    print(f"  Medium confidence (10-50K): {len(af0_med):>6,} variants")
    print(f"  Low confidence (AN<10K):    {len(af0_low):>6,} variants")
    print(f"  Total AF=0:                 {len(af0_high) + len(af0_med) + len(af0_low):>6,} variants")
    
    # Ultra-rare
    ultra_rare = unique_variants[(unique_variants['AF'] > 0) & (unique_variants['AF'] < 1e-5)]
    print(f"\nUltra-rare (0 < AF < 1e-5):   {len(ultra_rare):>6,} variants")
    
    # Very rare
    very_rare = unique_variants[(unique_variants['AF'] >= 1e-5) & (unique_variants['AF'] < 1e-4)]
    print(f"Very rare (1e-5 ≤ AF < 1e-4): {len(very_rare):>6,} variants")
    
    # Rare
    rare = unique_variants[(unique_variants['AF'] >= 1e-4) & (unique_variants['AF'] < 1e-3)]
    print(f"Rare (1e-4 ≤ AF < 1e-3):      {len(rare):>6,} variants")
    
    # Common
    common = unique_variants[unique_variants['AF'] >= 1e-3]
    print(f"Common (AF ≥ 1e-3):           {len(common):>6,} variants")
    
    # Coverage distribution
    print(f"\n{'='*70}")
    print(f"COVERAGE DISTRIBUTION")
    print(f"{'='*70}\n")
    
    print(f"Mean AN:   {unique_variants['AN'].mean():>10,.0f}")
    print(f"Median AN: {unique_variants['AN'].median():>10,.0f}")
    print(f"Min AN:    {unique_variants['AN'].min():>10,.0f}")
    print(f"Max AN:    {unique_variants['AN'].max():>10,.0f}")
    
    high_cov_all = unique_variants[unique_variants['AN'] >= 50000]
    print(f"\nVariants with AN≥50K: {len(high_cov_all):,} ({100*len(high_cov_all)/len(unique_variants):.1f}%)")
    
    # Path length context
    print(f"\n{'='*70}")
    print(f"PATH LENGTH CONTEXT")
    print(f"{'='*70}\n")
    
    for length in sorted(unique_variants['total_steps'].unique()):
        variants_in_length = unique_variants[unique_variants['total_steps'] == length]
        af0_in_length = variants_in_length[(variants_in_length['AF'] == 0) & 
                                           (variants_in_length['AN'] >= 50000)]
        
        print(f"{length}-step paths: {len(variants_in_length):>6,} observed variants " +
              f"({len(af0_in_length):>5,} AF=0 high-conf)")
    
    # High-priority candidates for AlphaGenome interpretation
    print(f"\n{'='*70}")
    print(f"HIGH-PRIORITY CANDIDATES")
    print(f"{'='*70}\n")
    
    high_priority = af0_high.copy()
    print(f"AF=0 with high coverage (AN≥50K): {len(high_priority):,} variants")
    print(f"\nThese variants represent:")
    print(f"  - Part of motif activation pathway")
    print(f"  - Never observed despite good sampling")
    print(f"  - Strong negative selection candidate")
    print(f"  - Ready for AlphaGenome impact scoring")
    
    # Save high-priority variants
    output_file = Path("results/gnomad_intersection/AP1/high_priority_variants.tsv")
    high_priority_full = observed[observed['variant_id'].isin(high_priority['variant_id'])].copy()
    high_priority_full.to_csv(output_file, sep='\t', index=False)
    
    print(f"\nSaved {len(high_priority_full):,} rows ({len(high_priority):,} unique variants) to:")
    print(f"  {output_file}")
    
    # Also save all observed variants with metadata for downstream analysis
    output_file_all = Path("results/gnomad_intersection/AP1/all_observed_variants.tsv")
    observed.to_csv(output_file_all, sep='\t', index=False)
    
    print(f"\nSaved all {len(observed):,} observed variant instances to:")
    print(f"  {output_file_all}")
    
    # Summary for paper/interpretation
    print(f"\n{'='*70}")
    print(f"SUMMARY")
    print(f"{'='*70}\n")
    
    print(f"Total unique variants analyzed: {len(unique_variants):,}")
    print(f"Strong constraint candidates (AF=0, AN≥50K): {len(af0_high):,} ({100*len(af0_high)/len(unique_variants):.1f}%)")
    print(f"Ultra-rare variants (AF<1e-5): {len(ultra_rare):,} ({100*len(ultra_rare)/len(unique_variants):.1f}%)")
    print(f"\nInterpretation:")
    print(f"  Each variant scored independently by AlphaGenome")
    print(f"  X-axis: -log10(AF) measures constraint")
    print(f"  Y-axis: ΔAlphaGenome measures predicted impact")
    print(f"  No path-completeness assumptions required")
    
    return 0

if __name__ == '__main__':
    sys.exit(main())
