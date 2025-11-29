#!/usr/bin/env python3
"""
Prepare unique variants for AlphaGenome scoring.

Takes paths_with_gnomad.tsv and creates a deduplicated list of 
variants that actually matched gnomAD (AF > 0).

This avoids scoring the same variant multiple times across different
mutation paths.
"""

import pandas as pd
import sys
from pathlib import Path

def main():
    base_dir = Path(__file__).parent.parent
    input_path = base_dir / "results" / "gnomad_intersection" / "AP1" / "paths_with_gnomad.tsv"
    output_path = base_dir / "results" / "gnomad_intersection" / "AP1" / "unique_variants.tsv"
    
    print(f"Loading: {input_path}")
    df = pd.read_csv(input_path, sep='\t', low_memory=False)
    print(f"Total mutation steps: {len(df):,}")
    
    # Handle VCF "." notation for missing values
    # bcftools outputs "." for missing INFO fields, which causes pandas to read as object type
    numeric_cols = ['AF', 'AC', 'AN', 'nhomalt']
    for col in numeric_cols:
        if col in df.columns:
            df[col] = df[col].replace('.', pd.NA)
            df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0.0)
    
    # Report on coverage confidence metrics
    if 'is_missing' in df.columns and 'coverage_confidence' in df.columns:
        num_missing = df['is_missing'].sum()
        num_observed = len(df) - num_missing
        print(f"\nCoverage breakdown:")
        print(f"  Observed in gnomAD: {num_observed:,} ({100*num_observed/len(df):.1f}%)")
        print(f"  Missing (AF=0, unknown coverage): {num_missing:,} ({100*num_missing/len(df):.1f}%)")
        
        if num_observed > 0:
            conf_counts = df[df['is_missing'] == 0]['coverage_confidence'].value_counts()
            print(f"\n  Coverage confidence (observed variants only):")
            for conf in ['high', 'medium', 'low', 'zero_an']:
                count = conf_counts.get(conf, 0)
                pct = 100 * count / num_observed if num_observed > 0 else 0
                print(f"    {conf}: {count:,} ({pct:.1f}%)")
    
    # Filter strategy:
    # KEEP: Variants observed in gnomAD (is_missing=0), including AF=0 with high AN
    # SKIP: Variants missing from gnomAD (is_missing=1, unknown coverage)
    # Rationale: AF=0 with high AN represents true constraint and should be scored
    if 'is_missing' in df.columns:
        df_matched = df[df['is_missing'] == 0].copy()
        print(f"\nSteps observed in gnomAD (including AF=0 with known AN): {len(df_matched):,}")
        
        # Report on AF=0 variants being kept
        af_zero_kept = (df_matched['AF'] == 0).sum()
        print(f"  Including AF=0 variants (high-confidence constraint): {af_zero_kept:,}")
    else:
        # Fallback: old behavior if columns don't exist
        df_matched = df[df['AF'] > 0].copy()
        print(f"\nSteps with gnomAD match (AF > 0): {len(df_matched):,}")
    
    # Create variant ID and deduplicate
    df_matched['variant_id'] = (
        df_matched['chr'].astype(str) + '_' +
        df_matched['genomic_position'].astype(str) + '_' +
        df_matched['ref_base'].astype(str) + '_' +
        df_matched['alt_base'].astype(str)
    )
    
    # Keep first occurrence of each unique variant (preserves metadata)
    df_unique = df_matched.drop_duplicates(subset='variant_id', keep='first')
    print(f"Unique variants: {len(df_unique):,}")
    
    # Report on AN distribution for unique variants
    if 'AN' in df_unique.columns:
        an_values = df_unique['AN']
        print(f"\nAN (allele number) distribution for unique variants:")
        print(f"  Mean: {an_values.mean():,.0f}")
        print(f"  Median: {an_values.median():,.0f}")
        print(f"  Min: {an_values.min():,.0f}")
        print(f"  Max: {an_values.max():,.0f}")
        high_conf = (an_values >= 50000).sum()
        print(f"  High confidence (AN≥50K): {high_conf:,} ({100*high_conf/len(df_unique):.1f}%)")
    
    # Final summary by AF category
    print(f"\nFinal unique variants by constraint category:")
    af_zero_unique = (df_unique['AF'] == 0).sum()
    af_nonzero_unique = (df_unique['AF'] > 0).sum()
    print(f"  AF=0 (constraint candidates): {af_zero_unique:,} ({100*af_zero_unique/len(df_unique):.1f}%)")
    print(f"  AF>0 (existing variation): {af_nonzero_unique:,} ({100*af_nonzero_unique/len(df_unique):.1f}%)")
    
    if 'coverage_confidence' in df_unique.columns:
        af_zero_high_conf = ((df_unique['AF'] == 0) & (df_unique['coverage_confidence'] == 'high')).sum()
        print(f"    → AF=0 with high confidence (AN≥50K): {af_zero_high_conf:,}")
    
    # Save
    df_unique.to_csv(output_path, sep='\t', index=False)
    print(f"\n✓ Saved: {output_path}")
    print(f"\nEstimated scoring time at 1.37 variants/sec:")
    minutes = len(df_unique) / 1.37 / 60
    print(f"  {minutes:.1f} minutes (~{minutes/60:.1f} hours)")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
