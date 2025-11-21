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
    
    # Filter to variants with gnomAD matches (AF > 0)
    df_matched = df[df['AF'] > 0].copy()
    print(f"Steps with gnomAD match: {len(df_matched):,}")
    
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
    
    # Save
    df_unique.to_csv(output_path, sep='\t', index=False)
    print(f"\nSaved: {output_path}")
    print(f"\nEstimated scoring time at 1.37 variants/sec:")
    minutes = len(df_unique) / 1.37 / 60
    print(f"  {minutes:.1f} minutes (~{minutes/60:.1f} hours)")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
