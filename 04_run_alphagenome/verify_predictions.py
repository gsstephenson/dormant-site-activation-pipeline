#!/usr/bin/env python3
"""
Verify AlphaGenome predictions data quality and completeness.
"""

import pandas as pd
import pyarrow.parquet as pq
from pathlib import Path
import sys

def verify_predictions(parquet_path):
    """Comprehensive verification of AlphaGenome predictions."""
    
    print("="*80)
    print("ALPHAGENOME PREDICTIONS VERIFICATION")
    print("="*80)
    
    # Read parquet metadata first (no data loading)
    parquet_file = pq.ParquetFile(parquet_path)
    
    print(f"\n1. FILE METADATA:")
    print(f"   File: {parquet_path}")
    print(f"   Size: {parquet_path.stat().st_size / 1e9:.2f} GB")
    print(f"   Rows: {parquet_file.metadata.num_rows:,}")
    print(f"   Columns: {parquet_file.schema.names}")
    
    # Read actual data in chunks for memory efficiency
    print(f"\n2. LOADING DATA (may take a moment)...")
    df = pd.read_parquet(parquet_path)
    
    print(f"\n3. DATA DIMENSIONS:")
    print(f"   Total rows: {len(df):,}")
    print(f"   Total columns: {len(df.columns)}")
    print(f"   Memory usage: {df.memory_usage(deep=True).sum() / 1e9:.2f} GB")
    
    print(f"\n4. VARIANT COVERAGE:")
    unique_variants = df['variant_id_str'].nunique()
    print(f"   Unique variants: {unique_variants:,}")
    avg_tracks = len(df) / unique_variants
    print(f"   Average tracks per variant: {avg_tracks:,.0f}")
    print(f"   Expected tracks per variant: ~89,000")
    coverage_pct = (avg_tracks / 89000) * 100
    print(f"   Coverage: {coverage_pct:.1f}%")
    
    if coverage_pct < 90:
        print(f"   ⚠️  WARNING: Only {coverage_pct:.1f}% of expected tracks!")
    else:
        print(f"   ✓ Good coverage")
    
    print(f"\n5. OUTPUT_TYPE DISTRIBUTION:")
    if 'output_type' in df.columns:
        output_type_counts = df.groupby('output_type').agg({
            'variant_id_str': 'count',
            'quantile_score': 'mean'
        }).rename(columns={'variant_id_str': 'total_rows'})
        
        output_type_counts['tracks_per_variant'] = output_type_counts['total_rows'] / unique_variants
        output_type_counts = output_type_counts.sort_values('total_rows', ascending=False)
        
        print(f"   Total output_types: {df['output_type'].nunique()}")
        print(f"\n   Breakdown:")
        print(f"   {'Output Type':<20} {'Total Rows':>15} {'Tracks/Variant':>15} {'Mean Score':>12}")
        print(f"   {'-'*20} {'-'*15} {'-'*15} {'-'*12}")
        for output_type, row in output_type_counts.iterrows():
            print(f"   {output_type:<20} {int(row['total_rows']):>15,} {int(row['tracks_per_variant']):>15,} {row['quantile_score']:>12.4f}")
        
        # Expected output types
        expected_types = {
            'ATAC', 'CAGE', 'CHIP_HISTONE', 'CHIP_TF', 'CONTACT_MAPS', 
            'DNASE', 'PROCAP', 'RNA_SEQ', 'SPLICE_JUNCTIONS', 
            'SPLICE_SITES', 'SPLICE_SITE_USAGE'
        }
        actual_types = set(df['output_type'].unique())
        missing_types = expected_types - actual_types
        
        if missing_types:
            print(f"\n   ⚠️  MISSING OUTPUT TYPES: {missing_types}")
        else:
            print(f"\n   ✓ All expected output_types present")
    else:
        print(f"   ⚠️  ERROR: No 'output_type' column found!")
    
    print(f"\n6. CELL TYPE/BIOSAMPLE COVERAGE:")
    if 'biosample_name' in df.columns:
        n_biosamples = df['biosample_name'].nunique()
        print(f"   Unique biosamples: {n_biosamples}")
        top_biosamples = df['biosample_name'].value_counts().head(10)
        print(f"\n   Top 10 biosamples by track count:")
        for biosample, count in top_biosamples.items():
            pct = (count / len(df)) * 100
            print(f"     {biosample}: {count:,} tracks ({pct:.1f}%)")
    
    print(f"\n7. SCORE QUALITY:")
    print(f"   Quantile score statistics:")
    print(df['quantile_score'].describe())
    
    # Check for missing/null scores
    null_scores = df['quantile_score'].isnull().sum()
    if null_scores > 0:
        print(f"   ⚠️  WARNING: {null_scores:,} null scores ({null_scores/len(df)*100:.2f}%)")
    else:
        print(f"   ✓ No null scores")
    
    print(f"\n8. SAMPLE VARIANTS:")
    sample_variants = df['variant_id_str'].unique()[:5]
    for var in sample_variants:
        var_df = df[df['variant_id_str'] == var]
        print(f"\n   {var}:")
        print(f"     Total tracks: {len(var_df):,}")
        if 'output_type' in df.columns:
            type_dist = var_df['output_type'].value_counts()
            print(f"     Output types: {len(type_dist)}")
            for ot, count in type_dist.head(5).items():
                print(f"       - {ot}: {count} tracks")
    
    print(f"\n9. METADATA COMPLETENESS:")
    metadata_cols = ['gnomad_AF', 'gnomad_AC', 'gnomad_AN', 'path_id', 'step_num']
    for col in metadata_cols:
        if col in df.columns:
            non_null = df[col].notna().sum()
            pct = (non_null / len(df)) * 100
            status = "✓" if pct > 95 else "⚠️"
            print(f"   {status} {col}: {non_null:,}/{len(df):,} ({pct:.1f}%)")
        else:
            print(f"   ⚠️  {col}: MISSING")
    
    print(f"\n10. EXPECTED VS ACTUAL TRACKS:")
    expected_tracks_per_variant = 89000
    actual_tracks_per_variant = len(df) / unique_variants
    difference = expected_tracks_per_variant - actual_tracks_per_variant
    
    print(f"   Expected: {expected_tracks_per_variant:,} tracks/variant")
    print(f"   Actual:   {int(actual_tracks_per_variant):,} tracks/variant")
    print(f"   Missing:  {int(difference):,} tracks/variant ({difference/expected_tracks_per_variant*100:.1f}%)")
    
    if difference > 10000:
        print(f"\n   ⚠️  SIGNIFICANT MISSING DATA DETECTED")
        print(f"   Possible causes:")
        print(f"     - Some scorers may have failed silently")
        print(f"     - API may have returned incomplete results")
        print(f"     - Certain cell types/tracks may be unavailable")
        print(f"     - Context window size may have limited track availability")
    
    print(f"\n{'='*80}")
    print(f"VERDICT:")
    print(f"{'='*80}")
    
    if coverage_pct >= 90 and not missing_types and null_scores == 0:
        print(f"✓ Data quality EXCELLENT - proceed to Module 05")
    elif coverage_pct >= 70:
        print(f"⚠️  Data quality ACCEPTABLE but incomplete - can proceed with caution")
        print(f"   Consider: Which analyses require the missing tracks?")
    else:
        print(f"❌ Data quality POOR - consider re-scoring or investigating")
    
    print(f"{'='*80}\n")
    
    return df

if __name__ == "__main__":
    parquet_path = Path("/mnt/work_1/gest9386/CU_Boulder/rotations/LAYER/dormant_site_activation_pipeline/results/alphagenome/AP1/predictions.parquet")
    
    if not parquet_path.exists():
        print(f"ERROR: File not found: {parquet_path}")
        sys.exit(1)
    
    df = verify_predictions(parquet_path)
