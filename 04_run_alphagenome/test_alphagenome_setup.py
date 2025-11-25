#!/usr/bin/env python3
"""
Quick test of AlphaGenome API setup.

Tests:
1. Import AlphaGenome modules
2. Initialize client with API key
3. Score a single test variant
"""

import os
import sys

def test_imports():
    """Test AlphaGenome module imports."""
    print("Testing AlphaGenome imports...")
    try:
        from alphagenome.data import genome
        from alphagenome.models import dna_client, variant_scorers
        print("✓ Successfully imported AlphaGenome modules")
        return True
    except ImportError as e:
        print(f"✗ Failed to import AlphaGenome: {e}")
        print("\nMake sure you're in the alphagenome conda environment:")
        print("  conda activate alphagenome-env")
        return False


def test_api_key():
    """Test API key availability."""
    print("\nTesting API key...")
    api_key = os.environ.get('ALPHA_GENOME_KEY')
    if api_key:
        print(f"✓ API key found (length: {len(api_key)})")
        return api_key
    else:
        print("✗ ALPHA_GENOME_KEY environment variable not set")
        print("\nSet your API key:")
        print("  export ALPHA_GENOME_KEY='your_key_here'")
        return None


def test_client(api_key):
    """Test AlphaGenome client initialization."""
    print("\nTesting client initialization...")
    try:
        from alphagenome.models import dna_client
        client = dna_client.create(api_key=api_key)
        print("✓ Successfully initialized AlphaGenome client")
        print(f"  Context window size: {dna_client.SEQUENCE_LENGTH_1MB:,} bp (1MB)")
        return client
    except Exception as e:
        print(f"✗ Failed to initialize client: {e}")
        return None


def test_variant_scoring(client):
    """Test scoring a single variant."""
    print("\nTesting variant scoring...")
    try:
        from alphagenome.data import genome
        from alphagenome.models import dna_client, variant_scorers
        
        # Create a test variant (rs334 - HBB sickle cell variant)
        variant = genome.Variant(
            chromosome="chr11",
            position=5227002,
            reference_bases="T",
            alternate_bases="A",
            name="rs334_test"
        )
        
        print(f"  Test variant: {variant.name} at chr11:5227002 T>A")
        
        # Get 1MB interval
        interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
        print(f"  Interval: {interval.chromosome}:{interval.start}-{interval.end} ({interval.end-interval.start:,} bp)")
        
        # Score with ALL scorers to capture all features (expression, ATAC, DNase, histones)
        scorers = list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.values())
        print(f"  Scoring with {len(scorers)} scorer(s): {list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.keys())}...")
        
        scores = client.score_variant(
            interval=interval,
            variant=variant,
            variant_scorers=scorers,
            organism=dna_client.Organism.HOMO_SAPIENS
        )
        
        # Convert to dataframe
        df = variant_scorers.tidy_scores([scores])
        
        print(f"✓ Successfully scored variant")
        print(f"  Returned {len(df)} track scores")
        print(f"  Columns: {df.columns.tolist()}")
        print(f"  Unique output_type values: {df['output_type'].unique().tolist() if 'output_type' in df.columns else 'N/A'}")
        print(f"  Unique histone_mark values: {sorted(df['histone_mark'].dropna().unique().tolist()) if 'histone_mark' in df.columns else 'N/A'}")
        print(f"  Sample track names: {df['track_name'].head(5).tolist() if 'track_name' in df.columns else 'N/A'}")
        print(f"  Mean quantile score: {df['quantile_score'].mean():.3f}")
        print(f"  Score range: [{df['quantile_score'].min():.3f}, {df['quantile_score'].max():.3f}]")
        
        return True
        
    except Exception as e:
        print(f"✗ Failed to score variant: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    print("="*60)
    print("AlphaGenome API Setup Test")
    print("="*60)
    
    # Run tests
    if not test_imports():
        return 1
    
    api_key = test_api_key()
    if not api_key:
        return 1
    
    client = test_client(api_key)
    if not client:
        return 1
    
    if not test_variant_scoring(client):
        return 1
    
    print("\n" + "="*60)
    print("✓ All tests passed! AlphaGenome is ready to use.")
    print("="*60)
    print("\nYou can now run:")
    print("  python 04_run_alphagenome/run_alphagenome_scoring.py --config pipeline_config.yaml")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
