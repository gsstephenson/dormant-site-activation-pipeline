#!/usr/bin/env python3
"""
Analyze forbidden variant enrichment near GWAS loci and ClinVar positions.

This is the main analysis script that:
1. Tests proximity enrichment (are forbidden sites near GWAS hits?)
2. Tests exact overlap (should be ~0 for both GWAS and ClinVar)
3. Generates matched null distributions for statistical testing
"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path
import logging
from scipy import stats
from collections import Counter
import subprocess
import tempfile
import os

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser(description='Analyze forbidden variant enrichment')
    parser.add_argument('--forbidden', required=True, help='Path to forbidden_variants.tsv')
    parser.add_argument('--forbidden-scores', help='Path to predictions_summary_ap1.tsv')
    parser.add_argument('--gwas', required=True, help='Path to processed GWAS BED directory')
    parser.add_argument('--clinvar', required=True, help='Path to processed ClinVar BED directory')
    parser.add_argument('--genome-callable', help='Path to callable sites BED (for null distribution)')
    parser.add_argument('--output', required=True, help='Output directory')
    parser.add_argument('--n-permutations', type=int, default=1000, help='Number of permutations for null')
    return parser.parse_args()


def load_forbidden_variants(filepath: str) -> pd.DataFrame:
    """Load forbidden variants and create BED-compatible format."""
    logger.info(f"Loading forbidden variants from {filepath}")
    
    df = pd.read_csv(filepath, sep='\t')
    logger.info(f"Loaded {len(df):,} forbidden variants")
    
    # Ensure chr prefix
    df['chr'] = df['chr'].astype(str)
    df['chr'] = df['chr'].apply(lambda x: f'chr{x}' if not x.startswith('chr') else x)
    
    return df


def forbidden_to_bed(df: pd.DataFrame, output_path: str):
    """Convert forbidden variants to BED format."""
    
    bed_df = pd.DataFrame({
        'chrom': df['chr'],
        'start': df['pos'] - 1,  # BED is 0-based
        'end': df['pos'],
        'name': df['variant_id'],
    })
    
    bed_df = bed_df.sort_values(['chrom', 'start'])
    bed_df.to_csv(output_path, sep='\t', index=False, header=False)
    
    return output_path


def run_bedtools_intersect(a_bed: str, b_bed: str, output: str = None, count_only: bool = False) -> int:
    """Run bedtools intersect and return count or write output."""
    
    if count_only:
        # Just count overlaps
        cmd = f"bedtools intersect -a {a_bed} -b {b_bed} -u | wc -l"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        return int(result.stdout.strip())
    else:
        # Write full output
        cmd = f"bedtools intersect -a {a_bed} -b {b_bed} -wa -wb > {output}"
        subprocess.run(cmd, shell=True, check=True)
        return output


def run_bedtools_window(a_bed: str, b_bed: str, window: int, count_only: bool = True) -> int:
    """Run bedtools window to find variants within window of targets."""
    
    if count_only:
        cmd = f"bedtools window -a {a_bed} -b {b_bed} -w {window} | cut -f1-4 | sort -u | wc -l"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        return int(result.stdout.strip())
    else:
        raise NotImplementedError("Non-count mode not implemented")


def analyze_proximity_enrichment(forbidden_bed: str, gwas_bed: str, 
                                  windows: list = [50000, 100000, 500000]) -> dict:
    """Analyze enrichment of forbidden variants near GWAS loci at different windows."""
    
    results = {}
    
    # Get total forbidden variants
    cmd = f"wc -l < {forbidden_bed}"
    total_forbidden = int(subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout.strip())
    
    logger.info(f"Total forbidden variants: {total_forbidden:,}")
    
    for window in windows:
        # Count forbidden variants within window of GWAS loci
        n_near_gwas = run_bedtools_window(forbidden_bed, gwas_bed, window)
        
        pct_near = 100 * n_near_gwas / total_forbidden if total_forbidden > 0 else 0
        
        results[f'window_{window//1000}kb'] = {
            'window_size': window,
            'n_forbidden_near_gwas': n_near_gwas,
            'pct_forbidden_near_gwas': pct_near,
            'total_forbidden': total_forbidden,
        }
        
        logger.info(f"Window {window//1000}kb: {n_near_gwas:,} ({pct_near:.2f}%) forbidden variants near GWAS loci")
    
    return results


def analyze_exact_overlap(forbidden_bed: str, target_bed: str, target_name: str) -> dict:
    """Analyze exact position overlap between forbidden variants and target database."""
    
    # Count exact overlaps
    n_overlap = run_bedtools_intersect(forbidden_bed, target_bed, count_only=True)
    
    # Get totals
    cmd = f"wc -l < {forbidden_bed}"
    total_forbidden = int(subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout.strip())
    
    cmd = f"wc -l < {target_bed}"
    total_target = int(subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout.strip())
    
    pct_overlap = 100 * n_overlap / total_forbidden if total_forbidden > 0 else 0
    
    result = {
        'target': target_name,
        'n_forbidden_overlap': n_overlap,
        'pct_forbidden_overlap': pct_overlap,
        'total_forbidden': total_forbidden,
        'total_target': total_target,
    }
    
    logger.info(f"{target_name}: {n_overlap:,} exact overlaps ({pct_overlap:.4f}%)")
    
    return result


def analyze_by_trait_category(forbidden_bed: str, gwas_dir: str, 
                               window: int = 100000) -> pd.DataFrame:
    """Analyze enrichment stratified by GWAS trait category."""
    
    results = []
    
    # Find category-specific BED files
    gwas_path = Path(gwas_dir)
    category_files = list(gwas_path.glob('gwas_*.bed'))
    
    # Filter to category files (exclude windowed files)
    category_files = [f for f in category_files if 'kb' not in f.name and 
                      f.name not in ['gwas_significant.bed', 'gwas_unique_loci.bed']]
    
    for cat_file in category_files:
        category = cat_file.stem.replace('gwas_', '')
        
        n_near = run_bedtools_window(forbidden_bed, str(cat_file), window)
        
        # Get number of GWAS hits in this category
        cmd = f"wc -l < {cat_file}"
        n_gwas = int(subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout.strip())
        
        results.append({
            'category': category,
            'n_gwas_loci': n_gwas,
            'n_forbidden_nearby': n_near,
            'window_kb': window // 1000,
        })
        
        logger.info(f"  {category}: {n_near:,} forbidden near {n_gwas:,} GWAS loci")
    
    return pd.DataFrame(results)


def generate_genome_null(forbidden_df: pd.DataFrame, n_permutations: int = 1000,
                         gwas_bed: str = None, window: int = 100000) -> dict:
    """Generate null distribution by sampling random genomic positions."""
    
    logger.info(f"Generating null distribution with {n_permutations} permutations...")
    
    # Get chromosome sizes (approximate for hg38)
    chrom_sizes = {
        'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555,
        'chr5': 181538259, 'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636,
        'chr9': 138394717, 'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
        'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189, 'chr16': 90338345,
        'chr17': 83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167,
        'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895,
    }
    
    # Match chromosome distribution of forbidden variants
    chrom_counts = forbidden_df['chr'].value_counts().to_dict()
    n_total = len(forbidden_df)
    
    null_counts = []
    
    for i in range(n_permutations):
        if (i + 1) % 100 == 0:
            logger.info(f"  Permutation {i+1}/{n_permutations}")
        
        # Generate random positions matching chromosome distribution
        with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
            for chrom, count in chrom_counts.items():
                if chrom in chrom_sizes:
                    positions = np.random.randint(1, chrom_sizes[chrom], size=count)
                    for pos in positions:
                        f.write(f"{chrom}\t{pos-1}\t{pos}\trandom_{i}\n")
            
            null_bed = f.name
        
        # Count how many null positions are near GWAS
        n_near = run_bedtools_window(null_bed, gwas_bed, window)
        null_counts.append(n_near)
        
        # Cleanup
        os.unlink(null_bed)
    
    return {
        'null_counts': null_counts,
        'null_mean': np.mean(null_counts),
        'null_std': np.std(null_counts),
        'null_median': np.median(null_counts),
        'null_95th': np.percentile(null_counts, 95),
    }


def calculate_enrichment_pvalue(observed: int, null_counts: list) -> dict:
    """Calculate enrichment statistics and p-value."""
    
    null_mean = np.mean(null_counts)
    null_std = np.std(null_counts)
    
    # Enrichment ratio
    enrichment = observed / null_mean if null_mean > 0 else np.inf
    
    # Empirical p-value (one-tailed, testing for enrichment)
    n_greater = sum(1 for x in null_counts if x >= observed)
    pvalue_enriched = (n_greater + 1) / (len(null_counts) + 1)
    
    # Also test for depletion
    n_less = sum(1 for x in null_counts if x <= observed)
    pvalue_depleted = (n_less + 1) / (len(null_counts) + 1)
    
    # Z-score
    zscore = (observed - null_mean) / null_std if null_std > 0 else 0
    
    return {
        'observed': observed,
        'expected': null_mean,
        'enrichment_ratio': enrichment,
        'zscore': zscore,
        'pvalue_enriched': pvalue_enriched,
        'pvalue_depleted': pvalue_depleted,
        'null_std': null_std,
    }


def write_report(results: dict, output_path: str):
    """Write summary report."""
    
    with open(output_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("FORBIDDEN VARIANTS: GWAS & ClinVar Validation Report\n")
        f.write("=" * 80 + "\n\n")
        
        # Proximity enrichment
        f.write("PROXIMITY ENRICHMENT (Forbidden variants near GWAS loci)\n")
        f.write("-" * 60 + "\n")
        for window_key, data in results.get('proximity', {}).items():
            f.write(f"\n{window_key}:\n")
            f.write(f"  Forbidden near GWAS: {data['n_forbidden_near_gwas']:,} / {data['total_forbidden']:,}\n")
            f.write(f"  Percentage: {data['pct_forbidden_near_gwas']:.2f}%\n")
            if 'enrichment' in data:
                f.write(f"  Enrichment ratio: {data['enrichment']['enrichment_ratio']:.2f}x\n")
                f.write(f"  P-value (enriched): {data['enrichment']['pvalue_enriched']:.4f}\n")
                f.write(f"  P-value (depleted): {data['enrichment']['pvalue_depleted']:.4f}\n")
        
        # Exact overlaps
        f.write("\n\nEXACT POSITION OVERLAPS\n")
        f.write("-" * 60 + "\n")
        for overlap in results.get('exact_overlaps', []):
            f.write(f"\n{overlap['target']}:\n")
            f.write(f"  Overlaps: {overlap['n_forbidden_overlap']:,}\n")
            f.write(f"  Percentage: {overlap['pct_forbidden_overlap']:.4f}%\n")
        
        # Trait categories
        if 'trait_enrichment' in results:
            f.write("\n\nTRAIT CATEGORY ENRICHMENT (100kb window)\n")
            f.write("-" * 60 + "\n")
            for _, row in results['trait_enrichment'].iterrows():
                f.write(f"  {row['category']}: {row['n_forbidden_nearby']:,} forbidden near {row['n_gwas_loci']:,} loci\n")
        
        f.write("\n" + "=" * 80 + "\n")
        f.write("Analysis complete.\n")


def main():
    args = parse_args()
    
    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load forbidden variants
    forbidden_df = load_forbidden_variants(args.forbidden)
    
    # Create temporary BED file for forbidden variants
    forbidden_bed = str(output_dir / 'forbidden_variants.bed')
    forbidden_to_bed(forbidden_df, forbidden_bed)
    
    results = {}
    
    # --- 1. Proximity Enrichment Analysis ---
    logger.info("\n=== Proximity Enrichment Analysis ===")
    
    gwas_bed = Path(args.gwas) / 'gwas_unique_loci.bed'
    if gwas_bed.exists():
        proximity_results = analyze_proximity_enrichment(
            forbidden_bed, 
            str(gwas_bed),
            windows=[50000, 100000, 500000]
        )
        results['proximity'] = proximity_results
        
        # Generate null distribution for 100kb window
        if args.n_permutations > 0:
            logger.info("\nGenerating null distribution for statistical testing...")
            null_dist = generate_genome_null(
                forbidden_df, 
                n_permutations=args.n_permutations,
                gwas_bed=str(gwas_bed),
                window=100000
            )
            
            # Calculate enrichment statistics
            observed = proximity_results['window_100kb']['n_forbidden_near_gwas']
            enrichment_stats = calculate_enrichment_pvalue(observed, null_dist['null_counts'])
            results['proximity']['window_100kb']['enrichment'] = enrichment_stats
            results['proximity']['window_100kb']['null_distribution'] = null_dist
            
            logger.info(f"\nEnrichment analysis (100kb window):")
            logger.info(f"  Observed: {enrichment_stats['observed']:,}")
            logger.info(f"  Expected: {enrichment_stats['expected']:.1f}")
            logger.info(f"  Enrichment: {enrichment_stats['enrichment_ratio']:.2f}x")
            logger.info(f"  P-value (enriched): {enrichment_stats['pvalue_enriched']:.4f}")
            logger.info(f"  P-value (depleted): {enrichment_stats['pvalue_depleted']:.4f}")
    
    # --- 2. Exact Overlap Analysis ---
    logger.info("\n=== Exact Overlap Analysis ===")
    
    exact_overlaps = []
    
    # GWAS exact overlap
    if gwas_bed.exists():
        overlap = analyze_exact_overlap(forbidden_bed, str(gwas_bed), "GWAS_lead_SNPs")
        exact_overlaps.append(overlap)
    
    # ClinVar overlaps
    clinvar_files = [
        ('clinvar_all.bed', 'ClinVar_all'),
        ('clinvar_pathogenic.bed', 'ClinVar_pathogenic'),
        ('clinvar_pathogenic_snvs.bed', 'ClinVar_pathogenic_SNVs'),
    ]
    
    for filename, label in clinvar_files:
        clinvar_bed = Path(args.clinvar) / filename
        if clinvar_bed.exists():
            overlap = analyze_exact_overlap(forbidden_bed, str(clinvar_bed), label)
            exact_overlaps.append(overlap)
    
    results['exact_overlaps'] = exact_overlaps
    
    # --- 3. Trait Category Analysis ---
    logger.info("\n=== Trait Category Analysis ===")
    
    if Path(args.gwas).exists():
        trait_df = analyze_by_trait_category(forbidden_bed, args.gwas)
        results['trait_enrichment'] = trait_df
        trait_df.to_csv(output_dir / 'trait_enrichment.tsv', sep='\t', index=False)
    
    # --- Save Results ---
    
    # Proximity results
    proximity_df = pd.DataFrame([
        {**{'window': k}, **{kk: vv for kk, vv in v.items() if not isinstance(vv, dict)}}
        for k, v in results.get('proximity', {}).items()
    ])
    proximity_df.to_csv(output_dir / 'proximity_enrichment.tsv', sep='\t', index=False)
    
    # Exact overlaps
    overlaps_df = pd.DataFrame(exact_overlaps)
    overlaps_df.to_csv(output_dir / 'exact_overlaps.tsv', sep='\t', index=False)
    
    # Summary report
    write_report(results, output_dir / 'validation_report.txt')
    
    logger.info(f"\nResults written to {output_dir}")
    logger.info("Analysis complete!")


if __name__ == '__main__':
    main()
