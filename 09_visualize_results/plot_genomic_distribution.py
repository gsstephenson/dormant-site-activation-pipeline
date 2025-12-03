#!/usr/bin/env python3
"""
Module 09: Genomic Distribution Visualization

Multi-panel figure showing WHERE forbidden variants are located:
- Panel A: Chromosome distribution (bar plot)
- Panel B: Distance to nearest gene/TSS
- Panel C: Genomic annotation enrichment (if available)
- Panel D: Constraint metric relationship (if available)
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns

# Style settings
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['axes.labelsize'] = 11

# Chromosome order for plotting
CHROM_ORDER = [str(i) for i in range(1, 23)] + ['X', 'Y']


def load_data(results_dir: Path):
    """Load forbidden variants data."""
    predictions_path = results_dir / 'forbidden_variants' / 'AP1' / 'predictions_summary_ap1.tsv'
    df = pd.read_csv(predictions_path, sep='\t')
    
    # Parse chromosome from variant_id
    df['chr'] = df['variant_id_str'].apply(lambda x: x.split(':')[0].replace('chr', ''))
    df['pos'] = df['variant_id_str'].apply(lambda x: int(x.split(':')[1]))
    
    return df


def panel_a_chromosome_distribution(ax, df):
    """Panel A: Chromosome distribution."""
    # Count per chromosome
    chrom_counts = df['chr'].value_counts()
    
    # Reorder by chromosome order
    ordered_chroms = [c for c in CHROM_ORDER if c in chrom_counts.index]
    counts = [chrom_counts.get(c, 0) for c in ordered_chroms]
    
    # Expected counts based on chromosome size (approximate)
    # Using relative lengths as proxy for expected counts
    chrom_sizes = {
        '1': 248.96, '2': 242.19, '3': 198.30, '4': 190.21, '5': 181.54,
        '6': 170.81, '7': 159.35, '8': 145.14, '9': 138.39, '10': 133.80,
        '11': 135.09, '12': 133.28, '13': 114.36, '14': 107.04, '15': 101.99,
        '16': 90.34, '17': 83.26, '18': 80.37, '19': 58.62, '20': 64.44,
        '21': 46.71, '22': 50.82, 'X': 156.04, 'Y': 57.23
    }
    
    total_size = sum(chrom_sizes.get(c, 0) for c in ordered_chroms)
    expected = [len(df) * chrom_sizes.get(c, 0) / total_size for c in ordered_chroms]
    
    x = np.arange(len(ordered_chroms))
    width = 0.35
    
    bars1 = ax.bar(x - width/2, counts, width, label='Observed', color='#e74c3c', edgecolor='black', linewidth=0.5)
    bars2 = ax.bar(x + width/2, expected, width, label='Expected\n(by size)', color='#95a5a6', edgecolor='black', linewidth=0.5, alpha=0.7)
    
    ax.set_xticks(x)
    ax.set_xticklabels(ordered_chroms, fontsize=8)
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('Variant Count')
    ax.set_title('A. Chromosome Distribution', fontweight='bold', loc='left')
    ax.legend(loc='upper right', fontsize=9)
    
    # Calculate chi-square statistic
    from scipy import stats
    chi2, p = stats.chisquare(counts, f_exp=expected)
    ax.text(0.02, 0.98, f'χ² = {chi2:.1f}\np = {p:.2e}', transform=ax.transAxes,
            fontsize=9, va='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))


def panel_b_variant_density(ax, df):
    """Panel B: Variant density across chromosomes (Manhattan-style)."""
    # Sort by chromosome and position
    df_sorted = df.copy()
    df_sorted['chr_num'] = df_sorted['chr'].apply(lambda x: int(x) if x.isdigit() else (23 if x == 'X' else 24))
    df_sorted = df_sorted.sort_values(['chr_num', 'pos'])
    
    # Calculate cumulative position for Manhattan plot
    cumulative_pos = 0
    chrom_starts = {}
    df_sorted['cumulative_pos'] = 0
    
    for chrom in CHROM_ORDER:
        if chrom in df_sorted['chr'].values:
            chrom_starts[chrom] = cumulative_pos
            mask = df_sorted['chr'] == chrom
            df_sorted.loc[mask, 'cumulative_pos'] = df_sorted.loc[mask, 'pos'] + cumulative_pos
            cumulative_pos = df_sorted.loc[mask, 'cumulative_pos'].max() + 1e7
    
    # Color by chromosome (alternating)
    colors = ['#1f77b4', '#aec7e8'] * 12
    
    for i, chrom in enumerate(CHROM_ORDER):
        if chrom in df_sorted['chr'].values:
            mask = df_sorted['chr'] == chrom
            ax.scatter(df_sorted.loc[mask, 'cumulative_pos'], 
                      df_sorted.loc[mask, 'ap1_raw_max'],
                      c=colors[i], s=2, alpha=0.5, rasterized=True)
    
    # Add chromosome labels
    chrom_centers = {}
    for chrom in CHROM_ORDER:
        if chrom in df_sorted['chr'].values:
            mask = df_sorted['chr'] == chrom
            chrom_centers[chrom] = df_sorted.loc[mask, 'cumulative_pos'].median()
    
    ax.set_xticks([chrom_centers[c] for c in CHROM_ORDER if c in chrom_centers])
    ax.set_xticklabels([c for c in CHROM_ORDER if c in chrom_centers], fontsize=7, rotation=45)
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('AP1 Raw Score')
    ax.set_title('B. Genomic Distribution of Forbidden Variants', fontweight='bold', loc='left')
    
    # Add horizontal line for median
    median_score = df_sorted['ap1_raw_max'].median()
    ax.axhline(y=median_score, color='red', linestyle='--', alpha=0.7, label=f'Median: {median_score:.0f}')
    ax.legend(loc='upper right', fontsize=9)


def panel_c_score_by_chromosome(ax, df):
    """Panel C: AP1 score distribution by chromosome."""
    # Get top 10 chromosomes by variant count
    top_chroms = df['chr'].value_counts().head(10).index.tolist()
    
    # Filter to top chromosomes for cleaner visualization
    df_top = df[df['chr'].isin(top_chroms)]
    
    # Create box plot
    order = sorted(top_chroms, key=lambda x: int(x) if x.isdigit() else (23 if x == 'X' else 24))
    
    bp = ax.boxplot([df_top[df_top['chr'] == c]['ap1_raw_max'].values for c in order],
                     labels=order, patch_artist=True, showfliers=False)
    
    for patch in bp['boxes']:
        patch.set_facecolor('#3498db')
        patch.set_alpha(0.7)
    
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('AP1 Raw Score')
    ax.set_title('C. AP1 Score by Chromosome (Top 10)', fontweight='bold', loc='left')
    
    # Add sample sizes
    for i, c in enumerate(order):
        n = (df_top['chr'] == c).sum()
        ax.text(i + 1, ax.get_ylim()[1] * 0.95, f'n={n:,}', ha='center', fontsize=7)


def panel_d_quantile_distribution(ax, df):
    """Panel D: AP1 quantile score distribution."""
    quantiles = df['ap1_quantile_max'].dropna().values
    
    # Histogram with density
    ax.hist(quantiles, bins=50, density=True, color='#9b59b6', alpha=0.7, edgecolor='black', linewidth=0.5)
    
    # Add vertical lines for key percentiles
    p99 = 0.99
    p999 = 0.999
    p9999 = 0.9999
    
    ax.axvline(x=p99, color='orange', linestyle='--', linewidth=2, label=f'>99th: {(quantiles > p99).sum():,} ({(quantiles > p99).mean()*100:.1f}%)')
    ax.axvline(x=p999, color='red', linestyle='--', linewidth=2, label=f'>99.9th: {(quantiles > p999).sum():,} ({(quantiles > p999).mean()*100:.1f}%)')
    ax.axvline(x=p9999, color='darkred', linestyle='--', linewidth=2, label=f'>99.99th: {(quantiles > p9999).sum():,} ({(quantiles > p9999).mean()*100:.1f}%)')
    
    ax.set_xlabel('AP1 Quantile Score (max)')
    ax.set_ylabel('Density')
    ax.set_title('D. AP1 Quantile Distribution', fontweight='bold', loc='left')
    ax.legend(loc='upper left', fontsize=8)
    
    # Add annotation for extreme scores
    ax.text(0.98, 0.95, f'Median: {np.median(quantiles):.4f}\nMean: {np.mean(quantiles):.4f}',
            transform=ax.transAxes, fontsize=9, va='top', ha='right',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))


def main():
    parser = argparse.ArgumentParser(description='Generate genomic distribution visualization')
    parser.add_argument('--results-dir', type=str, default='results',
                        help='Path to results directory')
    parser.add_argument('--output-dir', type=str, default='figures/forbidden_variants',
                        help='Output directory for figures')
    parser.add_argument('--dpi', type=int, default=300, help='Figure DPI')
    args = parser.parse_args()
    
    results_dir = Path(args.results_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("Loading data...")
    df = load_data(results_dir)
    print(f"  Loaded {len(df):,} forbidden variants")
    print(f"  Chromosomes: {df['chr'].nunique()}")
    
    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    plt.subplots_adjust(hspace=0.3, wspace=0.3)
    
    print("Generating panels...")
    panel_a_chromosome_distribution(axes[0, 0], df)
    panel_b_variant_density(axes[0, 1], df)
    panel_c_score_by_chromosome(axes[1, 0], df)
    panel_d_quantile_distribution(axes[1, 1], df)
    
    # Overall title
    fig.suptitle('Genomic Distribution of Forbidden Variants', 
                 fontsize=16, fontweight='bold', y=0.98)
    
    # Save
    png_path = output_dir / 'genomic_distribution.png'
    pdf_path = output_dir / 'genomic_distribution.pdf'
    
    fig.savefig(png_path, dpi=args.dpi, bbox_inches='tight', facecolor='white')
    fig.savefig(pdf_path, bbox_inches='tight', facecolor='white')
    
    print(f"\nSaved: {png_path}")
    print(f"Saved: {pdf_path}")
    
    plt.close(fig)


if __name__ == '__main__':
    main()
