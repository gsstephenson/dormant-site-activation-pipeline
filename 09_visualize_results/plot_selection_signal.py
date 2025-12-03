#!/usr/bin/env python3
"""
Module 09: Selection Signal Strength Visualization

Multi-panel figure characterizing the INTENSITY of purifying selection:
- Panel A: AP1 quantile distribution (showing extreme scores)
- Panel B: Allele number (AN) distribution showing coverage quality
- Panel C: Selection signal by score quantile
- Panel D: Top candidates with highest predicted impact
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

# Style settings
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['axes.labelsize'] = 11


def load_data(results_dir: Path):
    """Load forbidden variants data."""
    predictions_path = results_dir / 'forbidden_variants' / 'AP1' / 'predictions_summary_ap1.tsv'
    df = pd.read_csv(predictions_path, sep='\t')
    return df


def panel_a_quantile_extremity(ax, df):
    """Panel A: Show how extreme the AP1 quantile scores are."""
    quantiles = df['ap1_quantile_max'].dropna().values
    
    # Calculate percentages above various thresholds
    thresholds = [0.9, 0.95, 0.99, 0.999, 0.9999]
    percentages = [(quantiles > t).mean() * 100 for t in thresholds]
    counts = [(quantiles > t).sum() for t in thresholds]
    
    # Bar plot
    x = np.arange(len(thresholds))
    bars = ax.bar(x, percentages, color='#8e44ad', edgecolor='black', linewidth=1)
    
    ax.set_xticks(x)
    ax.set_xticklabels([f'>{t}' for t in thresholds])
    ax.set_xlabel('Quantile Threshold')
    ax.set_ylabel('Percentage of Variants (%)')
    ax.set_title('A. Proportion with Extreme AP1 Scores', fontweight='bold', loc='left')
    
    # Add count labels
    for bar, pct, count in zip(bars, percentages, counts):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
               f'{pct:.1f}%\n(n={count:,})', ha='center', va='bottom', fontsize=9)
    
    ax.set_ylim(0, max(percentages) * 1.25)
    
    # Add interpretation
    ax.text(0.98, 0.95, f'{(quantiles > 0.999).mean()*100:.1f}% of forbidden\nvariants exceed 99.9th\npercentile of all variants',
            transform=ax.transAxes, fontsize=9, va='top', ha='right',
            bbox=dict(boxstyle='round', facecolor='#f5eef8', edgecolor='#8e44ad'))


def panel_b_allele_number(ax, df):
    """Panel B: Allele number distribution showing coverage quality."""
    an_values = df['AN'].dropna().values
    
    # Histogram
    ax.hist(an_values, bins=50, color='#2980b9', alpha=0.7, edgecolor='black', linewidth=0.5)
    
    # Add vertical lines for key values
    median_an = np.median(an_values)
    max_an = np.max(an_values)
    
    ax.axvline(x=median_an, color='red', linestyle='--', linewidth=2, 
               label=f'Median: {median_an:,.0f}')
    ax.axvline(x=max_an, color='darkred', linestyle=':', linewidth=2,
               label=f'Max: {max_an:,.0f}')
    
    ax.set_xlabel('Allele Number (AN)')
    ax.set_ylabel('Count')
    ax.set_title('B. Sequencing Coverage Quality', fontweight='bold', loc='left')
    ax.legend(loc='upper left', fontsize=9)
    
    # Calculate what fraction has high coverage
    high_cov_thresh = 140000  # ~70K individuals * 2 alleles
    high_cov_pct = (an_values > high_cov_thresh).mean() * 100
    
    ax.text(0.98, 0.95, f'{high_cov_pct:.1f}% have AN > {high_cov_thresh:,}\n(high coverage)',
            transform=ax.transAxes, fontsize=9, va='top', ha='right',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))


def panel_c_score_distribution(ax, df):
    """Panel C: AP1 raw score distribution with log scale."""
    scores = df['ap1_raw_max'].dropna().values
    
    # Log-transform for better visualization
    log_scores = np.log10(scores + 1)
    
    # Histogram
    ax.hist(log_scores, bins=50, color='#e74c3c', alpha=0.7, edgecolor='black', linewidth=0.5)
    
    ax.set_xlabel('log₁₀(AP1 Raw Score + 1)')
    ax.set_ylabel('Count')
    ax.set_title('C. AP1 Score Distribution (log scale)', fontweight='bold', loc='left')
    
    # Add statistics
    ax.axvline(x=np.log10(np.median(scores) + 1), color='blue', linestyle='--', linewidth=2,
               label=f'Median: {np.median(scores):,.0f}')
    ax.axvline(x=np.log10(np.mean(scores) + 1), color='green', linestyle='--', linewidth=2,
               label=f'Mean: {np.mean(scores):,.0f}')
    
    ax.legend(loc='upper right', fontsize=9)
    
    # Add percentile annotations
    p95 = np.percentile(scores, 95)
    p99 = np.percentile(scores, 99)
    ax.text(0.98, 0.6, f'95th percentile: {p95:,.0f}\n99th percentile: {p99:,.0f}',
            transform=ax.transAxes, fontsize=9, va='top', ha='right',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))


def panel_d_top_candidates(ax, df):
    """Panel D: Top candidates table."""
    ax.axis('off')
    
    # Get top 10 by AP1 score
    top10 = df.nlargest(10, 'ap1_raw_max')[['variant_id_str', 'ap1_raw_max', 'ap1_best_tf', 'ap1_best_biosample', 'AN']]
    
    # Create table
    cell_text = []
    for _, row in top10.iterrows():
        vid = row['variant_id_str']
        # Shorten variant ID
        parts = vid.split(':')
        short_vid = f"{parts[0]}:{parts[1][:8]}..."
        cell_text.append([
            short_vid,
            f"{row['ap1_raw_max']:,.0f}",
            row['ap1_best_tf'],
            row['ap1_best_biosample'],
            f"{row['AN']:,.0f}"
        ])
    
    table = ax.table(
        cellText=cell_text,
        colLabels=['Variant', 'AP1 Score', 'Best TF', 'Best Tissue', 'AN'],
        loc='center',
        cellLoc='center',
        colColours=['#e8e8e8'] * 5
    )
    
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1.2, 1.5)
    
    ax.set_title('D. Top 10 Highest-Impact Forbidden Variants', fontweight='bold', loc='center', pad=20)


def main():
    parser = argparse.ArgumentParser(description='Generate selection signal strength visualization')
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
    
    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    plt.subplots_adjust(hspace=0.35, wspace=0.3)
    
    print("Generating panels...")
    panel_a_quantile_extremity(axes[0, 0], df)
    panel_b_allele_number(axes[0, 1], df)
    panel_c_score_distribution(axes[1, 0], df)
    panel_d_top_candidates(axes[1, 1], df)
    
    # Overall title
    fig.suptitle('Selection Signal Strength: Evidence for Strong Purifying Selection', 
                 fontsize=16, fontweight='bold', y=0.98)
    
    # Save
    png_path = output_dir / 'selection_signal_strength.png'
    pdf_path = output_dir / 'selection_signal_strength.pdf'
    
    fig.savefig(png_path, dpi=args.dpi, bbox_inches='tight', facecolor='white')
    fig.savefig(pdf_path, bbox_inches='tight', facecolor='white')
    
    print(f"\nSaved: {png_path}")
    print(f"Saved: {pdf_path}")
    
    plt.close(fig)


if __name__ == '__main__':
    main()
