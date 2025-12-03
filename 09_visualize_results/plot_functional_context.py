#!/usr/bin/env python3
"""
Module 09: Functional Context Visualization

Multi-panel figure showing functional properties of forbidden variants:
- Panel A: Enhancer activity distribution
- Panel B: Accessibility (ATAC/DNase) distribution
- Panel C: Multi-feature scatter (accessibility vs enhancer, colored by AP1)
- Panel D: High-priority candidates with multi-track evidence
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize
from matplotlib import cm
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


def panel_a_enhancer_distribution(ax, df):
    """Panel A: Enhancer activity distribution."""
    enhancer_scores = df['enhancer_quantile_max'].dropna().values
    
    # Histogram with KDE
    ax.hist(enhancer_scores, bins=50, density=True, color='#27ae60', alpha=0.7, 
            edgecolor='black', linewidth=0.5, label='Histogram')
    
    # Add vertical lines for percentiles
    p50 = np.median(enhancer_scores)
    p95 = np.percentile(enhancer_scores, 95)
    
    ax.axvline(x=p50, color='blue', linestyle='--', linewidth=2, label=f'Median: {p50:.3f}')
    ax.axvline(x=p95, color='red', linestyle='--', linewidth=2, label=f'95th: {p95:.3f}')
    
    ax.set_xlabel('Enhancer Quantile Score (max)')
    ax.set_ylabel('Density')
    ax.set_title('A. Enhancer Activity (H3K27ac)', fontweight='bold', loc='left')
    ax.legend(loc='upper left', fontsize=9)
    
    # Add high-enhancer percentage
    high_enh_pct = (enhancer_scores > 0.95).mean() * 100
    ax.text(0.98, 0.95, f'{high_enh_pct:.1f}% in top 5%\nof enhancer activity',
            transform=ax.transAxes, fontsize=9, va='top', ha='right',
            bbox=dict(boxstyle='round', facecolor='#e8f8f5', edgecolor='#27ae60'))


def panel_b_accessibility_distribution(ax, df):
    """Panel B: Accessibility distribution."""
    access_scores = df['accessibility_quantile_max'].dropna().values
    
    # Histogram
    ax.hist(access_scores, bins=50, density=True, color='#3498db', alpha=0.7,
            edgecolor='black', linewidth=0.5)
    
    # Add vertical lines for percentiles
    p50 = np.median(access_scores)
    p95 = np.percentile(access_scores, 95)
    
    ax.axvline(x=p50, color='blue', linestyle='--', linewidth=2, label=f'Median: {p50:.3f}')
    ax.axvline(x=p95, color='red', linestyle='--', linewidth=2, label=f'95th: {p95:.3f}')
    
    ax.set_xlabel('Accessibility Quantile Score (max)')
    ax.set_ylabel('Density')
    ax.set_title('B. Chromatin Accessibility (ATAC/DNase)', fontweight='bold', loc='left')
    ax.legend(loc='upper left', fontsize=9)
    
    # Add high-accessibility percentage
    high_acc_pct = (access_scores > 0.95).mean() * 100
    ax.text(0.98, 0.95, f'{high_acc_pct:.1f}% in top 5%\nof accessibility',
            transform=ax.transAxes, fontsize=9, va='top', ha='right',
            bbox=dict(boxstyle='round', facecolor='#ebf5fb', edgecolor='#3498db'))


def panel_c_multifeature_scatter(ax, df):
    """Panel C: Accessibility vs Enhancer, colored by AP1 score."""
    # Sample for performance
    df_sample = df.sample(n=min(10000, len(df)), random_state=42)
    
    x = df_sample['accessibility_quantile_max'].values
    y = df_sample['enhancer_quantile_max'].values
    c = df_sample['ap1_quantile_max'].values
    
    # Remove NaN
    mask = ~(np.isnan(x) | np.isnan(y) | np.isnan(c))
    x, y, c = x[mask], y[mask], c[mask]
    
    # Scatter with color
    scatter = ax.scatter(x, y, c=c, cmap='YlOrRd', alpha=0.5, s=10, 
                         norm=Normalize(vmin=0.9, vmax=1.0), rasterized=True)
    
    plt.colorbar(scatter, ax=ax, label='AP1 Quantile Score', shrink=0.8)
    
    ax.set_xlabel('Accessibility Quantile')
    ax.set_ylabel('Enhancer Quantile')
    ax.set_title('C. Multi-Feature Landscape', fontweight='bold', loc='left')
    
    # Mark high-impact quadrant
    ax.axvline(x=0.95, color='gray', linestyle='--', alpha=0.5)
    ax.axhline(y=0.95, color='gray', linestyle='--', alpha=0.5)
    
    # Count in each quadrant
    high_both = ((x > 0.95) & (y > 0.95)).sum()
    ax.text(0.97, 0.97, f'High Acc +\nHigh Enh:\nn={high_both}',
            transform=ax.transAxes, fontsize=8, ha='right', va='top',
            bbox=dict(boxstyle='round', facecolor='#ffcccc', alpha=0.9, edgecolor='red'))


def panel_d_high_priority(ax, df):
    """Panel D: High-priority candidates with multi-track evidence."""
    ax.axis('off')
    
    # Define high-priority as: high AP1 + high enhancer + high accessibility
    df_filtered = df.dropna(subset=['ap1_quantile_max', 'enhancer_quantile_max', 'accessibility_quantile_max'])
    
    high_priority = df_filtered[
        (df_filtered['ap1_quantile_max'] > 0.99) &
        (df_filtered['enhancer_quantile_max'] > 0.95) &
        (df_filtered['accessibility_quantile_max'] > 0.95)
    ].nlargest(10, 'ap1_raw_max')
    
    if len(high_priority) == 0:
        ax.text(0.5, 0.5, 'No variants meet all criteria:\n- AP1 quantile > 0.99\n- Enhancer quantile > 0.95\n- Accessibility quantile > 0.95',
                ha='center', va='center', fontsize=10, transform=ax.transAxes)
        ax.set_title('D. High-Priority Multi-Track Candidates', fontweight='bold', loc='center', pad=20)
        return
    
    # Create table
    cell_text = []
    for _, row in high_priority.iterrows():
        vid = row['variant_id_str']
        parts = vid.split(':')
        short_vid = f"{parts[0]}:{parts[1][:8]}..."
        cell_text.append([
            short_vid,
            f"{row['ap1_quantile_max']:.4f}",
            f"{row['enhancer_quantile_max']:.4f}",
            f"{row['accessibility_quantile_max']:.4f}",
            row['ap1_best_tf']
        ])
    
    table = ax.table(
        cellText=cell_text,
        colLabels=['Variant', 'AP1 Q', 'Enh Q', 'Acc Q', 'Best TF'],
        loc='center',
        cellLoc='center',
        colColours=['#e8e8e8'] * 5
    )
    
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1.2, 1.5)
    
    ax.set_title(f'D. Top {len(high_priority)} Multi-Track High-Priority Candidates', 
                 fontweight='bold', loc='center', pad=20)
    
    # Add criteria note
    ax.text(0.5, -0.05, 'Criteria: AP1 > 99th, Enhancer > 95th, Accessibility > 95th percentile',
            ha='center', va='top', fontsize=8, style='italic', transform=ax.transAxes)


def main():
    parser = argparse.ArgumentParser(description='Generate functional context visualization')
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
    panel_a_enhancer_distribution(axes[0, 0], df)
    panel_b_accessibility_distribution(axes[0, 1], df)
    panel_c_multifeature_scatter(axes[1, 0], df)
    panel_d_high_priority(axes[1, 1], df)
    
    # Overall title
    fig.suptitle('Functional Context of Forbidden Variants', 
                 fontsize=16, fontweight='bold', y=0.98)
    
    # Save
    png_path = output_dir / 'functional_context.png'
    pdf_path = output_dir / 'functional_context.pdf'
    
    fig.savefig(png_path, dpi=args.dpi, bbox_inches='tight', facecolor='white')
    fig.savefig(pdf_path, bbox_inches='tight', facecolor='white')
    
    print(f"\nSaved: {png_path}")
    print(f"Saved: {pdf_path}")
    
    plt.close(fig)


if __name__ == '__main__':
    main()
