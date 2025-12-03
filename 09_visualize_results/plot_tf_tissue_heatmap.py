#!/usr/bin/env python3
"""
Module 09: TF-Tissue Heatmap Visualization (Improved)

Two-panel heatmap figure showing TF-tissue combinations:
- Panel A: Median AP1 score by TF-tissue (filtered n≥10)
- Panel B: Variant counts by TF-tissue (filtered n≥10)

Improvements over original:
- Filter to only cells with n≥10 for reliable statistics
- Add hierarchical clustering for rows and columns
- Focus on top 10 TFs and tissues for cleaner visualization
- Add marginal totals
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LogNorm
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist

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


def create_pivot_tables(df, min_count=10, top_n_tfs=10, top_n_tissues=10):
    """Create pivot tables for median score and counts with filtering."""
    # Get top TFs and tissues by count
    top_tfs = df['ap1_best_tf'].value_counts().head(top_n_tfs).index.tolist()
    top_tissues = df['ap1_best_biosample'].value_counts().head(top_n_tissues).index.tolist()
    
    # Filter data to top TFs and tissues
    df_filtered = df[df['ap1_best_tf'].isin(top_tfs) & df['ap1_best_biosample'].isin(top_tissues)]
    
    # Create count pivot
    count_pivot = df_filtered.pivot_table(
        index='ap1_best_tf',
        columns='ap1_best_biosample',
        values='ap1_raw_max',
        aggfunc='count',
        fill_value=0
    )
    
    # Create median score pivot (only where count >= min_count)
    median_pivot = df_filtered.pivot_table(
        index='ap1_best_tf',
        columns='ap1_best_biosample',
        values='ap1_raw_max',
        aggfunc='median',
        fill_value=np.nan
    )
    
    # Mask low-count cells
    median_pivot_masked = median_pivot.copy()
    median_pivot_masked[count_pivot < min_count] = np.nan
    
    return median_pivot_masked, count_pivot, top_tfs, top_tissues


def cluster_and_reorder(pivot_df):
    """Apply hierarchical clustering to reorder rows and columns."""
    # Handle NaN for clustering
    filled = pivot_df.fillna(0)
    
    # Cluster rows
    if len(filled) > 1:
        row_linkage = linkage(pdist(filled.values), method='ward')
        row_order = dendrogram(row_linkage, no_plot=True)['leaves']
        row_labels = [filled.index[i] for i in row_order]
    else:
        row_labels = list(filled.index)
    
    # Cluster columns
    if len(filled.columns) > 1:
        col_linkage = linkage(pdist(filled.values.T), method='ward')
        col_order = dendrogram(col_linkage, no_plot=True)['leaves']
        col_labels = [filled.columns[i] for i in col_order]
    else:
        col_labels = list(filled.columns)
    
    return row_labels, col_labels


def main():
    parser = argparse.ArgumentParser(description='Generate TF-tissue heatmap')
    parser.add_argument('--results-dir', type=str, default='results',
                        help='Path to results directory')
    parser.add_argument('--output-dir', type=str, default='figures/forbidden_variants',
                        help='Output directory for figures')
    parser.add_argument('--min-count', type=int, default=10,
                        help='Minimum count for reliable statistics')
    parser.add_argument('--top-n-tfs', type=int, default=10,
                        help='Number of top TFs to show')
    parser.add_argument('--top-n-tissues', type=int, default=10,
                        help='Number of top tissues to show')
    parser.add_argument('--dpi', type=int, default=300, help='Figure DPI')
    args = parser.parse_args()
    
    results_dir = Path(args.results_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("Loading data...")
    df = load_data(results_dir)
    print(f"  Loaded {len(df):,} forbidden variants")
    
    print(f"Creating pivot tables (min_count={args.min_count})...")
    median_pivot, count_pivot, top_tfs, top_tissues = create_pivot_tables(
        df, min_count=args.min_count, top_n_tfs=args.top_n_tfs, top_n_tissues=args.top_n_tissues
    )
    
    # Cluster for better visualization
    print("Applying hierarchical clustering...")
    row_order, col_order = cluster_and_reorder(median_pivot)
    median_pivot = median_pivot.loc[row_order, col_order]
    count_pivot = count_pivot.loc[row_order, col_order]
    
    # Create figure
    fig, axes = plt.subplots(1, 2, figsize=(16, 8))
    plt.subplots_adjust(wspace=0.3)
    
    # Panel A: Median scores heatmap
    print("Generating Panel A: Median scores...")
    
    # Create mask for low-count cells
    mask = count_pivot < args.min_count
    
    # Draw heatmap
    sns.heatmap(median_pivot, ax=axes[0], cmap='YlOrRd', 
                annot=True, fmt='.0f', 
                mask=mask,
                cbar_kws={'label': 'Median AP1 Score', 'shrink': 0.8},
                linewidths=0.5, linecolor='white')
    
    axes[0].set_title(f'A. Median AP1 Score by TF-Tissue\n(cells with n≥{args.min_count} only)', 
                      fontweight='bold', loc='left')
    axes[0].set_xlabel('Tissue/Cell Type')
    axes[0].set_ylabel('Transcription Factor')
    
    # Rotate x labels
    axes[0].set_xticklabels(axes[0].get_xticklabels(), rotation=45, ha='right')
    
    # Add hatching for masked cells
    for i in range(len(row_order)):
        for j in range(len(col_order)):
            if mask.iloc[i, j]:
                axes[0].add_patch(plt.Rectangle((j, i), 1, 1, fill=False, 
                                                 hatch='///', edgecolor='gray', alpha=0.3))
    
    # Panel B: Counts heatmap
    print("Generating Panel B: Variant counts...")
    
    # Use log scale for counts
    count_for_plot = count_pivot.replace(0, np.nan)
    
    sns.heatmap(count_for_plot, ax=axes[1], cmap='Blues',
                annot=True, fmt='.0f',
                cbar_kws={'label': 'Variant Count', 'shrink': 0.8},
                linewidths=0.5, linecolor='white',
                norm=LogNorm(vmin=1, vmax=count_pivot.max().max()))
    
    axes[1].set_title('B. Variant Count by TF-Tissue\n(log color scale)', 
                      fontweight='bold', loc='left')
    axes[1].set_xlabel('Tissue/Cell Type')
    axes[1].set_ylabel('Transcription Factor')
    
    # Rotate x labels
    axes[1].set_xticklabels(axes[1].get_xticklabels(), rotation=45, ha='right')
    
    # Add marginal totals as annotations
    # Row totals
    row_totals = count_pivot.sum(axis=1)
    for i, (tf, total) in enumerate(row_totals.items()):
        axes[1].text(len(col_order) + 0.3, i + 0.5, f'Σ={total:,}', 
                     va='center', fontsize=8, color='darkblue')
    
    # Column totals
    col_totals = count_pivot.sum(axis=0)
    for j, (tissue, total) in enumerate(col_totals.items()):
        axes[1].text(j + 0.5, len(row_order) + 0.3, f'{total:,}', 
                     ha='center', va='top', fontsize=7, rotation=90, color='darkblue')
    
    # Overall title
    fig.suptitle(f'TF-Tissue Specificity of Forbidden Variants (Top {args.top_n_tfs} TFs × {args.top_n_tissues} Tissues)', 
                 fontsize=14, fontweight='bold', y=1.02)
    
    # Add legend for hatching
    hatch_patch = mpatches.Patch(facecolor='white', hatch='///', edgecolor='gray', 
                                  label=f'n < {args.min_count} (unreliable)')
    fig.legend(handles=[hatch_patch], loc='lower center', ncol=1, fontsize=9,
               bbox_to_anchor=(0.25, -0.02))
    
    # Save
    png_path = output_dir / 'tf_tissue_heatmap.png'
    pdf_path = output_dir / 'tf_tissue_heatmap.pdf'
    
    fig.savefig(png_path, dpi=args.dpi, bbox_inches='tight', facecolor='white')
    fig.savefig(pdf_path, bbox_inches='tight', facecolor='white')
    
    print(f"\nSaved: {png_path}")
    print(f"Saved: {pdf_path}")
    
    # Print summary statistics
    print("\n--- Summary ---")
    print(f"Total TF-tissue combinations: {len(row_order) * len(col_order)}")
    print(f"Cells with n≥{args.min_count}: {(~mask).sum().sum()}")
    print(f"Cells with n<{args.min_count} (masked): {mask.sum().sum()}")
    
    plt.close(fig)


if __name__ == '__main__':
    main()
