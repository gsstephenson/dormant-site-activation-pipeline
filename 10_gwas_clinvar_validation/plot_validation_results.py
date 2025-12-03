#!/usr/bin/env python3
"""
Generate visualization figures for GWAS/ClinVar validation results.

Creates multi-panel figures showing:
1. Proximity enrichment by window size
2. Trait category enrichment
3. Null distribution with observed value
4. Summary of validation findings
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Set style
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")


def parse_args():
    parser = argparse.ArgumentParser(description='Plot validation results')
    parser.add_argument('--input', required=True, help='Input directory with analysis results')
    parser.add_argument('--output', required=True, help='Output directory for figures')
    return parser.parse_args()


def plot_proximity_enrichment(proximity_df: pd.DataFrame, output_path: str):
    """Plot proximity enrichment by window size."""
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Panel A: Count by window size
    ax = axes[0]
    windows = proximity_df['window'].str.extract(r'(\d+)kb').astype(int).values.flatten()
    counts = proximity_df['n_forbidden_near_gwas'].values
    
    ax.bar(range(len(windows)), counts, color='steelblue', edgecolor='black')
    ax.set_xticks(range(len(windows)))
    ax.set_xticklabels([f'{w}kb' for w in windows])
    ax.set_xlabel('Window Size')
    ax.set_ylabel('Forbidden Variants Near GWAS')
    ax.set_title('A. Proximity to GWAS Loci')
    
    # Add count labels
    for i, (w, c) in enumerate(zip(windows, counts)):
        ax.text(i, c + 50, f'{c:,}', ha='center', fontsize=10)
    
    # Panel B: Percentage by window size
    ax = axes[1]
    pcts = proximity_df['pct_forbidden_near_gwas'].values
    
    ax.bar(range(len(windows)), pcts, color='coral', edgecolor='black')
    ax.set_xticks(range(len(windows)))
    ax.set_xticklabels([f'{w}kb' for w in windows])
    ax.set_xlabel('Window Size')
    ax.set_ylabel('% Forbidden Variants Near GWAS')
    ax.set_title('B. Percentage Near GWAS Loci')
    
    # Add percentage labels
    for i, (w, p) in enumerate(zip(windows, pcts)):
        ax.text(i, p + 0.5, f'{p:.1f}%', ha='center', fontsize=10)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_path.replace('.png', '.pdf'), bbox_inches='tight')
    plt.close()
    
    logger.info(f"Saved proximity enrichment figure to {output_path}")


def plot_trait_enrichment(trait_df: pd.DataFrame, output_path: str):
    """Plot enrichment by trait category."""
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Sort by number of forbidden variants nearby
    trait_df = trait_df.sort_values('n_forbidden_nearby', ascending=True)
    
    # Panel A: Horizontal bar chart
    ax = axes[0]
    colors = sns.color_palette("husl", len(trait_df))
    
    ax.barh(range(len(trait_df)), trait_df['n_forbidden_nearby'], color=colors, edgecolor='black')
    ax.set_yticks(range(len(trait_df)))
    ax.set_yticklabels(trait_df['category'])
    ax.set_xlabel('Forbidden Variants Within 100kb')
    ax.set_title('A. Forbidden Variants Near GWAS by Trait Category')
    
    # Panel B: Normalized by GWAS loci count
    ax = axes[1]
    
    # Variants per 1000 GWAS loci
    trait_df['rate_per_1000'] = 1000 * trait_df['n_forbidden_nearby'] / trait_df['n_gwas_loci']
    
    ax.barh(range(len(trait_df)), trait_df['rate_per_1000'], color=colors, edgecolor='black')
    ax.set_yticks(range(len(trait_df)))
    ax.set_yticklabels(trait_df['category'])
    ax.set_xlabel('Forbidden Variants per 1000 GWAS Loci')
    ax.set_title('B. Rate Normalized by GWAS Count')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_path.replace('.png', '.pdf'), bbox_inches='tight')
    plt.close()
    
    logger.info(f"Saved trait enrichment figure to {output_path}")


def plot_exact_overlaps(overlaps_df: pd.DataFrame, output_path: str):
    """Plot exact overlap results."""
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Create bar chart
    x = range(len(overlaps_df))
    colors = ['#e74c3c' if 'pathogenic' in t.lower() else '#3498db' 
              for t in overlaps_df['target']]
    
    bars = ax.bar(x, overlaps_df['n_forbidden_overlap'], color=colors, edgecolor='black')
    
    ax.set_xticks(x)
    ax.set_xticklabels(overlaps_df['target'], rotation=45, ha='right')
    ax.set_ylabel('Number of Exact Overlaps')
    ax.set_title('Exact Position Overlap: Forbidden Variants vs Disease Databases')
    
    # Add count labels
    for i, (bar, count) in enumerate(zip(bars, overlaps_df['n_forbidden_overlap'])):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1, 
                str(count), ha='center', fontsize=12, fontweight='bold')
    
    # Add annotation about expected result
    ax.text(0.98, 0.95, 
            'Expected: ~0 overlaps\n(forbidden variants don\'t exist in population)',
            transform=ax.transAxes, fontsize=10, verticalalignment='top',
            horizontalalignment='right', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_path.replace('.png', '.pdf'), bbox_inches='tight')
    plt.close()
    
    logger.info(f"Saved exact overlaps figure to {output_path}")


def plot_validation_summary(proximity_df: pd.DataFrame, overlaps_df: pd.DataFrame,
                            trait_df: pd.DataFrame, output_path: str):
    """Create 4-panel summary figure."""
    
    fig = plt.figure(figsize=(16, 12))
    
    # Panel A: Proximity by window
    ax1 = fig.add_subplot(2, 2, 1)
    windows = proximity_df['window'].str.extract(r'(\d+)kb').astype(int).values.flatten()
    pcts = proximity_df['pct_forbidden_near_gwas'].values
    
    ax1.bar(range(len(windows)), pcts, color='steelblue', edgecolor='black')
    ax1.set_xticks(range(len(windows)))
    ax1.set_xticklabels([f'{w}kb' for w in windows])
    ax1.set_xlabel('Window Size')
    ax1.set_ylabel('% Near GWAS')
    ax1.set_title('A. Proximity to GWAS Loci')
    
    for i, p in enumerate(pcts):
        ax1.text(i, p + 0.3, f'{p:.1f}%', ha='center', fontsize=9)
    
    # Panel B: Exact overlaps
    ax2 = fig.add_subplot(2, 2, 2)
    x = range(len(overlaps_df))
    colors = ['#e74c3c' if 'ClinVar' in t else '#3498db' for t in overlaps_df['target']]
    
    ax2.bar(x, overlaps_df['n_forbidden_overlap'], color=colors, edgecolor='black')
    ax2.set_xticks(x)
    ax2.set_xticklabels(overlaps_df['target'], rotation=45, ha='right', fontsize=8)
    ax2.set_ylabel('Exact Overlaps')
    ax2.set_title('B. Exact Position Overlaps')
    
    # Panel C: Top trait categories
    ax3 = fig.add_subplot(2, 2, 3)
    
    if trait_df is not None and len(trait_df) > 0:
        top_traits = trait_df.nlargest(8, 'n_forbidden_nearby')
        ax3.barh(range(len(top_traits)), top_traits['n_forbidden_nearby'], 
                 color=sns.color_palette("husl", len(top_traits)), edgecolor='black')
        ax3.set_yticks(range(len(top_traits)))
        ax3.set_yticklabels(top_traits['category'])
        ax3.set_xlabel('Forbidden Variants Nearby')
        ax3.set_title('C. Top Trait Categories (100kb window)')
    else:
        ax3.text(0.5, 0.5, 'No trait data available', ha='center', va='center')
        ax3.set_title('C. Trait Categories')
    
    # Panel D: Key findings text box
    ax4 = fig.add_subplot(2, 2, 4)
    ax4.axis('off')
    
    # Calculate summary statistics
    total_forbidden = proximity_df['total_forbidden'].iloc[0] if len(proximity_df) > 0 else 'N/A'
    pct_100kb = proximity_df[proximity_df['window'] == 'window_100kb']['pct_forbidden_near_gwas'].values
    pct_100kb = f'{pct_100kb[0]:.1f}%' if len(pct_100kb) > 0 else 'N/A'
    
    clinvar_overlap = overlaps_df[overlaps_df['target'].str.contains('ClinVar')]['n_forbidden_overlap'].sum()
    gwas_overlap = overlaps_df[overlaps_df['target'].str.contains('GWAS')]['n_forbidden_overlap'].sum()
    
    findings_text = f"""
    KEY FINDINGS
    ════════════════════════════════════════
    
    Total Forbidden Variants: {total_forbidden:,}
    
    PROXIMITY ANALYSIS:
    • {pct_100kb} within 100kb of GWAS loci
    
    EXACT OVERLAPS:
    • GWAS lead SNPs: {gwas_overlap}
    • ClinVar variants: {clinvar_overlap}
    
    INTERPRETATION:
    • Low exact overlap expected (forbidden variants
      don't exist in human populations)
    • Proximity enrichment suggests forbidden sites
      cluster near disease-relevant regulatory regions
    
    ════════════════════════════════════════
    """
    
    ax4.text(0.1, 0.9, findings_text, transform=ax4.transAxes, fontsize=11,
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.3))
    ax4.set_title('D. Summary of Findings')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_path.replace('.png', '.pdf'), bbox_inches='tight')
    plt.close()
    
    logger.info(f"Saved validation summary figure to {output_path}")


def main():
    args = parse_args()
    
    input_dir = Path(args.input)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load data
    proximity_df = None
    overlaps_df = None
    trait_df = None
    
    if (input_dir / 'proximity_enrichment.tsv').exists():
        proximity_df = pd.read_csv(input_dir / 'proximity_enrichment.tsv', sep='\t')
        plot_proximity_enrichment(proximity_df, str(output_dir / 'proximity_enrichment.png'))
    
    if (input_dir / 'exact_overlaps.tsv').exists():
        overlaps_df = pd.read_csv(input_dir / 'exact_overlaps.tsv', sep='\t')
        plot_exact_overlaps(overlaps_df, str(output_dir / 'exact_overlaps.png'))
    
    if (input_dir / 'trait_enrichment.tsv').exists():
        trait_df = pd.read_csv(input_dir / 'trait_enrichment.tsv', sep='\t')
        plot_trait_enrichment(trait_df, str(output_dir / 'trait_enrichment.png'))
    
    # Create summary figure
    if proximity_df is not None and overlaps_df is not None:
        plot_validation_summary(proximity_df, overlaps_df, trait_df,
                               str(output_dir / 'validation_summary.png'))
    
    logger.info("Visualization complete!")


if __name__ == '__main__':
    main()
