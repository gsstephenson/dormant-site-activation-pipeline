#!/usr/bin/env python3
"""
Module 05: Visualization - AP1 Activation Landscape Plot

Creates publication-quality 2D scatter plots showing:
- X-axis: Population accessibility (-log10(AF) × Hamming distance)
- Y-axis: AP1-specific functional impact (max quantile score)

Multiple panel options:
1. Main landscape with density contours
2. Quadrant-annotated version
3. Comparison of AP1-specific vs global max
4. Multi-panel with enhancer marks

Author: George Stephenson
Date: November 2025
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from pathlib import Path
import argparse
import logging
from scipy import stats

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Style settings for publication
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 10,
    'axes.labelsize': 12,
    'axes.titlesize': 14,
    'legend.fontsize': 9,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})


def load_landscape_data(filepath: Path) -> pd.DataFrame:
    """Load activation landscape data."""
    logger.info(f"Loading landscape data from {filepath}")
    df = pd.read_csv(filepath, sep='\t')
    logger.info(f"Loaded {len(df):,} variants")
    return df


def plot_main_landscape(
    data: pd.DataFrame,
    output_path: Path,
    title: str = "AP1 Dormant Site Activation Landscape"
):
    """
    Create the main activation landscape plot.
    
    This is the key figure showing which dormant sites are both:
    - Accessible through human variation (low X)
    - Functionally impactful when activated (high Y)
    
    Y-axis now uses LOG-TRANSFORMED RAW SCORES (not quantile scores) to show
    the true effect magnitude and preserve the selection signal.
    """
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Filter to variants with valid data
    plot_data = data.dropna(subset=['x_accessibility', 'y_ap1_impact'])
    logger.info(f"Plotting {len(plot_data):,} variants with complete data")
    
    if len(plot_data) == 0:
        logger.warning("No data to plot!")
        return
    
    # Color by allele frequency (common = blue, rare = red)
    # Clip to -6 to 0 range for better color visibility (data mostly in -2 to -6 range)
    af_colors = np.log10(plot_data['AF_final'].clip(lower=1e-6))
    
    # Main scatter plot
    # Y-axis is now log10(raw_score), typically ranging from ~0.5 to ~4.5
    scatter = ax.scatter(
        plot_data['x_accessibility'],
        plot_data['y_ap1_impact'],
        c=af_colors,
        cmap='RdYlBu',  # Red (rare) to Blue (common)
        alpha=0.6,
        s=15,
        edgecolors='none',
        vmin=-6,
        vmax=0
    )
    
    # Y-axis: log10(raw_score) ranges from ~0.5 to ~4.5
    # Set axis to show full range with some padding
    y_min_data = plot_data['y_ap1_impact'].min()
    y_max_data = plot_data['y_ap1_impact'].max()
    ax.set_ylim(0, y_max_data * 1.1)
    
    # Add reference lines at key thresholds (75th percentile in raw score space)
    y_75th = plot_data['y_ap1_impact'].quantile(0.75)
    ax.axhline(y=y_75th, color='green', linestyle=':', linewidth=1, alpha=0.7)
    ax.text(ax.get_xlim()[0] + 0.5, y_75th + 0.05, '75th %ile (strong effect)', fontsize=8, color='green')
    
    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax, label='log₁₀(Allele Frequency)', shrink=0.8)
    
    # Add quadrant lines at medians
    x_med = plot_data['x_accessibility'].median()
    y_med = plot_data['y_ap1_impact'].median()
    
    ax.axhline(y=y_med, color='gray', linestyle='--', linewidth=1, alpha=0.7)
    ax.axvline(x=x_med, color='gray', linestyle='--', linewidth=1, alpha=0.7)
    
    # Highlight high-priority quadrant
    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    
    # High priority = low X (accessible) + high Y (impactful)
    rect = Rectangle(
        (x_min, y_med), 
        x_med - x_min, 
        y_max - y_med,
        linewidth=2, 
        edgecolor='green', 
        facecolor='green',
        alpha=0.1,
        label='High Priority Zone'
    )
    ax.add_patch(rect)
    
    # Annotate quadrants
    ax.text(x_min + 0.05*(x_max-x_min), y_max - 0.05*(y_max-y_min), 
            'HIGH PRIORITY\n(Accessible + High Impact)', 
            fontsize=9, fontweight='bold', color='darkgreen',
            verticalalignment='top')
    
    ax.text(x_max - 0.05*(x_max-x_min), y_max - 0.05*(y_max-y_min),
            'High Impact\n(Hard to Access)',
            fontsize=9, color='gray', ha='right', va='top')
    
    ax.text(x_min + 0.05*(x_max-x_min), y_min + 0.05*(y_max-y_min),
            'Low Priority\n(Low Impact)',
            fontsize=9, color='gray', va='bottom')
    
    # Labels and title
    ax.set_xlabel('Population Accessibility Score\n[-log₁₀(AF) × Hamming Distance]\n← More Accessible | Harder to Access →', 
                  fontsize=11)
    ax.set_ylabel('AP1-Family TF Binding Impact\n[log₁₀(Raw Effect Score)]',
                  fontsize=11)
    ax.set_title(title, fontsize=14, fontweight='bold')
    
    # Add summary stats
    n_high_priority = len(plot_data[
        (plot_data['x_accessibility'] < x_med) & 
        (plot_data['y_ap1_impact'] > y_med)
    ])
    
    stats_text = (
        f"n = {len(plot_data):,} variants\n"
        f"High Priority: {n_high_priority:,} ({100*n_high_priority/len(plot_data):.1f}%)\n"
        f"Median X: {x_med:.1f}\n"
        f"Median Y: {y_med:.3f}"
    )
    ax.text(0.98, 0.02, stats_text, transform=ax.transAxes, fontsize=9,
            verticalalignment='bottom', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    logger.info(f"Saved main landscape to {output_path}")


def plot_comparison_panels(
    data: pd.DataFrame,
    output_path: Path
):
    """
    Create multi-panel comparison figure showing:
    1. AP1-specific Y-axis using RAW SCORES (our approach - shows selection signal)
    2. Global max Y-axis (naive approach - quantile scores)
    3. Enhancer marks validation (quantile scores)
    4. Correlation: raw AP1 scores vs enhancer marks
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    plot_data = data.dropna(subset=['x_accessibility', 'y_ap1_impact'])
    
    # Common settings
    x = plot_data['x_accessibility']
    
    # Panel A: AP1-specific using RAW SCORES (PRIMARY - our approach)
    ax = axes[0, 0]
    y = plot_data['y_ap1_impact']  # This is now log10(raw_score)
    y_min, y_max = y.min(), y.max()
    scatter = ax.scatter(x, y, c=y, cmap='Reds', alpha=0.5, s=10, vmin=y_min, vmax=y_max)
    ax.set_xlabel('Accessibility Score')
    ax.set_ylabel('AP1-Family TF Impact [log₁₀(Raw Score)]')
    ax.set_title('A) PRIMARY: AP1-Family TF Binding (Raw Scores)\n(JUND, JUN, FOS, FOSL1/2, ATF3, BATF)', fontweight='bold')
    plt.colorbar(scatter, ax=ax, label='log₁₀(Raw Score)')
    
    # Add 75th percentile reference line and quadrant lines
    y_75 = y.quantile(0.75)
    ax.axhline(y=y_75, color='green', linestyle=':', linewidth=1, alpha=0.7)
    ax.axhline(y=y.median(), color='gray', linestyle='--', alpha=0.5)
    ax.axvline(x=x.median(), color='gray', linestyle='--', alpha=0.5)
    
    # Panel B: Global max (for comparison - still quantile scores)
    ax = axes[0, 1]
    y_global = plot_data['global_max_score']
    scatter = ax.scatter(x, y_global, c=y_global, cmap='Purples', alpha=0.5, s=10, vmin=0, vmax=1)
    ax.set_xlabel('Accessibility Score')
    ax.set_ylabel('Global Max Impact [Quantile]')
    ax.set_title('B) COMPARISON: Global Max (All Tracks)\n(Less biologically specific)', fontweight='bold', color='gray')
    ax.set_ylim(0, 1.05)  # Full 0-1 scale for fair comparison
    plt.colorbar(scatter, ax=ax, label='Quantile Score')
    ax.axhline(y=0.9, color='green', linestyle=':', linewidth=1, alpha=0.7)
    ax.axhline(y=y_global.median(), color='gray', linestyle='--', alpha=0.5)
    ax.axvline(x=x.median(), color='gray', linestyle='--', alpha=0.5)
    
    # Panel C: Enhancer marks validation (quantile scores - for comparison)
    ax = axes[1, 0]
    y_enhancer = plot_data['enhancer_max_score'].fillna(0)
    scatter = ax.scatter(x, y_enhancer, c=y_enhancer, cmap='Greens', alpha=0.5, s=10, vmin=0, vmax=1)
    ax.set_xlabel('Accessibility Score')
    ax.set_ylabel('Enhancer Mark Impact [Quantile]')
    ax.set_title('C) VALIDATION: Enhancer Marks\n(H3K27ac, H3K4me1)', fontweight='bold')
    ax.set_ylim(0, 1.05)  # Full 0-1 scale for consistency
    plt.colorbar(scatter, ax=ax, label='Quantile Score')
    ax.axhline(y=0.9, color='green', linestyle=':', linewidth=1, alpha=0.7)
    ax.axhline(y=y_enhancer.median(), color='gray', linestyle='--', alpha=0.5)
    ax.axvline(x=x.median(), color='gray', linestyle='--', alpha=0.5)
    
    # Panel D: AP1 Raw Score vs Enhancer Quantile correlation
    ax = axes[1, 1]
    valid_both = plot_data.dropna(subset=['y_ap1_impact', 'enhancer_max_score'])
    if len(valid_both) > 10:
        scatter = ax.scatter(
            valid_both['y_ap1_impact'],  # log10(raw_score)
            valid_both['enhancer_max_score'],  # quantile
            c=valid_both['x_accessibility'],
            cmap='viridis_r',
            alpha=0.5, 
            s=10
        )
        plt.colorbar(scatter, ax=ax, label='Accessibility Score')
        
        # Add correlation
        r, p = stats.pearsonr(valid_both['y_ap1_impact'], valid_both['enhancer_max_score'])
        ax.text(0.05, 0.95, f'r = {r:.3f}\np = {p:.2e}', 
                transform=ax.transAxes, fontsize=10,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    ax.set_xlabel('AP1-Family TF Impact [log₁₀(Raw Score)]')
    ax.set_ylabel('Enhancer Mark Impact [Quantile]')
    ax.set_ylim(0, 1.05)
    ax.set_title('D) CONCORDANCE: AP1 vs Enhancer\n(Validates functional activation)', fontweight='bold')
    
    plt.suptitle('AP1 Dormant Site Activation Landscape - Multi-Modal Analysis', 
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    logger.info(f"Saved comparison panels to {output_path}")


def plot_tf_breakdown(data: pd.DataFrame, output_path: Path):
    """
    Show which AP1-family TFs drive the signal.
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Panel A: Distribution of best TF
    ax = axes[0]
    tf_counts = data['ap1_best_tf'].value_counts()
    colors = plt.cm.Set3(np.linspace(0, 1, len(tf_counts)))
    bars = ax.barh(tf_counts.index, tf_counts.values, color=colors)
    ax.set_xlabel('Number of Variants')
    ax.set_ylabel('Best AP1-Family TF')
    ax.set_title('A) Which TF Shows Strongest Signal?', fontweight='bold')
    
    # Add counts on bars
    for bar, count in zip(bars, tf_counts.values):
        ax.text(bar.get_width() + 50, bar.get_y() + bar.get_height()/2,
                f'{count:,}', va='center', fontsize=9)
    
    # Panel B: Score distribution by TF - USE RAW SCORES (log-transformed)
    ax = axes[1]
    tfs_to_plot = tf_counts.head(8).index.tolist()
    tf_data = data[data['ap1_best_tf'].isin(tfs_to_plot)]
    
    positions = range(len(tfs_to_plot))
    # Use y_ap1_impact which is now log10(raw_score)
    bp_data = [tf_data[tf_data['ap1_best_tf'] == tf]['y_ap1_impact'].dropna() 
               for tf in tfs_to_plot]
    
    bp = ax.boxplot(bp_data, positions=positions, patch_artist=True)
    for patch, color in zip(bp['boxes'], colors[:len(tfs_to_plot)]):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    ax.set_xticks(positions)
    ax.set_xticklabels(tfs_to_plot, rotation=45, ha='right')
    ax.set_ylabel('AP1 Impact [log₁₀(Raw Score)]')
    # Use data-driven y-limits for log-transformed raw scores
    y_75 = data['y_ap1_impact'].quantile(0.75)
    ax.axhline(y=y_75, color='green', linestyle=':', linewidth=1, alpha=0.7)
    ax.set_title('B) Effect Size Distribution by TF', fontweight='bold')
    
    plt.suptitle('AP1-Family Transcription Factor Analysis', fontsize=12, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    logger.info(f"Saved TF breakdown to {output_path}")


def plot_high_priority_detail(data: pd.DataFrame, output_path: Path):
    """
    Detailed view of high-priority candidates.
    """
    # Get high-priority variants
    x_med = data['x_accessibility'].median()
    y_med = data['y_ap1_impact'].median()
    
    high_priority = data[
        (data['x_accessibility'] < x_med) &
        (data['y_ap1_impact'] > y_med)
    ].copy()
    
    if len(high_priority) == 0:
        logger.warning("No high-priority variants to plot")
        return
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Panel A: Zoom on high-priority quadrant - use actual data range for log(raw_score)
    ax = axes[0, 0]
    scatter = ax.scatter(
        high_priority['x_accessibility'],
        high_priority['y_ap1_impact'],
        c=high_priority['AF_final'].apply(lambda x: np.log10(max(x, 1e-10))),
        cmap='RdYlBu',
        alpha=0.7,
        s=30
    )
    ax.set_xlabel('Accessibility Score')
    ax.set_ylabel('AP1 Impact [log₁₀(Raw Score)]')
    # Use actual data range for log-transformed raw scores
    y_min = high_priority['y_ap1_impact'].min() * 0.95
    y_max = high_priority['y_ap1_impact'].max() * 1.05
    ax.set_ylim(y_min, y_max)
    ax.set_title(f'A) High-Priority Candidates (n={len(high_priority):,})', fontweight='bold')
    
    # Add annotation showing 75th percentile
    y_75 = high_priority['y_ap1_impact'].quantile(0.75)
    ax.axhline(y=y_75, color='green', linestyle='--', alpha=0.5, linewidth=1)
    ax.text(ax.get_xlim()[1] * 0.95, y_75 + 0.02, '75th %ile', ha='right', fontsize=8, color='green')
    
    # Panel B: AF distribution of high-priority - use -log10(AF) so higher = rarer
    ax = axes[0, 1]
    af_neg_log = -np.log10(high_priority['AF_final'].clip(lower=1e-10))
    ax.hist(af_neg_log, bins=50, color='steelblue', edgecolor='white', alpha=0.7)
    ax.axvline(x=-np.log10(0.01), color='red', linestyle='--', label='1% MAF')
    ax.axvline(x=-np.log10(0.001), color='orange', linestyle='--', label='0.1% MAF')
    ax.set_xlabel('-log₁₀(Allele Frequency)\n← Common | Rare →')
    ax.set_ylabel('Count')
    ax.set_title('B) AF Distribution of High-Priority Variants', fontweight='bold')
    ax.legend()
    
    # Panel C: Hamming distance distribution
    ax = axes[1, 0]
    ax.hist(high_priority['hamming_distance'].dropna(), bins=range(1, 6), 
            color='coral', edgecolor='white', alpha=0.7, align='left')
    ax.set_xlabel('Hamming Distance (Mutations Needed)')
    ax.set_ylabel('Count')
    ax.set_title('C) Mutations Required for Activation', fontweight='bold')
    ax.set_xticks([1, 2, 3, 4])
    
    # Panel D: Top candidates table
    ax = axes[1, 1]
    ax.axis('off')
    
    top10 = high_priority.nlargest(10, 'y_ap1_impact')[
        ['variant_id_str', 'y_ap1_impact', 'ap1_best_tf', 'AF_final', 'hamming_distance']
    ].copy()
    top10['AF_final'] = top10['AF_final'].apply(lambda x: f'{x:.2e}')
    top10['y_ap1_impact'] = top10['y_ap1_impact'].apply(lambda x: f'{x:.3f}')
    top10.columns = ['Variant', 'AP1 Score', 'Best TF', 'AF', 'Steps']
    
    table = ax.table(
        cellText=top10.values,
        colLabels=top10.columns,
        loc='center',
        cellLoc='center'
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1.2, 1.5)
    ax.set_title('D) Top 10 High-Priority Candidates', fontweight='bold', y=0.95)
    
    plt.suptitle('High-Priority Dormant AP1 Site Candidates', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    logger.info(f"Saved high-priority detail to {output_path}")


def main():
    parser = argparse.ArgumentParser(description='Generate AP1 Activation Landscape Plots')
    parser.add_argument(
        '--input',
        type=Path,
        default=Path('results/landscape/AP1/AP1_activation_landscape.tsv'),
        help='Input landscape TSV file'
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=Path('figures/landscape'),
        help='Output directory for figures'
    )
    
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load data
    data = load_landscape_data(args.input)
    
    # Generate all plots
    plot_main_landscape(
        data, 
        args.output_dir / 'AP1_activation_landscape_main.png',
        title='AP1 Dormant Site Activation Landscape'
    )
    
    plot_comparison_panels(
        data,
        args.output_dir / 'AP1_activation_landscape_comparison.png'
    )
    
    plot_tf_breakdown(
        data,
        args.output_dir / 'AP1_tf_breakdown.png'
    )
    
    plot_high_priority_detail(
        data,
        args.output_dir / 'AP1_high_priority_candidates.png'
    )
    
    logger.info("=" * 60)
    logger.info("ALL PLOTS GENERATED SUCCESSFULLY")
    logger.info("=" * 60)


if __name__ == '__main__':
    main()
