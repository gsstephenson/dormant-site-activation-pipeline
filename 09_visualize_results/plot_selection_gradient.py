#!/usr/bin/env python3
"""
Module 09: Selection Gradient Visualization (Improved)

Multi-panel figure showing purifying selection evidence:
- Panel A: Survival rate by Hamming distance (with inset for H=1)
- Panel B: AP1 score distributions with complete statistical comparisons
- Panel C: Fold depletion with log scale
- Panel D: Key findings summary (scannable bullet format)

Improvements over original:
- Panel A: Added inset zoom for H=1 visibility (was invisible at scale)
- Panel B: Fixed y-axis truncation, added all pairwise comparisons
- Panel D: Converted to scannable bullet format with visual hierarchy
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

# Style settings
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10

# Color palette for Hamming distances
COLORS = {
    1: '#d62728',   # Red - forbidden
    2: '#ff7f0e',   # Orange - intermediate
    3: '#2ca02c',   # Green - observed
}


def load_data(results_dir: Path):
    """Load all required data files."""
    # Constraint by Hamming distance
    constraint_path = results_dir / 'purifying_selection' / 'AP1' / 'constraint_by_hamming.tsv'
    constraint_df = pd.read_csv(constraint_path, sep='\t')
    
    # Forbidden variants with AP1 scores
    forbidden_path = results_dir / 'forbidden_variants' / 'AP1' / 'predictions_summary_ap1.tsv'
    forbidden_df = pd.read_csv(forbidden_path, sep='\t')
    forbidden_df['hamming'] = 1  # All forbidden are H=1
    
    # Observed variants (activation landscape)
    landscape_path = results_dir / 'landscape' / 'AP1' / 'AP1_activation_landscape.tsv'
    landscape_df = pd.read_csv(landscape_path, sep='\t')
    
    return constraint_df, forbidden_df, landscape_df


def panel_a_survival_rate(ax, constraint_df):
    """Panel A: Survival rate with inset for H=1 visibility."""
    # Calculate survival rates
    survival_data = []
    for _, row in constraint_df.iterrows():
        h = int(row['hamming'])
        survival_rate = row['observed'] / row['total_possible'] * 100
        survival_data.append({
            'hamming': h,
            'survival_rate': survival_rate,
            'observed': row['observed'],
            'total': row['total_possible']
        })
    
    df = pd.DataFrame(survival_data)
    
    # Main bar plot (H=2 and H=3 visible, H=1 tiny)
    x = [0, 1, 2]
    bars = ax.bar(x, df['survival_rate'], color=[COLORS[1], COLORS[2], COLORS[3]], 
                  edgecolor='black', linewidth=1.5, width=0.6)
    
    ax.set_xticks(x)
    ax.set_xticklabels(['H=1\n(Forbidden)', 'H=2', 'H=3\n(Observed)'])
    ax.set_ylabel('Survival Rate (%)')
    ax.set_title('A. Survival Rate by Hamming Distance', fontweight='bold', loc='left')
    
    # Add value labels on bars
    for i, (bar, row) in enumerate(zip(bars, df.itertuples())):
        height = bar.get_height()
        if i == 0:  # H=1 is too small, add arrow annotation
            ax.annotate(f'{row.survival_rate:.4f}%\n(n={row.observed:,})',
                       xy=(bar.get_x() + bar.get_width()/2, height),
                       xytext=(0.3, 0.15),
                       textcoords='axes fraction',
                       fontsize=9,
                       ha='center',
                       arrowprops=dict(arrowstyle='->', color='#d62728', lw=1.5),
                       bbox=dict(boxstyle='round,pad=0.3', facecolor='#ffcccc', edgecolor='#d62728'))
        else:
            ax.text(bar.get_x() + bar.get_width()/2, height + 0.005,
                   f'{row.survival_rate:.2f}%', ha='center', va='bottom', fontsize=9)
    
    # Add expected reference line (if random)
    # Expected would be ~uniform, but we show actual H=3 as reference
    ax.axhline(y=df[df['hamming']==3]['survival_rate'].values[0], 
               color='gray', linestyle='--', alpha=0.5, label='H=3 baseline')
    
    ax.set_ylim(0, df['survival_rate'].max() * 1.25)
    
    # Add inset axes for H=1 zoom
    inset = ax.inset_axes([0.55, 0.5, 0.4, 0.4])
    inset.bar([0], [df[df['hamming']==1]['survival_rate'].values[0]], 
              color=COLORS[1], edgecolor='black', linewidth=1)
    inset.set_ylabel('Rate (%)', fontsize=8)
    inset.set_xticks([0])
    inset.set_xticklabels(['H=1'], fontsize=8)
    inset.set_title('Zoom: H=1', fontsize=9, fontweight='bold')
    inset.tick_params(axis='both', labelsize=8)
    
    # Connect inset to main plot
    ax.indicate_inset_zoom(inset, edgecolor='gray', alpha=0.5)


def panel_b_ap1_distributions(ax, forbidden_df, landscape_df):
    """Panel B: AP1 score distributions with complete statistical comparisons."""
    # Prepare data
    forbidden_scores = forbidden_df['ap1_raw_max'].dropna().values
    
    # Split landscape by hamming distance
    h2_scores = landscape_df[landscape_df['hamming_distance'] == 2]['ap1_max_raw_score'].dropna().values
    h3_scores = landscape_df[landscape_df['hamming_distance'] == 3]['ap1_max_raw_score'].dropna().values
    
    # Create combined dataframe for violin plot
    plot_data = []
    for score in forbidden_scores[:5000]:  # Sample for visualization
        plot_data.append({'Hamming': 'H=1\n(Forbidden)', 'AP1 Score': score, 'hamming': 1})
    for score in h2_scores[:5000]:
        plot_data.append({'Hamming': 'H=2', 'AP1 Score': score, 'hamming': 2})
    for score in h3_scores[:5000]:
        plot_data.append({'Hamming': 'H=3\n(Observed)', 'AP1 Score': score, 'hamming': 3})
    
    plot_df = pd.DataFrame(plot_data)
    
    # Create violin plot with box plots overlaid
    palette = {f'H=1\n(Forbidden)': COLORS[1], 'H=2': COLORS[2], f'H=3\n(Observed)': COLORS[3]}
    
    parts = ax.violinplot([forbidden_scores[:5000], h2_scores[:5000], h3_scores[:5000]],
                          positions=[0, 1, 2], showmedians=False, showextrema=False)
    
    for i, (pc, h) in enumerate(zip(parts['bodies'], [1, 2, 3])):
        pc.set_facecolor(COLORS[h])
        pc.set_edgecolor('black')
        pc.set_alpha(0.7)
    
    # Add box plots
    bp = ax.boxplot([forbidden_scores[:5000], h2_scores[:5000], h3_scores[:5000]],
                    positions=[0, 1, 2], widths=0.15, patch_artist=True,
                    showfliers=False, medianprops=dict(color='black', linewidth=2))
    
    for patch, h in zip(bp['boxes'], [1, 2, 3]):
        patch.set_facecolor('white')
        patch.set_edgecolor(COLORS[h])
        patch.set_linewidth(1.5)
    
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(['H=1\n(Forbidden)', 'H=2', 'H=3\n(Observed)'])
    ax.set_ylabel('AP1 Raw Score (max)')
    ax.set_title('B. AP1 Score Distribution by Variant Class', fontweight='bold', loc='left')
    
    # Auto-scale y-axis with padding
    all_scores = np.concatenate([forbidden_scores[:5000], h2_scores[:5000], h3_scores[:5000]])
    y_max = np.percentile(all_scores, 99.5) * 1.2
    ax.set_ylim(0, y_max)
    
    # Statistical comparisons - all pairwise
    def add_stat_annotation(ax, x1, x2, y, h, text):
        ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.2, c='black')
        ax.text((x1+x2)/2, y+h, text, ha='center', va='bottom', fontsize=8)
    
    # Perform tests
    stat_12, p_12 = stats.mannwhitneyu(forbidden_scores, h2_scores, alternative='two-sided')
    stat_13, p_13 = stats.mannwhitneyu(forbidden_scores, h3_scores, alternative='two-sided')
    stat_23, p_23 = stats.mannwhitneyu(h2_scores, h3_scores, alternative='two-sided')
    
    y_base = y_max * 0.75
    h_step = y_max * 0.08
    
    # H=1 vs H=2
    p_text = f'p={p_12:.2e}' if p_12 < 0.001 else f'p={p_12:.3f}'
    add_stat_annotation(ax, 0, 1, y_base, h_step, f'ns\n{p_text}' if p_12 > 0.05 else p_text)
    
    # H=2 vs H=3
    y_23 = y_base + h_step * 2
    p_text = f'p={p_23:.2e}' if p_23 < 0.001 else f'p={p_23:.3f}'
    add_stat_annotation(ax, 1, 2, y_23, h_step, f'ns\n{p_text}' if p_23 > 0.05 else p_text)
    
    # H=1 vs H=3
    y_13 = y_base + h_step * 4
    p_text = f'p={p_13:.2e}' if p_13 < 0.001 else f'p={p_13:.3f}'
    add_stat_annotation(ax, 0, 2, y_13, h_step, f'ns\n{p_text}' if p_13 > 0.05 else p_text)
    
    # Add median annotations
    for i, (scores, label) in enumerate([(forbidden_scores, 'H=1'), 
                                          (h2_scores, 'H=2'), 
                                          (h3_scores, 'H=3')]):
        med = np.median(scores)
        ax.text(i, med + y_max*0.02, f'med={med:.0f}', ha='center', fontsize=8, 
                color=COLORS[i+1], fontweight='bold')


def panel_c_fold_depletion(ax, constraint_df):
    """Panel C: Fold depletion bar chart (log scale)."""
    x = [0, 1, 2]
    fold_deps = constraint_df['fold_depletion'].values
    
    bars = ax.bar(x, fold_deps, color=[COLORS[1], COLORS[2], COLORS[3]], 
                  edgecolor='black', linewidth=1.5, width=0.6)
    
    ax.set_yscale('log')
    ax.set_xticks(x)
    ax.set_xticklabels(['H=1', 'H=2', 'H=3'])
    ax.set_ylabel('Fold Depletion (log scale)')
    ax.set_title('C. Selection Intensity', fontweight='bold', loc='left')
    
    # Add value labels
    for bar, row in zip(bars, constraint_df.itertuples()):
        height = bar.get_height()
        p_val = row.binom_p
        p_str = f'p<10⁻³⁰' if p_val < 1e-30 else f'p={p_val:.1e}'
        ax.text(bar.get_x() + bar.get_width()/2, height * 1.2,
               f'{height:.1f}×\n{p_str}', ha='center', va='bottom', fontsize=9)
    
    ax.axhline(y=1, color='gray', linestyle='--', alpha=0.7, label='No selection')
    ax.set_ylim(0.5, 200)
    ax.legend(loc='upper right', fontsize=9)


def panel_d_key_findings(ax):
    """Panel D: Key findings in scannable bullet format."""
    ax.axis('off')
    
    # Create styled text box
    box = FancyBboxPatch((0.02, 0.02), 0.96, 0.96,
                         boxstyle="round,pad=0.02,rounding_size=0.02",
                         facecolor='#f8f9fa', edgecolor='#2c3e50', linewidth=2)
    ax.add_patch(box)
    
    # Title
    ax.text(0.5, 0.92, 'KEY FINDING', fontsize=13, fontweight='bold',
            ha='center', va='top', color='#2c3e50')
    
    ax.text(0.5, 0.82, 'The Selection Paradox', fontsize=11, fontweight='bold',
            ha='center', va='top', color='#34495e', style='italic')
    
    # Bullet points with visual hierarchy (using text symbols for font compatibility)
    bullets = [
        ('▸', 'FORBIDDEN: Only 1 of 31,174 possible H=1 variants observed', '#d62728'),
        ('▸', 'SIMILAR IMPACT: H=1 AP1 scores ≈ H=2/H=3 scores', '#2980b9'),
        ('▸', 'TARGET: Selection purges site COMPLETION,\n      not binding STRENGTH', '#27ae60'),
        ('▸', 'INTENSITY: 76.5× depletion (p < 10⁻³⁰)', '#8e44ad'),
    ]
    
    y_pos = 0.68
    for icon, text, color in bullets:
        ax.text(0.08, y_pos, icon, fontsize=14, ha='left', va='top')
        ax.text(0.16, y_pos, text, fontsize=9, ha='left', va='top', 
                color=color, fontweight='bold',
                transform=ax.transAxes)
        y_pos -= 0.16
    
    # Bottom insight
    ax.text(0.5, 0.08, 
            'Interpretation: Human genomes are actively maintained\n'
            'one mutation away from functional AP1 sites.',
            fontsize=9, ha='center', va='bottom', style='italic', color='#7f8c8d')


def main():
    parser = argparse.ArgumentParser(description='Generate selection gradient visualization')
    parser.add_argument('--results-dir', type=str, default='results',
                        help='Path to results directory')
    parser.add_argument('--output-dir', type=str, default='figures/selection_gradient',
                        help='Output directory for figures')
    parser.add_argument('--dpi', type=int, default=300, help='Figure DPI')
    args = parser.parse_args()
    
    results_dir = Path(args.results_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("Loading data...")
    constraint_df, forbidden_df, landscape_df = load_data(results_dir)
    
    print(f"  Constraint data: {len(constraint_df)} Hamming distances")
    print(f"  Forbidden variants: {len(forbidden_df):,}")
    print(f"  Landscape variants: {len(landscape_df):,}")
    
    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    plt.subplots_adjust(hspace=0.3, wspace=0.25)
    
    print("Generating panels...")
    panel_a_survival_rate(axes[0, 0], constraint_df)
    panel_b_ap1_distributions(axes[0, 1], forbidden_df, landscape_df)
    panel_c_fold_depletion(axes[1, 0], constraint_df)
    panel_d_key_findings(axes[1, 1])
    
    # Overall title
    fig.suptitle('Purifying Selection Against AP1 Site Completion', 
                 fontsize=16, fontweight='bold', y=0.98)
    
    # Save
    png_path = output_dir / 'selection_gradient.png'
    pdf_path = output_dir / 'selection_gradient.pdf'
    
    fig.savefig(png_path, dpi=args.dpi, bbox_inches='tight', facecolor='white')
    fig.savefig(pdf_path, bbox_inches='tight', facecolor='white')
    
    print(f"\nSaved: {png_path}")
    print(f"Saved: {pdf_path}")
    
    plt.close(fig)


if __name__ == '__main__':
    main()
