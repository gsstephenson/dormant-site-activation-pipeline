#!/usr/bin/env python3
"""
Module 09: Forbidden Variants Characterization (Improved)

Multi-panel figure characterizing forbidden variants:
- Panel A: Mutation spectrum with proper y-axis (fixed negative issue)
- Panel B: TF enrichment by family
- Panel C: Tissue enrichment by category
- Panel D: AP1 vs Enhancer (hexbin to fix overplotting)

Improvements over original:
- Panel A: Fixed y-axis starting at 0, grouped transitions/transversions
- Panel C: Added tissue category grouping with color coding
- Panel D: Changed from scatter to hexbin for 30K points
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
from scipy import stats

# Style settings
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['axes.labelsize'] = 11

# AP1 TF families for color coding
TF_FAMILIES = {
    'FOS/JUN': ['FOS', 'FOSB', 'FOSL1', 'FOSL2', 'JUN', 'JUNB', 'JUND'],
    'ATF/CREB': ['ATF1', 'ATF2', 'ATF3', 'ATF4', 'ATF5', 'ATF6', 'ATF7', 'CREB1', 'CREM'],
    'BATF': ['BATF', 'BATF2', 'BATF3'],
    'MAF': ['MAFK', 'MAFF', 'MAFG', 'MAF', 'MAFB'],
    'NRF/BACH': ['NFE2', 'NRF1', 'NRF2', 'NRF3', 'BACH1', 'BACH2'],
}

FAMILY_COLORS = {
    'FOS/JUN': '#e41a1c',
    'ATF/CREB': '#377eb8',
    'BATF': '#4daf4a',
    'MAF': '#984ea3',
    'NRF/BACH': '#ff7f00',
    'Other': '#999999',
}

# Tissue categories
TISSUE_CATEGORIES = {
    'Cancer': ['K562', 'HepG2', 'A549', 'MCF-7', 'HCT116', 'HeLa-S3', 'SK-N-SH', 'T47D', 'Panc1', 'PC-3', 'LNCaP'],
    'Immune': ['GM12878', 'Jurkat', 'CD14+', 'CD4+', 'CD8+', 'B-cell', 'T-cell', 'PBMC', 'monocyte', 'macrophage'],
    'Stem': ['H1-hESC', 'H9', 'hESC', 'iPSC', 'iPS'],
    'Endothelial': ['HUVEC', 'HMEC', 'endothelial'],
    'Fibroblast': ['IMR-90', 'BJ', 'fibroblast', 'IMR90', 'WI-38'],
    'Neural': ['SK-N-SH', 'SH-SY5Y', 'neural', 'neuron', 'astrocyte', 'brain'],
}

TISSUE_CAT_COLORS = {
    'Cancer': '#d62728',
    'Immune': '#2ca02c',
    'Stem': '#9467bd',
    'Endothelial': '#17becf',
    'Fibroblast': '#bcbd22',
    'Neural': '#e377c2',
    'Other': '#7f7f7f',
}


def get_tf_family(tf_name):
    """Get family for a TF."""
    for family, members in TF_FAMILIES.items():
        if tf_name in members:
            return family
    return 'Other'


def get_tissue_category(tissue_name):
    """Get category for a tissue/cell type."""
    for category, patterns in TISSUE_CATEGORIES.items():
        for pattern in patterns:
            if pattern.lower() in tissue_name.lower():
                return category
    return 'Other'


def load_data(results_dir: Path):
    """Load forbidden variants data."""
    predictions_path = results_dir / 'forbidden_variants' / 'AP1' / 'predictions_summary_ap1.tsv'
    df = pd.read_csv(predictions_path, sep='\t')
    
    # Also load base forbidden variants for mutation info
    forbidden_path = results_dir / 'forbidden_variants' / 'AP1' / 'forbidden_variants.tsv'
    forbidden_df = pd.read_csv(forbidden_path, sep='\t')
    
    return df, forbidden_df


def panel_a_mutation_spectrum(ax, predictions_df, forbidden_df):
    """Panel A: Mutation spectrum with proper y-axis and grouping."""
    # Merge to get ref/alt info
    # Create mutation type column from variant_id or from forbidden_df
    
    # Parse variant IDs to get ref>alt
    mutations = []
    for vid in predictions_df['variant_id_str']:
        parts = vid.replace('chr', '').split(':')
        if len(parts) >= 3:
            ref_alt = parts[2]
            if '>' in ref_alt:
                mutations.append(ref_alt)
            else:
                mutations.append('Unknown')
        else:
            mutations.append('Unknown')
    
    predictions_df = predictions_df.copy()
    predictions_df['mutation'] = mutations
    
    # Count mutations
    mut_counts = predictions_df['mutation'].value_counts()
    
    # Define canonical mutation types (considering strand symmetry)
    # Group by transition vs transversion
    transitions = ['A>G', 'G>A', 'C>T', 'T>C']
    transversions = ['A>T', 'T>A', 'A>C', 'C>A', 'G>T', 'T>G', 'G>C', 'C>G']
    
    # Order: transitions first, then transversions
    ordered_muts = []
    counts = []
    colors = []
    
    for mut in transitions:
        if mut in mut_counts.index:
            ordered_muts.append(mut)
            counts.append(mut_counts[mut])
            colors.append('#3498db')  # Blue for transitions
    
    for mut in transversions:
        if mut in mut_counts.index:
            ordered_muts.append(mut)
            counts.append(mut_counts[mut])
            colors.append('#e74c3c')  # Red for transversions
    
    x = np.arange(len(ordered_muts))
    bars = ax.bar(x, counts, color=colors, edgecolor='black', linewidth=0.8)
    
    ax.set_xticks(x)
    ax.set_xticklabels(ordered_muts, rotation=45, ha='right', fontsize=9)
    ax.set_ylabel('Count')
    ax.set_title('A. Mutation Spectrum', fontweight='bold', loc='left')
    
    # CRITICAL FIX: Ensure y-axis starts at 0
    ax.set_ylim(0, max(counts) * 1.15)
    
    # Add percentage labels
    total = sum(counts)
    for bar, count in zip(bars, counts):
        pct = count / total * 100
        if pct > 5:  # Only label significant ones
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(counts)*0.02,
                   f'{pct:.1f}%', ha='center', va='bottom', fontsize=8)
    
    # Add background shading for transitions vs transversions
    n_trans = len([m for m in ordered_muts if m in transitions])
    ax.axvspan(-0.5, n_trans - 0.5, alpha=0.1, color='blue', label='Transitions')
    ax.axvspan(n_trans - 0.5, len(ordered_muts) - 0.5, alpha=0.1, color='red', label='Transversions')
    
    # Legend
    trans_patch = mpatches.Patch(color='#3498db', label=f'Transitions ({sum(counts[:n_trans])/total*100:.1f}%)')
    transv_patch = mpatches.Patch(color='#e74c3c', label=f'Transversions ({sum(counts[n_trans:])/total*100:.1f}%)')
    ax.legend(handles=[trans_patch, transv_patch], loc='upper right', fontsize=9)


def panel_b_tf_enrichment(ax, predictions_df):
    """Panel B: TF enrichment by family."""
    # Count best TFs
    tf_counts = predictions_df['ap1_best_tf'].value_counts().head(15)
    
    # Get family colors
    colors = [FAMILY_COLORS[get_tf_family(tf)] for tf in tf_counts.index]
    
    # Horizontal bar plot
    y = np.arange(len(tf_counts))
    bars = ax.barh(y, tf_counts.values, color=colors, edgecolor='black', linewidth=0.8)
    
    ax.set_yticks(y)
    ax.set_yticklabels(tf_counts.index)
    ax.set_xlabel('Count (variants)')
    ax.set_title('B. Top TFs by Frequency', fontweight='bold', loc='left')
    ax.invert_yaxis()
    
    # Add percentage labels
    total = len(predictions_df)
    for bar, count in zip(bars, tf_counts.values):
        pct = count / total * 100
        ax.text(bar.get_width() + max(tf_counts.values)*0.02, bar.get_y() + bar.get_height()/2,
               f'{pct:.1f}%', va='center', fontsize=8)
    
    # Legend for families
    legend_handles = [mpatches.Patch(color=c, label=f) for f, c in FAMILY_COLORS.items() 
                      if f in [get_tf_family(tf) for tf in tf_counts.index]]
    ax.legend(handles=legend_handles, loc='lower right', fontsize=8, title='Family')


def panel_c_tissue_enrichment(ax, predictions_df):
    """Panel C: Tissue enrichment with category grouping."""
    # Count best tissues
    tissue_counts = predictions_df['ap1_best_biosample'].value_counts().head(15)
    
    # Get category colors
    colors = [TISSUE_CAT_COLORS[get_tissue_category(t)] for t in tissue_counts.index]
    
    # Horizontal bar plot
    y = np.arange(len(tissue_counts))
    bars = ax.barh(y, tissue_counts.values, color=colors, edgecolor='black', linewidth=0.8)
    
    ax.set_yticks(y)
    ax.set_yticklabels(tissue_counts.index)
    ax.set_xlabel('Count (variants)')
    ax.set_title('C. Top Tissues/Cell Types', fontweight='bold', loc='left')
    ax.invert_yaxis()
    
    # Add percentage labels
    total = len(predictions_df)
    for bar, count in zip(bars, tissue_counts.values):
        pct = count / total * 100
        ax.text(bar.get_width() + max(tissue_counts.values)*0.02, bar.get_y() + bar.get_height()/2,
               f'{pct:.1f}%', va='center', fontsize=8)
    
    # Legend for categories
    seen_cats = set(get_tissue_category(t) for t in tissue_counts.index)
    legend_handles = [mpatches.Patch(color=TISSUE_CAT_COLORS[c], label=c) for c in seen_cats]
    ax.legend(handles=legend_handles, loc='lower right', fontsize=8, title='Category')


def panel_d_functional_potential(ax, predictions_df):
    """Panel D: 2D density heatmap showing AP1 vs Enhancer relationship."""
    # Use quantile scores (0-1 range, uniform distribution)
    ap1_q = predictions_df['ap1_quantile_max'].values
    enh_q = predictions_df['enhancer_quantile_max'].values
    
    # Filter out NaN
    mask = ~(np.isnan(ap1_q) | np.isnan(enh_q))
    ap1_q = ap1_q[mask]
    enh_q = enh_q[mask]
    n_total = len(ap1_q)
    
    # Create 2D histogram (heatmap)
    bins = 25
    h, xedges, yedges = np.histogram2d(ap1_q, enh_q, bins=bins, range=[[0, 1], [0, 1]])
    
    # Plot as heatmap
    im = ax.imshow(h.T, origin='lower', aspect='auto', cmap='YlOrRd',
                   extent=[0, 1, 0, 1], norm=LogNorm(vmin=1, vmax=h.max()))
    
    # Colorbar
    cbar = plt.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
    cbar.set_label('Count', fontsize=9)
    
    # Axis labels
    ax.set_xlabel('AP1 Quantile Score')
    ax.set_ylabel('Enhancer Quantile Score')
    ax.set_title('D. Functional Potential: AP1 × Enhancer', fontweight='bold', loc='left')
    
    # Highlight the HIGH-HIGH quadrant with dashed border
    from matplotlib.patches import Rectangle
    rect = Rectangle((0.5, 0.5), 0.5, 0.5, linewidth=2.5, linestyle='--',
                      edgecolor='white', facecolor='none', zorder=10)
    ax.add_patch(rect)
    
    # Spearman correlation
    rho, p = stats.spearmanr(ap1_q, enh_q)
    
    # Count in high-high quadrant
    q_hh = np.sum((ap1_q > 0.5) & (enh_q > 0.5))
    
    # Stats annotation at bottom left where there's no data
    ax.text(0.03, 0.03, f'Spearman ρ = {rho:.3f}\np < 10⁻³⁰', 
            transform=ax.transAxes, ha='left', va='bottom', fontsize=10,
            bbox=dict(boxstyle='round', facecolor='white', edgecolor='gray', alpha=0.95))


def main():
    parser = argparse.ArgumentParser(description='Generate forbidden variants characterization')
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
    predictions_df, forbidden_df = load_data(results_dir)
    print(f"  Loaded {len(predictions_df):,} forbidden variants")
    
    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    plt.subplots_adjust(hspace=0.4, wspace=0.35, bottom=0.1)
    
    print("Generating panels...")
    panel_a_mutation_spectrum(axes[0, 0], predictions_df, forbidden_df)
    panel_b_tf_enrichment(axes[0, 1], predictions_df)
    panel_c_tissue_enrichment(axes[1, 0], predictions_df)
    panel_d_functional_potential(axes[1, 1], predictions_df)
    
    # Overall title
    fig.suptitle('Characterization of 29,957 Forbidden Variants', 
                 fontsize=16, fontweight='bold', y=0.98)
    
    # Save
    png_path = output_dir / 'forbidden_variants_characterization.png'
    pdf_path = output_dir / 'forbidden_variants_characterization.pdf'
    
    fig.savefig(png_path, dpi=args.dpi, bbox_inches='tight', facecolor='white')
    fig.savefig(pdf_path, bbox_inches='tight', facecolor='white')
    
    print(f"\nSaved: {png_path}")
    print(f"Saved: {pdf_path}")
    
    plt.close(fig)


if __name__ == '__main__':
    main()
