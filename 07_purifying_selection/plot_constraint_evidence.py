#!/usr/bin/env python3
"""
Module 07: Plot Constraint Evidence for Dormant AP1 Sites

Generates publication-quality figures showing:
1. Coverage distribution across dormant site positions
2. Observed vs Expected by Hamming distance
3. Fold depletion (selection intensity)
4. Statistical summary panel

Author: George Stephenson
Date: December 2025
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

# Style settings
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'figure.dpi': 150,
    'savefig.dpi': 300,
})


def load_constraint_data(input_dir: Path, motif_family: str = "AP1") -> pd.DataFrame:
    """Load constraint statistics from TSV file."""
    # Try with motif family prefix first, then without
    constraint_file = input_dir / f"{motif_family}_constraint_by_hamming.tsv"
    if not constraint_file.exists():
        constraint_file = input_dir / "constraint_by_hamming.tsv"
    return pd.read_csv(constraint_file, sep='\t')


def load_coverage_data(input_dir: Path, motif_family: str = "AP1", hamming: int = 1) -> pd.DataFrame:
    """Load per-position coverage data if available."""
    coverage_file = input_dir / f"{motif_family}_h{hamming}_coverage.tsv"
    if coverage_file.exists():
        return pd.read_csv(coverage_file, sep='\t')
    return None


def plot_constraint_evidence(
    constraint_df: pd.DataFrame,
    output_dir: Path,
    motif_family: str = "AP1"
) -> None:
    """Generate main constraint evidence figure."""
    
    fig = plt.figure(figsize=(14, 10))
    
    hamming_vals = constraint_df['hamming'].tolist()
    colors = ['#d62728', '#ff7f0e', '#2ca02c'][:len(hamming_vals)]  # Red, orange, green
    
    # ========================================================================
    # Panel A: Survival Rate by Hamming Distance
    # ========================================================================
    ax1 = fig.add_subplot(2, 2, 1)
    
    survival_rates = [
        row['observed'] / row['total_possible'] * 100 
        for _, row in constraint_df.iterrows()
    ]
    
    bars = ax1.bar(hamming_vals, survival_rates, color=colors, edgecolor='black', linewidth=1.5)
    
    ax1.set_xlabel('Hamming Distance\n(Mutations from AP1 Consensus)', fontsize=11)
    ax1.set_ylabel('Variant Survival Rate (%)', fontsize=11)
    ax1.set_title('A) Variants Closer to AP1 Have Lower Survival Rates\n(Coverage-Validated)', 
                  fontsize=12, fontweight='bold')
    ax1.set_xticks(hamming_vals)
    ax1.set_xticklabels([f'H={h}' for h in hamming_vals])
    
    # Add annotations
    for i, (h, rate) in enumerate(zip(hamming_vals, survival_rates)):
        row = constraint_df[constraint_df['hamming'] == h].iloc[0]
        ax1.text(h, rate * 1.1, f'{rate:.4f}%', ha='center', va='bottom', 
                fontsize=10, fontweight='bold')
        ax1.text(h, rate * 0.5, f'{row["observed"]:,} / {row["total_possible"]:,}', 
                ha='center', va='center', fontsize=8, color='white')
    
    ax1.set_ylim(0, max(survival_rates) * 1.4)
    
    # ========================================================================
    # Panel B: Fold Depletion (Selection Intensity)
    # ========================================================================
    ax2 = fig.add_subplot(2, 2, 2)
    
    fold_depletions = constraint_df['fold_depletion'].tolist()
    
    bars = ax2.bar(hamming_vals, fold_depletions, color=colors, edgecolor='black', linewidth=1.5)
    ax2.axhline(y=1, color='gray', linestyle='--', linewidth=1.5, label='No selection')
    
    ax2.set_xlabel('Hamming Distance', fontsize=11)
    ax2.set_ylabel('Fold Depletion\n(Expected / Observed)', fontsize=11)
    ax2.set_title('B) Selection Intensity by Distance from Consensus', 
                  fontsize=12, fontweight='bold')
    ax2.set_xticks(hamming_vals)
    ax2.set_xticklabels([f'H={h}' for h in hamming_vals])
    
    for h, fd in zip(hamming_vals, fold_depletions):
        if not np.isinf(fd):
            ax2.text(h, fd * 1.05, f'{fd:.1f}×', ha='center', va='bottom', 
                    fontsize=14, fontweight='bold')
    
    ax2.set_ylim(0, max([f for f in fold_depletions if not np.isinf(f)]) * 1.3)
    ax2.legend(loc='upper right')
    
    # ========================================================================
    # Panel C: Coverage Quality
    # ========================================================================
    ax3 = fig.add_subplot(2, 2, 3)
    
    categories = ['High\n(AN≥100K)', 'Medium\n(50K-100K)', 'Low\n(<50K)', 'Missing']
    
    x = np.arange(len(hamming_vals))
    width = 0.2
    
    for i, (_, row) in enumerate(constraint_df.iterrows()):
        total = row['total_possible']
        values = [
            row['high_conf'] / total * 100,
            row['med_conf'] / total * 100,
            row['low_conf'] / total * 100,
            row['missing'] / total * 100
        ]
        
        bottom = 0
        for j, (val, cat) in enumerate(zip(values, categories)):
            color = ['#2ca02c', '#ffbb78', '#ff9896', '#d62728'][j]
            ax3.bar(i, val, bottom=bottom, color=color, 
                   label=cat if i == 0 else '', edgecolor='white', linewidth=0.5)
            bottom += val
    
    ax3.set_xlabel('Hamming Distance', fontsize=11)
    ax3.set_ylabel('Percentage of Positions', fontsize=11)
    ax3.set_title('C) Coverage Quality of Dormant Site Positions', 
                  fontsize=12, fontweight='bold')
    ax3.set_xticks(range(len(hamming_vals)))
    ax3.set_xticklabels([f'H={h}' for h in hamming_vals])
    ax3.legend(loc='upper right', fontsize=8, title='Coverage')
    ax3.set_ylim(0, 100)
    
    # ========================================================================
    # Panel D: Statistical Summary
    # ========================================================================
    ax4 = fig.add_subplot(2, 2, 4)
    ax4.axis('off')
    
    h1_row = constraint_df[constraint_df['hamming'] == 1]
    if len(h1_row) > 0:
        h1 = h1_row.iloc[0]
    else:
        h1 = constraint_df.iloc[0]
    
    summary_text = f"""
STATISTICAL SUMMARY: COVERAGE-VALIDATED CONSTRAINT
═══════════════════════════════════════════════════════════════

DATASET:
  gnomAD v4.1 genomes: 807,003 individuals
  Maximum AN: 1,614,006 (all samples called)

HAMMING DISTANCE = 1 (Closest to AP1 Consensus):
  Total possible mutations:     {h1['total_possible']:>12,}
  High-confidence (AN≥100K):    {h1['high_conf']:>12,} ({h1['high_conf']/h1['total_possible']*100:.1f}%)
  
  Observed in gnomAD:           {h1['observed']:>12,}
  Expected (uniform rate):      {h1['expected']:>12,.1f}
  
  Fold depletion:               {h1['fold_depletion']:>12.1f}×
  Binomial p-value:             {h1['binom_p']:>12.2e}

INTERPRETATION:
  Of the ~{h1['total_possible']:,} possible single-nucleotide mutations
  that would create functional AP1 binding sites:
  
  • {h1['high_conf']:,} are at well-covered positions
  • Only {h1['observed']} variant(s) observed in 807K individuals
  • This is {h1['fold_depletion']:.0f}× fewer than expected under neutrality
  
  → STRONG EVIDENCE of purifying selection against
    de novo AP1 site creation
"""
    
    ax4.text(0.02, 0.98, summary_text, transform=ax4.transAxes,
             fontsize=9, verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='lightyellow', edgecolor='gray', alpha=0.9))
    
    ax4.set_title('D) Statistical Summary', fontsize=12, fontweight='bold', y=1.0)
    
    # ========================================================================
    # Main title
    # ========================================================================
    fig.suptitle(
        f'Purifying Selection Against {motif_family} Site Activation\n'
        'Coverage-Validated Analysis Using gnomAD v4.1 Allele Numbers',
        fontsize=14, fontweight='bold', y=1.02
    )
    
    plt.tight_layout()
    
    output_file = output_dir / f"{motif_family}_constraint_evidence.png"
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    plt.close()
    
    print(f"✓ Saved: {output_file}")


def plot_coverage_distribution(
    coverage_df: pd.DataFrame,
    output_dir: Path,
    motif_family: str = "AP1",
    hamming: int = 1
) -> None:
    """Plot distribution of AN values for a specific Hamming distance."""
    
    if coverage_df is None:
        print(f"No coverage data available for H={hamming}")
        return
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Panel A: Histogram
    ax1 = axes[0]
    an_values = coverage_df[coverage_df['AN'] > 0]['AN']
    
    ax1.hist(an_values, bins=50, color='steelblue', edgecolor='black', alpha=0.7)
    ax1.axvline(x=100000, color='red', linestyle='--', linewidth=2, label='High-conf threshold')
    ax1.axvline(x=50000, color='orange', linestyle='--', linewidth=2, label='Medium-conf threshold')
    
    ax1.set_xlabel('Allele Number (AN)', fontsize=11)
    ax1.set_ylabel('Number of Positions', fontsize=11)
    ax1.set_title(f'A) Distribution of Coverage at H={hamming} Positions', 
                  fontsize=12, fontweight='bold')
    ax1.legend()
    
    # Panel B: Pie chart of confidence levels
    ax2 = axes[1]
    
    conf_counts = coverage_df['confidence'].value_counts()
    labels = ['High (AN≥100K)', 'Medium (50K-100K)', 'Low (<50K)', 'Missing']
    sizes = [
        conf_counts.get('high', 0),
        conf_counts.get('medium', 0),
        conf_counts.get('low', 0),
        conf_counts.get('missing', 0)
    ]
    colors_pie = ['#2ca02c', '#ffbb78', '#ff9896', '#d62728']
    explode = (0.05, 0, 0, 0.1)
    
    ax2.pie(sizes, explode=explode, labels=labels, colors=colors_pie,
            autopct='%1.1f%%', startangle=90)
    ax2.set_title(f'B) Coverage Quality for H={hamming} Positions',
                  fontsize=12, fontweight='bold')
    
    fig.suptitle(f'{motif_family} Dormant Sites: Coverage Analysis',
                 fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    
    output_file = output_dir / f"{motif_family}_h{hamming}_coverage_distribution.png"
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    plt.close()
    
    print(f"✓ Saved: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Plot constraint evidence for dormant AP1 sites"
    )
    
    parser.add_argument(
        '--input', '-i',
        type=Path,
        required=True,
        help='Input directory containing constraint analysis results'
    )
    
    parser.add_argument(
        '--output', '-o',
        type=Path,
        required=True,
        help='Output directory for figures'
    )
    
    parser.add_argument(
        '--motif-family',
        type=str,
        default='AP1',
        help='Motif family name (default: AP1)'
    )
    
    args = parser.parse_args()
    
    # Create output directory
    args.output.mkdir(parents=True, exist_ok=True)
    
    # Load constraint data
    constraint_df = load_constraint_data(args.input, args.motif_family)
    print(f"Loaded constraint data for {len(constraint_df)} Hamming distances")
    
    # Generate main figure
    plot_constraint_evidence(constraint_df, args.output, args.motif_family)
    
    # Generate coverage distribution plots if data available
    for hamming in constraint_df['hamming'].unique():
        coverage_df = load_coverage_data(args.input, args.motif_family, int(hamming))
        if coverage_df is not None:
            plot_coverage_distribution(coverage_df, args.output, args.motif_family, int(hamming))
    
    print("\nAll figures generated successfully!")
    return 0


if __name__ == '__main__':
    sys.exit(main())
