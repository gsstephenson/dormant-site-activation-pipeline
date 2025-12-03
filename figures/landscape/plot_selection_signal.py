#!/usr/bin/env python3
"""
Plot purifying selection signal: AF vs AlphaGenome effect size

Shows that variants with higher predicted AP1 binding activation 
are rarer in the population (lower AF) - evidence of purifying selection.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Load data
df = pd.read_csv('results/landscape/AP1/AP1_activation_landscape.tsv', sep='\t')
print(f"Loaded {len(df):,} variants")

# Use raw scores (not log-transformed) for effect size
# And AC (allele count) which is more intuitive than AF
effect = df['ap1_max_raw_score'].values
af = df['AF_final'].values
ac = df['AC_final'].values

# Filter to variants with valid data
mask = (effect > 0) & np.isfinite(effect) & np.isfinite(af)
effect = effect[mask]
af = af[mask]
ac = ac[mask]

print(f"Plotting {len(effect):,} variants with valid data")

# Create figure with multiple panels
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# =============================================================================
# Panel A: Raw scatter - AF vs Effect Size
# =============================================================================
ax = axes[0, 0]

# Use log scales for both axes
log_af = np.log10(np.clip(af, 1e-10, 1))  # Clip to avoid log(0)
log_effect = np.log10(effect)

# Color by Hamming distance
hamming = df['hamming_distance'].values[mask]
colors = {1: '#2ca02c', 2: '#ff7f0e', 3: '#d62728'}
c = [colors.get(int(h), '#d62728') for h in hamming]

scatter = ax.scatter(log_af, log_effect, c=c, alpha=0.4, s=15, edgecolors='none')

# Add regression line
slope, intercept, r, p, se = stats.linregress(log_af, log_effect)
x_line = np.linspace(log_af.min(), log_af.max(), 100)
y_line = slope * x_line + intercept
ax.plot(x_line, y_line, 'k--', linewidth=2, label=f'r = {r:.3f}, p = {p:.2e}')

ax.set_xlabel('Allele Frequency [log₁₀(AF)]', fontsize=12)
ax.set_ylabel('AP1 Binding Effect [log₁₀(raw score)]', fontsize=12)
ax.set_title('A) Purifying Selection: Rarer Variants Have Stronger Effects', fontsize=12, fontweight='bold')
ax.legend(loc='upper left')

# Add interpretation
if slope < 0:
    ax.text(0.95, 0.05, 'NEGATIVE correlation:\nStronger effects → Rarer variants\n= Purifying selection', 
            transform=ax.transAxes, fontsize=9, ha='right', va='bottom',
            bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))

# =============================================================================
# Panel B: Binned means - clearer signal
# =============================================================================
ax = axes[0, 1]

# Bin by AF
af_bins = [-10, -6, -5, -4, -3, -2, -1, 0]
af_labels = ['<10⁻⁶', '10⁻⁶-10⁻⁵', '10⁻⁵-10⁻⁴', '10⁻⁴-10⁻³', '10⁻³-10⁻²', '10⁻²-10⁻¹', '>10⁻¹']
bin_idx = np.digitize(log_af, af_bins)

means = []
stds = []
counts = []
for i in range(1, len(af_bins)):
    mask_bin = bin_idx == i
    if mask_bin.sum() > 0:
        means.append(np.mean(log_effect[mask_bin]))
        stds.append(np.std(log_effect[mask_bin]) / np.sqrt(mask_bin.sum()))  # SEM
        counts.append(mask_bin.sum())
    else:
        means.append(np.nan)
        stds.append(np.nan)
        counts.append(0)

x_pos = range(len(af_labels))
bars = ax.bar(x_pos, means, yerr=stds, capsize=5, color='steelblue', edgecolor='black', alpha=0.7)
ax.set_xticks(x_pos)
ax.set_xticklabels(af_labels, rotation=45, ha='right')
ax.set_xlabel('Allele Frequency Bin', fontsize=12)
ax.set_ylabel('Mean AP1 Effect [log₁₀(raw score)]', fontsize=12)
ax.set_title('B) Mean Effect Size by AF Bin', fontsize=12, fontweight='bold')

# Add counts on bars
for i, (bar, count) in enumerate(zip(bars, counts)):
    if count > 0:
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02, 
                f'n={count}', ha='center', va='bottom', fontsize=8)

# =============================================================================
# Panel C: Effect size distribution by rarity category
# =============================================================================
ax = axes[1, 0]

# Create rarity categories
def categorize_rarity(af_val):
    if af_val == 0:
        return 'Novel (AC=0)'
    elif af_val < 1e-5:
        return 'Ultra-rare\n(AF<10⁻⁵)'
    elif af_val < 1e-3:
        return 'Rare\n(10⁻⁵-10⁻³)'
    else:
        return 'Common\n(AF>10⁻³)'

categories = [categorize_rarity(a) for a in af]
cat_order = ['Novel (AC=0)', 'Ultra-rare\n(AF<10⁻⁵)', 'Rare\n(10⁻⁵-10⁻³)', 'Common\n(AF>10⁻³)']

# Box plot
data_by_cat = [log_effect[np.array(categories) == cat] for cat in cat_order]
bp = ax.boxplot(data_by_cat, labels=cat_order, patch_artist=True)

colors_box = ['purple', 'red', 'orange', 'blue']
for patch, color in zip(bp['boxes'], colors_box):
    patch.set_facecolor(color)
    patch.set_alpha(0.6)

ax.set_xlabel('Variant Rarity Category', fontsize=12)
ax.set_ylabel('AP1 Effect [log₁₀(raw score)]', fontsize=12)
ax.set_title('C) Effect Size Distribution by Rarity', fontsize=12, fontweight='bold')

# Add median values
for i, cat in enumerate(cat_order):
    data = log_effect[np.array(categories) == cat]
    if len(data) > 0:
        median = np.median(data)
        ax.text(i+1, median + 0.05, f'{median:.2f}', ha='center', fontsize=9, fontweight='bold')

# =============================================================================
# Panel D: Hamming distance effect
# =============================================================================
ax = axes[1, 1]

# Compare effect sizes across Hamming distances
hamming_cats = [1, 2, 3]
data_by_hamming = []
for h in hamming_cats:
    mask_h = hamming == h
    if mask_h.sum() > 0:
        data_by_hamming.append(log_effect[mask_h])
    else:
        data_by_hamming.append([])

bp = ax.boxplot(data_by_hamming, labels=['1 mutation\n(n=1)', '2 mutations\n(n=677)', '3 mutations\n(n=6,359)'], 
                patch_artist=True)

colors_h = ['#2ca02c', '#ff7f0e', '#d62728']
for patch, color in zip(bp['boxes'], colors_h):
    patch.set_facecolor(color)
    patch.set_alpha(0.6)

ax.set_xlabel('Hamming Distance (Mutations to Consensus)', fontsize=12)
ax.set_ylabel('AP1 Effect [log₁₀(raw score)]', fontsize=12)
ax.set_title('D) Effect Size by Mutational Distance', fontsize=12, fontweight='bold')

# Add note about selection
ax.text(0.02, 0.98, 
        'Note: Hamming=1 sites are under\nextreme purifying selection\n(only 1 of 31,174 survives in gnomAD)',
        transform=ax.transAxes, fontsize=9, va='top',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

# =============================================================================
# Final layout
# =============================================================================
plt.suptitle('Evidence for Purifying Selection Against AP1 Site Activation\n'
             'Variants with Stronger Predicted Effects are Rarer in the Population', 
             fontsize=14, fontweight='bold', y=1.02)

plt.tight_layout()
plt.savefig('figures/landscape/AP1_selection_signal.png', dpi=300, bbox_inches='tight')
plt.close()

print("\n✓ Saved: figures/landscape/AP1_selection_signal.png")

# Print summary statistics
print("\n=== SELECTION SIGNAL SUMMARY ===")
print(f"Correlation (log AF vs log effect): r = {r:.4f}, p = {p:.2e}")
print(f"Slope: {slope:.4f} (negative = purifying selection)")
print(f"\nInterpretation:")
if slope < 0 and p < 0.05:
    print("  ✓ SIGNIFICANT PURIFYING SELECTION DETECTED")
    print("  ✓ Variants that would create stronger AP1 binding are rarer")
    print("  ✓ This suggests selection acts against dormant site activation")
