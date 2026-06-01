"""
Create Visualizations for E:I Ratio vs rCMRGlc Analysis
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from scipy import stats

print('='*80)
print('CREATING E:I vs rCMRGlc VISUALIZATIONS')
print('='*80)

# Paths
PROJECT_DIR = Path(r"C:\Users\michaelproulx\MetabolicBrain")
DATA_DIR = PROJECT_DIR / "cellxgene_integration" / "ei_analysis" / "data" / "cortical"
OUTPUT_DIR = DATA_DIR / "visualizations"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Load data
print('\n[1/5] Loading integrated data...')
integrated = pd.read_csv(DATA_DIR / 'ei_metabolism_integrated.csv')
fine_grained = pd.read_csv(DATA_DIR / 'ei_ratios_whb_cortical_fine_grained.csv')
print(f'  Lobes: {len(integrated)}')
print(f'  Fine-grained regions: {len(fine_grained)}')

# Filter valid data for correlation
valid = integrated[integrated['EI_Ratio_aggregated'].notna() & integrated['rCMRGlc_Mean'].notna()]

# Set plot style
plt.style.use('seaborn-v0_8-darkgrid')
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 11

# 1. Scatter plot: E:I vs rCMRGlc with regression
print('\n[2/5] Creating scatter plot with regression...')
fig, ax = plt.subplots(figsize=(10, 8))

ei = valid['EI_Ratio_aggregated'].values
met = valid['rCMRGlc_Mean'].values
lobes = valid['Lobe'].values

# Calculate regression
slope, intercept, r_value, p_value, std_err = stats.linregress(ei, met)
line_x = np.array([ei.min() - 0.2, ei.max() + 0.2])
line_y = slope * line_x + intercept

# Plot points with labels
colors = plt.cm.Set2(range(len(ei)))
for i, (x, y, lobe, color) in enumerate(zip(ei, met, lobes, colors)):
    ax.scatter(x, y, s=200, alpha=0.7, color=color, edgecolors='black', linewidth=1.5)
    # Clean lobe name for label
    label = lobe.replace('_', ' ').replace(' lobe', '').replace(' cortex', '')
    ax.annotate(label, (x, y), xytext=(5, 5), textcoords='offset points',
                fontsize=10, fontweight='bold')

# Regression line
ax.plot(line_x, line_y, 'r--', linewidth=2, alpha=0.8,
        label=f'y = {slope:.2f}x + {intercept:.1f}\nR² = {r_value**2:.3f}, p = {p_value:.3f}')

ax.set_xlabel('E:I Ratio (Excitatory/Inhibitory)', fontsize=13, fontweight='bold')
ax.set_ylabel('rCMRGlc (μmol/100g/min)', fontsize=13, fontweight='bold')
ax.set_title('E:I Ratio vs Regional Cerebral Glucose Metabolism\nWhole Human Brain Atlas (n=6 lobes)',
             fontsize=14, fontweight='bold', pad=20)
ax.legend(loc='upper left', fontsize=11, framealpha=0.9)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'ei_vs_rcmrglc_scatter.png', dpi=300, bbox_inches='tight')
print(f'  Saved: ei_vs_rcmrglc_scatter.png')
plt.close()

# 2. Bar chart: E:I ratios by lobe
print('\n[3/5] Creating E:I ratios bar chart...')
fig, ax = plt.subplots(figsize=(12, 7))

lobes_sorted = integrated.sort_values('EI_Ratio_aggregated', ascending=False)
x_pos = np.arange(len(lobes_sorted))
bars = ax.bar(x_pos, lobes_sorted['EI_Ratio_aggregated'],
              color=plt.cm.viridis(np.linspace(0.2, 0.9, len(lobes_sorted))),
              edgecolor='black', linewidth=1.5, alpha=0.8)

# Add value labels on bars
for i, (idx, row) in enumerate(lobes_sorted.iterrows()):
    height = row['EI_Ratio_aggregated']
    ax.text(i, height + 0.1, f"{height:.2f}", ha='center', va='bottom',
            fontsize=11, fontweight='bold')

# Clean lobe names for x-axis
labels = [lobe.replace('_', ' ') for lobe in lobes_sorted['Lobe']]
ax.set_xticks(x_pos)
ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=11)
ax.set_ylabel('E:I Ratio', fontsize=13, fontweight='bold')
ax.set_title('Excitatory:Inhibitory Neuron Ratios Across Cortical Lobes\nWhole Human Brain Atlas',
             fontsize=14, fontweight='bold', pad=20)
ax.axhline(y=3.0, color='red', linestyle='--', linewidth=2, alpha=0.5, label='Typical cortical E:I (3:1)')
ax.legend(fontsize=11)
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'ei_ratios_by_lobe.png', dpi=300, bbox_inches='tight')
print(f'  Saved: ei_ratios_by_lobe.png')
plt.close()

# 3. Dual-axis plot: E:I and rCMRGlc comparison
print('\n[4/5] Creating dual-axis comparison...')
fig, ax1 = plt.subplots(figsize=(14, 8))

lobes_sorted = valid.sort_values('EI_Ratio_aggregated')
x_pos = np.arange(len(lobes_sorted))
labels = [lobe.replace('_', ' ') for lobe in lobes_sorted['Lobe']]

# E:I ratio bars
color1 = '#2E86AB'
bars1 = ax1.bar(x_pos - 0.2, lobes_sorted['EI_Ratio_aggregated'],
                width=0.4, label='E:I Ratio', color=color1, alpha=0.8, edgecolor='black')
ax1.set_xlabel('Cortical Lobe', fontsize=13, fontweight='bold')
ax1.set_ylabel('E:I Ratio', fontsize=13, fontweight='bold', color=color1)
ax1.tick_params(axis='y', labelcolor=color1)
ax1.set_xticks(x_pos)
ax1.set_xticklabels(labels, rotation=45, ha='right', fontsize=11)

# rCMRGlc bars
ax2 = ax1.twinx()
color2 = '#A23B72'
bars2 = ax2.bar(x_pos + 0.2, lobes_sorted['rCMRGlc_Mean'],
                width=0.4, label='rCMRGlc', color=color2, alpha=0.8, edgecolor='black')
ax2.set_ylabel('rCMRGlc (μmol/100g/min)', fontsize=13, fontweight='bold', color=color2)
ax2.tick_params(axis='y', labelcolor=color2)

# Add value labels
for i, (idx, row) in enumerate(lobes_sorted.iterrows()):
    ax1.text(i - 0.2, row['EI_Ratio_aggregated'] + 0.1, f"{row['EI_Ratio_aggregated']:.2f}",
             ha='center', fontsize=9, fontweight='bold', color=color1)
    ax2.text(i + 0.2, row['rCMRGlc_Mean'] + 1, f"{row['rCMRGlc_Mean']:.1f}",
             ha='center', fontsize=9, fontweight='bold', color=color2)

ax1.set_title('E:I Ratio vs rCMRGlc: Side-by-Side Comparison\nWhole Human Brain Atlas',
              fontsize=14, fontweight='bold', pad=20)

# Combined legend
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left', fontsize=11)

ax1.grid(axis='y', alpha=0.3)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'ei_rcmrglc_comparison.png', dpi=300, bbox_inches='tight')
print(f'  Saved: ei_rcmrglc_comparison.png')
plt.close()

# 4. Heatmap-style visualization of fine-grained regions
print('\n[5/5] Creating fine-grained region heatmap...')

# Filter regions with valid E:I ratios
fine_valid = fine_grained[fine_grained['EI_Ratio'].notna()].copy()
fine_valid = fine_valid.sort_values('EI_Ratio', ascending=False).head(20)  # Top 20

fig, ax = plt.subplots(figsize=(12, 10))

y_pos = np.arange(len(fine_valid))
colors = plt.cm.RdYlGn_r(np.linspace(0.2, 0.9, len(fine_valid)))
bars = ax.barh(y_pos, fine_valid['EI_Ratio'], color=colors,
               edgecolor='black', linewidth=1, alpha=0.85)

# Add value labels
for i, (idx, row) in enumerate(fine_valid.iterrows()):
    ax.text(row['EI_Ratio'] + 0.1, i, f"{row['EI_Ratio']:.2f}  (n={row['Total']:,})",
            va='center', fontsize=9, fontweight='bold')

ax.set_yticks(y_pos)
ax.set_yticklabels(fine_valid['Region'], fontsize=10)
ax.set_xlabel('E:I Ratio', fontsize=13, fontweight='bold')
ax.set_title('Top 20 Cortical Regions by E:I Ratio\nWhole Human Brain Atlas (Fine-Grained)',
             fontsize=14, fontweight='bold', pad=20)
ax.axvline(x=3.0, color='red', linestyle='--', linewidth=2, alpha=0.5,
           label='Typical cortical E:I (3:1)')
ax.legend(fontsize=11)
ax.grid(axis='x', alpha=0.3)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'ei_ratios_fine_grained_top20.png', dpi=300, bbox_inches='tight')
print(f'  Saved: ei_ratios_fine_grained_top20.png')
plt.close()

# Create summary report
print('\nCreating summary statistics report...')
summary = f"""
E:I RATIO vs rCMRGlc - VISUALIZATION SUMMARY
{'='*80}

Generated: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}

DATASETS:
- Lobes analyzed: {len(integrated)}
- Fine-grained regions: {len(fine_grained)}
- Regions with valid E:I ratios: {len(fine_valid)}

KEY STATISTICS:
- Pearson correlation: r = {stats.pearsonr(ei, met)[0]:.3f}, p = {stats.pearsonr(ei, met)[1]:.4f}
- Linear regression: rCMRGlc = {slope:.2f} * EI_Ratio + {intercept:.1f}
- R-squared: {r_value**2:.3f}

E:I RATIOS BY LOBE (sorted by ratio):
"""
for _, row in integrated.sort_values('EI_Ratio_aggregated', ascending=False).iterrows():
    lobe_name = row['Lobe'].replace('_', ' ')
    summary += f"\n  {lobe_name:30s} E:I = {row['EI_Ratio_aggregated']:.2f}  rCMRGlc = {row['rCMRGlc_Mean']:.1f}"

summary += f"""

VISUALIZATIONS CREATED:
1. ei_vs_rcmrglc_scatter.png - Scatter plot with regression line
2. ei_ratios_by_lobe.png - Bar chart of E:I ratios
3. ei_rcmrglc_comparison.png - Dual-axis comparison
4. ei_ratios_fine_grained_top20.png - Top 20 regions by E:I ratio

All visualizations saved to: {OUTPUT_DIR}

{'='*80}
"""

with open(OUTPUT_DIR / 'visualization_summary.txt', 'w') as f:
    f.write(summary)

print(summary)

print('\n' + '='*80)
print('VISUALIZATION COMPLETE!')
print(f'Output directory: {OUTPUT_DIR}')
print('='*80)
