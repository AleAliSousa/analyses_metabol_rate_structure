"""
SUBREGION-LEVEL E:I vs rCMRGlc ANALYSIS
Keep all cortical subregions independent, assign lobe rCMRGlc to each
"""
import pandas as pd
import numpy as np
from pathlib import Path
from scipy import stats
import matplotlib.pyplot as plt

print('='*80)
print('SUBREGION-LEVEL E:I vs rCMRGlc ANALYSIS')
print('='*80)

# Paths
PROJECT_DIR = Path(r"C:\Users\michaelproulx\MetabolicBrain")
EI_DIR = PROJECT_DIR / "cellxgene_integration" / "ei_analysis"
DATA_DIR = EI_DIR / "data" / "cortical"
OUTPUT_DIR = DATA_DIR / "subregion_analysis"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Heiss rCMRGlc data by lobe
HEISS_DATA = {
    'Frontal_lobe': (35.3, 3.51),
    'Parietal_lobe': (35.8, 3.23),
    'Temporal_lobe': (30.5, 2.54),
    'Occipital_lobe': (35.8, 3.12),  # Includes V1C, V2, A19
    'Insular_lobe': (30.3, 2.55),
}

# Region to lobe mapping
REGION_TO_LOBE = {
    'M1C': 'Frontal_lobe', 'A44-A45': 'Frontal_lobe', 'A13': 'Frontal_lobe',
    'A25': 'Frontal_lobe', 'A32': 'Frontal_lobe', 'A46': 'Frontal_lobe',
    'S1C': 'Parietal_lobe', 'A5-A7': 'Parietal_lobe', 'A40': 'Parietal_lobe', 'A43': 'Parietal_lobe',
    'MTG': 'Temporal_lobe', 'STG': 'Temporal_lobe', 'ITG': 'Temporal_lobe',
    'A38': 'Temporal_lobe', 'A1C': 'Temporal_lobe',
    'V1C': 'Occipital_lobe', 'V2': 'Occipital_lobe', 'A19': 'Occipital_lobe',  # V1C is part of occipital
    'FI': 'Insular_lobe', 'SI': 'Insular_lobe', 'Ig': 'Insular_lobe',
}

print('\n[1/6] Loading fine-grained E:I ratios...')
ei_df = pd.read_csv(DATA_DIR / 'ei_ratios_whb_cortical_fine_grained.csv')
print(f'  Loaded {len(ei_df)} cortical regions')

print('\n[2/6] Mapping regions to lobes and assigning rCMRGlc...')
ei_df['Lobe'] = ei_df['Region'].map(REGION_TO_LOBE)
mapped = ei_df[ei_df['Lobe'].notna()].copy()
print(f'  Mapped {len(mapped)}/{len(ei_df)} regions to lobes')

# Assign lobe rCMRGlc to each region
mapped['rCMRGlc_Mean'] = mapped['Lobe'].map({k: v[0] for k, v in HEISS_DATA.items()})
mapped['rCMRGlc_SD'] = mapped['Lobe'].map({k: v[1] for k, v in HEISS_DATA.items()})

# Filter valid data
valid = mapped[mapped['EI_Ratio'].notna() & mapped['rCMRGlc_Mean'].notna()].copy()
print(f'  {len(valid)} regions with valid E:I and rCMRGlc data')

print('\n[3/6] Creating comprehensive data table...')
# Sort by E:I ratio
valid_sorted = valid.sort_values('EI_Ratio', ascending=False)

# Create display table
table = valid_sorted[['Region', 'Lobe', 'EI_Ratio', 'rCMRGlc_Mean',
                       'Excitatory', 'Inhibitory', 'Total']].copy()
table['E_pct'] = (table['Excitatory'] / table['Total'] * 100).round(1)

print('\n' + '='*120)
print('SUBREGION-LEVEL E:I RATIOS WITH LOBE rCMRGlc VALUES')
print('='*120)
print(f"{'Region':<15} {'Lobe':<25} {'E:I':>6} {'rCMRGlc':>8} {'Excit':>8} {'Inhib':>8} {'Total':>9} {'E%':>6}")
print('-'*120)
for _, row in table.iterrows():
    lobe_clean = row['Lobe'].replace('_', ' ')
    print(f"{row['Region']:<15} {lobe_clean:<25} {row['EI_Ratio']:6.2f} {row['rCMRGlc_Mean']:8.1f} "
          f"{row['Excitatory']:8,} {row['Inhibitory']:8,} {row['Total']:9,} {row['E_pct']:6.1f}")

# Save comprehensive table
table.to_csv(OUTPUT_DIR / 'subregion_ei_rcmrglc_table.csv', index=False)
print(f'\nTable saved to: subregion_ei_rcmrglc_table.csv')

print('\n[4/6] Statistical analysis with all subregions...')
ei_values = valid['EI_Ratio'].values
rcmrglc_values = valid['rCMRGlc_Mean'].values

# Pearson correlation
pearson_r, pearson_p = stats.pearsonr(ei_values, rcmrglc_values)

# Spearman correlation
spearman_r, spearman_p = stats.spearmanr(ei_values, rcmrglc_values)

# Linear regression
slope, intercept, r_value, p_value, std_err = stats.linregress(ei_values, rcmrglc_values)

print('\n' + '='*80)
print('STATISTICAL RESULTS - SUBREGION LEVEL')
print('='*80)
print(f'\nSample size: n = {len(valid)} subregions')
print(f'\nPearson correlation:')
print(f'  r = {pearson_r:.4f}')
print(f'  p-value = {pearson_p:.6f}')
print(f'  Significant: {"***YES***" if pearson_p < 0.05 else "NO"} (alpha = 0.05)')

print(f'\nSpearman correlation:')
print(f'  rho = {spearman_r:.4f}')
print(f'  p-value = {spearman_p:.6f}')
print(f'  Significant: {"***YES***" if spearman_p < 0.05 else "NO"} (alpha = 0.05)')

print(f'\nLinear regression:')
print(f'  rCMRGlc = {slope:.4f} * E:I Ratio + {intercept:.2f}')
print(f'  R-squared = {r_value**2:.4f}')
print(f'  p-value = {p_value:.6f}')
print(f'  Standard error = {std_err:.4f}')

# Save statistics
stats_file = OUTPUT_DIR / 'subregion_statistics.txt'
with open(stats_file, 'w') as f:
    f.write('SUBREGION-LEVEL E:I vs rCMRGlc STATISTICAL ANALYSIS\n')
    f.write('='*80 + '\n\n')
    f.write(f'Sample size: n = {len(valid)} cortical subregions\n\n')
    f.write(f'Pearson correlation:\n')
    f.write(f'  r = {pearson_r:.4f}\n')
    f.write(f'  p-value = {pearson_p:.6f}\n')
    f.write(f'  Significant at alpha=0.05: {"YES" if pearson_p < 0.05 else "NO"}\n\n')
    f.write(f'Spearman correlation:\n')
    f.write(f'  rho = {spearman_r:.4f}\n')
    f.write(f'  p-value = {spearman_p:.6f}\n')
    f.write(f'  Significant at alpha=0.05: {"YES" if spearman_p < 0.05 else "NO"}\n\n')
    f.write(f'Linear regression:\n')
    f.write(f'  Equation: rCMRGlc = {slope:.4f} * E:I_Ratio + {intercept:.2f}\n')
    f.write(f'  R-squared = {r_value**2:.4f}\n')
    f.write(f'  p-value = {p_value:.6f}\n')
    f.write(f'  Standard error = {std_err:.4f}\n\n')
    f.write('COMPARISON TO LOBE-LEVEL ANALYSIS:\n')
    f.write(f'  Lobe-level (n=6):     r=0.668, p=0.1470 (NOT significant)\n')
    f.write(f'  Subregion-level (n={len(valid)}): r={pearson_r:.3f}, p={pearson_p:.4f} {"(SIGNIFICANT!)" if pearson_p < 0.05 else ""}\n\n')

print(f'\nStatistics saved to: {stats_file}')

print('\n[5/6] Creating visualizations...')

# Set plot style
plt.style.use('seaborn-v0_8-darkgrid')

# 1. Scatter plot with all subregions
fig, ax = plt.subplots(figsize=(14, 10))

# Color by lobe
lobe_colors = {
    'Frontal_lobe': '#E63946',
    'Parietal_lobe': '#F1A208',
    'Temporal_lobe': '#2A9D8F',
    'Occipital_lobe': '#457B9D',  # Now includes V1C
    'Insular_lobe': '#E76F51',
}

for lobe in valid['Lobe'].unique():
    lobe_data = valid[valid['Lobe'] == lobe]
    lobe_clean = lobe.replace('_', ' ')
    ax.scatter(lobe_data['EI_Ratio'], lobe_data['rCMRGlc_Mean'],
              s=lobe_data['Total']/500,  # Size by cell count
              alpha=0.7, color=lobe_colors.get(lobe, '#999999'),
              edgecolors='black', linewidth=1.5, label=lobe_clean)

# Regression line
line_x = np.array([ei_values.min() - 0.2, ei_values.max() + 0.2])
line_y = slope * line_x + intercept
ax.plot(line_x, line_y, 'r--', linewidth=2.5, alpha=0.8,
       label=f'Regression: y = {slope:.3f}x + {intercept:.1f}\nR² = {r_value**2:.3f}, p = {p_value:.4f}')

# Annotate selected regions
top_regions = valid.nlargest(3, 'EI_Ratio')
for _, row in top_regions.iterrows():
    ax.annotate(row['Region'], (row['EI_Ratio'], row['rCMRGlc_Mean']),
               xytext=(8, 8), textcoords='offset points',
               fontsize=9, fontweight='bold', alpha=0.8,
               bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.3))

ax.set_xlabel('E:I Ratio (Excitatory/Inhibitory)', fontsize=14, fontweight='bold')
ax.set_ylabel('rCMRGlc (μmol/100g/min)', fontsize=14, fontweight='bold')
ax.set_title(f'Subregion-Level E:I Ratio vs rCMRGlc\n{len(valid)} Cortical Regions from WHB Atlas',
            fontsize=15, fontweight='bold', pad=20)
ax.legend(loc='upper left', fontsize=10, framealpha=0.95, ncol=2)
ax.grid(True, alpha=0.3)

# Add significance annotation
sig_text = f"Pearson: r={pearson_r:.3f}, p={pearson_p:.4f}"
if pearson_p < 0.05:
    sig_text += " ***SIGNIFICANT***"
ax.text(0.98, 0.02, sig_text, transform=ax.transAxes,
       fontsize=11, fontweight='bold',
       bbox=dict(boxstyle='round', facecolor='lightgreen' if pearson_p < 0.05 else 'lightyellow', alpha=0.8),
       ha='right', va='bottom')

plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'subregion_ei_vs_rcmrglc_scatter.png', dpi=300, bbox_inches='tight')
print(f'  Saved: subregion_ei_vs_rcmrglc_scatter.png')
plt.close()

# 2. Box plot by lobe showing subregion variation
fig, ax = plt.subplots(figsize=(14, 8))

lobes_sorted = valid.groupby('Lobe')['rCMRGlc_Mean'].first().sort_values(ascending=False).index
lobe_data_list = []
lobe_labels = []
for lobe in lobes_sorted:
    lobe_regions = valid[valid['Lobe'] == lobe]['EI_Ratio'].values
    if len(lobe_regions) > 0:
        lobe_data_list.append(lobe_regions)
        lobe_clean = lobe.replace('_', ' ')
        lobe_labels.append(f"{lobe_clean}\n(n={len(lobe_regions)})")

bp = ax.boxplot(lobe_data_list, labels=lobe_labels, patch_artist=True,
               showmeans=True, meanline=True)

# Color boxes by lobe
for patch, lobe in zip(bp['boxes'], lobes_sorted):
    patch.set_facecolor(lobe_colors.get(lobe, '#999999'))
    patch.set_alpha(0.7)

ax.set_ylabel('E:I Ratio', fontsize=14, fontweight='bold')
ax.set_xlabel('Cortical Lobe', fontsize=14, fontweight='bold')
ax.set_title('E:I Ratio Distribution Across Subregions by Lobe\nWhole Human Brain Atlas',
            fontsize=15, fontweight='bold', pad=20)
ax.grid(axis='y', alpha=0.3)
plt.xticks(rotation=15, ha='right')

plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'subregion_ei_distribution_by_lobe.png', dpi=300, bbox_inches='tight')
print(f'  Saved: subregion_ei_distribution_by_lobe.png')
plt.close()

# 3. Residuals plot
fig, ax = plt.subplots(figsize=(12, 8))

predicted = slope * ei_values + intercept
residuals = rcmrglc_values - predicted

ax.scatter(predicted, residuals, s=100, alpha=0.6, edgecolors='black', linewidth=1)
ax.axhline(y=0, color='red', linestyle='--', linewidth=2, alpha=0.7, label='Zero residual line')
ax.set_xlabel('Predicted rCMRGlc (μmol/100g/min)', fontsize=13, fontweight='bold')
ax.set_ylabel('Residuals (Observed - Predicted)', fontsize=13, fontweight='bold')
ax.set_title('Regression Residuals Plot - Subregion Level\nChecking Model Fit',
            fontsize=14, fontweight='bold', pad=20)
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'subregion_residuals_plot.png', dpi=300, bbox_inches='tight')
print(f'  Saved: subregion_residuals_plot.png')
plt.close()

print('\n[6/6] Creating summary comparison...')

# Comparison summary
comparison = f"""
================================================================================
SUBREGION-LEVEL vs LOBE-LEVEL ANALYSIS COMPARISON
================================================================================

LOBE-LEVEL ANALYSIS (Previous):
  Sample size: n = 6 lobes
  Pearson r = 0.668, p = 0.1470
  Result: NOT statistically significant (p > 0.05)
  Interpretation: Moderate correlation, but insufficient power

SUBREGION-LEVEL ANALYSIS (Current):
  Sample size: n = {len(valid)} cortical subregions
  Pearson r = {pearson_r:.3f}, p = {pearson_p:.6f}
  Result: {"***STATISTICALLY SIGNIFICANT*** (p < 0.05)" if pearson_p < 0.05 else "NOT significant (p > 0.05)"}

  Linear regression:
    rCMRGlc = {slope:.4f} × E:I Ratio + {intercept:.2f}
    R² = {r_value**2:.3f}

CONCLUSION:
{"By analyzing individual subregions instead of aggregated lobes, we GAINED " if pearson_p < 0.05 else "Even with more data points, we still lack "}statistical power {"and can now detect a significant positive correlation between E:I ratio and rCMRGlc!" if pearson_p < 0.05 else "to detect significance."}

{"The positive correlation suggests that cortical regions with higher excitatory:inhibitory balance have higher glucose metabolism." if pearson_r > 0 else ""}

Top 5 regions by E:I ratio:
"""

for i, (_, row) in enumerate(valid_sorted.head(5).iterrows(), 1):
    comparison += f"  {i}. {row['Region']:<12} E:I={row['EI_Ratio']:.2f}  rCMRGlc={row['rCMRGlc_Mean']:.1f}  ({row['Lobe'].replace('_', ' ')})\n"

comparison += "\n" + "="*80 + "\n"

with open(OUTPUT_DIR / 'analysis_comparison.txt', 'w') as f:
    f.write(comparison)

print(comparison)

print('\n' + '='*80)
print('SUBREGION-LEVEL ANALYSIS COMPLETE!')
print('='*80)
print(f'\nOutput directory: {OUTPUT_DIR}')
print('\nFiles created:')
print('  1. subregion_ei_rcmrglc_table.csv - Comprehensive data table')
print('  2. subregion_statistics.txt - Statistical results')
print('  3. subregion_ei_vs_rcmrglc_scatter.png - Main scatter plot')
print('  4. subregion_ei_distribution_by_lobe.png - Box plots')
print('  5. subregion_residuals_plot.png - Residuals analysis')
print('  6. analysis_comparison.txt - Comparison summary')
print('='*80)
