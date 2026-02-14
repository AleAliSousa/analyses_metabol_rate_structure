"""
PER-LOBE ANALYSIS & METABOLIC CLUSTERING EXPLORATION
"""
import pandas as pd
import numpy as np
from pathlib import Path
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

print('='*80)
print('PER-LOBE ANALYSIS & METABOLIC CLUSTERING')
print('='*80)

# Paths
PROJECT_DIR = Path(r"C:\Users\michaelproulx\MetabolicBrain")
EI_DIR = PROJECT_DIR / "cellxgene_integration" / "ei_analysis"
DATA_DIR = EI_DIR / "data" / "cortical"
OUTPUT_DIR = DATA_DIR / "per_lobe_analysis"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Load data
print('\n[1/4] Loading subregion data...')
subregion_data = pd.read_csv(DATA_DIR / 'subregion_analysis' / 'subregion_ei_rcmrglc_table.csv')
print(f'  Loaded {len(subregion_data)} subregions across {subregion_data["Lobe"].nunique()} lobes')

# Part 1: Per-lobe scatter plots
print('\n[2/4] Creating per-lobe scatter plots...')

lobes = subregion_data['Lobe'].unique()
lobe_stats = []

fig, axes = plt.subplots(2, 3, figsize=(18, 12))
axes = axes.flatten()

for idx, lobe in enumerate(sorted(lobes)):
    ax = axes[idx]
    lobe_data = subregion_data[subregion_data['Lobe'] == lobe]

    # Get E:I values for this lobe
    ei_values = lobe_data['EI_Ratio'].values
    n_regions = len(lobe_data)

    # Stats
    mean_ei = ei_values.mean()
    std_ei = ei_values.std()
    min_ei = ei_values.min()
    max_ei = ei_values.max()
    cv_ei = (std_ei / mean_ei * 100) if mean_ei > 0 else 0  # Coefficient of variation

    rcmrglc = lobe_data['rCMRGlc_Mean'].iloc[0]

    # Scatter plot - E:I on x-axis, just jitter on y for visibility
    y_jitter = np.random.normal(0, 0.02, n_regions)

    ax.scatter(ei_values, y_jitter, s=lobe_data['Total']/400,
              alpha=0.7, edgecolors='black', linewidth=1.5, color='steelblue')

    # Add region labels
    for _, row in lobe_data.iterrows():
        ax.annotate(row['Region'], (row['EI_Ratio'], np.random.normal(0, 0.02)),
                   fontsize=9, ha='center', alpha=0.8)

    # Vertical line at mean
    ax.axvline(mean_ei, color='red', linestyle='--', linewidth=2, alpha=0.7,
              label=f'Mean E:I = {mean_ei:.2f}')

    # Shade the range
    ax.axvspan(min_ei, max_ei, alpha=0.1, color='gray')

    ax.set_ylim(-0.15, 0.15)
    ax.set_yticks([])
    ax.set_xlabel('E:I Ratio', fontsize=12, fontweight='bold')
    ax.set_title(f'{lobe.replace("_", " ")}\nrCMRGlc = {rcmrglc:.1f}  |  n = {n_regions} regions',
                fontsize=13, fontweight='bold', pad=10)
    ax.legend(loc='upper right', fontsize=10)
    ax.grid(axis='x', alpha=0.3)

    # Add stats box
    stats_text = f'Range: {min_ei:.2f} - {max_ei:.2f}\nSD: {std_ei:.2f}\nCV: {cv_ei:.1f}%'
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
           fontsize=9, verticalalignment='top',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # Store stats
    lobe_stats.append({
        'Lobe': lobe,
        'rCMRGlc': rcmrglc,
        'n_regions': n_regions,
        'Mean_EI': mean_ei,
        'SD_EI': std_ei,
        'CV_EI': cv_ei,
        'Min_EI': min_ei,
        'Max_EI': max_ei,
        'Range_EI': max_ei - min_ei
    })

# Hide extra subplot
if len(lobes) < 6:
    axes[5].set_visible(False)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'per_lobe_ei_distribution.png', dpi=300, bbox_inches='tight')
print(f'  Saved: per_lobe_ei_distribution.png')
plt.close()

# Create stats table
lobe_stats_df = pd.DataFrame(lobe_stats).sort_values('rCMRGlc', ascending=False)
lobe_stats_df.to_csv(OUTPUT_DIR / 'per_lobe_statistics.csv', index=False)

print('\n' + '='*80)
print('PER-LOBE E:I STATISTICS')
print('='*80)
for _, row in lobe_stats_df.iterrows():
    lobe_clean = row['Lobe'].replace('_', ' ')
    print(f"\n{lobe_clean}:")
    print(f"  rCMRGlc: {row['rCMRGlc']:.1f}")
    print(f"  Regions: {row['n_regions']}")
    print(f"  E:I Mean: {row['Mean_EI']:.2f} ± {row['SD_EI']:.2f}")
    print(f"  E:I Range: {row['Min_EI']:.2f} - {row['Max_EI']:.2f} (span: {row['Range_EI']:.2f})")
    print(f"  Coefficient of Variation: {row['CV_EI']:.1f}%")

# Part 2: Metabolic clustering analysis
print('\n[3/4] Analyzing metabolic clustering...')

# Define metabolic groups
HIGH_METABOLISM = [35.3, 35.8]
LOW_METABOLISM = [30.3, 30.5]

subregion_data['Metabolic_Group'] = subregion_data['rCMRGlc_Mean'].apply(
    lambda x: 'High (35+)' if x in HIGH_METABOLISM else 'Low (30-31)')

high_met = subregion_data[subregion_data['Metabolic_Group'] == 'High (35+)']
low_met = subregion_data[subregion_data['Metabolic_Group'] == 'Low (30-31)']

print(f"\nMetabolic Groups:")
print(f"  High metabolism (rCMRGlc ~35): {len(high_met)} subregions")
print(f"    Lobes: {', '.join(high_met['Lobe'].unique())}")
print(f"    Mean E:I: {high_met['EI_Ratio'].mean():.2f} ± {high_met['EI_Ratio'].std():.2f}")
print(f"    E:I Range: {high_met['EI_Ratio'].min():.2f} - {high_met['EI_Ratio'].max():.2f}")

print(f"\n  Low metabolism (rCMRGlc ~30): {len(low_met)} subregions")
print(f"    Lobes: {', '.join(low_met['Lobe'].unique())}")
print(f"    Mean E:I: {low_met['EI_Ratio'].mean():.2f} ± {low_met['EI_Ratio'].std():.2f}")
print(f"    E:I Range: {low_met['EI_Ratio'].min():.2f} - {low_met['EI_Ratio'].max():.2f}")

# Statistical test: Do high-metabolism regions have different E:I than low-metabolism?
t_stat, p_val = stats.ttest_ind(high_met['EI_Ratio'], low_met['EI_Ratio'])
mann_u_stat, mann_p = stats.mannwhitneyu(high_met['EI_Ratio'], low_met['EI_Ratio'])

print(f"\nStatistical Comparison:")
print(f"  T-test: t = {t_stat:.3f}, p = {p_val:.4f}")
print(f"  Mann-Whitney U: U = {mann_u_stat:.1f}, p = {mann_p:.4f}")
print(f"  Result: {'SIGNIFICANT' if p_val < 0.05 else 'NOT significant'} difference in E:I between groups")

# Create visualization comparing metabolic groups
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# Plot 1: Box plot comparison
ax = axes[0]
data_to_plot = [high_met['EI_Ratio'], low_met['EI_Ratio']]
bp = ax.boxplot(data_to_plot, labels=['High Metabolism\n(rCMRGlc ~35)', 'Low Metabolism\n(rCMRGlc ~30)'],
                patch_artist=True, showmeans=True)
bp['boxes'][0].set_facecolor('coral')
bp['boxes'][1].set_facecolor('lightblue')
ax.set_ylabel('E:I Ratio', fontsize=13, fontweight='bold')
ax.set_title(f'E:I Ratio by Metabolic Group\np = {p_val:.4f}', fontsize=14, fontweight='bold')
ax.grid(axis='y', alpha=0.3)

# Add n values
ax.text(1, ax.get_ylim()[1]*0.95, f'n={len(high_met)}', ha='center', fontsize=11, fontweight='bold')
ax.text(2, ax.get_ylim()[1]*0.95, f'n={len(low_met)}', ha='center', fontsize=11, fontweight='bold')

# Plot 2: Violin plot with individual points
ax = axes[1]
parts = ax.violinplot([high_met['EI_Ratio'], low_met['EI_Ratio']], positions=[1, 2],
                       showmeans=True, showmedians=True)
for pc, color in zip(parts['bodies'], ['coral', 'lightblue']):
    pc.set_facecolor(color)
    pc.set_alpha(0.7)

# Overlay points
ax.scatter(np.ones(len(high_met)) + np.random.normal(0, 0.04, len(high_met)),
          high_met['EI_Ratio'], alpha=0.5, s=50, color='darkred', edgecolors='black', linewidth=0.5)
ax.scatter(2*np.ones(len(low_met)) + np.random.normal(0, 0.04, len(low_met)),
          low_met['EI_Ratio'], alpha=0.5, s=50, color='darkblue', edgecolors='black', linewidth=0.5)

ax.set_xticks([1, 2])
ax.set_xticklabels(['High Metabolism\n(rCMRGlc ~35)', 'Low Metabolism\n(rCMRGlc ~30)'])
ax.set_ylabel('E:I Ratio', fontsize=13, fontweight='bold')
ax.set_title('E:I Distribution by Metabolic Group', fontsize=14, fontweight='bold')
ax.grid(axis='y', alpha=0.3)

# Plot 3: Histogram comparison
ax = axes[2]
ax.hist(high_met['EI_Ratio'], bins=8, alpha=0.6, label=f'High Met (n={len(high_met)})',
       color='coral', edgecolor='black', linewidth=1.5)
ax.hist(low_met['EI_Ratio'], bins=8, alpha=0.6, label=f'Low Met (n={len(low_met)})',
       color='lightblue', edgecolor='black', linewidth=1.5)
ax.set_xlabel('E:I Ratio', fontsize=13, fontweight='bold')
ax.set_ylabel('Frequency', fontsize=13, fontweight='bold')
ax.set_title('E:I Distribution Overlap', fontsize=14, fontweight='bold')
ax.legend(fontsize=11)
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'metabolic_group_comparison.png', dpi=300, bbox_inches='tight')
print(f'  Saved: metabolic_group_comparison.png')
plt.close()

# Part 3: Alternative analyses
print('\n[4/4] Running alternative analyses...')

# Analysis 1: Excitatory percentage instead of ratio
subregion_data['E_percentage'] = subregion_data['E_pct']

fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# E% vs rCMRGlc
ax = axes[0]
colors = {'Frontal_lobe': '#E63946', 'Parietal_lobe': '#F1A208', 'Temporal_lobe': '#2A9D8F',
          'Occipital_lobe': '#457B9D', 'Insular_lobe': '#E76F51'}
for lobe in subregion_data['Lobe'].unique():
    lobe_data = subregion_data[subregion_data['Lobe'] == lobe]
    ax.scatter(lobe_data['E_percentage'], lobe_data['rCMRGlc_Mean'],
              s=lobe_data['Total']/500, alpha=0.7, color=colors.get(lobe, '#999'),
              edgecolors='black', linewidth=1.5, label=lobe.replace('_', ' '))

r_epct, p_epct = stats.pearsonr(subregion_data['E_percentage'], subregion_data['rCMRGlc_Mean'])
ax.set_xlabel('Excitatory Percentage (%)', fontsize=13, fontweight='bold')
ax.set_ylabel('rCMRGlc (μmol/100g/min)', fontsize=13, fontweight='bold')
ax.set_title(f'Excitatory % vs rCMRGlc\nr = {r_epct:.3f}, p = {p_epct:.4f}', fontsize=14, fontweight='bold')
ax.legend(fontsize=9, loc='best')
ax.grid(True, alpha=0.3)

# Within-group E:I variance
ax = axes[1]
lobe_order = lobe_stats_df.sort_values('rCMRGlc', ascending=False)['Lobe'].values
y_pos = np.arange(len(lobe_order))
cv_values = [lobe_stats_df[lobe_stats_df['Lobe']==lobe]['CV_EI'].values[0] for lobe in lobe_order]
colors_ordered = [colors.get(lobe, '#999') for lobe in lobe_order]

bars = ax.barh(y_pos, cv_values, color=colors_ordered, alpha=0.7, edgecolor='black', linewidth=1.5)
ax.set_yticks(y_pos)
ax.set_yticklabels([lobe.replace('_', ' ') for lobe in lobe_order])
ax.set_xlabel('Coefficient of Variation (%)', fontsize=13, fontweight='bold')
ax.set_title('Within-Lobe E:I Variability\n(Higher = more heterogeneous)', fontsize=14, fontweight='bold')
ax.grid(axis='x', alpha=0.3)

# Add rCMRGlc on right y-axis
ax2 = ax.twinx()
ax2.set_yticks(y_pos)
rcmrglc_labels = [f"{lobe_stats_df[lobe_stats_df['Lobe']==lobe]['rCMRGlc'].values[0]:.1f}"
                  for lobe in lobe_order]
ax2.set_yticklabels(rcmrglc_labels)
ax2.set_ylabel('rCMRGlc', fontsize=12, fontweight='bold', rotation=270, labelpad=20)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'alternative_analyses.png', dpi=300, bbox_inches='tight')
print(f'  Saved: alternative_analyses.png')
plt.close()

print(f"\nAlternative Analysis Results:")
print(f"  Excitatory % vs rCMRGlc: r = {r_epct:.3f}, p = {p_epct:.4f}")

# Analysis 2: Cell density patterns
print(f"\nCell Density Analysis:")
for group_name, group_data in [('High Metabolism', high_met), ('Low Metabolism', low_met)]:
    avg_total = group_data['Total'].mean()
    avg_exc = group_data['Excitatory'].mean()
    avg_inh = group_data['Inhibitory'].mean()
    print(f"\n  {group_name}:")
    print(f"    Avg cells per region: {avg_total:,.0f}")
    print(f"    Avg excitatory: {avg_exc:,.0f}")
    print(f"    Avg inhibitory: {avg_inh:,.0f}")

# Save summary report
with open(OUTPUT_DIR / 'clustering_analysis_summary.txt', 'w') as f:
    f.write('METABOLIC CLUSTERING ANALYSIS\n')
    f.write('='*80 + '\n\n')
    f.write('OBSERVATION: Two metabolic clusters exist\n')
    f.write('  High metabolism: rCMRGlc ~35 (Frontal, Parietal, Occipital)\n')
    f.write('  Low metabolism: rCMRGlc ~30 (Temporal, Insular)\n\n')
    f.write(f'High metabolism regions (n={len(high_met)}):\n')
    f.write(f'  Mean E:I: {high_met["EI_Ratio"].mean():.2f} ± {high_met["EI_Ratio"].std():.2f}\n')
    f.write(f'  Range: {high_met["EI_Ratio"].min():.2f} - {high_met["EI_Ratio"].max():.2f}\n\n')
    f.write(f'Low metabolism regions (n={len(low_met)}):\n')
    f.write(f'  Mean E:I: {low_met["EI_Ratio"].mean():.2f} ± {low_met["EI_Ratio"].std():.2f}\n')
    f.write(f'  Range: {low_met["EI_Ratio"].min():.2f} - {low_met["EI_Ratio"].max():.2f}\n\n')
    f.write(f'Statistical comparison:\n')
    f.write(f'  T-test: p = {p_val:.4f}\n')
    f.write(f'  Result: {"SIGNIFICANT" if p_val < 0.05 else "NOT significant"}\n\n')
    f.write(f'Alternative correlation (Excitatory % vs rCMRGlc):\n')
    f.write(f'  r = {r_epct:.3f}, p = {p_epct:.4f}\n')

print('\n' + '='*80)
print('ANALYSIS COMPLETE!')
print('='*80)
print(f'\nOutput directory: {OUTPUT_DIR}')
print('\nFiles created:')
print('  1. per_lobe_ei_distribution.png - Individual lobe E:I distributions')
print('  2. per_lobe_statistics.csv - Per-lobe statistics table')
print('  3. metabolic_group_comparison.png - High vs Low metabolism comparison')
print('  4. alternative_analyses.png - E% and within-lobe variability')
print('  5. clustering_analysis_summary.txt - Summary report')
print('='*80)
