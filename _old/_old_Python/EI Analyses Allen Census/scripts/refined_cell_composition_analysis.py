"""
REFINED INTERNEURON SUBTYPE ANALYSIS
Using supercluster-level annotations for interneuron subtypes
"""
import pandas as pd
import numpy as np
from pathlib import Path
from scipy import stats
import matplotlib.pyplot as plt

print('='*80)
print('REFINED INTERNEURON SUBTYPE ANALYSIS')
print('Using supercluster annotations')
print('='*80)

# Paths
PROJECT_DIR = Path(r"C:\Users\michaelproulx\MetabolicBrain")
EI_DIR = PROJECT_DIR / "cellxgene_integration" / "ei_analysis"
DATA_DIR = PROJECT_DIR / "cellxgene_integration" / "data" / "cortical"
OUTPUT_DIR = EI_DIR / "data" / "cortical" / "followup_analyses"

# Load data
print('\n[1/5] Loading data...')
subregion_data = pd.read_csv(OUTPUT_DIR / 'subregions_with_functional_classification.csv')
annotations = pd.read_csv(EI_DIR / 'WHB_cluster_annotation_term.csv')
membership = pd.read_csv(EI_DIR / 'WHB_cluster_to_cluster_annotation_membership.csv')

print('\n[2/5] Analyzing supercluster annotations...')
supercluster_annot = annotations[annotations['cluster_annotation_term_set_name'] == 'supercluster'].copy()

print(f'\nSupercluster types ({len(supercluster_annot)}):')
for name in sorted(supercluster_annot['name'].unique()):
    count = len(supercluster_annot[supercluster_annot['name'] == name])
    print(f'  - {name} ({count})')

# Define interneuron developmental origin mapping
# MGE (Medial Ganglionic Eminence) -> PV and SST interneurons
# CGE (Caudal Ganglionic Eminence) -> VIP and some others
# LAMP5-LHX6 -> LAMP5 interneurons and Chandelier cells

def classify_supercluster(sc_name):
    if pd.isna(sc_name):
        return 'Unknown'
    name = str(sc_name).lower()

    if 'mge interneuron' in name:
        return 'MGE_Interneuron'  # PV + SST dominant
    elif 'cge interneuron' in name:
        return 'CGE_Interneuron'  # VIP dominant
    elif 'lamp5' in name or 'chandelier' in name:
        return 'LAMP5_Chandelier'
    elif 'intratelencephalic' in name:
        return 'Excitatory_IT'
    elif 'corticothalamic' in name or 'near-projecting' in name:
        return 'Excitatory_CT_NP'
    elif 'ca1' in name or 'ca4' in name or 'hippocampal' in name:
        return 'Hippocampal'
    elif 'astrocyte' in name:
        return 'Astrocyte'
    elif 'oligodendrocyte' in name:
        return 'Oligodendrocyte'
    elif 'microglia' in name:
        return 'Microglia'
    else:
        return 'Other'

supercluster_annot['Cell_Type_Category'] = supercluster_annot['name'].apply(classify_supercluster)

print('\n\nCell type categories:')
print(supercluster_annot['Cell_Type_Category'].value_counts())

# Build mapping
print('\n[3/5] Building cluster to supercluster mapping...')
cluster_to_supercluster = {}
for _, row in supercluster_annot.iterrows():
    members = membership[membership['cluster_annotation_term_label'] == row['label']]
    for clust_alias in members['cluster_alias'].values:
        cluster_to_supercluster[clust_alias] = row['Cell_Type_Category']

print(f'Mapped {len(cluster_to_supercluster)} clusters')

# Load cells
print('\n[4/5] Loading cell metadata and counting cell types per region...')
print('(Loading 672 MB file...)')
cell_df = pd.read_csv(DATA_DIR / 'WHB-10Xv3_cell_metadata.csv')
print(f'Loaded {len(cell_df):,} cells')

# Map to cell types
cell_df['cell_type_category'] = cell_df['cluster_alias'].map(cluster_to_supercluster).fillna('Unknown')

# Filter cortical
cell_df['region_clean'] = cell_df['region_of_interest_label'].str.replace('Human ', '')
cortical = cell_df[cell_df['anatomical_division_label'] == 'Cerebral cortex'].copy()
print(f'Cortical cells: {len(cortical):,}')

# Count by region
print('\nCalculating cell type composition by region...')
region_composition = []

for region in subregion_data['Region'].unique():
    region_cells = cortical[cortical['region_clean'] == region]

    if len(region_cells) == 0:
        continue

    total = len(region_cells)
    mge = (region_cells['cell_type_category'] == 'MGE_Interneuron').sum()
    cge = (region_cells['cell_type_category'] == 'CGE_Interneuron').sum()
    lamp5 = (region_cells['cell_type_category'] == 'LAMP5_Chandelier').sum()
    excit_it = (region_cells['cell_type_category'] == 'Excitatory_IT').sum()
    excit_ct = (region_cells['cell_type_category'] == 'Excitatory_CT_NP').sum()
    astro = (region_cells['cell_type_category'] == 'Astrocyte').sum()
    oligo = (region_cells['cell_type_category'] == 'Oligodendrocyte').sum()
    micro = (region_cells['cell_type_category'] == 'Microglia').sum()

    total_interneuron = mge + cge + lamp5
    total_excitatory = excit_it + excit_ct
    total_neuron = total_interneuron + total_excitatory
    total_glia = astro + oligo + micro

    region_composition.append({
        'Region': region,
        'Total_Cells': total,
        'MGE_Interneuron': mge,
        'CGE_Interneuron': cge,
        'LAMP5_Chandelier': lamp5,
        'Total_Interneuron': total_interneuron,
        'Excitatory_IT': excit_it,
        'Excitatory_CT_NP': excit_ct,
        'Total_Excitatory': total_excitatory,
        'Astrocyte': astro,
        'Oligodendrocyte': oligo,
        'Microglia': micro,
        'Total_Glia': total_glia,
        'Total_Neuron': total_neuron,
        # Percentages
        'MGE_pct': mge/total*100 if total > 0 else 0,
        'CGE_pct': cge/total*100 if total > 0 else 0,
        'LAMP5_pct': lamp5/total*100 if total > 0 else 0,
        'Interneuron_pct': total_interneuron/total*100 if total > 0 else 0,
        'Excitatory_pct': total_excitatory/total*100 if total > 0 else 0,
        'Glia_pct': total_glia/total*100 if total > 0 else 0,
        'Neuron_pct': total_neuron/total*100 if total > 0 else 0,
    })

comp_df = pd.DataFrame(region_composition)

# Merge with metabolic data
merged = subregion_data.merge(comp_df, on='Region', how='left')
merged.to_csv(OUTPUT_DIR / 'subregions_with_cell_composition.csv', index=False)

print('\n[5/5] Analyzing correlations with metabolism...')

print('\nCell type composition summary:')
print(f'  MGE Interneurons: {comp_df["MGE_pct"].mean():.2f}% (SD={comp_df["MGE_pct"].std():.2f})')
print(f'  CGE Interneurons: {comp_df["CGE_pct"].mean():.2f}% (SD={comp_df["CGE_pct"].std():.2f})')
print(f'  LAMP5/Chandelier: {comp_df["LAMP5_pct"].mean():.2f}% (SD={comp_df["LAMP5_pct"].std():.2f})')
print(f'  Total Interneurons: {comp_df["Interneuron_pct"].mean():.2f}% (SD={comp_df["Interneuron_pct"].std():.2f})')
print(f'  Excitatory neurons: {comp_df["Excitatory_pct"].mean():.2f}% (SD={comp_df["Excitatory_pct"].std():.2f})')
print(f'  Glia: {comp_df["Glia_pct"].mean():.2f}% (SD={comp_df["Glia_pct"].std():.2f})')

print('\nCorrelations with rCMRGlc:')
corr_results = []
for var in ['MGE_pct', 'CGE_pct', 'LAMP5_pct', 'Interneuron_pct', 'Excitatory_pct', 'Glia_pct']:
    valid = merged[[var, 'rCMRGlc_Mean']].dropna()
    if len(valid) > 3:
        r, p = stats.pearsonr(valid[var], valid['rCMRGlc_Mean'])
        sig = "***" if p < 0.05 else "**" if p < 0.01 else "*" if p < 0.1 else ""
        print(f'  {var:20s}: r={r:6.3f}, p={p:.4f} {sig}')
        corr_results.append({'Variable': var, 'r': r, 'p': p})

# Create comprehensive visualization
fig, axes = plt.subplots(3, 3, figsize=(18, 16))

# Plot 1-3: Interneuron subtypes vs metabolism
for idx, (var, title) in enumerate([
    ('MGE_pct', 'MGE Interneurons (PV+SST dominant)'),
    ('CGE_pct', 'CGE Interneurons (VIP dominant)'),
    ('LAMP5_pct', 'LAMP5/Chandelier Interneurons')
]):
    ax = axes[0, idx]
    valid = merged[[var, 'rCMRGlc_Mean', 'Region']].dropna()

    ax.scatter(valid[var], valid['rCMRGlc_Mean'],
              s=200, alpha=0.7, edgecolors='black', linewidth=1.5, c=valid['rCMRGlc_Mean'],
              cmap='viridis')

    for _, row in valid.iterrows():
        ax.annotate(row['Region'], (row[var], row['rCMRGlc_Mean']),
                   fontsize=7, alpha=0.6)

    if len(valid) > 3:
        r, p = stats.pearsonr(valid[var], valid['rCMRGlc_Mean'])
        ax.set_title(f'{title}\nr={r:.3f}, p={p:.3f}', fontsize=11, fontweight='bold')

        # Add regression line if significant
        if p < 0.1:
            z = np.polyfit(valid[var], valid['rCMRGlc_Mean'], 1)
            p_fit = np.poly1d(z)
            x_line = np.linspace(valid[var].min(), valid[var].max(), 100)
            ax.plot(x_line, p_fit(x_line), 'r--', linewidth=2, alpha=0.7)

    ax.set_xlabel(f'% of cells', fontsize=10, fontweight='bold')
    ax.set_ylabel('rCMRGlc (μmol/100g/min)', fontsize=10, fontweight='bold')
    ax.grid(True, alpha=0.3)

# Plot 4-6: Neuron/Glia composition
for idx, (var, title) in enumerate([
    ('Excitatory_pct', 'Excitatory Neurons'),
    ('Interneuron_pct', 'Total Interneurons'),
    ('Glia_pct', 'Glia (Astro+Oligo+Micro)')
]):
    ax = axes[1, idx]
    valid = merged[[var, 'rCMRGlc_Mean', 'Region']].dropna()

    ax.scatter(valid[var], valid['rCMRGlc_Mean'],
              s=200, alpha=0.7, edgecolors='black', linewidth=1.5, c=valid['rCMRGlc_Mean'],
              cmap='plasma')

    for _, row in valid.iterrows():
        ax.annotate(row['Region'], (row[var], row['rCMRGlc_Mean']),
                   fontsize=7, alpha=0.6)

    if len(valid) > 3:
        r, p = stats.pearsonr(valid[var], valid['rCMRGlc_Mean'])
        ax.set_title(f'{title}\nr={r:.3f}, p={p:.3f}', fontsize=11, fontweight='bold')

        if p < 0.1:
            z = np.polyfit(valid[var], valid['rCMRGlc_Mean'], 1)
            p_fit = np.poly1d(z)
            x_line = np.linspace(valid[var].min(), valid[var].max(), 100)
            ax.plot(x_line, p_fit(x_line), 'r--', linewidth=2, alpha=0.7)

    ax.set_xlabel(f'% of cells', fontsize=10, fontweight='bold')
    ax.set_ylabel('rCMRGlc (μmol/100g/min)', fontsize=10, fontweight='bold')
    ax.grid(True, alpha=0.3)

# Plot 7: Stacked composition by region
ax = axes[2, 0]
valid_regions = merged[['Region', 'Excitatory_pct', 'Interneuron_pct', 'Glia_pct', 'rCMRGlc_Mean']].dropna()
valid_regions = valid_regions.sort_values('rCMRGlc_Mean', ascending=False)

x_pos = np.arange(len(valid_regions))
exc_vals = valid_regions['Excitatory_pct'].values
inh_vals = valid_regions['Interneuron_pct'].values
glia_vals = valid_regions['Glia_pct'].values

ax.bar(x_pos, exc_vals, label='Excitatory', color='#E63946', alpha=0.8)
ax.bar(x_pos, inh_vals, bottom=exc_vals, label='Interneuron', color='#F1A208', alpha=0.8)
ax.bar(x_pos, glia_vals, bottom=exc_vals+inh_vals, label='Glia', color='#457B9D', alpha=0.8)

ax.set_xticks(x_pos)
ax.set_xticklabels(valid_regions['Region'].values, rotation=90, fontsize=8)
ax.set_ylabel('% of cells', fontsize=11, fontweight='bold')
ax.set_title('Cell Composition by Region\n(sorted by rCMRGlc)', fontsize=11, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(axis='y', alpha=0.3)

# Plot 8: Interneuron subtype breakdown
ax = axes[2, 1]
valid_inh = merged[['Region', 'MGE_pct', 'CGE_pct', 'LAMP5_pct', 'rCMRGlc_Mean']].dropna()
valid_inh = valid_inh.sort_values('rCMRGlc_Mean', ascending=False)

x_pos = np.arange(len(valid_inh))
mge_vals = valid_inh['MGE_pct'].values
cge_vals = valid_inh['CGE_pct'].values
lamp5_vals = valid_inh['LAMP5_pct'].values

ax.bar(x_pos, mge_vals, label='MGE (PV+SST)', color='#2A9D8F', alpha=0.8)
ax.bar(x_pos, cge_vals, bottom=mge_vals, label='CGE (VIP)', color='#E76F51', alpha=0.8)
ax.bar(x_pos, lamp5_vals, bottom=mge_vals+cge_vals, label='LAMP5/Chan', color='#9D4EDD', alpha=0.8)

ax.set_xticks(x_pos)
ax.set_xticklabels(valid_inh['Region'].values, rotation=90, fontsize=8)
ax.set_ylabel('% of cells', fontsize=11, fontweight='bold')
ax.set_title('Interneuron Subtypes by Region\n(sorted by rCMRGlc)', fontsize=11, fontweight='bold')
ax.legend(fontsize=9)
ax.grid(axis='y', alpha=0.3)

# Plot 9: Correlation heatmap
ax = axes[2, 2]
corr_df = pd.DataFrame(corr_results)
if len(corr_df) > 0:
    vars_sorted = corr_df.sort_values('r', ascending=False)
    y_pos = np.arange(len(vars_sorted))

    colors = ['green' if p < 0.05 else 'orange' if p < 0.1 else 'gray'
              for p in vars_sorted['p']]

    ax.barh(y_pos, vars_sorted['r'], color=colors, alpha=0.7, edgecolor='black', linewidth=1.5)
    ax.set_yticks(y_pos)
    ax.set_yticklabels([v.replace('_pct', '').replace('_', ' ') for v in vars_sorted['Variable']], fontsize=10)
    ax.set_xlabel('Pearson r', fontsize=11, fontweight='bold')
    ax.set_title('Correlations with rCMRGlc\nGreen=p<0.05, Orange=p<0.1', fontsize=11, fontweight='bold')
    ax.axvline(0, color='black', linewidth=1)
    ax.grid(axis='x', alpha=0.3)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'analysis2_refined_cell_composition.png', dpi=300, bbox_inches='tight')
print('\n  Saved: analysis2_refined_cell_composition.png')
plt.close()

print('\n' + '='*80)
print('REFINED CELL COMPOSITION ANALYSIS COMPLETE!')
print('='*80)
print(f'\nOutput: {OUTPUT_DIR}')
print('\nKey files:')
print('  - subregions_with_cell_composition.csv')
print('  - analysis2_refined_cell_composition.png')
print('='*80)
