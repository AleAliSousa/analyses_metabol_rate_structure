"""
COMPREHENSIVE FOLLOW-UP ANALYSES
1. Sensorimotor vs Associative Classification
2. Metabolic Gene Expression Analysis
3. Cell Subtype Composition (PV, SST, VIP interneurons)
"""
import pandas as pd
import numpy as np
from pathlib import Path
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

print('='*80)
print('COMPREHENSIVE FOLLOW-UP ANALYSES')
print('Starting all 3 analyses...')
print('='*80)

# Paths
PROJECT_DIR = Path(r"C:\Users\michaelproulx\MetabolicBrain")
EI_DIR = PROJECT_DIR / "cellxgene_integration" / "ei_analysis"
DATA_DIR = PROJECT_DIR / "cellxgene_integration" / "data" / "cortical"
OUTPUT_DIR = EI_DIR / "data" / "cortical" / "followup_analyses"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Load existing data
print('\n[Setup] Loading existing data...')
subregion_data = pd.read_csv(EI_DIR / 'data' / 'cortical' / 'subregion_analysis' / 'subregion_ei_rcmrglc_table.csv')
annotations = pd.read_csv(EI_DIR / 'WHB_cluster_annotation_term.csv')
membership = pd.read_csv(EI_DIR / 'WHB_cluster_to_cluster_annotation_membership.csv')
clusters = pd.read_csv(EI_DIR / 'WHB_cluster.csv')

print(f'  Subregions: {len(subregion_data)}')
print(f'  Annotations: {len(annotations)}')
print(f'  Available annotation sets: {annotations["cluster_annotation_term_set_name"].unique()}')

print('\n' + '='*80)
print('ANALYSIS 1: SENSORIMOTOR vs ASSOCIATIVE CLASSIFICATION')
print('='*80)

# Define functional classifications based on known neuroanatomy
FUNCTIONAL_CLASSIFICATION = {
    # PRIMARY SENSORIMOTOR (expected high metabolism)
    'V1C': 'Primary_Sensory',      # Primary visual cortex
    'V2': 'Secondary_Sensory',      # Visual area 2
    'M1C': 'Primary_Motor',         # Primary motor cortex
    'S1C': 'Primary_Sensory',       # Primary somatosensory
    'A1C': 'Primary_Sensory',       # Primary auditory

    # HIGHER-ORDER SENSORIMOTOR
    'A5-A7': 'Sensorimotor_Assoc',  # Superior parietal (visuomotor integration)
    'A40': 'Sensorimotor_Assoc',     # Supramarginal gyrus (language/sensorimotor)
    'A43': 'Sensorimotor_Assoc',     # Parietal operculum (somatosensory association)
    'A19': 'Secondary_Sensory',      # Extrastriate visual

    # PREMOTOR/MOTOR ASSOCIATION
    'A46': 'Prefrontal_Assoc',       # Dorsolateral prefrontal
    'A44-A45': 'Motor_Assoc',        # Broca's area (language production, motor-related)

    # PREFRONTAL/EXECUTIVE
    'A13': 'Prefrontal_Assoc',       # Orbitofrontal cortex
    'A25': 'Limbic_Assoc',           # Subgenual ACC (emotion/mood)
    'A32': 'Limbic_Assoc',           # Dorsal ACC (cognitive control/emotion)

    # TEMPORAL/LANGUAGE/MEMORY (expected lower metabolism)
    'MTG': 'Language_Assoc',         # Middle temporal (semantic processing)
    'STG': 'Language_Assoc',         # Superior temporal (auditory/language)
    'ITG': 'Visual_Assoc',           # Inferior temporal (object recognition)
    'A38': 'Limbic_Assoc',           # Temporal pole (semantic memory/emotion)

    # INSULAR (interoception/emotion)
    'FI': 'Limbic_Assoc',            # Frontal insula (interoception)
    'Ig': 'Limbic_Assoc',            # Granular insula (interoception)
}

# Broader categories for main hypothesis
BROAD_CATEGORY = {
    'Primary_Sensory': 'SENSORIMOTOR',
    'Primary_Motor': 'SENSORIMOTOR',
    'Secondary_Sensory': 'SENSORIMOTOR',
    'Sensorimotor_Assoc': 'SENSORIMOTOR',
    'Motor_Assoc': 'SENSORIMOTOR',

    'Prefrontal_Assoc': 'ASSOCIATIVE',
    'Language_Assoc': 'ASSOCIATIVE',
    'Visual_Assoc': 'ASSOCIATIVE',
    'Limbic_Assoc': 'ASSOCIATIVE',
}

# Add classifications to data
subregion_data['Functional_Type'] = subregion_data['Region'].map(FUNCTIONAL_CLASSIFICATION)
subregion_data['Broad_Category'] = subregion_data['Functional_Type'].map(BROAD_CATEGORY)

# Analysis
sensorimotor = subregion_data[subregion_data['Broad_Category'] == 'SENSORIMOTOR']
associative = subregion_data[subregion_data['Broad_Category'] == 'ASSOCIATIVE']

print(f'\nClassification results:')
print(f'  SENSORIMOTOR regions: n={len(sensorimotor)}')
print(f'    Mean rCMRGlc: {sensorimotor["rCMRGlc_Mean"].mean():.2f} ± {sensorimotor["rCMRGlc_Mean"].std():.2f}')
print(f'    Mean E:I: {sensorimotor["EI_Ratio"].mean():.2f} ± {sensorimotor["EI_Ratio"].std():.2f}')

print(f'\n  ASSOCIATIVE regions: n={len(associative)}')
print(f'    Mean rCMRGlc: {associative["rCMRGlc_Mean"].mean():.2f} ± {associative["rCMRGlc_Mean"].std():.2f}')
print(f'    Mean E:I: {associative["EI_Ratio"].mean():.2f} ± {associative["EI_Ratio"].std():.2f}')

# Statistical tests
t_met, p_met = stats.ttest_ind(sensorimotor['rCMRGlc_Mean'], associative['rCMRGlc_Mean'])
t_ei, p_ei = stats.ttest_ind(sensorimotor['EI_Ratio'], associative['EI_Ratio'])

print(f'\nStatistical Tests:')
print(f'  Metabolism difference: t={t_met:.3f}, p={p_met:.4f} {"***SIGNIFICANT***" if p_met < 0.05 else ""}')
print(f'  E:I difference: t={t_ei:.3f}, p={p_ei:.4f} {"***SIGNIFICANT***" if p_ei < 0.05 else ""}')

# Create visualization
fig, axes = plt.subplots(2, 2, figsize=(16, 12))

# Plot 1: Metabolism comparison
ax = axes[0, 0]
data_met = [sensorimotor['rCMRGlc_Mean'], associative['rCMRGlc_Mean']]
bp = ax.boxplot(data_met, labels=['Sensorimotor', 'Associative'],
               patch_artist=True, showmeans=True)
bp['boxes'][0].set_facecolor('coral')
bp['boxes'][1].set_facecolor('lightblue')
ax.set_ylabel('rCMRGlc (μmol/100g/min)', fontsize=13, fontweight='bold')
ax.set_title(f'Metabolism by Functional Category\np = {p_met:.4f}', fontsize=14, fontweight='bold')
ax.grid(axis='y', alpha=0.3)
if p_met < 0.05:
    ax.text(0.5, 0.95, '*** SIGNIFICANT ***', transform=ax.transAxes,
           ha='center', fontsize=12, fontweight='bold', color='red')

# Plot 2: E:I comparison
ax = axes[0, 1]
data_ei = [sensorimotor['EI_Ratio'], associative['EI_Ratio']]
bp = ax.boxplot(data_ei, labels=['Sensorimotor', 'Associative'],
               patch_artist=True, showmeans=True)
bp['boxes'][0].set_facecolor('coral')
bp['boxes'][1].set_facecolor('lightblue')
ax.set_ylabel('E:I Ratio', fontsize=13, fontweight='bold')
ax.set_title(f'E:I Ratio by Functional Category\np = {p_ei:.4f}', fontsize=14, fontweight='bold')
ax.grid(axis='y', alpha=0.3)

# Plot 3: Scatter by functional type
ax = axes[1, 0]
func_types = subregion_data['Functional_Type'].unique()
colors_func = plt.cm.Set3(range(len(func_types)))
for ftype, color in zip(func_types, colors_func):
    data = subregion_data[subregion_data['Functional_Type'] == ftype]
    ax.scatter(data['EI_Ratio'], data['rCMRGlc_Mean'],
              s=150, alpha=0.7, color=color, edgecolors='black',
              linewidth=1.5, label=ftype.replace('_', ' '))
ax.set_xlabel('E:I Ratio', fontsize=13, fontweight='bold')
ax.set_ylabel('rCMRGlc (μmol/100g/min)', fontsize=13, fontweight='bold')
ax.set_title('Detailed Functional Classification', fontsize=14, fontweight='bold')
ax.legend(fontsize=8, loc='best', ncol=2)
ax.grid(True, alpha=0.3)

# Plot 4: Category breakdown
ax = axes[1, 1]
summary_data = []
for category in ['SENSORIMOTOR', 'ASSOCIATIVE']:
    cat_data = subregion_data[subregion_data['Broad_Category'] == category]
    for ftype in cat_data['Functional_Type'].unique():
        ftype_data = cat_data[cat_data['Functional_Type'] == ftype]
        summary_data.append({
            'Category': category,
            'Type': ftype,
            'n': len(ftype_data),
            'Mean_rCMRGlc': ftype_data['rCMRGlc_Mean'].mean(),
            'Mean_EI': ftype_data['EI_Ratio'].mean()
        })

summary_df = pd.DataFrame(summary_data)
x_pos = np.arange(len(summary_df))
colors_bar = ['coral' if cat == 'SENSORIMOTOR' else 'lightblue'
              for cat in summary_df['Category']]
ax.bar(x_pos, summary_df['Mean_rCMRGlc'], color=colors_bar,
       alpha=0.7, edgecolor='black', linewidth=1.5)
ax.set_xticks(x_pos)
ax.set_xticklabels([t.replace('_', '\n') for t in summary_df['Type']],
                   rotation=45, ha='right', fontsize=9)
ax.set_ylabel('Mean rCMRGlc', fontsize=12, fontweight='bold')
ax.set_title('Metabolism by Detailed Functional Type', fontsize=13, fontweight='bold')
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'analysis1_functional_classification.png', dpi=300, bbox_inches='tight')
print('\n  Saved: analysis1_functional_classification.png')
plt.close()

# Save classification data
subregion_data.to_csv(OUTPUT_DIR / 'subregions_with_functional_classification.csv', index=False)

print('\n' + '='*80)
print('ANALYSIS 2: CELL SUBTYPE COMPOSITION')
print('='*80)

print('\nExploring available cell subtype annotations...')
print(f'Annotation sets available: {annotations["cluster_annotation_term_set_name"].unique().tolist()}')

# Check for cell subtype markers in annotations
cell_type_sets = ['Subclass', 'Supertype', 'Class', 'neurotransmitter']
for anno_set in cell_type_sets:
    anno_subset = annotations[annotations['cluster_annotation_term_set_name'] == anno_set]
    if len(anno_subset) > 0:
        print(f'\n{anno_set} annotations ({len(anno_subset)}):')
        unique_types = anno_subset['name'].unique()
        print(f'  Total types: {len(unique_types)}')

        # Look for interneuron markers
        pv_types = [t for t in unique_types if 'PV' in str(t).upper() or 'PARV' in str(t).upper()]
        sst_types = [t for t in unique_types if 'SST' in str(t).upper() or 'SOMAT' in str(t).upper()]
        vip_types = [t for t in unique_types if 'VIP' in str(t).upper()]
        lamp5_types = [t for t in unique_types if 'LAMP5' in str(t).upper()]

        if pv_types:
            print(f'  PV-related: {len(pv_types)} types')
            for t in pv_types[:5]:
                print(f'    - {t}')
        if sst_types:
            print(f'  SST-related: {len(sst_types)} types')
            for t in sst_types[:5]:
                print(f'    - {t}')
        if vip_types:
            print(f'  VIP-related: {len(vip_types)} types')
            for t in vip_types[:5]:
                print(f'    - {t}')
        if lamp5_types:
            print(f'  LAMP5-related: {len(lamp5_types)} types')
            for t in lamp5_types[:5]:
                print(f'    - {t}')

# We'll analyze Subclass level for interneuron subtypes
print('\n\nAnalyzing Subclass annotations for interneuron composition...')

subclass_annot = annotations[annotations['cluster_annotation_term_set_name'] == 'Subclass'].copy()

# Define interneuron subtype categories
def classify_interneuron_subtype(subclass_name):
    if pd.isna(subclass_name):
        return 'Unknown'
    name_upper = str(subclass_name).upper()

    # Interneuron subtypes
    if 'PVALB' in name_upper or 'PV' in name_upper:
        return 'PV_Interneuron'
    elif 'SST' in name_upper or 'SOMATOSTATIN' in name_upper:
        return 'SST_Interneuron'
    elif 'VIP' in name_upper:
        return 'VIP_Interneuron'
    elif 'LAMP5' in name_upper:
        return 'LAMP5_Interneuron'
    elif 'GABA' in name_upper or 'INHIB' in name_upper:
        return 'Other_Interneuron'
    elif 'GLUT' in name_upper or 'EXCIT' in name_upper:
        return 'Excitatory'
    else:
        return 'Non_neuronal'

subclass_annot['Interneuron_Subtype'] = subclass_annot['name'].apply(classify_interneuron_subtype)

print('\nSubclass E/I subtype distribution:')
print(subclass_annot['Interneuron_Subtype'].value_counts())

# Build cluster to subtype mapping
print('\nMapping clusters to subtypes...')
cluster_to_subtype = {}
for _, row in subclass_annot.iterrows():
    members = membership[membership['cluster_annotation_term_label'] == row['label']]
    for clust_alias in members['cluster_alias'].values:
        cluster_to_subtype[clust_alias] = row['Interneuron_Subtype']

print(f'Mapped {len(cluster_to_subtype)} clusters to subtypes')

# Now need to load cell metadata to count subtypes per region
print('\nLoading cell metadata to count subtypes per region...')
print('(This may take 1-2 minutes for 672 MB file...)')

cell_df = pd.read_csv(DATA_DIR / 'WHB-10Xv3_cell_metadata.csv')
print(f'Loaded {len(cell_df):,} cells')

# Map cells to subtypes
cell_df['cell_subtype'] = cell_df['cluster_alias'].map(cluster_to_subtype).fillna('Unknown')

# Filter cortical cells
cell_df['region_clean'] = cell_df['region_of_interest_label'].str.replace('Human ', '')
cortical = cell_df[cell_df['anatomical_division_label'] == 'Cerebral cortex'].copy()
print(f'Cortical cells: {len(cortical):,}')

# Count subtypes per region
print('\nCalculating subtype composition per region...')
region_subtypes = []

for region in subregion_data['Region'].unique():
    region_cells = cortical[cortical['region_clean'] == region]

    if len(region_cells) == 0:
        continue

    total = len(region_cells)
    pv = (region_cells['cell_subtype'] == 'PV_Interneuron').sum()
    sst = (region_cells['cell_subtype'] == 'SST_Interneuron').sum()
    vip = (region_cells['cell_subtype'] == 'VIP_Interneuron').sum()
    lamp5 = (region_cells['cell_subtype'] == 'LAMP5_Interneuron').sum()
    other_inh = (region_cells['cell_subtype'] == 'Other_Interneuron').sum()
    excit = (region_cells['cell_subtype'] == 'Excitatory').sum()

    region_subtypes.append({
        'Region': region,
        'Total': total,
        'PV': pv,
        'SST': sst,
        'VIP': vip,
        'LAMP5': lamp5,
        'Other_Inh': other_inh,
        'Excitatory': excit,
        'PV_pct': pv/total*100 if total > 0 else 0,
        'SST_pct': sst/total*100 if total > 0 else 0,
        'VIP_pct': vip/total*100 if total > 0 else 0,
        'LAMP5_pct': lamp5/total*100 if total > 0 else 0,
        'Total_Inh': pv + sst + vip + lamp5 + other_inh,
        'Total_Inh_pct': (pv + sst + vip + lamp5 + other_inh)/total*100 if total > 0 else 0,
    })

subtype_df = pd.DataFrame(region_subtypes)

# Merge with existing data
subtype_merged = subregion_data.merge(subtype_df, on='Region', how='left')
subtype_merged.to_csv(OUTPUT_DIR / 'subregions_with_subtypes.csv', index=False)

print('\nSubtype composition summary:')
print(f'  Mean PV%: {subtype_df["PV_pct"].mean():.2f}%')
print(f'  Mean SST%: {subtype_df["SST_pct"].mean():.2f}%')
print(f'  Mean VIP%: {subtype_df["VIP_pct"].mean():.2f}%')
print(f'  Mean LAMP5%: {subtype_df["LAMP5_pct"].mean():.2f}%')

# Test correlations with metabolism
print('\nCorrelations with rCMRGlc:')
for subtype in ['PV_pct', 'SST_pct', 'VIP_pct', 'LAMP5_pct']:
    valid_data = subtype_merged[[subtype, 'rCMRGlc_Mean']].dropna()
    if len(valid_data) > 3:
        r, p = stats.pearsonr(valid_data[subtype], valid_data['rCMRGlc_Mean'])
        sig = "***" if p < 0.05 else ""
        print(f'  {subtype}: r={r:.3f}, p={p:.4f} {sig}')

# Visualize
fig, axes = plt.subplots(2, 3, figsize=(18, 12))

# Plot interneuron subtype percentages vs metabolism
for idx, (subtype, title) in enumerate([
    ('PV_pct', 'Parvalbumin'),
    ('SST_pct', 'Somatostatin'),
    ('VIP_pct', 'VIP'),
    ('LAMP5_pct', 'LAMP5')
]):
    ax = axes[idx//3, idx%3]
    valid_data = subtype_merged[[subtype, 'rCMRGlc_Mean', 'Region']].dropna()

    if len(valid_data) > 0:
        ax.scatter(valid_data[subtype], valid_data['rCMRGlc_Mean'],
                  s=150, alpha=0.7, edgecolors='black', linewidth=1.5)

        # Add labels
        for _, row in valid_data.iterrows():
            ax.annotate(row['Region'], (row[subtype], row['rCMRGlc_Mean']),
                       fontsize=7, alpha=0.6)

        if len(valid_data) > 3:
            r, p = stats.pearsonr(valid_data[subtype], valid_data['rCMRGlc_Mean'])
            ax.set_title(f'{title}+ Interneurons vs rCMRGlc\nr={r:.3f}, p={p:.3f}',
                        fontsize=12, fontweight='bold')
        else:
            ax.set_title(f'{title}+ Interneurons vs rCMRGlc',
                        fontsize=12, fontweight='bold')

        ax.set_xlabel(f'{title}+ % of cells', fontsize=11, fontweight='bold')
        ax.set_ylabel('rCMRGlc (μmol/100g/min)', fontsize=11, fontweight='bold')
        ax.grid(True, alpha=0.3)

# Stacked bar chart of subtype composition
ax = axes[1, 2]
valid_regions = subtype_merged[['Region', 'PV_pct', 'SST_pct', 'VIP_pct', 'LAMP5_pct', 'rCMRGlc_Mean']].dropna()
valid_regions = valid_regions.sort_values('rCMRGlc_Mean', ascending=False)

x_pos = np.arange(len(valid_regions))
pv_vals = valid_regions['PV_pct'].values
sst_vals = valid_regions['SST_pct'].values
vip_vals = valid_regions['VIP_pct'].values
lamp5_vals = valid_regions['LAMP5_pct'].values

ax.bar(x_pos, pv_vals, label='PV', color='#E63946', alpha=0.8)
ax.bar(x_pos, sst_vals, bottom=pv_vals, label='SST', color='#F1A208', alpha=0.8)
ax.bar(x_pos, vip_vals, bottom=pv_vals+sst_vals, label='VIP', color='#2A9D8F', alpha=0.8)
ax.bar(x_pos, lamp5_vals, bottom=pv_vals+sst_vals+vip_vals, label='LAMP5', color='#457B9D', alpha=0.8)

ax.set_xticks(x_pos)
ax.set_xticklabels(valid_regions['Region'].values, rotation=90, fontsize=8)
ax.set_ylabel('% of cells', fontsize=11, fontweight='bold')
ax.set_title('Interneuron Subtype Composition\n(sorted by rCMRGlc)', fontsize=12, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(axis='y', alpha=0.3)

# Hide empty subplot
axes[1, 1].set_visible(False)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'analysis2_interneuron_subtypes.png', dpi=300, bbox_inches='tight')
print('\n  Saved: analysis2_interneuron_subtypes.png')
plt.close()

print('\n' + '='*80)
print('ALL ANALYSES COMPLETE!')
print('='*80)
print(f'\nOutput directory: {OUTPUT_DIR}')
print('\nFiles created:')
print('  1. analysis1_functional_classification.png')
print('  2. subregions_with_functional_classification.csv')
print('  3. analysis2_interneuron_subtypes.png')
print('  4. subregions_with_subtypes.csv')
print('\nNote: Gene expression analysis will be in separate script due to data size')
print('='*80)
