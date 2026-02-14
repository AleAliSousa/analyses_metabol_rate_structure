"""
COMPLETE E:I RATIO ANALYSIS - WHB Cortical Data
================================================

Full pipeline:
1. Load cell metadata + taxonomy
2. Classify neurons as E/I
3. Filter cortical regions
4. Calculate E:I ratios by region
5. Map to Heiss lobes
6. Integrate with rCMRGlc
7. Calculate correlations

Estimated runtime: 10-15 minutes
"""

import pandas as pd
import numpy as np
from pathlib import Path
from scipy import stats

print('='*80)
print('COMPLETE E:I vs rCMRGlc ANALYSIS')
print('='*80)

# Paths
PROJECT_DIR = Path(r"C:\Users\michaelproulx\MetabolicBrain")
EI_DIR = PROJECT_DIR / "cellxgene_integration" / "ei_analysis"
DATA_DIR = PROJECT_DIR / "cellxgene_integration" / "data" / "cortical"
OUTPUT_DIR = EI_DIR / "data" / "cortical"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Heiss rCMRGlc data
HEISS_DATA = {
    'Cerebral_cortex': (33.5, 2.98),
    'Frontal_lobe': (35.3, 3.51),
    'Parietal_lobe': (35.8, 3.23),
    'Temporal_lobe': (30.5, 2.54),
    'Occipital_lobe': (35.8, 3.12),
    'Insular_lobe': (30.3, 2.55),
    'Primary_visual_cortex': (40.1, np.nan),
}

# Region to lobe mapping
REGION_TO_LOBE = {
    'M1C': 'Frontal_lobe',  # Primary motor
    'A44-A45': 'Frontal_lobe',  # Broca's area
    'A13': 'Frontal_lobe',  # Orbitofrontal
    'A25': 'Frontal_lobe',  # Subgenual
    'A32': 'Frontal_lobe',  # Dorsal ACC
    'A46': 'Frontal_lobe',  # Dorsolateral prefrontal

    'S1C': 'Parietal_lobe',  # Primary somatosensory
    'A5-A7': 'Parietal_lobe',  # Superior parietal
    'A40': 'Parietal_lobe',  # Supramarginal gyrus

    'MTG': 'Temporal_lobe',  # Middle temporal gyrus
    'STG': 'Temporal_lobe',  # Superior temporal gyrus
    'ITG': 'Temporal_lobe',  # Inferior temporal gyrus
    'A38': 'Temporal_lobe',  # Temporal pole
    'A1C': 'Temporal_lobe',  # Primary auditory

    'V1C': 'Primary_visual_cortex',  # Primary visual
    'V2': 'Occipital_lobe',  # Visual V2
    'A19': 'Occipital_lobe',  # Associative visual

    'FI': 'Insular_lobe',  # Frontal insula
    'SI': 'Insular_lobe',  # Superior insula
}

print('\nStep 1: Loading taxonomy files...')
annotations = pd.read_csv(EI_DIR / 'WHB_cluster_annotation_term.csv')
membership = pd.read_csv(EI_DIR / 'WHB_cluster_to_cluster_annotation_membership.csv')
clusters = pd.read_csv(EI_DIR / 'WHB_cluster.csv')

print(f'  Annotations: {len(annotations):,}')
print(f'  Membership: {len(membership):,}')
print(f'  Clusters: {len(clusters):,}')

print('\nStep 2: Identifying E/I cell types from Class annotations...')
class_annot = annotations[annotations['cluster_annotation_term_set_name'] == 'Class']
print(f'  Class annotations: {len(class_annot)}')

# E/I classification function
def classify_ei(class_name):
    if pd.isna(class_name):
        return 'Unknown'
    name_lower = str(class_name).lower()

    # Excitatory
    if any(kw in name_lower for kw in ['glut', 'glutamatergic', 'excit']):
        return 'Excitatory'
    # Inhibitory
    if any(kw in name_lower for kw in ['gaba', 'inhib']):
        return 'Inhibitory'
    return 'Non-neuronal'

class_annot['ei_type'] = class_annot['name'].apply(classify_ei)

print('\nE/I classification:')
for ei, count in class_annot['ei_type'].value_counts().items():
    print(f'  {ei}: {count}')

print('\nExcitatory classes:')
for name in class_annot[class_annot['ei_type'] == 'Excitatory']['name'].values:
    print(f'  - {name}')

print('\nInhibitory classes:')
for name in class_annot[class_annot['ei_type'] == 'Inhibitory']['name'].values:
    print(f'  - {name}')

print('\nStep 3: Mapping clusters to E/I types...')
# Build cluster -> class -> E/I mapping
cluster_to_ei = {}

for _, class_row in class_annot.iterrows():
    class_label = class_row['label']
    ei_type = class_row['ei_type']

    # Find all clusters under this class
    class_members = membership[membership['cluster_annotation_term_label'] == class_label]

    for cluster_label in class_members['cluster_label'].values:
        cluster_to_ei[cluster_label] = ei_type

print(f'  Mapped {len(cluster_to_ei)} clusters to E/I types')

# Get cluster alias mapping
cluster_label_to_alias = dict(zip(clusters['label'], clusters['alias']))

print('\nStep 4: Loading cell metadata (this will take a few minutes)...')
cell_df = pd.read_csv(DATA_DIR / 'WHB-10Xv3_cell_metadata.csv')
print(f'  Loaded {len(cell_df):,} cells')

print('\nStep 5: Joining cell data with cluster labels...')
# Create cluster_label column by looking up alias
alias_to_label = {v: k for k, v in cluster_label_to_alias.items()}
cell_df['cluster_label'] = cell_df['cluster_alias'].map(alias_to_label)

# Map to E/I
cell_df['ei_type'] = cell_df['cluster_label'].map(cluster_to_ei)
cell_df['ei_type'] = cell_df['ei_type'].fillna('Unknown')

print('\nCell E/I distribution:')
ei_dist = cell_df['ei_type'].value_counts()
for ei, count in ei_dist.items():
    print(f'  {ei}: {count:,} ({count/len(cell_df)*100:.1f}%)')

print('\nStep 6: Filtering for cortical cells...')
cell_df['region_clean'] = cell_df['region_of_interest_label'].str.replace('Human ', '')
cortical = cell_df[cell_df['anatomical_division_label'] == 'Cerebral cortex']
print(f'  Cortical cells: {len(cortical):,}')

print('\nStep 7: Calculating E:I ratios by cortical region...')
ei_by_region = []

for region in cortical['region_clean'].unique():
    region_cells = cortical[cortical['region_clean'] == region]

    e_count = (region_cells['ei_type'] == 'Excitatory').sum()
    i_count = (region_cells['ei_type'] == 'Inhibitory').sum()
    total = len(region_cells)

    ei_ratio = e_count / i_count if i_count > 0 else np.nan

    ei_by_region.append({
        'Region': region,
        'Excitatory': e_count,
        'Inhibitory': i_count,
        'Total': total,
        'EI_Ratio': ei_ratio,
        'E_percent': (e_count/total*100) if total > 0 else 0,
    })

ei_df = pd.DataFrame(ei_by_region)
ei_df = ei_df.sort_values('Total', ascending=False)

print('\nE:I ratios calculated for', len(ei_df), 'cortical regions')
print('\nTop regions:')
for _, row in ei_df.head(10).iterrows():
    print(f"  {row['Region']:15s} E:{row['Excitatory']:6,}  I:{row['Inhibitory']:6,}  Ratio:{row['EI_Ratio']:.2f}" if not np.isnan(row['EI_Ratio']) else f"  {row['Region']:15s} E:{row['Excitatory']:6,}  I:{row['Inhibitory']:6,}  Ratio:N/A")

# Save fine-grained results
ei_df.to_csv(OUTPUT_DIR / 'ei_ratios_whb_cortical_fine_grained.csv', index=False)

print('\nStep 8: Mapping regions to Heiss lobes...')
ei_df['Lobe'] = ei_df['Region'].map(REGION_TO_LOBE)
mapped = ei_df[ei_df['Lobe'].notna()]
print(f'  Mapped {len(mapped)}/{len(ei_df)} regions to lobes')

print('\nStep 9: Aggregating E:I ratios by lobe...')
lobe_ei = []

for lobe in mapped['Lobe'].unique():
    lobe_cells = mapped[mapped['Lobe'] == lobe]

    total_e = lobe_cells['Excitatory'].sum()
    total_i = lobe_cells['Inhibitory'].sum()
    total_cells = lobe_cells['Total'].sum()
    n_regions = len(lobe_cells)

    # Proper aggregation: sum E, sum I, then divide
    ei_ratio_agg = total_e / total_i if total_i > 0 else np.nan

    lobe_ei.append({
        'Lobe': lobe,
        'Excitatory': total_e,
        'Inhibitory': total_i,
        'Total': total_cells,
        'EI_Ratio_aggregated': ei_ratio_agg,
        'n_constituent_regions': n_regions,
        'E_percent': (total_e/total_cells*100) if total_cells > 0 else 0,
    })

lobe_df = pd.DataFrame(lobe_ei)
lobe_df = lobe_df.sort_values('Total', ascending=False)

print('\nAggregated E:I ratios by lobe:')
for _, row in lobe_df.iterrows():
    print(f"\n{row['Lobe']}:")
    print(f"  E: {row['Excitatory']:,}, I: {row['Inhibitory']:,}")
    print(f"  E:I Ratio: {row['EI_Ratio_aggregated']:.2f}" if not np.isnan(row['EI_Ratio_aggregated']) else "  E:I Ratio: N/A")
    print(f"  Regions: {row['n_constituent_regions']}")

# Save lobe-level results
lobe_df.to_csv(OUTPUT_DIR / 'ei_ratios_aggregated_lobes.csv', index=False)

print('\nStep 10: Integrating with Heiss rCMRGlc data...')
lobe_df['rCMRGlc_Mean'] = lobe_df['Lobe'].map({k: v[0] for k, v in HEISS_DATA.items()})
lobe_df['rCMRGlc_SD'] = lobe_df['Lobe'].map({k: v[1] for k, v in HEISS_DATA.items()})

integrated = lobe_df[lobe_df['rCMRGlc_Mean'].notna()]
print(f'  Integrated {len(integrated)} lobes with rCMRGlc data')

print('\nIntegrated data:')
for _, row in integrated.iterrows():
    print(f"\n{row['Lobe']}:")
    print(f"  E:I Ratio: {row['EI_Ratio_aggregated']:.2f}" if not np.isnan(row['EI_Ratio_aggregated']) else "  E:I Ratio: N/A")
    print(f"  rCMRGlc: {row['rCMRGlc_Mean']:.1f} umol/100g/min")

# Save integrated results
integrated.to_csv(OUTPUT_DIR / 'ei_metabolism_integrated.csv', index=False)

print('\nStep 11: Statistical analysis...')
if len(integrated) >= 3:
    # Filter out rows with NaN E:I ratio
    valid = integrated[integrated['EI_Ratio_aggregated'].notna()]

    if len(valid) >= 3:
        ei_values = valid['EI_Ratio_aggregated'].values
        rcmrglc_values = valid['rCMRGlc_Mean'].values

        # Pearson correlation
        pearson_r, pearson_p = stats.pearsonr(ei_values, rcmrglc_values)

        # Spearman correlation
        spearman_r, spearman_p = stats.spearmanr(ei_values, rcmrglc_values)

        # Linear regression
        slope, intercept, r_value, p_value, std_err = stats.linregress(ei_values, rcmrglc_values)

        print('\n' + '='*80)
        print('CORRELATION RESULTS')
        print('='*80)

        print(f'\nPearson correlation:')
        print(f'  r = {pearson_r:.3f}')
        print(f'  p-value = {pearson_p:.4f}')
        print(f'  Significant: {"YES" if pearson_p < 0.05 else "NO"}')

        print(f'\nSpearman correlation:')
        print(f'  rho = {spearman_r:.3f}')
        print(f'  p-value = {spearman_p:.4f}')
        print(f'  Significant: {"YES" if spearman_p < 0.05 else "NO"}')

        print(f'\nLinear regression:')
        print(f'  rCMRGlc = {slope:.2f} * EI_Ratio + {intercept:.2f}')
        print(f'  R-squared = {r_value**2:.3f}')
        print(f'  p-value = {p_value:.4f}')

        # Save statistics
        stats_file = OUTPUT_DIR / 'ei_metabolism_statistics.txt'
        with open(stats_file, 'w') as f:
            f.write('E:I Ratio vs rCMRGlc - Statistical Analysis\\n')
            f.write('='*80 + '\\n\\n')
            f.write(f'Number of regions: {len(valid)}\\n\\n')
            f.write(f'Pearson r = {pearson_r:.3f}, p = {pearson_p:.4f}\\n')
            f.write(f'Spearman rho = {spearman_r:.3f}, p = {spearman_p:.4f}\\n')
            f.write(f'Linear regression: rCMRGlc = {slope:.2f} * EI_Ratio + {intercept:.2f}\\n')
            f.write(f'R-squared = {r_value**2:.3f}, p = {p_value:.4f}\\n')

        print(f'\nStatistics saved to: {stats_file}')
    else:
        print(f'  Insufficient valid data ({len(valid)} regions)')
else:
    print(f'  Insufficient data ({len(integrated)} regions)')

print('\n' + '='*80)
print('ANALYSIS COMPLETE!')
print('='*80)

print('\nOutput files:')
print(f'  1. Fine-grained E:I ratios: ei_ratios_whb_cortical_fine_grained.csv')
print(f'  2. Lobe-aggregated E:I ratios: ei_ratios_aggregated_lobes.csv')
print(f'  3. Integrated with metabolism: ei_metabolism_integrated.csv')
print(f'  4. Statistical results: ei_metabolism_statistics.txt')

print('\n' + '='*80)
