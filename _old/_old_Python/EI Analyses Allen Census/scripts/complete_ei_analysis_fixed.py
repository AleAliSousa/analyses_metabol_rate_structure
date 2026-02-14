"""
COMPLETE E:I RATIO ANALYSIS - WHB Cortical Data (FIXED)
========================================================

Uses neurotransmitter-level annotations for E/I classification.
Runs full pipeline from cell metadata to final correlation results.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from scipy import stats
import time

start_time = time.time()

print('='*80)
print('COMPLETE E:I vs rCMRGlc ANALYSIS - WHB Dataset')
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
    'M1C': 'Frontal_lobe',
    'A44-A45': 'Frontal_lobe',
    'A13': 'Frontal_lobe',
    'A25': 'Frontal_lobe',
    'A32': 'Frontal_lobe',
    'A46': 'Frontal_lobe',
    'S1C': 'Parietal_lobe',
    'A5-A7': 'Parietal_lobe',
    'A40': 'Parietal_lobe',
    'A43': 'Parietal_lobe',
    'MTG': 'Temporal_lobe',
    'STG': 'Temporal_lobe',
    'ITG': 'Temporal_lobe',
    'A38': 'Temporal_lobe',
    'A1C': 'Temporal_lobe',
    'V1C': 'Primary_visual_cortex',
    'V2': 'Occipital_lobe',
    'A19': 'Occipital_lobe',
    'FI': 'Insular_lobe',
    'SI': 'Insular_lobe',
    'Ig': 'Insular_lobe',
}

print('\nStep 1: Loading taxonomy files...')
annotations = pd.read_csv(EI_DIR / 'WHB_cluster_annotation_term.csv')
membership = pd.read_csv(EI_DIR / 'WHB_cluster_to_cluster_annotation_membership.csv')
clusters = pd.read_csv(EI_DIR / 'WHB_cluster.csv')

print(f'  Annotations: {len(annotations):,}')
print(f'  Membership: {len(membership):,}')
print(f'  Clusters: {len(clusters):,}')

print('\nStep 2: Identifying E/I from neurotransmitter annotations...')
nt_annot = annotations[annotations['cluster_annotation_term_set_name'] == 'neurotransmitter'].copy()
print(f'  Neurotransmitter annotations: {len(nt_annot)}')

# E/I classification based on neurotransmitter
def classify_ei_nt(nt_name):
    if pd.isna(nt_name):
        return 'Unknown'
    name_str = str(nt_name).upper()

    # Excitatory: VGLUT (vesicular glutamate transporters)
    if 'VGLUT' in name_str:
        return 'Excitatory'
    # Inhibitory: GABA, GLY (glycine)
    if 'GABA' in name_str or name_str == 'GLY':
        return 'Inhibitory'
    # Other neurotransmitters
    return 'Other'

nt_annot['ei_type'] = nt_annot['name'].apply(classify_ei_nt)

print('\nNeurotransmitter E/I classification:')
for ei, count in nt_annot['ei_type'].value_counts().items():
    print(f'  {ei}: {count}')

print('\nExcitatory neurotransmitters:')
for name in nt_annot[nt_annot['ei_type'] == 'Excitatory']['name'].values:
    print(f'  - {name}')

print('\nInhibitory neurotransmitters:')
for name in nt_annot[nt_annot['ei_type'] == 'Inhibitory']['name'].values:
    print(f'  - {name}')

print('\nStep 3: Building cluster -> E/I mapping...')
cluster_to_ei = {}

for _, nt_row in nt_annot.iterrows():
    nt_label = nt_row['label']
    ei_type = nt_row['ei_type']

    # Find all clusters with this neurotransmitter
    nt_members = membership[membership['cluster_annotation_term_label'] == nt_label]

    for cluster_label in nt_members['cluster_label'].values:
        cluster_to_ei[cluster_label] = ei_type

print(f'  Mapped {len(cluster_to_ei)} clusters to E/I types')

# Get cluster alias mapping
cluster_label_to_alias = dict(zip(clusters['label'], clusters['cluster_alias']))
cluster_alias_to_label = {v: k for k, v in cluster_label_to_alias.items()}

print(f'  Cluster aliases: {len(cluster_alias_to_label)}')

print('\nStep 4: Loading cell metadata...')
print('  (This is a large file - will take 2-3 minutes)')
cell_df = pd.read_csv(DATA_DIR / 'WHB-10Xv3_cell_metadata.csv')
print(f'  Loaded {len(cell_df):,} cells')

print('\nStep 5: Mapping cells to E/I types...')
# Map cluster alias -> cluster label -> E/I type
cell_df['cluster_label'] = cell_df['cluster_alias'].map(cluster_alias_to_label)
cell_df['ei_type'] = cell_df['cluster_label'].map(cluster_to_ei)
cell_df['ei_type'] = cell_df['ei_type'].fillna('Unknown')

print('\nCell E/I distribution:')
ei_dist = cell_df['ei_type'].value_counts()
total_cells = len(cell_df)
for ei, count in ei_dist.items():
    print(f'  {ei}: {count:,} ({count/total_cells*100:.1f}%)')

print('\nStep 6: Filtering for cortical cells...')
cell_df['region_clean'] = cell_df['region_of_interest_label'].str.replace('Human ', '')
cortical = cell_df[cell_df['anatomical_division_label'] == 'Cerebral cortex'].copy()
print(f'  Cortical cells: {len(cortical):,} ({len(cortical)/total_cells*100:.1f}%)')

cortical_ei = cortical['ei_type'].value_counts()
print('\nCortical cell E/I distribution:')
for ei, count in cortical_ei.items():
    print(f'  {ei}: {count:,} ({count/len(cortical)*100:.1f}%)')

print('\nStep 7: Calculating E:I ratios by cortical region...')
ei_by_region = []

for region in sorted(cortical['region_clean'].unique()):
    region_cells = cortical[cortical['region_clean'] == region]

    e_count = (region_cells['ei_type'] == 'Excitatory').sum()
    i_count = (region_cells['ei_type'] == 'Inhibitory').sum()
    other_count = (region_cells['ei_type'] == 'Other').sum()
    unknown_count = (region_cells['ei_type'] == 'Unknown').sum()
    total = len(region_cells)

    ei_ratio = e_count / i_count if i_count > 0 else np.nan

    ei_by_region.append({
        'Region': region,
        'Excitatory': e_count,
        'Inhibitory': i_count,
        'Other': other_count,
        'Unknown': unknown_count,
        'Total': total,
        'EI_Ratio': ei_ratio,
        'E_percent': (e_count/total*100) if total > 0 else 0,
        'I_percent': (i_count/total*100) if total > 0 else 0,
    })

ei_df = pd.DataFrame(ei_by_region)
ei_df = ei_df.sort_values('Total', ascending=False)

print(f'\nE:I ratios calculated for {len(ei_df)} cortical regions')
print('\nTop 20 cortical regions:')
print('-'*100)
for _, row in ei_df.head(20).iterrows():
    ratio_str = f"{row['EI_Ratio']:.2f}" if not np.isnan(row['EI_Ratio']) else "N/A"
    print(f"{row['Region']:15s} Total:{row['Total']:7,}  E:{row['Excitatory']:6,}  I:{row['Inhibitory']:6,}  E:I={ratio_str:6s}  E%:{row['E_percent']:5.1f}")

# Save fine-grained results
fine_output = OUTPUT_DIR / 'ei_ratios_whb_cortical_fine_grained.csv'
ei_df.to_csv(fine_output, index=False)
print(f'\nFine-grained results saved: {fine_output}')

print('\nStep 8: Mapping cortical regions to Heiss lobes...')
ei_df['Lobe'] = ei_df['Region'].map(REGION_TO_LOBE)
mapped = ei_df[ei_df['Lobe'].notna()]
unmapped = ei_df[ei_df['Lobe'].isna()]

print(f'  Mapped: {len(mapped)} regions')
print(f'  Unmapped: {len(unmapped)} regions')

if len(unmapped) > 0 and len(unmapped) <= 10:
    print('\n  Unmapped regions:')
    for region in unmapped['Region'].values:
        print(f'    - {region}')

print('\nStep 9: Aggregating E:I ratios by lobe...')
lobe_ei = []

for lobe in sorted(mapped['Lobe'].unique()):
    lobe_cells = mapped[mapped['Lobe'] == lobe]

    # Sum across all regions in this lobe
    total_e = lobe_cells['Excitatory'].sum()
    total_i = lobe_cells['Inhibitory'].sum()
    total_cells = lobe_cells['Total'].sum()
    n_regions = len(lobe_cells)

    # PROPER aggregation: sum E, sum I, THEN divide
    ei_ratio_agg = total_e / total_i if total_i > 0 else np.nan

    # Standard deviation of E:I ratios across constituent regions
    valid_ratios = lobe_cells['EI_Ratio'].dropna()
    ei_ratio_sd = valid_ratios.std() if len(valid_ratios) > 1 else 0

    lobe_ei.append({
        'Lobe': lobe,
        'Excitatory': total_e,
        'Inhibitory': total_i,
        'Total': total_cells,
        'EI_Ratio_aggregated': ei_ratio_agg,
        'EI_Ratio_SD': ei_ratio_sd,
        'n_constituent_regions': n_regions,
        'E_percent': (total_e/total_cells*100) if total_cells > 0 else 0,
        'I_percent': (total_i/total_cells*100) if total_cells > 0 else 0,
    })

lobe_df = pd.DataFrame(lobe_ei)
lobe_df = lobe_df.sort_values('Total', ascending=False)

print('\nAggregated E:I ratios by lobe:')
print('='*100)
for _, row in lobe_df.iterrows():
    ratio_str = f"{row['EI_Ratio_aggregated']:.2f}" if not np.isnan(row['EI_Ratio_aggregated']) else "N/A"
    print(f"\n{row['Lobe']}:")
    print(f"  Total cells: {row['Total']:,}")
    print(f"  Excitatory: {row['Excitatory']:,} ({row['E_percent']:.1f}%)")
    print(f"  Inhibitory: {row['Inhibitory']:,} ({row['I_percent']:.1f}%)")
    print(f"  E:I Ratio: {ratio_str}")
    print(f"  Constituent regions: {row['n_constituent_regions']}")

# Save lobe-level results
lobe_output = OUTPUT_DIR / 'ei_ratios_aggregated_lobes.csv'
lobe_df.to_csv(lobe_output, index=False)
print(f'\nLobe-aggregated results saved: {lobe_output}')

print('\nStep 10: Integrating with Heiss rCMRGlc data...')
lobe_df['rCMRGlc_Mean'] = lobe_df['Lobe'].map({k: v[0] for k, v in HEISS_DATA.items()})
lobe_df['rCMRGlc_SD'] = lobe_df['Lobe'].map({k: v[1] for k, v in HEISS_DATA.items()})

integrated = lobe_df[lobe_df['rCMRGlc_Mean'].notna()].copy()
print(f'  Integrated {len(integrated)} lobes with rCMRGlc data')

print('\n' + '='*100)
print('INTEGRATED E:I RATIOS WITH METABOLIC DATA')
print('='*100)

for _, row in integrated.iterrows():
    ratio_str = f"{row['EI_Ratio_aggregated']:.2f}" if not np.isnan(row['EI_Ratio_aggregated']) else "N/A"
    print(f"\n{row['Lobe']}:")
    print(f"  E:I Ratio: {ratio_str}")
    print(f"  rCMRGlc: {row['rCMRGlc_Mean']:.1f} +/- {row['rCMRGlc_SD']:.1f} umol/100g/min" if not np.isnan(row['rCMRGlc_SD']) else f"  rCMRGlc: {row['rCMRGlc_Mean']:.1f} umol/100g/min")

# Save integrated results
integrated_output = OUTPUT_DIR / 'ei_metabolism_integrated.csv'
integrated.to_csv(integrated_output, index=False)
print(f'\nIntegrated results saved: {integrated_output}')

print('\nStep 11: Statistical analysis - E:I vs rCMRGlc correlation...')

valid = integrated[integrated['EI_Ratio_aggregated'].notna()].copy()

if len(valid) >= 3:
    ei_values = valid['EI_Ratio_aggregated'].values
    rcmrglc_values = valid['rCMRGlc_Mean'].values

    # Pearson correlation
    pearson_r, pearson_p = stats.pearsonr(ei_values, rcmrglc_values)

    # Spearman correlation
    spearman_r, spearman_p = stats.spearmanr(ei_values, rcmrglc_values)

    # Linear regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(ei_values, rcmrglc_values)

    print('\n' + '='*100)
    print('HYPOTHESIS TEST: E:I RATIO vs rCMRGlc CORRELATION')
    print('='*100)

    print(f'\nSample size: {len(valid)} lobes')
    print(f'Lobes analyzed: {", ".join(valid["Lobe"].values)}')

    print(f'\nPearson correlation:')
    print(f'  r = {pearson_r:.4f}')
    print(f'  p-value = {pearson_p:.4f}')
    print(f'  Significant at p<0.05: {"YES" if pearson_p < 0.05 else "NO"}')
    print(f'  Interpretation: {"Positive" if pearson_r > 0 else "Negative"} correlation')

    print(f'\nSpearman rank correlation:')
    print(f'  rho = {spearman_r:.4f}')
    print(f'  p-value = {spearman_p:.4f}')
    print(f'  Significant at p<0.05: {"YES" if spearman_p < 0.05 else "NO"}')

    print(f'\nLinear regression:')
    print(f'  Model: rCMRGlc = {slope:.3f} * EI_Ratio + {intercept:.3f}')
    print(f'  R-squared = {r_value**2:.4f}')
    print(f'  p-value = {p_value:.4f}')
    print(f'  Standard error = {std_err:.3f}')

    print('\n' + '='*100)
    print('BIOLOGICAL INTERPRETATION')
    print('='*100)

    if pearson_p < 0.05:
        if pearson_r > 0:
            print('\nRESULT: Significant POSITIVE correlation detected!')
            print('Higher E:I ratio is associated with HIGHER regional glucose metabolism.')
            print('This suggests that excitatory neurotransmission may drive metabolic demand.')
        else:
            print('\nRESULT: Significant NEGATIVE correlation detected!')
            print('Higher E:I ratio is associated with LOWER regional glucose metabolism.')
            print('This is contrary to typical expectations.')
    else:
        print('\nRESULT: No significant linear correlation detected.')
        print('E:I ratio and rCMRGlc may not be linearly related, or sample size is insufficient.')

    print(f'\nNote: Literature reports cortical E:I ratios typically 3:1 to 5:1')
    print(f'Observed range in this data: {ei_values.min():.2f} to {ei_values.max():.2f}')

    # Save statistics
    stats_file = OUTPUT_DIR / 'ei_metabolism_statistics.txt'
    with open(stats_file, 'w') as f:
        f.write('E:I Ratio vs rCMRGlc - Statistical Analysis\n')
        f.write('='*100 + '\n\n')
        f.write('Dataset: WHB-10Xv3 (Whole Human Brain)\n')
        f.write(f'Total cells analyzed: {len(cell_df):,}\n')
        f.write(f'Cortical cells: {len(cortical):,}\n')
        f.write(f'Number of lobes: {len(valid)}\n')
        f.write(f'Lobes: {", ".join(valid["Lobe"].values)}\n\n')

        f.write('Correlation Analysis:\n')
        f.write('-'*100 + '\n')
        f.write(f'Pearson r = {pearson_r:.4f}, p = {pearson_p:.4f}\n')
        f.write(f'Spearman rho = {spearman_r:.4f}, p = {spearman_p:.4f}\n\n')

        f.write('Linear Regression:\n')
        f.write('-'*100 + '\n')
        f.write(f'rCMRGlc = {slope:.3f} * EI_Ratio + {intercept:.3f}\n')
        f.write(f'R-squared = {r_value**2:.4f}\n')
        f.write(f'p-value = {p_value:.4f}\n')
        f.write(f'Standard error = {std_err:.3f}\n\n')

        f.write('Data by Lobe:\n')
        f.write('-'*100 + '\n')
        for _, row in valid.iterrows():
            f.write(f"{row['Lobe']:30s} E:I={row['EI_Ratio_aggregated']:.2f}  rCMRGlc={row['rCMRGlc_Mean']:.1f}\n")

    print(f'\nStatistics saved: {stats_file}')
else:
    print(f'\nInsufficient data for correlation analysis ({len(valid)} lobes with valid E:I ratios)')

elapsed = time.time() - start_time
print('\n' + '='*100)
print('ANALYSIS COMPLETE!')
print('='*100)

print(f'\nTotal processing time: {elapsed/60:.1f} minutes')

print('\nOutput files created:')
print(f'  1. {fine_output.name}')
print(f'  2. {lobe_output.name}')
print(f'  3. {integrated_output.name}')
if len(valid) >= 3:
    print(f'  4. {stats_file.name}')

print('\nAll files saved in:')
print(f'  {OUTPUT_DIR}')

print('\n' + '='*100)
print('READY FOR REVIEW')
print('='*100)
