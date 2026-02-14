import pandas as pd
import numpy as np
from pathlib import Path
from scipy import stats
import time

start_time = time.time()
print('='*80)
print('E:I vs rCMRGlc ANALYSIS - WHB Cortical Data')
print('Starting:', time.strftime('%Y-%m-%d %H:%M:%S'))
print('='*80)

PROJECT_DIR = Path(r"C:\Users\michaelproulx\MetabolicBrain")
EI_DIR = PROJECT_DIR / "cellxgene_integration" / "ei_analysis"
DATA_DIR = PROJECT_DIR / "cellxgene_integration" / "data" / "cortical"
OUTPUT_DIR = EI_DIR / "data" / "cortical"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

HEISS_DATA = {
    'Frontal_lobe': (35.3, 3.51),
    'Parietal_lobe': (35.8, 3.23),
    'Temporal_lobe': (30.5, 2.54),
    'Occipital_lobe': (35.8, 3.12),  # Includes V1C, V2, A19
    'Insular_lobe': (30.3, 2.55),
}

REGION_TO_LOBE = {
    'M1C': 'Frontal_lobe', 'A44-A45': 'Frontal_lobe', 'A13': 'Frontal_lobe',
    'A25': 'Frontal_lobe', 'A32': 'Frontal_lobe', 'A46': 'Frontal_lobe',
    'S1C': 'Parietal_lobe', 'A5-A7': 'Parietal_lobe', 'A40': 'Parietal_lobe', 'A43': 'Parietal_lobe',
    'MTG': 'Temporal_lobe', 'STG': 'Temporal_lobe', 'ITG': 'Temporal_lobe',
    'A38': 'Temporal_lobe', 'A1C': 'Temporal_lobe',
    'V1C': 'Occipital_lobe', 'V2': 'Occipital_lobe', 'A19': 'Occipital_lobe',  # V1C is part of occipital
    'FI': 'Insular_lobe', 'SI': 'Insular_lobe', 'Ig': 'Insular_lobe',
}

print('\n[1/11] Loading taxonomy...')
annotations = pd.read_csv(EI_DIR / 'WHB_cluster_annotation_term.csv')
membership = pd.read_csv(EI_DIR / 'WHB_cluster_to_cluster_annotation_membership.csv')
clusters = pd.read_csv(EI_DIR / 'WHB_cluster.csv')
print(f'  Loaded: {len(annotations):,} annotations, {len(membership):,} memberships, {len(clusters):,} clusters')

print('\n[2/11] Classifying E/I from neurotransmitters...')
nt_annot = annotations[annotations['cluster_annotation_term_set_name'] == 'neurotransmitter'].copy()

def classify_ei(nt):
    if pd.isna(nt): return 'Unknown'
    n = str(nt).upper()
    if 'VGLUT' in n: return 'Excitatory'
    if 'GABA' in n or n == 'GLY': return 'Inhibitory'
    return 'Other'

nt_annot['ei_type'] = nt_annot['name'].apply(classify_ei)
print(f'  Excitatory: {(nt_annot["ei_type"]=="Excitatory").sum()}')
print(f'  Inhibitory: {(nt_annot["ei_type"]=="Inhibitory").sum()}')

print('\n[3/11] Building cluster->E/I map...')
cluster_to_ei = {}
for _, row in nt_annot.iterrows():
    members = membership[membership['cluster_annotation_term_label'] == row['label']]
    # Membership uses 'cluster_alias', not 'cluster_label'
    for clust_alias in members['cluster_alias'].values:
        cluster_to_ei[clust_alias] = row['ei_type']
print(f'  Mapped {len(cluster_to_ei)} clusters')

print('\n[4/11] Loading cell metadata (672 MB - will take 2-3 min)...')
cell_df = pd.read_csv(DATA_DIR / 'WHB-10Xv3_cell_metadata.csv')
print(f'  Loaded {len(cell_df):,} cells')

print('\n[5/11] Mapping cells to E/I...')
cell_df['ei_type'] = cell_df['cluster_alias'].map(cluster_to_ei).fillna('Unknown')
print(f'  Excitatory: {(cell_df["ei_type"]=="Excitatory").sum():,} ({(cell_df["ei_type"]=="Excitatory").sum()/len(cell_df)*100:.1f}%)')
print(f'  Inhibitory: {(cell_df["ei_type"]=="Inhibitory").sum():,} ({(cell_df["ei_type"]=="Inhibitory").sum()/len(cell_df)*100:.1f}%)')

print('\n[6/11] Filtering cortical cells...')
cell_df['region_clean'] = cell_df['region_of_interest_label'].str.replace('Human ', '')
cortical = cell_df[cell_df['anatomical_division_label'] == 'Cerebral cortex'].copy()
print(f'  Cortical: {len(cortical):,} cells ({len(cortical)/len(cell_df)*100:.1f}%)')

print('\n[7/11] Calculating E:I by region...')
regions = []
for region in sorted(cortical['region_clean'].unique()):
    rc = cortical[cortical['region_clean'] == region]
    e = (rc['ei_type'] == 'Excitatory').sum()
    i = (rc['ei_type'] == 'Inhibitory').sum()
    regions.append({
        'Region': region, 'Excitatory': e, 'Inhibitory': i, 'Total': len(rc),
        'EI_Ratio': e/i if i > 0 else np.nan,
        'E_percent': (e/len(rc)*100) if len(rc) > 0 else 0
    })

ei_df = pd.DataFrame(regions).sort_values('Total', ascending=False)
print(f'  Calculated E:I for {len(ei_df)} cortical regions')
print('\n  Top 10 regions:')
for _, r in ei_df.head(10).iterrows():
    ratio = f"{r['EI_Ratio']:.2f}" if not np.isnan(r['EI_Ratio']) else "N/A"
    print(f"    {r['Region']:12s} Total:{r['Total']:6,} E:{r['Excitatory']:5,} I:{r['Inhibitory']:5,} Ratio:{ratio:>5s}")

ei_df.to_csv(OUTPUT_DIR / 'ei_ratios_whb_cortical_fine_grained.csv', index=False)

print('\n[8/11] Mapping to Heiss lobes...')
ei_df['Lobe'] = ei_df['Region'].map(REGION_TO_LOBE)
mapped = ei_df[ei_df['Lobe'].notna()]
print(f'  Mapped {len(mapped)}/{len(ei_df)} regions')

print('\n[9/11] Aggregating by lobe...')
lobes = []
for lobe in sorted(mapped['Lobe'].unique()):
    lc = mapped[mapped['Lobe'] == lobe]
    e_sum = lc['Excitatory'].sum()
    i_sum = lc['Inhibitory'].sum()
    lobes.append({
        'Lobe': lobe, 'Excitatory': e_sum, 'Inhibitory': i_sum,
        'Total': lc['Total'].sum(), 'EI_Ratio_aggregated': e_sum/i_sum if i_sum > 0 else np.nan,
        'n_constituent_regions': len(lc),
        'E_percent': (e_sum/lc['Total'].sum()*100) if lc['Total'].sum() > 0 else 0
    })

lobe_df = pd.DataFrame(lobes).sort_values('Total', ascending=False)
print(f'  {len(lobe_df)} lobes')
for _, r in lobe_df.iterrows():
    ratio = f"{r['EI_Ratio_aggregated']:.2f}" if not np.isnan(r['EI_Ratio_aggregated']) else "N/A"
    print(f"    {r['Lobe']:30s} E:I={ratio:>5s} ({r['n_constituent_regions']} regions)")

lobe_df.to_csv(OUTPUT_DIR / 'ei_ratios_aggregated_lobes.csv', index=False)

print('\n[10/11] Integrating with Heiss rCMRGlc...')
lobe_df['rCMRGlc_Mean'] = lobe_df['Lobe'].map({k: v[0] for k, v in HEISS_DATA.items()})
lobe_df['rCMRGlc_SD'] = lobe_df['Lobe'].map({k: v[1] for k, v in HEISS_DATA.items()})
integrated = lobe_df[lobe_df['rCMRGlc_Mean'].notna()].copy()
print(f'  {len(integrated)} lobes with rCMRGlc data')

print('\n  Integrated data:')
for _, r in integrated.iterrows():
    print(f"    {r['Lobe']:30s} E:I={r['EI_Ratio_aggregated']:.2f}  rCMRGlc={r['rCMRGlc_Mean']:.1f}")

integrated.to_csv(OUTPUT_DIR / 'ei_metabolism_integrated.csv', index=False)

print('\n[11/11] Statistical analysis...')
valid = integrated[integrated['EI_Ratio_aggregated'].notna()]

if len(valid) >= 3:
    ei = valid['EI_Ratio_aggregated'].values
    met = valid['rCMRGlc_Mean'].values

    pr, pp = stats.pearsonr(ei, met)
    sr, sp = stats.spearmanr(ei, met)
    slope, intercept, rv, pv, se = stats.linregress(ei, met)

    print(f'\n  Pearson:  r={pr:6.3f}, p={pp:.4f} {"***SIGNIFICANT***" if pp < 0.05 else ""}')
    print(f'  Spearman: rho={sr:6.3f}, p={sp:.4f}')
    print(f'  Regression: rCMRGlc = {slope:.3f}*EI + {intercept:.2f}, R²={rv**2:.3f}, p={pv:.4f}')

    with open(OUTPUT_DIR / 'ei_metabolism_statistics.txt', 'w') as f:
        f.write(f'E:I Ratio vs rCMRGlc Analysis\n{"="*80}\n\n')
        f.write(f'Regions: {len(valid)}\n')
        f.write(f'Pearson r={pr:.4f}, p={pp:.4f}\n')
        f.write(f'Spearman rho={sr:.4f}, p={sp:.4f}\n')
        f.write(f'Regression: rCMRGlc = {slope:.3f}*EI + {intercept:.2f}\n')
        f.write(f'R²={rv**2:.4f}, p={pv:.4f}\n\n')
        for _, r in valid.iterrows():
            f.write(f"{r['Lobe']:30s} E:I={r['EI_Ratio_aggregated']:.2f}  rCMRGlc={r['rCMRGlc_Mean']:.1f}\n")
else:
    print(f'  Insufficient data ({len(valid)} regions)')

elapsed = time.time() - start_time
print(f'\n{"="*80}')
print(f'COMPLETE! Time: {elapsed/60:.1f} minutes')
print(f'Output: {OUTPUT_DIR}')
print(f'{"="*80}')
