"""
Verify WHB Cortical Data and Prepare for E:I Analysis
"""
import pandas as pd
from pathlib import Path

data_file = Path(r'C:\Users\michaelproulx\MetabolicBrain\cellxgene_integration\data\cortical\WHB-10Xv3_cell_metadata.csv')

print('='*80)
print('STEP 1: LOADING AND VERIFYING WHB DATA')
print('='*80)

# Load sample to check structure
print(f'\nLoading sample rows from {data_file.name}...')
df_sample = pd.read_csv(data_file, nrows=10000)

print(f'Sample: {len(df_sample):,} rows, {len(df_sample.columns)} columns')

# Show all columns
print('\nAll columns:')
for i, col in enumerate(df_sample.columns, 1):
    print(f'  {i:2d}. {col}')

# Find region and cell type columns
region_cols = [col for col in df_sample.columns if any(kw in col.lower() for kw in ['region', 'structure', 'parcel', 'division'])]
celltype_cols = [col for col in df_sample.columns if any(kw in col.lower() for kw in ['class', 'subclass', 'type', 'cluster', 'supertype'])]

print(f'\nRegion columns: {region_cols}')
print(f'Cell type columns: {celltype_cols}')

if region_cols:
    region_col = region_cols[0]
    print(f'\n{region_col} - Unique values (top 30):')
    for i, (val, count) in enumerate(df_sample[region_col].value_counts().head(30).items(), 1):
        print(f'  {i:2d}. {val}: {count:,}')

    # Check for cortical
    cortical_kw = ['ctx', 'cortex', 'frontal', 'temporal', 'parietal', 'occipital', 'motor', 'visual', 'isocortex']
    is_cortical = df_sample[region_col].astype(str).str.contains('|'.join(cortical_kw), case=False, na=False)
    print(f'\nCortical cells: {is_cortical.sum():,} ({is_cortical.sum()/len(df_sample)*100:.1f}%)')

if celltype_cols:
    ct_col = celltype_cols[0]
    print(f'\n{ct_col} - Unique values (top 20):')
    for i, (val, count) in enumerate(df_sample[ct_col].value_counts().head(20).items(), 1):
        print(f'  {i:2d}. {val}: {count:,}')

print('\n' + '='*80)
print('VERIFICATION COMPLETE - Data structure identified')
print('='*80)
