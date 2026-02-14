"""
Check existing metadata for cortical data
"""
import pandas as pd
from pathlib import Path

# Check existing metadata
metadata_path = Path(r'C:\Users\michaelproulx\MetabolicBrain\cellxgene_integration\data\cell_metadata.csv')

print('=' * 80)
print('CHECKING EXISTING METADATA')
print('=' * 80)

print(f'\nFile: {metadata_path}')
print(f'Size: {metadata_path.stat().st_size / (1024**2):.1f} MB')

# Read first few rows to check structure
print('\nLoading sample rows...')
df = pd.read_csv(metadata_path, nrows=5000)

print(f'Sample size: {len(df):,} rows')
print(f'Total columns: {len(df.columns)}')

# Show columns
print('\nColumns:')
for i, col in enumerate(df.columns, 1):
    print(f'  {i:2d}. {col}')

# Check for region column
print('\n' + '=' * 80)
print('REGION COLUMNS')
print('=' * 80)

region_cols = [col for col in df.columns if any(kw in col.lower() for kw in ['region', 'structure', 'tissue', 'parcellation'])]
print(f'\nFound {len(region_cols)} region-related columns:')
for col in region_cols:
    print(f'  - {col}')

if region_cols:
    region_col = region_cols[0]
    print(f'\nSample regions from "{region_col}":')
    unique_regions = df[region_col].unique()[:20]
    for i, region in enumerate(unique_regions, 1):
        print(f'  {i:2d}. {region}')

    if len(df[region_col].unique()) > 20:
        print(f'  ... and {len(df[region_col].unique()) - 20} more')

    # Check if cortical
    cortical_keywords = ['cortex', 'frontal', 'temporal', 'parietal', 'occipital', 'visual', 'motor', 'ctx']
    is_cortical = df[region_col].astype(str).str.contains('|'.join(cortical_keywords), case=False, na=False)
    n_cortical = is_cortical.sum()

    print(f'\nCortical cells in sample: {n_cortical:,} ({n_cortical/len(df)*100:.1f}%)')

    if n_cortical > 0:
        print('\nCortical regions found:')
        cortical_regions = df[is_cortical][region_col].unique()[:10]
        for region in cortical_regions:
            print(f'  - {region}')

# Check for cell type columns
print('\n' + '=' * 80)
print('CELL TYPE COLUMNS')
print('=' * 80)

celltype_cols = [col for col in df.columns if any(kw in col.lower() for kw in ['cell_type', 'celltype', 'class', 'cluster', 'supertype'])]
print(f'\nFound {len(celltype_cols)} cell type-related columns:')
for col in celltype_cols[:10]:
    print(f'  - {col}')

if celltype_cols:
    ct_col = celltype_cols[0]
    print(f'\nSample cell types from "{ct_col}":')
    unique_types = df[ct_col].unique()[:15]
    for i, ctype in enumerate(unique_types, 1):
        print(f'  {i:2d}. {ctype}')

print('\n' + '=' * 80)
print('ASSESSMENT COMPLETE')
print('=' * 80)
