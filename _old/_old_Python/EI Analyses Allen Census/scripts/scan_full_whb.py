"""
Process FULL WHB Dataset to Find ALL Regions Including Cortical
"""
import pandas as pd
from pathlib import Path

whb_file = Path(r'C:\Users\michaelproulx\MetabolicBrain\cellxgene_integration\data\cortical\WHB-10Xv3_cell_metadata.csv')

print('='*80)
print('PROCESSING FULL WHB DATASET')
print('='*80)

print(f'\nFile: {whb_file}')
print(f'Size: {whb_file.stat().st_size / (1024**2):.1f} MB')

print('\nScanning entire dataset for all unique regions...')
print('(This will take a few minutes)')

# Read in chunks to handle large file
chunk_size = 200000
all_regions = {}
total_rows = 0

for i, chunk in enumerate(pd.read_csv(whb_file, chunksize=chunk_size), 1):
    total_rows += len(chunk)

    if 'region_of_interest_label' in chunk.columns:
        region_counts = chunk['region_of_interest_label'].value_counts()

        for region, count in region_counts.items():
            if region in all_regions:
                all_regions[region] += count
            else:
                all_regions[region] = count

    if i % 5 == 0:
        print(f'  Processed {total_rows:,} rows, found {len(all_regions)} unique regions...')

print(f'\n  Total rows processed: {total_rows:,}')
print(f'  Total unique regions: {len(all_regions)}')

# Sort by count
sorted_regions = sorted(all_regions.items(), key=lambda x: x[1], reverse=True)

print('\n' + '='*80)
print('ALL REGIONS FOUND:')
print('='*80)

for i, (region, count) in enumerate(sorted_regions, 1):
    print(f'{i:3d}. {region:40s} {count:>10,} cells')

# Check for cortical
cortical_kw = ['ctx', 'cortex', 'frontal', 'temporal', 'parietal', 'occipital', 'motor', 'visual', 'isocortex', 'sensory']
cortical_regions = [(r, c) for r, c in sorted_regions if any(kw in str(r).lower() for kw in cortical_kw)]

print('\n' + '='*80)
if cortical_regions:
    print(f'CORTICAL REGIONS FOUND: {len(cortical_regions)}')
    print('='*80)
    for region, count in cortical_regions:
        print(f'  {region}: {count:,} cells')
else:
    print('NO CORTICAL REGIONS FOUND')
    print('='*80)
    print('\nDataset appears to be subcortical/brainstem only')
    print('Will need alternative data source for cortical analysis')

print('\n' + '='*80)
print('SCAN COMPLETE')
print('='*80)
