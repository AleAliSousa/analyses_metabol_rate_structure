"""
Calculate E:I Ratios from HMBA Basal Ganglia Data
==================================================

Uses existing downloaded data with full E/I annotations.
Demonstrates complete methodology that can be applied to cortical data later.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict

PROJECT_DIR = Path(r"C:\Users\michaelproulx\MetabolicBrain")
DATA_DIR = PROJECT_DIR / "cellxgene_integration" / "data"
OUTPUT_DIR = PROJECT_DIR / "cellxgene_integration" / "ei_analysis" / "data" / "cortical"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print('='*80)
print('E:I RATIO CALCULATION - HMBA BASAL GANGLIA DATA')
print('='*80)

# Load annotations
print('\nLoading cluster annotations...')
annotations = pd.read_csv(DATA_DIR / 'cluster_annotation_term.csv')
print(f'  Loaded {len(annotations):,} annotations')

# Load cell-to-cluster membership
print('Loading cell-to-cluster membership...')
membership = pd.read_csv(DATA_DIR / 'cell_to_cluster_membership.csv')
print(f'  Loaded {len(membership):,} cell-cluster mappings')

# Get Class level annotations (has E/I information)
class_annotations = annotations[annotations['cluster_annotation_term_set_name'] == 'Class'].copy()
print(f'\nClass annotations: {len(class_annotations)}')

# Create E/I classification
def classify_ei(class_name):
    if pd.isna(class_name):
        return 'Unknown'

    name_lower = str(class_name).lower()

    # Excitatory keywords
    if any(kw in name_lower for kw in ['glut', 'glutamatergic', 'f m glut', 'f glut']):
        return 'Excitatory'

    # Inhibitory keywords
    if any(kw in name_lower for kw in ['gaba', 'gabaergic', 'cn', 'inhibitory']):
        return 'Inhibitory'

    return 'Non-neuronal'

class_annotations['ei_type'] = class_annotations['name'].apply(classify_ei)

print('\nE/I Classification:')
for ei_type, count in class_annotations['ei_type'].value_counts().items():
    print(f'  {ei_type}: {count}')

print('\nExcitatory classes:')
excit = class_annotations[class_annotations['ei_type'] == 'Excitatory']
for name in excit['name'].values:
    print(f'  - {name}')

print('\nInhibitory classes:')
inhib = class_annotations[class_annotations['ei_type'] == 'Inhibitory']
for name in inhib['name'].values:
    print(f'  - {name}')

# Get cluster to class mapping
print('\nMapping clusters to classes...')
cluster_to_class = {}

for _, row in class_annotations.iterrows():
    label = row['label']
    ei_type = row['ei_type']

    # Find all clusters under this class
    children = annotations[annotations['parent_term_label'] == label]

    for _, child in children.iterrows():
        cluster_to_class[child['label']] = ei_type

print(f'Mapped {len(cluster_to_class)} clusters to E/I types')

# Map cells to E/I types
print('\nMapping cells to E/I types...')
membership['ei_type'] = membership['cluster_label'].map(cluster_to_class)
membership['ei_type'] = membership['ei_type'].fillna('Unknown')

print('\nCell E/I distribution:')
ei_counts = membership['ei_type'].value_counts()
total = len(membership)
for ei_type, count in ei_counts.items():
    print(f'  {ei_type}: {count:,} ({count/total*100:.1f}%)')

# Calculate E:I ratios by neighborhood
print('\n' + '='*80)
print('CALCULATING E:I RATIOS BY REGION')
print('='*80)

# Get neighborhood annotations
neighborhoods = annotations[annotations['cluster_annotation_term_set_name'] == 'Neighborhood']
print(f'\nNeighborhoods: {len(neighborhoods)}')
for name in neighborhoods['name'].values:
    print(f'  - {name}')

# Map clusters to neighborhoods
cluster_to_neighborhood = {}
for _, row in neighborhoods.iterrows():
    label = row['label']
    name = row['name']

    # Find all children recursively
    def get_all_children(parent_label):
        children = annotations[annotations['parent_term_label'] == parent_label]
        all_labels = list(children['label'].values)

        for child_label in list(all_labels):
            all_labels.extend(get_all_children(child_label))

        return all_labels

    child_labels = get_all_children(label)
    for child_label in child_labels:
        cluster_to_neighborhood[child_label] = name

membership['neighborhood'] = membership['cluster_label'].map(cluster_to_neighborhood)

# Calculate E:I by neighborhood
ei_by_region = []

for neighborhood in neighborhoods['name'].values:
    cells_in_region = membership[membership['neighborhood'] == neighborhood]

    if len(cells_in_region) == 0:
        continue

    e_count = (cells_in_region['ei_type'] == 'Excitatory').sum()
    i_count = (cells_in_region['ei_type'] == 'Inhibitory').sum()
    total_count = len(cells_in_region)

    if i_count > 0:
        ei_ratio = e_count / i_count
    else:
        ei_ratio = np.nan

    ei_by_region.append({
        'Region': neighborhood,
        'Excitatory': e_count,
        'Inhibitory': i_count,
        'Total': total_count,
        'EI_Ratio': ei_ratio,
        'E_percent': (e_count / total_count * 100) if total_count > 0 else 0,
        'I_percent': (i_count / total_count * 100) if total_count > 0 else 0
    })

ei_df = pd.DataFrame(ei_by_region)
ei_df = ei_df.sort_values('Total', ascending=False)

print('\nE:I Ratios by Region:')
print('='*80)
for _, row in ei_df.iterrows():
    print(f"\n{row['Region']}:")
    print(f"  Excitatory: {row['Excitatory']:,}")
    print(f"  Inhibitory: {row['Inhibitory']:,}")
    print(f"  E:I Ratio: {row['EI_Ratio']:.2f}" if not np.isnan(row['EI_Ratio']) else "  E:I Ratio: N/A")
    print(f"  Total cells: {row['Total']:,}")

# Save results
output_file = OUTPUT_DIR / 'ei_ratios_hmba_basal_ganglia.csv'
ei_df.to_csv(output_file, index=False)

print('\n' + '='*80)
print('RESULTS SAVED')
print('='*80)
print(f'\nOutput: {output_file}')
print(f'Size: {output_file.stat().st_size / 1024:.1f} KB')

print('\n' + '='*80)
print('E:I RATIO CALCULATION COMPLETE')
print('='*80)
