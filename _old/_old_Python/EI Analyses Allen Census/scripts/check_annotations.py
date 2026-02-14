"""Check annotation content"""
import pandas as pd
from pathlib import Path

data_dir = Path(r'C:\Users\michaelproulx\MetabolicBrain\cellxgene_integration\data')
annotations = pd.read_csv(data_dir / 'cluster_annotation_term.csv')

# Check annotation set names
print('=' * 80)
print('ANNOTATION SETS')
print('=' * 80)
sets = annotations['cluster_annotation_term_set_name'].unique()
for s in sets:
    count = len(annotations[annotations['cluster_annotation_term_set_name'] == s])
    print(f'  {s}: {count} annotations')

# Look at neighborhoods (regions)
print('\n' + '=' * 80)
print('NEIGHBORHOODS (Regions)')
print('=' * 80)
neigh = annotations[annotations['cluster_annotation_term_set_name'] == 'Neighborhood']
if len(neigh) > 0:
    for i, name in enumerate(neigh['name'].values, 1):
        print(f'  {i:2d}. {name}')

# Look at classes
print('\n' + '=' * 80)
print('CLASSES (Broad Cell Types)')
print('=' * 80)
classes = annotations[annotations['cluster_annotation_term_set_name'] == 'Class']
if len(classes) > 0:
    for i, name in enumerate(classes['name'].values, 1):
        print(f'  {i:2d}. {name}')

# Look at subclasses
print('\n' + '=' * 80)
print('SUBCLASSES (Detailed Cell Types) - First 30')
print('=' * 80)
subclass = annotations[annotations['cluster_annotation_term_set_name'] == 'Subclass']
if len(subclass) > 0:
    for i, name in enumerate(subclass['name'].values[:30], 1):
        print(f'  {i:2d}. {name}')
    if len(subclass) > 30:
        print(f'  ... and {len(subclass) - 30} more')
