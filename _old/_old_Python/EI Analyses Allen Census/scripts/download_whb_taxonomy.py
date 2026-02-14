"""
Download WHB Taxonomy Files for E/I Annotation
"""
import urllib.request
from pathlib import Path

S3_BASE = "https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com"
OUTPUT_DIR = Path(r"C:\Users\michaelproulx\MetabolicBrain\cellxgene_integration\data\cortical")

files_to_download = [
    'cluster_annotation_term.csv',
    'cluster_to_cluster_annotation_membership.csv',
]

print('='*80)
print('DOWNLOADING WHB TAXONOMY FILES')
print('='*80)

for filename in files_to_download:
    url = f"{S3_BASE}/metadata/WHB-taxonomy/metadata/{filename}"
    output_path = OUTPUT_DIR / f"WHB_{filename}"

    print(f'\nDownloading: {filename}')
    print(f'URL: {url}')

    try:
        urllib.request.urlretrieve(url, output_path)
        size_kb = output_path.stat().st_size / 1024
        print(f'SUCCESS: {size_kb:.1f} KB')
    except Exception as e:
        print(f'ERROR: {e}')

print('\n' + '='*80)
print('DOWNLOAD COMPLETE')
print('='*80)
