"""
Download Cortical Data from Allen Brain WHB (Whole Human Brain) Dataset
========================================================================

Uses the SAME method that successfully downloaded basal ganglia data.
"""

import urllib.request
from pathlib import Path

PROJECT_DIR = Path(r"C:\Users\michaelproulx\MetabolicBrain")
INTEGRATION_DIR = PROJECT_DIR / "cellxgene_integration"
DATA_DIR = INTEGRATION_DIR / "data" / "cortical"
DATA_DIR.mkdir(parents=True, exist_ok=True)

# S3 base URL
S3_BASE = "https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com"

# Datasets to try (from manifest)
DATASETS_TO_TRY = [
    {
        'name': 'WHB-10Xv3',
        'description': 'Whole Human Brain - 10x v3',
        'dates': ['20230830', '20240330'],
        'files': ['cell_metadata.csv']
    },
    {
        'name': 'ASAP-PMDBS-10X',
        'description': 'ASAP Postmortem Brain (whole brain)',
        'dates': ['metadata'],
        'files': ['cell_metadata.csv']
    }
]


def download_with_progress(url, output_path):
    """Download with progress indicator"""
    print(f"\nAttempting: {url}")
    print(f"Saving to: {output_path.name}")

    try:
        def report_hook(block_num, block_size, total_size):
            if total_size > 0 and block_num % 100 == 0:
                downloaded = block_num * block_size
                percent = min(100, downloaded * 100 / total_size)
                mb_downloaded = downloaded / (1024**2)
                mb_total = total_size / (1024**2)
                print(f"  Progress: {percent:.1f}% ({mb_downloaded:.1f}/{mb_total:.1f} MB)", end='\r')

        urllib.request.urlretrieve(url, output_path, reporthook=report_hook)
        print(f"\n  SUCCESS! Downloaded {output_path.stat().st_size / (1024**2):.1f} MB")
        return True

    except urllib.error.HTTPError as e:
        print(f"\n  HTTP Error {e.code}: {e.reason}")
        return False
    except Exception as e:
        print(f"\n  Error: {e}")
        return False


def main():
    print("=" * 80)
    print("DOWNLOADING WHOLE HUMAN BRAIN METADATA")
    print("Using proven S3 direct download method")
    print("=" * 80)

    downloaded_files = []

    for dataset in DATASETS_TO_TRY:
        print(f"\n\n{'='*80}")
        print(f"Trying: {dataset['name']} - {dataset['description']}")
        print("="*80)

        for date in dataset['dates']:
            for filename in dataset['files']:
                # Construct URL
                url = f"{S3_BASE}/metadata/{dataset['name']}/{date}/{filename}"
                output_path = DATA_DIR / f"{dataset['name']}_{filename}"

                success = download_with_progress(url, output_path)

                if success:
                    downloaded_files.append((dataset['name'], output_path))
                    print(f"\n  File saved: {output_path}")
                    break  # Move to next dataset after first success

            if downloaded_files and downloaded_files[-1][0] == dataset['name']:
                break  # Already got this dataset

    print("\n" + "=" * 80)
    print("DOWNLOAD SUMMARY")
    print("=" * 80)

    if downloaded_files:
        print(f"\nSuccessfully downloaded {len(downloaded_files)} datasets:")
        for name, path in downloaded_files:
            size_mb = path.stat().st_size / (1024**2)
            print(f"  {name}: {size_mb:.1f} MB - {path.name}")

        print("\nNext step: Filter for cortical cells and extract E:I ratios")
        return downloaded_files
    else:
        print("\nNo downloads succeeded.")
        print("Datasets tried:", [d['name'] for d in DATASETS_TO_TRY])
        return None


if __name__ == "__main__":
    result = main()
    if result:
        print("\nSUCCESS!")
    else:
        print("\nFAILED - check URLs or network connection")
