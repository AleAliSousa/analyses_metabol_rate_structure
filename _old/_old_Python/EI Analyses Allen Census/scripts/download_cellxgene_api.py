"""
Direct CELLxGENE API Access - Find and Download Cortical Datasets
==================================================================

Uses CELLxGENE REST API to programmatically find and download cortical datasets.
"""

import requests
import json
import urllib.request
from pathlib import Path
import pandas as pd

PROJECT_DIR = Path(r"C:\Users\michaelproulx\MetabolicBrain")
INTEGRATION_DIR = PROJECT_DIR / "cellxgene_integration"
DATA_DIR = INTEGRATION_DIR / "data" / "cortical"
DATA_DIR.mkdir(parents=True, exist_ok=True)

# CELLxGENE API endpoints
API_BASE = "https://api.cellxgene.cziscience.com/curation/v1"


def search_cortical_collections():
    """Search for cortical brain datasets via API"""
    print("=" * 80)
    print("SEARCHING CELLxGENE FOR CORTICAL DATASETS")
    print("=" * 80)

    try:
        # Get all collections
        response = requests.get(f"{API_BASE}/collections", timeout=30)
        response.raise_for_status()

        collections = response.json()
        print(f"\nTotal collections available: {len(collections)}")

        # Filter for brain/cortical collections
        cortical_keywords = [
            'brain', 'cortex', 'cortical', 'cerebral', 'frontal',
            'temporal', 'parietal', 'occipital', 'motor', 'visual'
        ]

        cortical_collections = []

        for collection in collections:
            name = collection.get('name', '').lower()
            description = str(collection.get('description', '')).lower()

            # Check if cortical
            is_cortical = any(kw in name or kw in description for kw in cortical_keywords)

            # Also check if human
            is_human = 'human' in name or 'human' in description or 'homo sapiens' in description

            if is_cortical and is_human:
                cortical_collections.append(collection)
                print(f"\nFound: {collection.get('name')}")
                print(f"  ID: {collection.get('id')}")
                print(f"  DOI: {collection.get('doi', 'N/A')}")

        return cortical_collections

    except requests.exceptions.RequestException as e:
        print(f"\nAPI Error: {e}")
        return None


def get_collection_datasets(collection_id):
    """Get datasets for a specific collection"""
    try:
        response = requests.get(f"{API_BASE}/collections/{collection_id}", timeout=30)
        response.raise_for_status()

        collection_data = response.json()
        datasets = collection_data.get('datasets', [])

        return datasets

    except requests.exceptions.RequestException as e:
        print(f"Error fetching collection {collection_id}: {e}")
        return []


def download_dataset_metadata(dataset):
    """Download metadata for a dataset"""
    dataset_id = dataset.get('id')
    dataset_title = dataset.get('title', 'unknown')

    print(f"\n  Dataset: {dataset_title}")
    print(f"  ID: {dataset_id}")
    print(f"  Cells: {dataset.get('cell_count', 'unknown'):,}" if dataset.get('cell_count') else "")

    # Get download links
    assets = dataset.get('assets', [])

    for asset in assets:
        filetype = asset.get('filetype')

        # Look for H5AD or RDS files with metadata
        if filetype in ['H5AD', 'RDS']:
            url = asset.get('url')
            filesize = asset.get('filesize', 0)

            print(f"    Found {filetype}: {filesize / (1024**2):.1f} MB")

            # Only download if manageable size (<2GB)
            if filesize > 0 and filesize < 2 * 1024**3:
                filename = f"{dataset_id}_{filetype.lower()}.{filetype.lower()}"
                output_path = DATA_DIR / filename

                print(f"    Downloading to: {output_path.name}...")

                try:
                    # Download with progress
                    def report_hook(block_num, block_size, total_size):
                        if total_size > 0 and block_num % 100 == 0:
                            downloaded = block_num * block_size
                            percent = min(100, downloaded * 100 / total_size)
                            print(f"      Progress: {percent:.1f}%", end='\r')

                    urllib.request.urlretrieve(url, output_path, reporthook=report_hook)
                    print(f"\n      SUCCESS: Downloaded {output_path.name}")
                    return output_path

                except Exception as e:
                    print(f"\n      Download failed: {e}")
            else:
                print(f"      Skipping (file too large: {filesize/(1024**3):.1f} GB)")

    return None


def main():
    print("\n" + "=" * 80)
    print("CELLxGENE API - AUTOMATED CORTICAL DATA DOWNLOAD")
    print("=" * 80)

    # Search for cortical collections
    collections = search_cortical_collections()

    if not collections:
        print("\nNo cortical collections found or API error occurred.")
        return

    print(f"\n\nFound {len(collections)} cortical/brain collections")

    # Process each collection
    downloaded_files = []

    for i, collection in enumerate(collections[:3], 1):  # Limit to first 3 to avoid overwhelming
        print(f"\n{'=' * 80}")
        print(f"Processing Collection {i}/{min(3, len(collections))}: {collection.get('name')}")
        print("=" * 80)

        collection_id = collection.get('id')
        datasets = get_collection_datasets(collection_id)

        print(f"\n  Datasets in this collection: {len(datasets)}")

        # Download metadata from first dataset with cortical data
        for dataset in datasets[:2]:  # Try first 2 datasets
            result = download_dataset_metadata(dataset)
            if result:
                downloaded_files.append(result)
                break  # Move to next collection after first successful download

    print("\n" + "=" * 80)
    print("DOWNLOAD SUMMARY")
    print("=" * 80)

    if downloaded_files:
        print(f"\nSuccessfully downloaded {len(downloaded_files)} files:")
        for f in downloaded_files:
            print(f"  - {f}")
    else:
        print("\nNo files downloaded. Possible reasons:")
        print("  - All files too large (>2GB)")
        print("  - Download URLs not accessible")
        print("  - API rate limiting")

    return downloaded_files


if __name__ == "__main__":
    result = main()
