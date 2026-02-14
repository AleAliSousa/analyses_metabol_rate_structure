"""
Assess Data Requirements for E:I Analysis
==========================================

This script checks available cortical datasets and estimates download sizes
WITHOUT downloading any data. Run this first to understand storage needs.

Author: Claude Code
Date: 2026-02-08
"""

import requests
import json
from pathlib import Path

PROJECT_DIR = Path(r"C:\Users\michaelproulx\MetabolicBrain")
INTEGRATION_DIR = PROJECT_DIR / "cellxgene_integration"

def assess_hca_cortical_data():
    """
    Check HCA Brain Atlas v1.0 cortical datasets via CELLxGENE API
    """
    print("\n" + "=" * 80)
    print("ASSESSING HCA BRAIN ATLAS CORTICAL DATA")
    print("=" * 80)

    # CELLxGENE API endpoint for collections
    api_base = "https://api.cellxgene.cziscience.com/curation/v1"

    print("\nSearching for Human Brain Cell Atlas datasets...")
    print("(This only queries metadata - no downloads)")

    try:
        # Search for collections
        response = requests.get(f"{api_base}/collections")

        if response.status_code != 200:
            print(f"  API request failed: {response.status_code}")
            return None

        collections = response.json()

        # Find brain-related collections
        brain_collections = []
        for collection in collections:
            name = collection.get('name', '').lower()
            description = collection.get('description', '').lower()

            if any(term in name or term in description for term in
                   ['brain', 'cortex', 'cortical', 'cerebral']):
                brain_collections.append(collection)

        print(f"\nFound {len(brain_collections)} brain-related collections")

        # Get details for each collection
        total_size = 0
        total_cells = 0
        cortical_datasets = []

        for collection in brain_collections:
            collection_id = collection.get('id')
            name = collection.get('name', 'Unknown')

            print(f"\n  Collection: {name}")
            print(f"  ID: {collection_id}")

            # Get datasets in this collection
            datasets_response = requests.get(
                f"{api_base}/collections/{collection_id}"
            )

            if datasets_response.status_code == 200:
                collection_details = datasets_response.json()
                datasets = collection_details.get('datasets', [])

                for dataset in datasets:
                    tissue = dataset.get('tissue', [])
                    n_cells = dataset.get('cell_count', 0)

                    # Check if cortical
                    tissue_str = str(tissue).lower()
                    is_cortical = any(term in tissue_str for term in
                                    ['cortex', 'frontal', 'temporal', 'parietal',
                                     'occipital', 'insula', 'visual'])

                    if is_cortical or 'brain' in tissue_str:
                        # Get file assets to estimate size
                        assets = dataset.get('assets', [])
                        dataset_size = 0

                        for asset in assets:
                            if asset.get('filetype') in ['H5AD', 'RDS']:
                                size_bytes = asset.get('filesize', 0)
                                dataset_size += size_bytes

                        if dataset_size > 0:
                            cortical_datasets.append({
                                'collection': name,
                                'dataset_id': dataset.get('id'),
                                'tissue': tissue,
                                'n_cells': n_cells,
                                'size_bytes': dataset_size,
                                'size_mb': dataset_size / (1024**2),
                                'size_gb': dataset_size / (1024**3)
                            })

                            total_size += dataset_size
                            total_cells += n_cells

        return {
            'cortical_datasets': cortical_datasets,
            'total_size_gb': total_size / (1024**3),
            'total_cells': total_cells,
            'n_datasets': len(cortical_datasets)
        }

    except Exception as e:
        print(f"\n  Error: {e}")
        return None


def assess_abc_atlas_data():
    """
    Check Allen Brain Cell Atlas data requirements
    """
    print("\n" + "=" * 80)
    print("ASSESSING ALLEN BRAIN CELL (ABC) ATLAS DATA")
    print("=" * 80)

    print("\nABC Atlas access methods:")
    print("1. Via AWS S3 (direct download)")
    print("2. Via abc_atlas_access Python package (metadata queries)")

    # Known ABC Atlas datasets from documentation
    abc_datasets = {
        'ASAP-PMDBS Human Postmortem': {
            'cells': 3_000_000,
            'donors': 220,
            'metadata_size_mb': 100,  # Estimated
            'expression_size_gb': 50,  # Estimated for full matrix
            'description': 'Whole brain coverage including cortex'
        },
        'WHB (Whole Human Brain)': {
            'cells': 4_000_000,
            'donors': 'multiple',
            'metadata_size_mb': 150,  # Estimated
            'expression_size_gb': 80,  # Estimated
            'description': 'Complete brain atlas with cortical regions'
        }
    }

    print("\nKnown ABC Atlas datasets:")
    total_metadata_mb = 0
    total_expression_gb = 0

    for name, info in abc_datasets.items():
        print(f"\n  {name}:")
        print(f"    Cells: {info['cells']:,}")
        print(f"    Donors: {info['donors']}")
        print(f"    Metadata (estimated): {info['metadata_size_mb']} MB")
        print(f"    Expression matrix (estimated): {info['expression_size_gb']} GB")
        print(f"    Description: {info['description']}")

        total_metadata_mb += info['metadata_size_mb']
        total_expression_gb += info['expression_size_gb']

    return {
        'datasets': abc_datasets,
        'total_metadata_mb': total_metadata_mb,
        'total_expression_gb': total_expression_gb
    }


def generate_recommendation(hca_info, abc_info):
    """
    Generate storage and processing recommendations
    """
    print("\n" + "=" * 80)
    print("STORAGE AND PROCESSING RECOMMENDATIONS")
    print("=" * 80)

    print("\n" + "─" * 80)
    print("METADATA-ONLY APPROACH (RECOMMENDED FOR LIMITED RESOURCES)")
    print("─" * 80)

    metadata_total_mb = 0

    if hca_info:
        # Assume metadata is ~10% of full H5AD size
        hca_metadata_mb = (hca_info['total_size_gb'] * 1024) * 0.1
        metadata_total_mb += hca_metadata_mb
        print(f"\nHCA Metadata only: ~{hca_metadata_mb:.0f} MB")
        print(f"  ({hca_info['total_cells']:,} cells from {hca_info['n_datasets']} datasets)")

    if abc_info:
        abc_metadata_mb = abc_info['total_metadata_mb']
        metadata_total_mb += abc_metadata_mb
        print(f"\nABC Metadata only: ~{abc_metadata_mb:.0f} MB")

    print(f"\n  TOTAL METADATA: ~{metadata_total_mb:.0f} MB (~{metadata_total_mb/1024:.1f} GB)")
    print("\n  This approach:")
    print("    ✓ Downloads only cell annotations and region labels")
    print("    ✓ Enables COUNT-BASED E:I ratio calculation")
    print("    ✓ Fast processing on laptop")
    print("    ✓ Sufficient for correlation analysis")

    print("\n" + "─" * 80)
    print("FULL EXPRESSION DATA (NOT RECOMMENDED FOR LAPTOP)")
    print("─" * 80)

    full_data_gb = 0
    if hca_info:
        full_data_gb += hca_info['total_size_gb']
    if abc_info:
        full_data_gb += abc_info['total_expression_gb']

    print(f"\n  TOTAL WITH EXPRESSION: ~{full_data_gb:.1f} GB")
    print("\n  This approach:")
    print("    ✗ Large downloads (may take hours)")
    print("    ✗ High memory requirements (16+ GB RAM)")
    print("    ✗ Slow processing on laptop")
    print("    + Enables EXPRESSION-BASED E:I validation")

    print("\n" + "─" * 80)
    print("RECOMMENDED STRATEGY")
    print("─" * 80)

    print("\nFor your laptop:")
    print("  1. Download METADATA ONLY (~0.5-1 GB total)")
    print("  2. Calculate E:I ratios using cell type annotations (count-based)")
    print("  3. This provides sufficient data for hypothesis testing")
    print("  4. Expression-based validation can be done later if needed")

    print(f"\n  Estimated disk space needed: {metadata_total_mb/1024:.1f} GB")
    print(f"  Estimated RAM needed: 4-8 GB")
    print(f"  Estimated processing time: 15-30 minutes")

    return {
        'metadata_only_gb': metadata_total_mb / 1024,
        'full_data_gb': full_data_gb,
        'recommendation': 'metadata_only'
    }


def main():
    """
    Main assessment function
    """
    print("\n" + "=" * 80)
    print("E:I RATIO ANALYSIS - DATA REQUIREMENTS ASSESSMENT")
    print("=" * 80)

    print("\nThis script checks available datasets and estimates storage needs.")
    print("No data will be downloaded during this assessment.")

    # Assess HCA data
    hca_info = assess_hca_cortical_data()

    # Assess ABC data
    abc_info = assess_abc_atlas_data()

    # Generate recommendations
    if hca_info or abc_info:
        recommendation = generate_recommendation(hca_info, abc_info)

        # Save summary
        summary_path = INTEGRATION_DIR / "ei_analysis" / "data_requirements_summary.txt"
        summary_path.parent.mkdir(parents=True, exist_ok=True)

        with open(summary_path, 'w') as f:
            f.write("E:I Ratio Analysis - Data Requirements Summary\n")
            f.write("=" * 80 + "\n\n")

            if hca_info:
                f.write(f"HCA Datasets: {hca_info['n_datasets']}\n")
                f.write(f"HCA Total Cells: {hca_info['total_cells']:,}\n")
                f.write(f"HCA Metadata Size: ~{(hca_info['total_size_gb'] * 0.1):.1f} GB\n")
                f.write(f"HCA Full Size: ~{hca_info['total_size_gb']:.1f} GB\n\n")

            if abc_info:
                f.write(f"ABC Metadata Size: ~{abc_info['total_metadata_mb']/1024:.1f} GB\n")
                f.write(f"ABC Full Size: ~{abc_info['total_expression_gb']:.1f} GB\n\n")

            f.write(f"RECOMMENDED DOWNLOAD: {recommendation['metadata_only_gb']:.1f} GB (metadata only)\n")
            f.write(f"FULL DOWNLOAD (not recommended): {recommendation['full_data_gb']:.1f} GB\n")

        print(f"\n\nSummary saved to: {summary_path}")

    else:
        print("\nCould not assess data requirements. Check internet connection.")

    print("\n" + "=" * 80)
    print("ASSESSMENT COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    main()
