"""
Phase 1: Data Access and Exploration for HCA Brain v1.0
========================================================

This script provides tools to:
1. Guide downloading Human Cell Atlas Brain v1.0 data from CELLxGENE
2. Explore the downloaded H5AD file structure
3. Extract metadata about brain regions and cell types

Author: Claude Code
Date: 2026-01-30
"""

import os
import pandas as pd
import numpy as np
import anndata as ad
from pathlib import Path

# Configuration
PROJECT_DIR = Path(r"C:\Sandbox\michaelproulx\3pAgentBox\user\current\MetabolicBrain")
INTEGRATION_DIR = PROJECT_DIR / "cellxgene_integration"
DATA_DIR = INTEGRATION_DIR / "data"

# Create directories
DATA_DIR.mkdir(parents=True, exist_ok=True)

def print_download_instructions():
    """
    Print instructions for downloading HCA Brain v1.0 data from CELLxGENE
    """
    print("=" * 80)
    print("HUMAN CELL ATLAS BRAIN v1.0 - DATA DOWNLOAD INSTRUCTIONS")
    print("=" * 80)
    print("\nRecommended Data Source: HCA Brain v1.0 via CELLxGENE")
    print("\nSTEPS:")
    print("\n1. Visit the HCA Data Portal:")
    print("   https://data.humancellatlas.org/hca-bio-networks/nervous-system/atlases/brain-v1-0")
    print("\n2. OR visit CELLxGENE directly:")
    print("   https://cellxgene.cziscience.com/")
    print("   Search for 'Human Brain Cell Atlas' or 'brain v1.0'")
    print("\n3. Look for datasets with:")
    print("   - Organism: Homo sapiens (Human)")
    print("   - Tissue: Brain (multiple regions)")
    print("   - Assay: snRNA-seq or scRNA-seq")
    print("   - Adult samples (to match PET metabolic data)")
    print("\n4. Download the dataset(s) in H5AD format")
    print(f"   Save to: {DATA_DIR}")
    print("\n5. Preferred dataset characteristics:")
    print("   - 3+ million cells")
    print("   - Coverage of major brain regions (cortex, hippocampus, thalamus, etc.)")
    print("   - Cell type annotations (neurons, astrocytes, oligodendrocytes, microglia)")
    print("\nALTERNATIVE: If direct download is not available, the script will guide you")
    print("through using the CELLxGENE API or web interface.")
    print("\n" + "=" * 80)
    print()


def explore_h5ad_file(h5ad_path):
    """
    Explore the structure of a downloaded H5AD file

    Args:
        h5ad_path (str): Path to the H5AD file

    Returns:
        dict: Dictionary containing exploration results
    """
    print(f"\nExploring H5AD file: {h5ad_path}")
    print("=" * 80)

    try:
        # Load the AnnData object
        adata = ad.read_h5ad(h5ad_path)

        results = {
            'n_cells': adata.n_obs,
            'n_genes': adata.n_vars,
            'obs_columns': list(adata.obs.columns),
            'var_columns': list(adata.var.columns),
        }

        # Basic information
        print(f"\nDataset Overview:")
        print(f"  Number of cells: {adata.n_obs:,}")
        print(f"  Number of genes: {adata.n_vars:,}")

        # Cell metadata (obs)
        print(f"\n  Cell Metadata Columns ({len(adata.obs.columns)}):")
        for col in adata.obs.columns[:20]:  # Show first 20
            n_unique = adata.obs[col].nunique()
            dtype = adata.obs[col].dtype
            print(f"    - {col}: {n_unique} unique values ({dtype})")

        if len(adata.obs.columns) > 20:
            print(f"    ... and {len(adata.obs.columns) - 20} more columns")

        # Gene metadata (var)
        print(f"\n  Gene Metadata Columns ({len(adata.var.columns)}):")
        for col in adata.var.columns[:10]:  # Show first 10
            print(f"    - {col}")

        if len(adata.var.columns) > 10:
            print(f"    ... and {len(adata.var.columns) - 10} more columns")

        # Look for region annotations
        print("\n  Searching for brain region annotations...")
        region_cols = [col for col in adata.obs.columns
                      if any(keyword in col.lower() for keyword in
                            ['region', 'tissue', 'area', 'structure', 'anatomy', 'location'])]

        if region_cols:
            print(f"  Found {len(region_cols)} potential region column(s):")
            for col in region_cols:
                unique_regions = adata.obs[col].unique()
                print(f"\n    {col}: {len(unique_regions)} unique regions")
                if len(unique_regions) <= 50:
                    for region in sorted(unique_regions)[:20]:
                        count = (adata.obs[col] == region).sum()
                        print(f"      - {region}: {count:,} cells")
                    if len(unique_regions) > 20:
                        print(f"      ... and {len(unique_regions) - 20} more regions")
            results['region_columns'] = region_cols
        else:
            print("  No obvious region columns found. Check obs columns manually.")
            results['region_columns'] = []

        # Look for cell type annotations
        print("\n  Searching for cell type annotations...")
        celltype_cols = [col for col in adata.obs.columns
                        if any(keyword in col.lower() for keyword in
                              ['cell_type', 'celltype', 'cluster', 'annotation', 'class'])]

        if celltype_cols:
            print(f"  Found {len(celltype_cols)} potential cell type column(s):")
            for col in celltype_cols:
                unique_types = adata.obs[col].unique()
                print(f"\n    {col}: {len(unique_types)} unique cell types")
                if len(unique_types) <= 30:
                    for celltype in sorted(unique_types)[:15]:
                        count = (adata.obs[col] == celltype).sum()
                        print(f"      - {celltype}: {count:,} cells")
                    if len(unique_types) > 15:
                        print(f"      ... and {len(unique_types) - 15} more cell types")
            results['celltype_columns'] = celltype_cols
        else:
            print("  No obvious cell type columns found. Check obs columns manually.")
            results['celltype_columns'] = []

        # Summary statistics
        print("\n" + "=" * 80)
        print("SUMMARY:")
        print(f"  Dataset contains {adata.n_obs:,} cells from {adata.n_vars:,} genes")
        print(f"  Region columns identified: {len(results.get('region_columns', []))}")
        print(f"  Cell type columns identified: {len(results.get('celltype_columns', []))}")
        print("=" * 80)

        return results

    except Exception as e:
        print(f"Error exploring H5AD file: {e}")
        return None


def list_available_h5ad_files():
    """
    List all H5AD files in the data directory
    """
    h5ad_files = list(DATA_DIR.glob("*.h5ad"))

    if not h5ad_files:
        print(f"\nNo H5AD files found in {DATA_DIR}")
        print("\nPlease download HCA Brain v1.0 data following the instructions above.")
        return []

    print(f"\nFound {len(h5ad_files)} H5AD file(s) in {DATA_DIR}:")
    for i, file_path in enumerate(h5ad_files, 1):
        file_size = file_path.stat().st_size / (1024**3)  # GB
        print(f"  {i}. {file_path.name} ({file_size:.2f} GB)")

    return h5ad_files


def create_metadata_summary(h5ad_path, output_path=None):
    """
    Create a summary CSV of metadata from the H5AD file

    Args:
        h5ad_path (str): Path to the H5AD file
        output_path (str): Path to save the summary CSV
    """
    if output_path is None:
        output_path = INTEGRATION_DIR / "metadata_summary.csv"

    print(f"\nCreating metadata summary from: {h5ad_path}")

    try:
        adata = ad.read_h5ad(h5ad_path)

        # Save obs metadata to CSV
        adata.obs.to_csv(output_path)
        print(f"Metadata saved to: {output_path}")
        print(f"  Shape: {adata.obs.shape[0]:,} cells × {adata.obs.shape[1]} metadata columns")

        return output_path

    except Exception as e:
        print(f"Error creating metadata summary: {e}")
        return None


def main():
    """
    Main function to run Phase 1: Data Access and Exploration
    """
    print("\n" + "=" * 80)
    print("PHASE 1: DATA ACCESS AND EXPLORATION")
    print("Human Cell Atlas Brain v1.0")
    print("=" * 80)

    # Step 1: Check for existing H5AD files
    h5ad_files = list_available_h5ad_files()

    if not h5ad_files:
        # No files found - provide download instructions
        print_download_instructions()
        print("\nAfter downloading, run this script again to explore the data.")
        return

    # Step 2: Explore each H5AD file
    print("\n" + "=" * 80)
    print("EXPLORING H5AD FILES")
    print("=" * 80)

    for h5ad_file in h5ad_files:
        results = explore_h5ad_file(h5ad_file)

        # Create metadata summary
        if results:
            create_metadata_summary(h5ad_file)

    print("\n" + "=" * 80)
    print("PHASE 1 COMPLETE")
    print("=" * 80)
    print("\nNext Steps:")
    print("1. Review the metadata summary to identify:")
    print("   - The correct column name for brain regions")
    print("   - The correct column name for cell types")
    print("2. Run Phase 2 (region mapping) to match HCA regions to your CSV regions")
    print("=" * 80)


if __name__ == "__main__":
    main()
