"""
Phase 3: Cell Count Extraction
================================

This script extracts cell counts from HCA Brain v1.0 data and aggregates them
by brain region and cell type.

Author: Claude Code
Date: 2026-01-30
"""

import pandas as pd
import numpy as np
import anndata as ad
from pathlib import Path
from collections import defaultdict

# Configuration
PROJECT_DIR = Path(r"C:\Sandbox\michaelproulx\3pAgentBox\user\current\MetabolicBrain")
INTEGRATION_DIR = PROJECT_DIR / "cellxgene_integration"

# Cell type classification
# Map detailed cell type annotations to broad categories
CELL_TYPE_CATEGORIES = {
    'neuron': [
        'neuron', 'neuronal', 'excitatory', 'inhibitory', 'glutamatergic',
        'gabaergic', 'dopaminergic', 'cholinergic', 'serotonergic',
        'pyramidal', 'granule', 'purkinje', 'interneuron'
    ],
    'astrocyte': [
        'astrocyte', 'astro', 'astroglial'
    ],
    'oligodendrocyte': [
        'oligodendrocyte', 'oligo', 'opc', 'oligodendrocyte precursor'
    ],
    'microglia': [
        'microglia', 'microglial'
    ],
    'other_glia': [
        'ependymal', 'choroid plexus', 'schwann', 'radial glia'
    ],
    'vascular': [
        'endothelial', 'pericyte', 'vascular', 'vsmc'
    ],
}


def classify_cell_type(cell_type_str, categories=CELL_TYPE_CATEGORIES):
    """
    Classify a cell type annotation into broad categories

    Args:
        cell_type_str (str): Cell type annotation
        categories (dict): Category mapping

    Returns:
        str: Classified category or 'unknown'
    """
    if pd.isna(cell_type_str):
        return 'unknown'

    cell_type_lower = str(cell_type_str).lower()

    for category, keywords in categories.items():
        for keyword in keywords:
            if keyword in cell_type_lower:
                return category

    return 'unknown'


def extract_cell_counts(h5ad_path, region_mapping_path, region_col=None, celltype_col=None):
    """
    Extract cell counts from H5AD file aggregated by region and cell type

    Args:
        h5ad_path (str): Path to H5AD file
        region_mapping_path (str): Path to region mapping CSV
        region_col (str): Column name for regions in H5AD metadata
        celltype_col (str): Column name for cell types in H5AD metadata

    Returns:
        pd.DataFrame: Cell count dataframe
    """
    print("\n" + "=" * 80)
    print("EXTRACTING CELL COUNTS")
    print("=" * 80)

    # Load data
    print(f"\nLoading H5AD file: {h5ad_path}")
    adata = ad.read_h5ad(h5ad_path)
    print(f"  Loaded {adata.n_obs:,} cells")

    # Load region mapping
    print(f"\nLoading region mapping: {region_mapping_path}")
    region_mapping = pd.read_csv(region_mapping_path)

    # Create mapping dictionary (HCA region -> CSV region)
    region_map = dict(zip(
        region_mapping['HCA_region_name'],
        region_mapping['CSV_region_name']
    ))

    # Auto-detect region and cell type columns if not provided
    if region_col is None:
        # Look for region-related columns
        potential_cols = [col for col in adata.obs.columns
                         if any(keyword in col.lower() for keyword in
                               ['region', 'tissue', 'area', 'structure'])]
        if potential_cols:
            region_col = potential_cols[0]
            print(f"\nAuto-detected region column: {region_col}")
        else:
            raise ValueError("Could not auto-detect region column. Please specify region_col parameter.")

    if celltype_col is None:
        # Look for cell type-related columns
        potential_cols = [col for col in adata.obs.columns
                         if any(keyword in col.lower() for keyword in
                               ['cell_type', 'celltype', 'cluster', 'annotation'])]
        if potential_cols:
            celltype_col = potential_cols[0]
            print(f"Auto-detected cell type column: {celltype_col}")
        else:
            raise ValueError("Could not auto-detect cell type column. Please specify celltype_col parameter.")

    # Add metadata columns
    print("\nProcessing cell metadata...")

    # Map HCA regions to CSV regions
    adata.obs['CSV_region'] = adata.obs[region_col].map(region_map)

    # Classify cell types
    adata.obs['cell_category'] = adata.obs[celltype_col].apply(classify_cell_type)

    # Count cells that were successfully mapped
    mapped_cells = adata.obs['CSV_region'].notna().sum()
    print(f"  Cells successfully mapped to CSV regions: {mapped_cells:,} ({mapped_cells/adata.n_obs*100:.1f}%)")

    # Count by cell category
    print("\nCell type classification:")
    category_counts = adata.obs['cell_category'].value_counts()
    for category, count in category_counts.items():
        print(f"  {category}: {count:,} cells ({count/adata.n_obs*100:.1f}%)")

    # Aggregate counts by region and cell type
    print("\nAggregating cell counts by region and cell type...")

    # Initialize results structure
    results = defaultdict(lambda: {
        'Neuron_N': 0,
        'Astro_N': 0,
        'Oligo_N': 0,
        'Microglia_N': 0,
        'Glia_N': 0,  # Total glia
        'Other_N': 0,
        'Total_N': 0,
        'cell_type_details': defaultdict(int)
    })

    # Group by region and cell category
    grouped = adata.obs.groupby(['CSV_region', 'cell_category']).size()

    for (region, category), count in grouped.items():
        if pd.isna(region):
            continue

        results[region]['Total_N'] += count
        results[region]['cell_type_details'][category] += count

        if category == 'neuron':
            results[region]['Neuron_N'] += count
        elif category == 'astrocyte':
            results[region]['Astro_N'] += count
            results[region]['Glia_N'] += count
        elif category == 'oligodendrocyte':
            results[region]['Oligo_N'] += count
            results[region]['Glia_N'] += count
        elif category == 'microglia':
            results[region]['Microglia_N'] += count
            results[region]['Glia_N'] += count
        elif category == 'other_glia':
            results[region]['Glia_N'] += count
        else:
            results[region]['Other_N'] += count

    # Convert to DataFrame
    count_data = []
    for region, counts in results.items():
        row = {
            'Region': region,
            'Neuron_N': counts['Neuron_N'],
            'Glia_N': counts['Glia_N'],
            'Astro_N': counts['Astro_N'],
            'Oligo_N': counts['Oligo_N'],
            'Microglia_N': counts['Microglia_N'],
            'Other_N': counts['Other_N'],
            'Total_N': counts['Total_N'],
        }
        count_data.append(row)

    counts_df = pd.DataFrame(count_data)
    counts_df = counts_df.sort_values('Total_N', ascending=False)

    # Summary
    print("\n" + "=" * 80)
    print(f"EXTRACTION SUMMARY")
    print("=" * 80)
    print(f"Regions with data: {len(counts_df)}")
    print(f"Total cells extracted: {counts_df['Total_N'].sum():,}")
    print(f"Total neurons: {counts_df['Neuron_N'].sum():,}")
    print(f"Total glia: {counts_df['Glia_N'].sum():,}")
    print(f"  - Astrocytes: {counts_df['Astro_N'].sum():,}")
    print(f"  - Oligodendrocytes: {counts_df['Oligo_N'].sum():,}")
    print(f"  - Microglia: {counts_df['Microglia_N'].sum():,}")

    print("\nTop 10 regions by cell count:")
    print(counts_df[['Region', 'Total_N', 'Neuron_N', 'Glia_N']].head(10).to_string(index=False))

    return counts_df


def calculate_standard_deviations(h5ad_path, region_mapping_path, counts_df,
                                 region_col=None, celltype_col=None, donor_col=None):
    """
    Calculate standard deviations across donors/samples if available

    Args:
        h5ad_path (str): Path to H5AD file
        region_mapping_path (str): Path to region mapping CSV
        counts_df (pd.DataFrame): Counts dataframe from extract_cell_counts
        region_col (str): Column name for regions
        celltype_col (str): Column name for cell types
        donor_col (str): Column name for donor/sample ID

    Returns:
        pd.DataFrame: Counts with standard deviations
    """
    print("\n" + "=" * 80)
    print("CALCULATING STANDARD DEVIATIONS")
    print("=" * 80)

    # Load data
    adata = ad.read_h5ad(h5ad_path)
    region_mapping = pd.read_csv(region_mapping_path)
    region_map = dict(zip(region_mapping['HCA_region_name'], region_mapping['CSV_region_name']))

    # Auto-detect donor column if not provided
    if donor_col is None:
        potential_cols = [col for col in adata.obs.columns
                         if any(keyword in col.lower() for keyword in
                               ['donor', 'sample', 'individual', 'subject'])]
        if potential_cols:
            donor_col = potential_cols[0]
            print(f"\nAuto-detected donor column: {donor_col}")
        else:
            print("\nNo donor column found. Cannot calculate standard deviations.")
            print("Returning counts without SDs.")
            return counts_df

    # Count unique donors
    n_donors = adata.obs[donor_col].nunique()
    print(f"Number of unique donors/samples: {n_donors}")

    if n_donors < 2:
        print("\nInsufficient donors for standard deviation calculation.")
        return counts_df

    # Add metadata
    if region_col is None or celltype_col is None:
        # Use same auto-detection as in extract_cell_counts
        if region_col is None:
            potential = [col for col in adata.obs.columns if 'region' in col.lower()]
            region_col = potential[0] if potential else None
        if celltype_col is None:
            potential = [col for col in adata.obs.columns if 'cell_type' in col.lower()]
            celltype_col = potential[0] if potential else None

    adata.obs['CSV_region'] = adata.obs[region_col].map(region_map)
    adata.obs['cell_category'] = adata.obs[celltype_col].apply(classify_cell_type)

    # Calculate counts per donor
    donor_counts = []

    for donor in adata.obs[donor_col].unique():
        donor_data = adata.obs[adata.obs[donor_col] == donor]

        for region in counts_df['Region'].unique():
            region_data = donor_data[donor_data['CSV_region'] == region]

            neuron_n = (region_data['cell_category'] == 'neuron').sum()
            astro_n = (region_data['cell_category'] == 'astrocyte').sum()
            oligo_n = (region_data['cell_category'] == 'oligodendrocyte').sum()
            micro_n = (region_data['cell_category'] == 'microglia').sum()
            glia_n = astro_n + oligo_n + micro_n

            donor_counts.append({
                'donor': donor,
                'Region': region,
                'Neuron_N': neuron_n,
                'Astro_N': astro_n,
                'Oligo_N': oligo_n,
                'Microglia_N': micro_n,
                'Glia_N': glia_n,
            })

    donor_df = pd.DataFrame(donor_counts)

    # Calculate SDs
    sd_df = donor_df.groupby('Region').agg({
        'Neuron_N': 'std',
        'Astro_N': 'std',
        'Oligo_N': 'std',
        'Microglia_N': 'std',
        'Glia_N': 'std',
    }).reset_index()

    sd_df.columns = ['Region', 'NeuronSD', 'AstroSD', 'OligoSD', 'MicroSD', 'GliaSD']

    # Merge with counts
    counts_with_sd = counts_df.merge(sd_df, on='Region', how='left')
    counts_with_sd['n_donors'] = n_donors

    print("\nStandard deviations calculated successfully.")

    return counts_with_sd


def save_extracted_counts(counts_df, output_path=None):
    """
    Save extracted cell counts to CSV

    Args:
        counts_df (pd.DataFrame): Counts dataframe
        output_path (str): Output path
    """
    if output_path is None:
        output_path = INTEGRATION_DIR / "extracted_cell_counts.csv"

    counts_df.to_csv(output_path, index=False)
    print(f"\nExtracted counts saved to: {output_path}")


def main():
    """
    Main function for Phase 3: Cell Count Extraction
    """
    print("\n" + "=" * 80)
    print("PHASE 3: CELL COUNT EXTRACTION")
    print("=" * 80)

    print("\nThis script should be run after Phase 1 and Phase 2 are complete.")
    print("\nExample usage:")
    print("  from phase3_cell_extraction import extract_cell_counts, calculate_standard_deviations")
    print("  ")
    print("  # Extract counts")
    print("  counts_df = extract_cell_counts(")
    print("      h5ad_path='data/brain_atlas.h5ad',")
    print("      region_mapping_path='region_mapping.csv'")
    print("  )")
    print("  ")
    print("  # Add standard deviations")
    print("  counts_with_sd = calculate_standard_deviations(")
    print("      h5ad_path='data/brain_atlas.h5ad',")
    print("      region_mapping_path='region_mapping.csv',")
    print("      counts_df=counts_df")
    print("  )")
    print("  ")
    print("  # Save results")
    print("  save_extracted_counts(counts_with_sd)")


if __name__ == "__main__":
    main()
