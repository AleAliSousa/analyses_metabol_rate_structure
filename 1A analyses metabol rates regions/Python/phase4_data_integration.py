"""
Phase 4: Data Integration
===========================

This script integrates extracted cell counts from HCA Brain v1.0 into the
user's existing brain_struct_metabol_anat.csv file.

Author: Claude Code
Date: 2026-01-30
"""

import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime

# Configuration
PROJECT_DIR = Path(r"C:\Sandbox\michaelproulx\3pAgentBox\user\current\MetabolicBrain")
INTEGRATION_DIR = PROJECT_DIR / "cellxgene_integration"
CSV_PATH = PROJECT_DIR / "R" / "brain_struct_metabol_anat.csv"


def load_data(original_csv_path, extracted_counts_path):
    """
    Load original CSV and extracted counts

    Args:
        original_csv_path (str): Path to original CSV
        extracted_counts_path (str): Path to extracted counts CSV

    Returns:
        tuple: (original_df, extracted_df)
    """
    print("\n" + "=" * 80)
    print("LOADING DATA")
    print("=" * 80)

    original_df = pd.read_csv(original_csv_path)
    print(f"\nOriginal CSV: {original_csv_path}")
    print(f"  Shape: {original_df.shape[0]} regions × {original_df.shape[1]} columns")

    extracted_df = pd.read_csv(extracted_counts_path)
    print(f"\nExtracted counts: {extracted_counts_path}")
    print(f"  Shape: {extracted_df.shape[0]} regions × {extracted_df.shape[1]} columns")

    return original_df, extracted_df


def analyze_missing_data(original_df):
    """
    Analyze missing data in the original CSV

    Args:
        original_df (pd.DataFrame): Original dataframe

    Returns:
        dict: Missing data statistics
    """
    print("\n" + "=" * 80)
    print("ANALYZING MISSING DATA")
    print("=" * 80)

    cell_count_cols = [
        'Neuron_N', 'NeuronSD',
        'Glia_N', 'GliaSD',
        'Astro_N', 'AstroSD',
        'Oligo_N', 'OligoSD',
        'Microglia_N', 'MicroSD'
    ]

    stats = {}

    print("\nMissing data by column:")
    for col in cell_count_cols:
        if col in original_df.columns:
            missing = original_df[col].isna().sum()
            total = len(original_df)
            pct = missing / total * 100
            stats[col] = {'missing': missing, 'total': total, 'pct': pct}
            print(f"  {col}: {missing}/{total} missing ({pct:.1f}%)")
        else:
            print(f"  {col}: Column not found")

    # Count regions with complete cell count data
    complete_rows = original_df[cell_count_cols].notna().all(axis=1).sum()
    print(f"\nRegions with complete cell count data: {complete_rows}/{len(original_df)}")

    return stats


def integrate_data(original_df, extracted_df, overwrite=False, fill_only_missing=True):
    """
    Integrate extracted data into original dataframe

    Args:
        original_df (pd.DataFrame): Original dataframe
        extracted_df (pd.DataFrame): Extracted counts dataframe
        overwrite (bool): Whether to overwrite existing values
        fill_only_missing (bool): Only fill missing values (if overwrite=False)

    Returns:
        tuple: (integrated_df, integration_log)
    """
    print("\n" + "=" * 80)
    print("INTEGRATING DATA")
    print("=" * 80)

    print(f"\nIntegration mode:")
    print(f"  Overwrite existing values: {overwrite}")
    print(f"  Fill only missing values: {fill_only_missing}")

    # Create a copy to avoid modifying original
    integrated_df = original_df.copy()

    # Columns to integrate
    count_columns = {
        'Neuron_N': 'Neuron_N',
        'NeuronSD': 'NeuronSD',
        'Glia_N': 'Glia_N',
        'GliaSD': 'GliaSD',
        'Astro_N': 'Astro_N',
        'AstroSD': 'AstroSD',
        'Oligo_N': 'Oligo_N',
        'OligoSD': 'OligoSD',
        'Microglia_N': 'Microglia_N',
        'MicroSD': 'MicroSD',
    }

    # Integration log
    log = []
    log.append(f"Integration started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    log.append(f"Overwrite mode: {overwrite}")
    log.append(f"Fill only missing: {fill_only_missing}")
    log.append("")

    # Track changes
    changes = {
        'filled': 0,
        'overwritten': 0,
        'unchanged': 0,
        'not_available': 0,
    }

    # Merge data
    for idx, row in integrated_df.iterrows():
        region = row['Region']

        # Find matching region in extracted data
        extracted_row = extracted_df[extracted_df['Region'] == region]

        if extracted_row.empty:
            log.append(f"Region '{region}': No extracted data available")
            continue

        extracted_row = extracted_row.iloc[0]

        # Process each column
        for orig_col, extr_col in count_columns.items():
            if orig_col not in integrated_df.columns:
                # Add column if it doesn't exist
                integrated_df[orig_col] = np.nan

            original_val = integrated_df.loc[idx, orig_col]
            extracted_val = extracted_row.get(extr_col, np.nan)

            # Determine action
            if pd.isna(extracted_val):
                # No extracted data available for this column
                changes['not_available'] += 1
                continue

            if pd.isna(original_val) or original_val == '':
                # Fill missing value
                integrated_df.loc[idx, orig_col] = extracted_val
                changes['filled'] += 1
                log.append(f"  {region} - {orig_col}: Filled {extracted_val}")

            elif overwrite:
                # Overwrite existing value
                old_val = original_val
                integrated_df.loc[idx, orig_col] = extracted_val
                changes['overwritten'] += 1
                log.append(f"  {region} - {orig_col}: Changed {old_val} -> {extracted_val}")

            else:
                # Keep original value
                changes['unchanged'] += 1

    # Summary
    print("\n" + "=" * 80)
    print("INTEGRATION SUMMARY")
    print("=" * 80)
    print(f"  Values filled (was missing): {changes['filled']}")
    print(f"  Values overwritten: {changes['overwritten']}")
    print(f"  Values unchanged: {changes['unchanged']}")
    print(f"  Values not available: {changes['not_available']}")

    log.append("")
    log.append("=" * 80)
    log.append("SUMMARY")
    log.append("=" * 80)
    log.append(f"Values filled: {changes['filled']}")
    log.append(f"Values overwritten: {changes['overwritten']}")
    log.append(f"Values unchanged: {changes['unchanged']}")
    log.append(f"Values not available: {changes['not_available']}")

    return integrated_df, log


def calculate_densities(integrated_df):
    """
    Calculate cell densities from counts and volumes

    Args:
        integrated_df (pd.DataFrame): Integrated dataframe

    Returns:
        pd.DataFrame: Dataframe with densities calculated
    """
    print("\n" + "=" * 80)
    print("CALCULATING DENSITIES")
    print("=" * 80)

    df = integrated_df.copy()

    # Check if Volume column exists
    if 'Volume' not in df.columns:
        print("Volume column not found. Cannot calculate densities.")
        return df

    # Calculate densities where volume is available
    density_calculations = {
        'NeurDensity': ('Neuron_N', 'Volume'),
        'GliaDensity': ('Glia_N', 'Volume'),
        'AstroDensity': ('Astro_N', 'Volume'),
        'OligoDensity': ('Oligo_N', 'Volume'),
        'MicroDensity': ('Microglia_N', 'Volume'),
    }

    calculated = 0

    for density_col, (count_col, volume_col) in density_calculations.items():
        if density_col not in df.columns:
            df[density_col] = np.nan

        # Calculate where both count and volume are available
        mask = df[count_col].notna() & df[volume_col].notna() & (df[volume_col] > 0)

        if mask.sum() > 0:
            df.loc[mask, density_col] = df.loc[mask, count_col] / df.loc[mask, volume_col]
            calculated += mask.sum()
            print(f"  {density_col}: Calculated for {mask.sum()} regions")

    # Calculate Glia/Neuron density ratio
    if 'GliaNeurDensityRatio' not in df.columns:
        df['GliaNeurDensityRatio'] = np.nan

    mask = df['GliaDensity'].notna() & df['NeurDensity'].notna() & (df['NeurDensity'] > 0)
    if mask.sum() > 0:
        df.loc[mask, 'GliaNeurDensityRatio'] = df.loc[mask, 'GliaDensity'] / df.loc[mask, 'NeurDensity']
        print(f"  GliaNeurDensityRatio: Calculated for {mask.sum()} regions")

    print(f"\nTotal density calculations: {calculated}")

    return df


def compare_before_after(original_df, integrated_df):
    """
    Compare data completeness before and after integration

    Args:
        original_df (pd.DataFrame): Original dataframe
        integrated_df (pd.DataFrame): Integrated dataframe

    Returns:
        pd.DataFrame: Comparison dataframe
    """
    print("\n" + "=" * 80)
    print("BEFORE/AFTER COMPARISON")
    print("=" * 80)

    cell_count_cols = [
        'Neuron_N', 'Glia_N', 'Astro_N', 'Oligo_N', 'Microglia_N',
        'NeurDensity', 'GliaDensity', 'AstroDensity', 'OligoDensity', 'MicroDensity'
    ]

    comparison = []

    for col in cell_count_cols:
        if col not in original_df.columns:
            original_complete = 0
        else:
            original_complete = original_df[col].notna().sum()

        if col not in integrated_df.columns:
            integrated_complete = 0
        else:
            integrated_complete = integrated_df[col].notna().sum()

        improvement = integrated_complete - original_complete

        comparison.append({
            'Column': col,
            'Before': original_complete,
            'After': integrated_complete,
            'Improvement': improvement,
            'Pct_Before': f"{original_complete/len(original_df)*100:.1f}%",
            'Pct_After': f"{integrated_complete/len(integrated_df)*100:.1f}%",
        })

    comparison_df = pd.DataFrame(comparison)

    print("\nData Completeness Comparison:")
    print(comparison_df.to_string(index=False))

    # Overall completeness
    total_before = comparison_df['Before'].sum()
    total_after = comparison_df['After'].sum()
    total_possible = len(cell_count_cols) * len(original_df)

    print(f"\nOverall Completeness:")
    print(f"  Before: {total_before}/{total_possible} cells filled ({total_before/total_possible*100:.1f}%)")
    print(f"  After:  {total_after}/{total_possible} cells filled ({total_after/total_possible*100:.1f}%)")
    print(f"  Improvement: +{total_after - total_before} cells")

    return comparison_df


def save_integrated_data(integrated_df, log, output_csv_path=None, output_log_path=None):
    """
    Save integrated dataframe and log

    Args:
        integrated_df (pd.DataFrame): Integrated dataframe
        log (list): Integration log
        output_csv_path (str): Output CSV path
        output_log_path (str): Output log path
    """
    if output_csv_path is None:
        output_csv_path = PROJECT_DIR / "R" / "brain_struct_metabol_anat_INTEGRATED.csv"

    if output_log_path is None:
        output_log_path = INTEGRATION_DIR / "integration_log.txt"

    # Save CSV
    integrated_df.to_csv(output_csv_path, index=False)
    print(f"\nIntegrated data saved to: {output_csv_path}")

    # Save log
    with open(output_log_path, 'w') as f:
        f.write('\n'.join(log))
    print(f"Integration log saved to: {output_log_path}")


def main():
    """
    Main function for Phase 4: Data Integration
    """
    print("\n" + "=" * 80)
    print("PHASE 4: DATA INTEGRATION")
    print("=" * 80)

    # Check for required files
    csv_path = CSV_PATH
    extracted_path = INTEGRATION_DIR / "extracted_cell_counts.csv"

    if not csv_path.exists():
        print(f"\nERROR: Original CSV not found at {csv_path}")
        return

    if not extracted_path.exists():
        print(f"\nERROR: Extracted counts not found at {extracted_path}")
        print("Please run Phase 3 first to extract cell counts.")
        return

    # Load data
    original_df, extracted_df = load_data(csv_path, extracted_path)

    # Analyze missing data
    analyze_missing_data(original_df)

    # Integrate data
    integrated_df, log = integrate_data(
        original_df,
        extracted_df,
        overwrite=False,
        fill_only_missing=True
    )

    # Calculate densities
    integrated_df = calculate_densities(integrated_df)

    # Compare before/after
    comparison_df = compare_before_after(original_df, integrated_df)

    # Save results
    save_integrated_data(integrated_df, log)

    print("\n" + "=" * 80)
    print("PHASE 4 COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    main()
