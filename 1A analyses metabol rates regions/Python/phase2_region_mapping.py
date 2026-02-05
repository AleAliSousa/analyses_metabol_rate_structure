"""
Phase 2: Region Mapping
========================

This script creates a mapping between HCA Brain v1.0 region names and
the user's CSV region names using flexible hierarchical matching.

Author: Claude Code
Date: 2026-01-30
"""

import pandas as pd
import numpy as np
from pathlib import Path
import re
from difflib import SequenceMatcher

# Configuration
PROJECT_DIR = Path(r"C:\Sandbox\michaelproulx\3pAgentBox\user\current\MetabolicBrain")
INTEGRATION_DIR = PROJECT_DIR / "cellxgene_integration"
CSV_PATH = PROJECT_DIR / "R" / "brain_struct_metabol_anat.csv"

# Region synonym dictionary for mapping
REGION_SYNONYMS = {
    # Amygdala variants
    'corpus_amygdaloideum': 'Corpus_amygdaloideum',
    'amygdala': 'Corpus_amygdaloideum',
    'amygdaloid': 'Corpus_amygdaloideum',

    # Hippocampus subregions
    'ca1': 'Hippocampus',
    'ca2': 'Hippocampus',
    'ca3': 'Hippocampus',
    'ca4': 'Hippocampus',
    'dentate_gyrus': 'Hippocampus',
    'dg': 'Hippocampus',
    'subiculum': 'Hippocampus',

    # Cortical regions
    'frontal_cortex': 'Frontal_lobe',
    'prefrontal_cortex': 'Frontal_lobe',
    'motor_cortex': 'Frontal_lobe',
    'parietal_cortex': 'Parietal_lobe',
    'temporal_cortex': 'Temporal_lobe',
    'occipital_cortex': 'Occipital_lobe',
    'visual_cortex': 'Occipital_lobe',
    'insula': 'Insular_lobe',

    # Striatum
    'caudate': 'Caudatum',
    'caudate_nucleus': 'Caudatum',
    'putamen': 'Putamen',
    'nucleus_accumbens': 'Nucleus_accumbens',
    'accumbens': 'Nucleus_accumbens',

    # Thalamus
    'lateral_geniculate': 'Corpus_geniculatum_laterale',
    'lgn': 'Corpus_geniculatum_laterale',
    'medial_geniculate': 'Corpus_geniculatum_mediale',
    'mgn': 'Corpus_geniculatum_mediale',

    # Midbrain
    'substantia_nigra': 'Substantia_nigra',
    'red_nucleus': 'Nucleus_ruber',
    'superior_colliculus': 'Colliculus_superior',
    'inferior_colliculus': 'Colliculus_inferior',

    # Cerebellum
    'cerebellar_cortex': 'Cerebellar_cortex',
    'dentate_nucleus': 'Nucleus_dentatus_cerebelli',

    # White matter
    'corpus_callosum': 'Centrum_semiovale',
    'internal_capsule': 'Capsula_interna',
}

# Hierarchical relationships (child -> parent)
HIERARCHICAL_MAP = {
    # Brodmann areas -> Lobes
    'BA1': 'Parietal_lobe',
    'BA2': 'Parietal_lobe',
    'BA3': 'Parietal_lobe',
    'BA4': 'Frontal_lobe',
    'BA5': 'Parietal_lobe',
    'BA6': 'Frontal_lobe',
    'BA7': 'Parietal_lobe',
    'BA8': 'Frontal_lobe',
    'BA9': 'Frontal_lobe',
    'BA10': 'Frontal_lobe',
    'BA11': 'Frontal_lobe',
    'BA17': 'Occipital_lobe',
    'BA18': 'Occipital_lobe',
    'BA19': 'Occipital_lobe',
    'BA20': 'Temporal_lobe',
    'BA21': 'Temporal_lobe',
    'BA22': 'Temporal_lobe',
    'BA37': 'Temporal_lobe',
    'BA38': 'Temporal_lobe',
    'BA39': 'Parietal_lobe',
    'BA40': 'Parietal_lobe',
    'BA41': 'Temporal_lobe',
    'BA42': 'Temporal_lobe',
    'BA44': 'Frontal_lobe',
    'BA45': 'Frontal_lobe',
    'BA46': 'Frontal_lobe',
    'BA47': 'Frontal_lobe',

    # Hippocampal subregions
    'CA1': 'Hippocampus',
    'CA2': 'Hippocampus',
    'CA3': 'Hippocampus',
    'CA4': 'Hippocampus',
    'dentate gyrus': 'Hippocampus',
    'DG': 'Hippocampus',

    # Striatal components
    'caudate': 'Striatum',
    'putamen': 'Striatum',
}


def normalize_region_name(name):
    """
    Normalize region name for comparison

    Args:
        name (str): Region name

    Returns:
        str: Normalized name
    """
    if pd.isna(name):
        return ""

    # Convert to lowercase
    name = str(name).lower()

    # Remove common prefixes/suffixes
    name = re.sub(r'^(left|right|bilateral)\s+', '', name)
    name = re.sub(r'\s+(left|right|bilateral)$', '', name)

    # Replace spaces and hyphens with underscores
    name = re.sub(r'[\s-]+', '_', name)

    # Remove special characters
    name = re.sub(r'[^\w_]', '', name)

    return name


def calculate_similarity(str1, str2):
    """
    Calculate similarity between two strings

    Args:
        str1 (str): First string
        str2 (str): Second string

    Returns:
        float: Similarity score (0-1)
    """
    return SequenceMatcher(None, str1, str2).ratio()


def map_region_flexible(hca_region, csv_regions, threshold=0.6):
    """
    Map an HCA region name to a CSV region name using flexible matching

    Args:
        hca_region (str): HCA region name
        csv_regions (list): List of CSV region names
        threshold (float): Minimum similarity threshold

    Returns:
        tuple: (matched_region, mapping_type, confidence, similarity_score)
    """
    hca_norm = normalize_region_name(hca_region)

    # 1. Check for exact match (after normalization)
    for csv_region in csv_regions:
        csv_norm = normalize_region_name(csv_region)
        if hca_norm == csv_norm:
            return (csv_region, 'direct', 'high', 1.0)

    # 2. Check synonym dictionary
    if hca_norm in REGION_SYNONYMS:
        matched = REGION_SYNONYMS[hca_norm]
        if matched in csv_regions:
            return (matched, 'synonym', 'high', 0.95)

    # 3. Check hierarchical mapping
    if hca_region in HIERARCHICAL_MAP:
        matched = HIERARCHICAL_MAP[hca_region]
        if matched in csv_regions:
            return (matched, 'hierarchical', 'high', 0.90)

    # 4. Fuzzy string matching
    best_match = None
    best_score = 0

    for csv_region in csv_regions:
        csv_norm = normalize_region_name(csv_region)

        # Calculate similarity
        similarity = calculate_similarity(hca_norm, csv_norm)

        # Also check if one is contained in the other
        containment_bonus = 0
        if hca_norm in csv_norm or csv_norm in hca_norm:
            containment_bonus = 0.2

        total_score = similarity + containment_bonus

        if total_score > best_score:
            best_score = total_score
            best_match = csv_region

    if best_score >= threshold:
        confidence = 'high' if best_score >= 0.85 else 'medium' if best_score >= 0.7 else 'low'
        return (best_match, 'fuzzy', confidence, best_score)

    # 5. No match found
    return (None, 'no_match', 'none', 0.0)


def create_region_mapping(hca_regions, csv_path):
    """
    Create mapping between HCA regions and CSV regions

    Args:
        hca_regions (list): List of HCA region names
        csv_path (str): Path to user's CSV file

    Returns:
        pd.DataFrame: Mapping dataframe
    """
    # Load CSV to get region names
    csv_df = pd.read_csv(csv_path)
    csv_regions = csv_df['Region'].tolist()

    print(f"\nCreating region mapping...")
    print(f"  HCA regions: {len(hca_regions)}")
    print(f"  CSV regions: {len(csv_regions)}")

    # Create mapping
    mappings = []

    for hca_region in hca_regions:
        matched_region, mapping_type, confidence, score = map_region_flexible(
            hca_region, csv_regions
        )

        mappings.append({
            'HCA_region_name': hca_region,
            'HCA_region_normalized': normalize_region_name(hca_region),
            'CSV_region_name': matched_region,
            'mapping_type': mapping_type,
            'confidence': confidence,
            'similarity_score': score
        })

    mapping_df = pd.DataFrame(mappings)

    # Summary statistics
    print("\nMapping Summary:")
    print(f"  Total HCA regions: {len(mapping_df)}")
    print(f"  Mapped regions: {mapping_df['CSV_region_name'].notna().sum()}")
    print(f"  Unmapped regions: {mapping_df['CSV_region_name'].isna().sum()}")

    print("\nBy Mapping Type:")
    print(mapping_df['mapping_type'].value_counts())

    print("\nBy Confidence:")
    print(mapping_df['confidence'].value_counts())

    return mapping_df


def review_mapping(mapping_df, show_all=False):
    """
    Review the region mapping and highlight potential issues

    Args:
        mapping_df (pd.DataFrame): Mapping dataframe
        show_all (bool): Show all mappings or just uncertain ones
    """
    print("\n" + "=" * 80)
    print("REGION MAPPING REVIEW")
    print("=" * 80)

    if show_all:
        print("\nAll Mappings:")
        for idx, row in mapping_df.iterrows():
            print(f"\n  {row['HCA_region_name']}")
            print(f"    -> {row['CSV_region_name']}")
            print(f"    Type: {row['mapping_type']}, Confidence: {row['confidence']}, Score: {row['similarity_score']:.2f}")
    else:
        # Show uncertain mappings
        uncertain = mapping_df[
            (mapping_df['confidence'].isin(['low', 'medium'])) |
            (mapping_df['CSV_region_name'].isna())
        ]

        if len(uncertain) > 0:
            print(f"\nUncertain or Unmapped Regions ({len(uncertain)}):")
            for idx, row in uncertain.iterrows():
                print(f"\n  {row['HCA_region_name']}")
                if pd.notna(row['CSV_region_name']):
                    print(f"    -> {row['CSV_region_name']} (confidence: {row['confidence']}, score: {row['similarity_score']:.2f})")
                else:
                    print(f"    -> NO MATCH FOUND")
        else:
            print("\nAll regions mapped with high confidence!")

        # Show high-confidence mappings count
        high_conf = mapping_df[mapping_df['confidence'] == 'high']
        print(f"\nHigh Confidence Mappings: {len(high_conf)}")


def save_mapping(mapping_df, output_path=None):
    """
    Save the region mapping to CSV

    Args:
        mapping_df (pd.DataFrame): Mapping dataframe
        output_path (str): Output path
    """
    if output_path is None:
        output_path = INTEGRATION_DIR / "region_mapping.csv"

    mapping_df.to_csv(output_path, index=False)
    print(f"\nRegion mapping saved to: {output_path}")


def main(hca_regions=None):
    """
    Main function for Phase 2: Region Mapping

    Args:
        hca_regions (list): List of HCA region names (if None, will prompt for input)
    """
    print("\n" + "=" * 80)
    print("PHASE 2: REGION MAPPING")
    print("=" * 80)

    # If no HCA regions provided, this is just a test/example
    if hca_regions is None:
        print("\nNo HCA regions provided.")
        print("This script should be run after Phase 1 with the actual HCA region names.")
        print("\nExample usage:")
        print("  from phase2_region_mapping import create_region_mapping")
        print("  hca_regions = ['hippocampus', 'CA1', 'frontal cortex', ...]")
        print("  mapping_df = create_region_mapping(hca_regions, csv_path)")
        return None

    # Create mapping
    mapping_df = create_region_mapping(hca_regions, CSV_PATH)

    # Review mapping
    review_mapping(mapping_df, show_all=False)

    # Save mapping
    save_mapping(mapping_df)

    print("\n" + "=" * 80)
    print("PHASE 2 COMPLETE")
    print("=" * 80)

    return mapping_df


if __name__ == "__main__":
    main()
