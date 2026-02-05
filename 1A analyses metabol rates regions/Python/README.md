# CELLxGENE Data Integration Pipeline

## Overview

This pipeline integrates cell count data from the **Human Cell Atlas (HCA) Brain v1.0** into your metabolic brain analysis dataset (`brain_struct_metabol_anat.csv`). The goal is to fill missing cell count data (neurons, glia, astrocytes, oligodendrocytes, microglia) to improve statistical power in your analysis of brain region metabolic rates.

## Data Source: HCA Brain v1.0

**Recommended Dataset**: Human Cell Atlas Brain v1.0

### Key Characteristics:
- **Scale**: 3+ million nuclei from adult human postmortem brains
- **Coverage**: ~100 dissections across major brain divisions
- **Donors**: 3 adult human donors
- **Cell Types**: 2.5 million neurons + 888,300 non-neuronal cells
- **Brain Regions**: 88 tissues including:
  - Cerebral cortex (multiple gyri and Brodmann areas)
  - Hippocampus (CA1-CA4, dentate gyrus)
  - Amygdala, Thalamus, Hypothalamus
  - Cerebellum, Midbrain, Pons, Medulla
  - Basal nuclei/ganglia (striatum, pallidum)
  - Spinal cord

### Access:
- **Primary Source**: [HCA Data Portal - Brain v1.0](https://data.humancellatlas.org/hca-bio-networks/nervous-system/atlases/brain-v1-0)
- **Alternative**: [CELLxGENE Data Portal](https://cellxgene.cziscience.com/)
- **Format**: H5AD (AnnData format for Python)
- **License**: CC BY 4.0

## Pipeline Phases

### Phase 1: Data Access and Exploration
**Script**: `phase1_data_access.py`

**Purpose**: Download and explore HCA Brain v1.0 data

**Steps**:
1. Download H5AD file(s) from HCA Data Portal or CELLxGENE
2. Place in `cellxgene_integration/data/` directory
3. Run Phase 1 to explore data structure
4. Identify region and cell type annotation columns

**Outputs**:
- Metadata summary CSV
- Console report of data structure

### Phase 2: Region Mapping
**Script**: `phase2_region_mapping.py`

**Purpose**: Map HCA region names to your CSV region names

**Mapping Strategy**:
1. **Direct matching**: Exact matches (e.g., Hippocampus → Hippocampus)
2. **Synonym mapping**: Known anatomical synonyms (e.g., Corpus amygdaloideum → Amygdala)
3. **Hierarchical aggregation**: Child regions → parent regions (e.g., CA1, CA2, CA3, CA4 → Hippocampus)
4. **Fuzzy matching**: String similarity for near-matches

**Outputs**:
- `region_mapping.csv`: Mapping table with confidence scores
- Requires manual review before Phase 3

### Phase 3: Cell Count Extraction
**Script**: `phase3_cell_extraction.py`

**Purpose**: Extract and aggregate cell counts by region and cell type

**Cell Type Classification**:
- **Neurons**: All neuronal cell types
- **Astrocytes**: Astrocytes and astroglia
- **Oligodendrocytes**: Oligodendrocytes and OPCs
- **Microglia**: Microglial cells
- **Glia (total)**: Sum of astrocytes, oligodendrocytes, microglia, and other glia

**Outputs**:
- `extracted_cell_counts.csv`: Cell counts by region with standard deviations

### Phase 4: Data Integration
**Script**: `phase4_data_integration.py`

**Purpose**: Merge extracted counts into your original CSV

**Integration Process**:
1. Load original CSV and extracted counts
2. Analyze missing data patterns
3. Fill missing values (default: only fill empty cells, don't overwrite)
4. Calculate cell densities (if volume data available)
5. Generate before/after comparison

**Outputs**:
- `brain_struct_metabol_anat_INTEGRATED.csv`: Your CSV with filled data
- `integration_log.txt`: Detailed log of all changes

## Installation

### Prerequisites
- Python 3.8+
- pip package manager

### Install Required Packages

```bash
cd C:\Sandbox\michaelproulx\3pAgentBox\user\current\MetabolicBrain\cellxgene_integration
python -m pip install pandas numpy anndata scanpy h5py
```

## Usage

### Option 1: Run All Phases Automatically

```bash
python run_integration.py
```

This will run all phases sequentially with checkpoints for manual review.

### Option 2: Run Individual Phases

```bash
# Phase 1: Data exploration
python run_integration.py --phase 1

# Phase 2: Region mapping (requires H5AD file)
python run_integration.py --phase 2 --h5ad "data/brain_atlas.h5ad"

# Phase 3: Cell count extraction
python run_integration.py --phase 3 --h5ad "data/brain_atlas.h5ad"

# Phase 4: Data integration
python run_integration.py --phase 4
```

### Option 3: Use Scripts Directly in Python

```python
# Phase 1
import phase1_data_access as phase1
phase1.main()

# Phase 2
import phase2_region_mapping as phase2
hca_regions = ['hippocampus', 'frontal cortex', ...]  # From Phase 1
mapping_df = phase2.create_region_mapping(hca_regions, csv_path)
phase2.save_mapping(mapping_df)

# Phase 3
import phase3_cell_extraction as phase3
counts_df = phase3.extract_cell_counts(
    h5ad_path='data/brain_atlas.h5ad',
    region_mapping_path='region_mapping.csv'
)
counts_with_sd = phase3.calculate_standard_deviations(
    h5ad_path='data/brain_atlas.h5ad',
    region_mapping_path='region_mapping.csv',
    counts_df=counts_df
)
phase3.save_extracted_counts(counts_with_sd)

# Phase 4
import phase4_data_integration as phase4
phase4.main()
```

## Data Download Instructions

### Method 1: HCA Data Portal

1. Visit: https://data.humancellatlas.org/hca-bio-networks/nervous-system/atlases/brain-v1-0
2. Click on "Access Data" or "Download"
3. Select datasets covering your brain regions of interest
4. Download in **H5AD format**
5. Save to: `cellxgene_integration/data/`

### Method 2: CELLxGENE

1. Visit: https://cellxgene.cziscience.com/
2. Search for: "Human Brain Cell Atlas" or "brain v1.0"
3. Filter by:
   - Organism: Homo sapiens
   - Tissue: Brain
   - Assay: snRNA-seq or scRNA-seq
4. Click "Download" and select **H5AD format**
5. Save to: `cellxgene_integration/data/`

## Directory Structure

```
cellxgene_integration/
├── data/                          # Place H5AD files here
│   └── (your_downloaded_file.h5ad)
├── phase1_data_access.py          # Phase 1 script
├── phase2_region_mapping.py       # Phase 2 script
├── phase3_cell_extraction.py      # Phase 3 script
├── phase4_data_integration.py     # Phase 4 script
├── run_integration.py             # Master orchestration script
├── README.md                      # This file
├── region_mapping.csv             # Output from Phase 2
├── extracted_cell_counts.csv      # Output from Phase 3
├── metadata_summary.csv           # Output from Phase 1
└── integration_log.txt            # Output from Phase 4
```

## Output Files

### Main Output
**`brain_struct_metabol_anat_INTEGRATED.csv`**
- Location: `MetabolicBrain/R/`
- Your original CSV with cell count data filled in
- Ready for re-analysis with `brain_struct_costs.Rmd`

### Supporting Files
- **`region_mapping.csv`**: HCA → CSV region name mappings
- **`extracted_cell_counts.csv`**: Raw extracted cell counts by region
- **`integration_log.txt`**: Detailed log of integration changes
- **`metadata_summary.csv`**: Full metadata from H5AD file

## Quality Control

### Manual Review Checkpoints

1. **After Phase 1**:
   - Verify H5AD file contains expected brain regions
   - Identify correct column names for regions and cell types

2. **After Phase 2**:
   - Review `region_mapping.csv`
   - Check for incorrectly mapped regions
   - Manually edit mapping if needed
   - Focus on low-confidence and unmapped regions

3. **After Phase 3**:
   - Review `extracted_cell_counts.csv`
   - Check for biological plausibility of counts
   - Verify cell type classifications

4. **After Phase 4**:
   - Compare before/after data completeness
   - Review `integration_log.txt` for changes
   - Verify no unexpected overwrites

## Expected Outcomes

### Optimistic Scenario:
- Fill 60-80% of missing cell count data
- Improve from ~7 to ~30+ complete observations
- Dramatically increase statistical power

### Realistic Scenario:
- Fill 40-60% of missing data
- Improve from ~7 to ~20+ complete observations
- Moderate increase in statistical power
- Some regions may not have exact matches

### Success Metrics:
- Reduce "observations deleted due to missingness" from 39 to <20
- Improved R² and p-values in regression models
- Maintained biological plausibility of values

## Re-running Analysis

After integration, re-run your analysis:

```r
# In R
library(readr)

# Load integrated data
brain_data <- read_csv("brain_struct_metabol_anat_INTEGRATED.csv")

# Run your existing analysis from brain_struct_costs.Rmd
# Compare results before/after integration
```

## Troubleshooting

### Issue: "No H5AD files found"
**Solution**: Download HCA Brain v1.0 data and place in `cellxgene_integration/data/`

### Issue: "Could not auto-detect region column"
**Solution**: Manually specify column name using `region_col` parameter

### Issue: "Many regions unmapped"
**Solution**:
1. Review HCA region names in metadata summary
2. Add custom mappings to `REGION_SYNONYMS` in `phase2_region_mapping.py`
3. Re-run Phase 2

### Issue: "Cell type classification looks wrong"
**Solution**:
1. Review HCA cell type annotations
2. Update `CELL_TYPE_CATEGORIES` in `phase3_cell_extraction.py`
3. Re-run Phase 3

### Issue: "Integration overwrote existing data"
**Solution**:
1. Check `integration_log.txt` to see what changed
2. Restore original CSV from backup
3. Re-run Phase 4 with `overwrite=False` (default)

## Data Sources and References

### Human Cell Atlas Brain v1.0
- Portal: https://data.humancellatlas.org/hca-bio-networks/nervous-system/atlases/brain-v1-0
- Collection: Part of CZ Initiative CELLxGENE
- License: CC BY 4.0
- DOI: (Check HCA portal for specific dataset DOI)

### Your Original Data Sources
- Cellular composition: Hanson et al. 2018; Karlsen & Pakkenberg 2011; Morgan et al. 2014
- Metabolic rates: Heiss et al. 2004 (high-resolution PET)

## Notes and Limitations

### Strengths:
- Human data (matches your PET metabolic data)
- Large sample size (3+ million cells)
- Comprehensive brain region coverage
- Multiple donors for standard deviation calculations

### Limitations:
- Single-cell RNA-seq counts are sample-based, not total population counts
- May not perfectly match stereological counting methods from original studies
- Region definitions may differ between atlases
- Volume data still needed for density calculations (not provided by scRNA-seq)

### Important Considerations:
- Extracted counts represent sampled cells, not absolute regional totals
- Use these to fill missing data gaps, not to replace existing literature-derived counts
- Standard deviations reflect donor variability, not measurement uncertainty
- Always validate biological plausibility of filled values

## Contact and Support

For issues with this pipeline:
1. Check troubleshooting section above
2. Review integration log for specific errors
3. Examine intermediate outputs (mapping, extracted counts)

## Version History

- **v1.0** (2026-01-30): Initial pipeline implementation
  - Phase 1: Data access and exploration
  - Phase 2: Region mapping with flexible matching
  - Phase 3: Cell count extraction with classification
  - Phase 4: Data integration with logging

## License

This pipeline is provided as-is for research purposes. The HCA Brain v1.0 data is licensed under CC BY 4.0. Please cite the appropriate sources when publishing results.
