# E:I Ratio Analysis - Data Requirements Estimate

## Summary for Limited-Resource Laptop

**RECOMMENDED APPROACH: Metadata-Only Download**
- **Total Download Size**: ~0.5-1.5 GB
- **Disk Space Needed**: ~2-3 GB (with processing overhead)
- **RAM Required**: 4-8 GB
- **Processing Time**: 15-45 minutes

---

## Data Sources

### 1. Human Brain Cell Atlas (HCA) v1.0 - Cortical Regions

**Available via CELLxGENE Platform**

**Cortical Coverage**:
- Frontal lobe: M1C, A44-45, A13-14, supplementary motor
- Temporal lobe: MTG, STG, ITG, A38, entorhinal
- Parietal lobe: S1C, A5-A7, A40
- Occipital lobe: V1C, V2, A19
- Insular cortex

**Data Specifications**:
- **Total cells**: ~3 million cortical neurons
- **Cell type annotations**: 461 clusters, 3313 subclusters
- **Metadata size**: ~300-600 MB (cell annotations, regions, types)
- **Full H5AD size**: ~15-30 GB (includes expression matrices)

**What We Need**:
- ✓ Cell metadata CSV: Region labels, cell type annotations, cluster assignments
- ✗ Expression matrices: NOT needed for count-based E:I ratios

### 2. Allen Brain Cell (ABC) Atlas - ASAP-PMDBS

**Available via Allen Institute Portal**

**Coverage**:
- Whole human brain including all cortical regions
- 220 donors, 3 million cells

**Data Specifications**:
- **Metadata size**: ~100-200 MB
- **Full expression data**: ~50-80 GB

**What We Need**:
- ✓ Cell metadata: Region, supercluster, cluster, subcluster annotations
- ✗ Expression matrices: NOT needed for primary analysis

---

## Download Strategy: Metadata-Only Approach

### Phase 1: HCA Cortical Metadata (~500 MB)

**Download Method**: CELLxGENE API or direct CSV download
- Cell metadata tables only
- Contains: cell_id, region, tissue, cell_type, cluster, supercluster

**Files**:
```
data/cortical/hca_cortical_metadata.csv        ~300-500 MB
```

### Phase 2: ABC Cortical Metadata (~100 MB)

**Download Method**: AbcProjectCache with metadata_only=True
- Cell annotations only
- Contains: cell_label, region_label, supertype, class, subclass

**Files**:
```
data/cortical/abc_cortical_metadata.csv        ~100-200 MB
```

### Total Metadata Download: ~0.5-1 GB

---

## E:I Ratio Calculation Methods

### Method 1: Count-Based (Primary - Works with Metadata Only)

**How it works**:
- Classify cells as Excitatory or Inhibitory based on annotations
- Excitatory keywords: "excitatory", "glutamatergic", "Glut", "IT", "ET", "CT"
- Inhibitory keywords: "inhibitory", "gabaergic", "GABA", "PVALB", "SST", "VIP"
- Count cells per region
- Calculate E:I ratio = (# Excitatory) / (# Inhibitory)

**Requirements**:
- ✓ Only needs metadata
- ✓ Fast processing (<30 min)
- ✓ Sufficient for correlation analysis

### Method 2: Expression-Based (Optional Validation - Requires Full Data)

**How it works**:
- Calculate marker gene expression: SLC17A7 (excitatory), GAD1/GAD2 (inhibitory)
- Requires full expression matrices

**Requirements**:
- ✗ Needs 30-80 GB downloads
- ✗ High RAM (16+ GB)
- ✗ Slow processing (hours)
- Not recommended for laptop

---

## Recommended Workflow

### Step 1: Assessment (This Script)
- Verify dataset availability
- Check disk space

### Step 2: Download Metadata Only (~30-60 min)
```bash
python download_hca_cortical_metadata.py --metadata-only
python download_abc_cortical_metadata.py --metadata-only
```

**Download size**: ~0.5-1 GB
**Disk space**: ~2 GB (with decompression)

### Step 3: Calculate E:I Ratios (~10-20 min)
```bash
python calculate_ei_ratios_count_based.py --source hca
python calculate_ei_ratios_count_based.py --source abc
```

**Processing**: In-memory, minimal disk writes
**Output**: Small CSV files (<10 MB)

### Step 4: Integration & Analysis (~5 min)
```bash
python aggregate_ei_ratios_to_lobes.py
python integrate_ei_metabolism.py
python statistical_analysis_ei_metabolism.py
```

**Output**: Final analysis CSV (~1 MB)

---

## Storage Requirements Summary

| Component | Size | Required? |
|-----------|------|-----------|
| HCA metadata | 300-500 MB | ✓ Yes |
| ABC metadata | 100-200 MB | ✓ Yes |
| Processing temp files | 500 MB | ✓ Yes |
| Output CSVs | 50 MB | ✓ Yes |
| **TOTAL REQUIRED** | **~1.5-2.5 GB** | |
| | | |
| HCA expression matrices | 15-30 GB | ✗ No |
| ABC expression matrices | 50-80 GB | ✗ No |

---

## Laptop Compatibility Check

**Minimum Requirements**:
- ✓ Disk space: 3 GB free
- ✓ RAM: 4 GB (8 GB recommended)
- ✓ Python with pandas, numpy
- ✓ Internet connection for download

**Expected Performance**:
- Download time: 30-90 minutes (depends on connection)
- Processing time: 15-45 minutes
- Total time: 1-2 hours

**If you have less than 3 GB disk space**:
- Option 1: Download HCA only (skip ABC) - reduces to ~1 GB
- Option 2: Download by lobe (Frontal first, then Temporal, etc.) - ~200 MB each

---

## Next Steps

1. **Check your available disk space**:
   - Required: 3 GB minimum
   - Recommended: 5 GB (for safety margin)

2. **Confirm approach**:
   - Metadata-only (recommended): 1-2 hours, ~2 GB
   - Single-source (HCA or ABC only): ~1 hour, ~1 GB
   - By-lobe incremental: Multiple sessions, ~200 MB each

3. **Proceed with download scripts** once disk space is verified

---

## Alternative: Cloud Processing

If laptop resources are insufficient:
- Use Google Colab (free, 12 GB RAM, 100 GB disk)
- Download data in cloud
- Run analysis there
- Download only final results (~10 MB)

Would you like me to create Colab-compatible notebooks instead?
