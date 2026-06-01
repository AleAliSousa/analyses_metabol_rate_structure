# Subregion-Level E:I Ratios with Lobe rCMRGlc Values

**Analysis Date:** 2026-02-09  
**Sample Size:** n = 20 cortical subregions  
**Data Source:** Allen Brain Cell Atlas WHB-10Xv3

---

## Complete Data Table

| Rank | Region | Lobe | E:I Ratio | rCMRGlc | Excitatory | Inhibitory | Total Cells | E% |
|------|--------|------|-----------|---------|------------|------------|-------------|-----|
| 1 | **V1C** | Primary visual cortex | **4.42** | **40.1** | 24,848 | 5,628 | 33,484 | 74.2 |
| 2 | Ig | Insular lobe | 3.13 | 30.3 | 26,495 | 8,464 | 38,937 | 68.0 |
| 3 | M1C | Frontal lobe | 3.05 | 35.3 | 82,188 | 26,987 | 116,576 | 70.5 |
| 4 | V2 | Occipital lobe | 2.94 | 35.8 | 21,967 | 7,480 | 32,638 | 67.3 |
| 5 | MTG | Temporal lobe | 2.87 | 30.5 | 62,380 | 21,705 | 107,301 | 58.1 |
| 6 | A1C | Temporal lobe | 2.68 | 30.5 | 21,823 | 8,132 | 32,790 | 66.6 |
| 7 | A13 | Frontal lobe | 2.67 | 35.3 | 23,407 | 8,778 | 34,919 | 67.0 |
| 8 | A19 | Occipital lobe | 2.50 | 35.8 | 18,819 | 7,541 | 29,674 | 63.4 |
| 9 | A43 | Parietal lobe | 2.49 | 35.8 | 27,370 | 10,978 | 42,834 | 63.9 |
| 10 | A5-A7 | Parietal lobe | 2.48 | 35.8 | 22,805 | 9,178 | 35,290 | 64.6 |
| 11 | A38 | Temporal lobe | 2.41 | 30.5 | 24,062 | 9,966 | 37,642 | 63.9 |
| 12 | S1C | Parietal lobe | 2.41 | 35.8 | 24,673 | 10,247 | 38,188 | 64.6 |
| 13 | FI | Insular lobe | 2.39 | 30.3 | 22,760 | 9,534 | 34,949 | 65.1 |
| 14 | A25 | Frontal lobe | 2.38 | 35.3 | 24,374 | 10,221 | 37,767 | 64.5 |
| 15 | STG | Temporal lobe | 2.35 | 30.5 | 18,972 | 8,070 | 40,565 | 46.8 |
| 16 | A44-A45 | Frontal lobe | 2.34 | 35.3 | 24,705 | 10,538 | 38,259 | 64.6 |
| 17 | A46 | Frontal lobe | 2.29 | 35.3 | 19,145 | 8,376 | 31,065 | 61.6 |
| 18 | A40 | Parietal lobe | 2.04 | 35.8 | 23,646 | 11,583 | 39,053 | 60.5 |
| 19 | ITG | Temporal lobe | 1.73 | 30.5 | 16,517 | 9,561 | 29,522 | 55.9 |
| 20 | A32 | Frontal lobe | 1.48 | 35.3 | 13,083 | 8,827 | 24,808 | 52.7 |

---

## Statistical Results

### Correlation Analysis
- **Pearson correlation:** r = 0.313, p = 0.180 (NOT significant)
- **Spearman correlation:** rho = 0.105, p = 0.659 (NOT significant)

### Linear Regression
- **Equation:** rCMRGlc = 1.51 × E:I Ratio + 30.14
- **R-squared:** 0.098 (9.8% variance explained)
- **P-value:** 0.180

---

## Comparison: Lobe-Level vs Subregion-Level

| Metric | Lobe-Level (n=6) | Subregion-Level (n=20) |
|--------|------------------|------------------------|
| Pearson r | 0.668 | 0.313 |
| P-value | 0.147 | 0.180 |
| R-squared | 0.446 | 0.098 |
| Significance | No | No |
| Correlation | Moderate | Weak |

**Key Finding:** The correlation is WEAKER at the subregion level, suggesting the lobe-level correlation was partially an aggregation artifact.

---

## Within-Lobe E:I Variation

| Lobe | n | E:I Range | rCMRGlc | Variation |
|------|---|-----------|---------|-----------|
| **Frontal** | 6 | **1.48 - 3.05** | 35.3 | 2.1× range |
| Temporal | 5 | 1.73 - 2.87 | 30.5 | 1.7× range |
| Parietal | 4 | 2.04 - 2.49 | 35.8 | 1.2× range |
| Insular | 2 | 2.39 - 3.13 | 30.3 | 1.3× range |
| Occipital | 2 | 2.50 - 2.94 | 35.8 | 1.2× range |
| Primary visual | 1 | 4.42 | 40.1 | N/A |

**Note:** Frontal lobe shows the widest E:I variation (1.48 to 3.05), but all subregions share the same rCMRGlc value (35.3).

---

## Key Observations

### ✅ Supporting Evidence
- **V1C (Primary visual cortex)** has BOTH the highest E:I ratio (4.42) AND highest metabolism (40.1)
- Positive correlation trend (r=0.31), consistent with hypothesis direction

### ⚠️ Limitations
1. **Same rCMRGlc for all subregions within a lobe** - masks true subregion-specific relationships
2. **Horizontal clustering** in scatter plot - limited ability to detect correlations
3. **High within-lobe E:I variation** - suggests need for subregion-specific metabolic data

### 📊 Statistical Power
- Even with n=20 (3× more than lobe-level), correlation remains non-significant
- Need region-specific (not lobe-averaged) rCMRGlc data for proper testing

---

## Recommendations

1. **Obtain finer-grained metabolic data** - Brodmann area or gyrus-level PET measurements
2. **Focus on primary cortices** - V1, M1, S1, A1 show most variation
3. **Layer-specific analysis** - E:I varies by cortical layer (2/3 vs 5/6)
4. **Metabolic gene expression** - Analyze glycolysis/oxidative metabolism markers in RNA-seq

---

## Files

- **Data:** `subregion_ei_rcmrglc_table.csv`
- **Statistics:** `subregion_statistics.txt`
- **Visualizations:** `subregion_ei_vs_rcmrglc_scatter.png`, `subregion_ei_distribution_by_lobe.png`
- **Location:** `data/cortical/subregion_analysis/`
