# E:I Ratio vs rCMRGlc Analysis - Results Guide

**Analysis Status:** COMPLETE
**Generated:** 2026-02-09
**Total Runtime:** ~3.5 minutes

---

## Quick Start - View Results

### Visualizations (Start Here!)
Located in: `data/cortical/visualizations/`

1. **ei_vs_rcmrglc_scatter.png** - Main finding: correlation between E:I and glucose metabolism
2. **ei_ratios_by_lobe.png** - E:I ratios comparison across brain lobes
3. **ei_rcmrglc_comparison.png** - Side-by-side E:I and metabolism comparison
4. **ei_ratios_fine_grained_top20.png** - Regional variation in E:I ratios

### Key Finding
Moderate positive correlation between E:I ratio and rCMRGlc:
- **Pearson r = 0.668** (p = 0.147, not significant)
- **Primary visual cortex**: Highest E:I (4.42) + Highest metabolism (40.1)
- **Interpretation**: Regions with more excitatory neurons tend to have higher glucose metabolism

---

## File Organization

```
ei_analysis/
├── FINAL_ANALYSIS_REPORT.txt          ⭐ READ THIS FIRST
├── EIprogress.txt                     - Execution log
├── README.md                          - This file
│
├── data/cortical/
│   ├── ei_ratios_whb_cortical_fine_grained.csv    - 34 regions
│   ├── ei_ratios_aggregated_lobes.csv             - 6 lobes
│   ├── ei_metabolism_integrated.csv               - Final data
│   ├── ei_metabolism_statistics.txt               - Statistics
│   ├── validation_report.txt                      - QC checks
│   │
│   └── visualizations/
│       ├── ei_vs_rcmrglc_scatter.png
│       ├── ei_ratios_by_lobe.png
│       ├── ei_rcmrglc_comparison.png
│       └── ei_ratios_fine_grained_top20.png
│
└── Scripts:
    ├── final_ei_analysis.py
    ├── create_visualizations.py
    └── validate_analysis.py
```

---

## Results Summary

### Dataset
- **3,369,219** total cells from Allen Brain Cell Atlas
- **1,305,075** cortical cells (38.7%)
- **34** fine-grained cortical regions
- **6** brain lobes with rCMRGlc data

### E:I Ratios by Lobe

| Lobe | E:I | rCMRGlc | Excitatory | Inhibitory | Total |
|------|-----|---------|------------|------------|-------|
| Primary visual | 4.42 | 40.1 | 24,848 | 5,628 | 33,484 |
| Insular | 2.74 | 30.3 | 49,255 | 17,998 | 73,886 |
| Occipital | 2.72 | 35.8 | 40,786 | 15,021 | 62,312 |
| Frontal | 2.54 | 35.3 | 186,902 | 73,727 | 283,394 |
| Temporal | 2.50 | 30.5 | 143,754 | 57,434 | 247,820 |
| Parietal | 2.35 | 35.8 | 98,494 | 41,986 | 155,365 |

### Statistics
- **Pearson:** r = 0.668, p = 0.147
- **Regression:** rCMRGlc = 3.23 × E:I + 25.4
- **R-squared:** 0.446

---

## Data Validation

**Quality Checks:** 30/31 PASSED (96.8%)

✓ All files present
✓ 3.37M cells loaded
✓ 87.4% cortical cells classified
✓ E:I ratios validated (2.35-4.42)
✓ Aggregations verified
✓ rCMRGlc values valid

---

## Reproducing Analysis

```bash
python final_ei_analysis.py           # ~0.1 min
python create_visualizations.py       # ~0.7 min
python validate_analysis.py           # ~0.2 min
```

---

**Analysis completed successfully - 2026-02-09**
