# Metabolic dataset by region coverage

Created 12 November 2026. Scope: the five papers audited in this project
(Heiss 2004, Heiss 1984, Heiss 1985, Borghammer 2012 and Andersen 2022).
This is not an exhaustive bibliography of Heiss-coauthored datasets.
Andersen is included for paired FDG/UCB-J convergence, independently of authorship.

The matrix contains 5 datasets × 70 distinct anatomical territories, with
96 reported dataset–territory cells. Territories include nested aggregates,
subregions and source-specific tissue compartments. Neither 70 nor 96 is an
independent regional sample size. Participants and regional observations are
counted separately in the dataset inventory. Available means include uncertain
anatomical matches; those counts are not regression eligibility counts.

| Code | Meaning |
|---|---|
| A | Regional mean available for a named atlas ROI or cortical aggregate |
| S | Regional mean available for a selected spatial sample |
| U | Regional mean available, but anatomical label or tissue boundary uncertain |
| P | Native means published, but checked transcription not yet available here |
| T | Region appears in statistical tests, without a regional mean in that table |
| - | No separate measurement in the audited tables |

An A is not a declaration of an exact match between papers. Establish pairwise
anatomical compatibility before joining outcome data. A missing cell is not
filled from a parent, child, neighboring structure or average of selected nuclei.
In particular, whole atlas thalamus, broad sampled thalamus and individual
thalamic nuclei remain separate; so do lentiform/putamen, unsplit pallidum/GPe/GPi,
and cerebellar cortex/gray matter/atlas cerebellum/central white matter.

## Sources and interpretation

The builder reads the checked 2004 intermediate and source tables in Evo-M1,
including Andersen region-label notes and errata. Heiss 1984 native coverage was
checked against Table 3, p.220 of https://doi.org/10.1038/jcbfm.1984.30.
That table's 22 territories are manually registered; their native rates have
not been transcribed by this inventory. The older comparison column in Heiss
2004 reuses 1984 data and is not counted as another dataset. It cannot transfer
2004 nuclear labels to the older broad territories.

Andersen Table S1 supplies means for 13 regions. HIP adds a fourteenth reported
territory through Table S2, which supplies age-comparison p-values, not regional
means. SN is only conventionally interpreted as substantia nigra in the source
label reference, so it receives U. The unresolved young-CBL discrepancy is
retained as a limitation, without altering published values.

Borghammer contributes eight subcortical structures and two selected WM
reference VOIs. All three normalizations belong to the same 11 controls.
Heiss 1985 contributes 23 regional rows plus a weighted whole-brain aggregate;
the n=3 nondreamer and n=1 dreamer subsets and repeated states remain identifiable.
Unknown cross-paper participant overlap remains unknown, rather than being
counted as independent replication. None of these coverage files refits primary
analyses or changes their inclusion policies.

## Files and regeneration

- `tables/anatomy/metabolic_dataset_region_coverage.csv`: literal dataset rows ×
  region columns, with 350 categorical coverage cells.
- `tables/anatomy/metabolic_dataset_region_evidence.csv`: all 350 cells in long
  form, including source labels, availability, anatomical scope and notes.
- `tables/anatomy/metabolic_dataset_inventory.csv`: cohort, measurement scale,
  dependencies, source and counts for each paper.
- `data_intermediate/metabolic_dataset_region_coverage.json`: workbook input,
  including SHA-256 fingerprints of the machine-readable source files.
- `outputs/01a1d166-5e50-7ef2-a6aa-ee24c6b4b692/metabolic_dataset_region_coverage.xlsx`:
  filterable workbook with regions as rows for readability, plus cohort details
  and evidence for the 96 reported cells.

Run `python3 checks/build_metabolic_coverage_matrix.py`; use `--evo-root` for a
different Evo-M1 location. Assertions check source label sets, unique keys,
regional counts and several critical non-substitutions.

To rebuild Excel, copy `checks/build_metabolic_coverage_workbook.mjs` to a
temporary directory with a `node_modules` symlink to the bundled artifact-tool
runtime. Run the bundled Node executable with the copied script, the absolute
JSON path and the absolute output directory as its two arguments. The builder
recalculates, inspects, renders all three worksheets and exports the workbook.
