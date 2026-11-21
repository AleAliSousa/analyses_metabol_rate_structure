# analyses_metabol_rate_structure_nonhuman

Mouse replication of the human s1b analyses (regional cell-type composition vs
regional glucose utilisation). See `NOTES.md` for the project note.

## Data
- **Cell types per structure:** Allen Brain Cell Atlas whole-mouse-brain 10x
  dataset (Yao et al. 2023, Nature 624:317), release 20231215: 4.04 M cells,
  29 dissected structures (ROIs), taxonomy class/subclass/supertype/cluster.
  Fetched by `mouse_s1b_1` (about 0.9 GB, git-ignored); only the aggregated
  counts `data_intermediate/allen_wmb10x_counts_roi_donor_celltype.csv` are kept.
- **Mouse rates:** `Evo-M1-Trait-Data/__merging_cerebral_metabolic_rate`
  (`Mus musculus`, `CMRgl`, umol/100 g/min; 22 regions from Jay 1985,
  Bouilleret 2000, Yamaguchi 1992, Itoh 2001, ...). Snapshotted by `mouse_s1b_0`.

## Pipeline (`scripts/run_mouse_s1b.sh`)
| Script | Human counterpart | Output |
|---|---|---|
| `mouse_s1b_0_sync_mouse_rates.R` | `00_sync_evo_m1_inputs.R` | `data_raw/Evo_M1_mouse_CMRgl_by_region.csv` |
| `mouse_s1b_1_extract_allen_wmb10x.R` | `s1b_1_extract_transcriptomic` | counts by ROI x donor x class/subclass |
| `mouse_s1b_2_validate_roi_map.R` | `s1b_2_mapping_*` | `mouse_s1b_roi_mapping_coverage.csv` |
| `mouse_s1b_3_n_neuronal_composition.R` | `s1b_3_n_*` | `tables/`, `figs/mouse_s1b/neuronal_*` |
| `mouse_s1b_3_nn_nonneuronal_composition.R` | `s1b_3_nn_*`, `s1b_6`, `s1b_7` | `nonneuronal_*` |
| `mouse_s1b_5_EI_ratio_cortex.R` | `s1b_5` | `EI_*` incl. S1 (SSp) vs other isocortex |

Crosswalks (all decisions live here, scripts only validate them):
`metadata/mouse_roi_to_cmrgl_region_map.csv`,
`metadata/mouse_class_to_s1b_category_map.csv`,
`metadata/mouse_subclass_to_nonneuronal_category_map.csv`.

## Caveats to review
1. **Data are whole-cell 10x scRNA-seq (Yao 2023), not snRNA-seq.** No Allen
   whole-brain snRNA atlas with per-structure dissection exists for mouse, so this
   is the closest equivalent; neuron-enriched sorting means neuron:non-neuron
   ratios are not trusted. Non-neuronal proportions are computed within
   non-neuronal cells for this reason.
2. Allen ROIs are coarse dissections and mouse rate regions are coarse
   (e.g. Sensorimotor_cortex pools S1 and M1, so SSp and MOp are averaged).
   `primary_included` marks exact/contained matches; partial matches are
   sensitivity only. See the map notes.
3. Few donors per ROI (SSp has 3). n = 12 regions with rates, so results are
   exploratory, as for the human s1b nulls.
4. Class categories and exclusions (MSN, neuromodulatory) are a proposal that
   mirrors the human grouping; please review the class map.
