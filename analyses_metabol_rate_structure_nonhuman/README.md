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

## Pipeline (`scripts/run_mouse_pipeline.sh`)
| Script | Human counterpart | Output |
|---|---|---|
| `mouse_s1b_0_sync_mouse_rates.R` | `00_sync_evo_m1_inputs.R` | `data_raw/Evo_M1_mouse_CMRgl_by_region.csv` |
| `mouse_s1b_1_extract_allen_wmb10x.R` | `s1b_1_extract_transcriptomic` | counts by ROI x donor x class/subclass |
| `mouse_s1b_2_validate_roi_map.R` | `s1b_2_mapping_*` | `mouse_s1b_roi_mapping_coverage.csv` |
| `mouse_s1b_3_n_neuronal_composition.R` | `s1b_3_n_*` | `tables/`, `figs/mouse_s1b/neuronal_*` |
| `mouse_s1b_3_nn_nonneuronal_composition.R` | `s1b_3_nn_*`, `s1b_6`, `s1b_7` | `nonneuronal_*` |
| `mouse_s1b_5_EI_ratio_cortex.R` | `s1b_5` | `EI_*` incl. S1 (SSp) vs other isocortex |

### s1a (stereology analogue) and s1c (synaptic density)
| Script | Human counterpart | Output |
|---|---|---|
| `mouse_s1a_0_sync_isotropic_fractionator.R` | `00_sync_evo_m1_inputs.R` | mouse Herculano-Houzel 2011 counts (validity check only) |
| `mouse_s1a_1_extract_allen_merfish_ccf.R` | `s1a_1` input assembly | counts by CCF structure x cell class; CCF volumes |
| `mouse_s1a_2_regional_cell_composition_vs_rcmr.R` | `s1a_1_stereology_cell_types` | `tables/mouse_s1a/`, `figs/mouse_s1a/` |
| `mouse_s1c_synaptic_density_vs_rcmr.R` | `s1c_synaptic_density_metabolic_rate` | `tables/mouse_s1c/` (needs input, see below) |

s1a uses the Allen whole-mouse-brain **MERFISH** atlas (in situ, 3.74 M cells, one donor,
CCF-registered), so neuron:glia balance is not altered by cell sorting as in the 10x data.
Densities are cells per mm^3 of CCF volume: a relative measure, use ranks only.
17 CMRgl regions are matched (`primary_all`); grey-matter-only and with-partial-region
sensitivity sets are also written.

**s1c is not runnable yet:** no mouse regional synaptic-density table exists in
Evo-M1-Trait-Data and the candidate sources (SV2A Mouse Synaptome Atlas, Edinburgh)
were unreachable from the cloud container. The script is tested (missing-input stop and a
synthetic dry run) and waits for `data_raw/mouse_synaptic_density_by_region.csv`; the
contract and candidate sources are in `metadata/mouse_synaptic_density_input_spec.md`.

Crosswalks (all decisions live here, scripts only validate them):
`metadata/mouse_roi_to_cmrgl_region_map.csv`,
`metadata/mouse_class_to_s1b_category_map.csv`,
`metadata/mouse_subclass_to_nonneuronal_category_map.csv`,
`metadata/mouse_ccf_structure_to_cmrgl_region_map.csv` (s1a).

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
5. **MERFISH neuron fractions disagree with the isotropic fractionator**
   (`tables/mouse_s1a/validation_vs_isotropic_fractionator.csv`): cortex 0.70 vs 0.53,
   cerebellum 0.72 vs 0.86. Segmentation of small dense cells (cerebellar granule
   cells) is the likely cause. Treat neuron:glia ratios as relative between regions,
   not absolute.
6. Only two mouse regions (cortex, cerebellum) can be compared with isotropic-fractionator
   counts, so those are not used for regression.
