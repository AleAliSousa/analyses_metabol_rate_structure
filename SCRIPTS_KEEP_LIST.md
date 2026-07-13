# Scripts keep list — outputs that reach the presentation / manuscript

**Generated:** 2026-07-11
**Presentation:** `Brain energetics/Energetic Constraints on Brain Organization in Human Paleoneurology_11072026.pptx` (69 slides)
**Manuscript:** `Brain energetics/MS Is human brain organization economical_04072026_TRACKED.docx`

## Method

Every embedded image in the PPTX (95 rasters) and DOCX (8 rasters) was hashed
(exact MD5 + perceptual dHash) and matched against all 319 figure files in
`figs/`. A match means the figure a script produces is literally embedded in the
deck. Producing scripts were then traced back through their input files to
include the data-prep steps they depend on.

**Manuscript note:** none of the 8 images in the tracked `.docx` are
script-generated figures — they are schematics/photographs. So the keep list
below is driven entirely by the presentation.

---

## Tier A — Keep: figure is confirmed embedded in the deck (pixel match)

| Script | Figure(s) produced | Slide(s) |
|---|---|---|
| `neocortex_grey_white.R` | `figs/traits/neocortex_gray_white.png` | 8 |
| `s1a_1_stereology_cell_types_30052026.R` | `figs/s1a/p_counts.jpg`, `figs/s1a/p_densities.jpg` | 58, 18 |
| `s1b_4_n_supercluster_rcmr_correlation_matrix_telencephalon_13062026.R` | `figs/s1b/scatter_rcmr_supercluster_{all_regions,telencephalon,non_telencephalon}.png` | 57, 60, 61 |
| `s2_stress_volume_01062026.R` | `figs/s2/mean_brain_region_volumes_boxplot.jpg`, `figs/s2/rCMRGlc_volume_change.jpg` | 30, 31 |
| `s4_endocranial.R` | `figs/s4/endocranial_region_budgets.png` (65), `endocranial_region_cost_pie_3species_volctx.png` (45), `rcmrglc_vs_fig3a_deviation.png` (43), `volume_timeline_logage_flipped.png` (46) | 43, 45, 46, 65 |
| `s4_endocranial_cerebellum.R` | `figs/s4/cerebellum_volume_timeline_logage.png` | 68 |

## Tier B — Keep: analysis + plot are in the deck, but the embedded image is an earlier render

These scripts are the current producers for topics the deck clearly presents,
but the picture on the slide is an older/cropped render, so it does not
pixel-match the current file. Keep them — they are core content.

| Script | Deck topic | Slide(s) |
|---|---|---|
| `s1b_5_n_EI_ratio_original_vs_jorstad_overlay_two_MSN_plots_raw_EI_only_16062026.R` | Excitatory:Inhibitory ratio (Study 1) | 24, 25 |
| `s1b_6_nn_type1_type2_astrocyte_compositional_rcmr_26052026.R` | Astrocyte composition vs rCMRglc (Study 1) | 21, 22 |
| `s3_predicValuesPGLS_MERGED_variant.R` | Study 3 evolutionary deviation engine | 33–39 |
| `s3_run_and_compare_configs.R` | Study 3 driver (sources the engine per config) | 33–39 |
| `s3_run_frontal_vermis_anthro.R` | Study 3 driver (frontal/vermis focus) | 33–39 |

## Tier C — Keep: upstream data-prep (no figure of their own, but required inputs)

| Script | Produces / does | Feeds |
|---|---|---|
| `0_Heiss_Stephan_and_table1_30052026.R` | `data_intermediate/Heiss_Stephan_data.csv` (shared rCMRglc × volume reference) | s1b_4/5/6, s3 |
| `0_bind_matano_1985a_to_stephan.R` | augments `Stephan_primates` volumes | s3, s4 |
| `s1b_1_n_extract_transcriptomic_neuronal_30052026.R` | `..._neuronal.rds` | s1b_4, s1b_5 |
| `s1b_1_nn_extract_transcriptomic_nonneuronal_30052026.R` | `..._nonneuronal.rds` | s1b_6 |
| `s1b_2_mapping_rcmrglc_transcriptomic_cells_anatomy_21052026.R` | maps cells → anatomy → `Heiss_Stephan_data.csv` | s1b_4/5/6 |
| `s1b_3_*` transcriptomic (neuronal/nonneuronal ± telencephalon) | add proportion columns to `Heiss_Stephan_data.csv` | s1b_4/5/6 |
| `R/plot_settings.R` | shared ggplot theme | sourced by most s1b/s3 scripts |

---

## Not required for the deck or manuscript (context)

Verified to have **no** output in either file:

- `v2_0_Heiss_Stephan_and_table1_05062026.R` — writes `Heiss_Stephan_data_v2.csv`, which no other script reads (orphan; the deck pipeline uses the v1 file).
- `s1a_2_stereology_proportions_30052026.R` — stereology proportion pies; not embedded.
- `s1b_5_n_EI_ratio_telencephalon_26052026.R` — superseded by the `..._raw_EI_only_16062026` version.
- QC / diagnostic / robustness (no slide): `s1b_2_check_dissection_roi.R`, `s1b_x_check_dissection_roi.R`, `s3_0_missingness_clade_diagnostic_04062026.R`, `s3_1_phylo_multiple_imputation_04062026.R`, `s3_compare_stephan_vs_merged.R`, `network_residual_autocorrelation_analysis.R`.

## Snapshot fixes (2026-07-11)

Figures that appeared in the deck but were saved by hand (no reproducible export
in the script) have been repaired so each writes both a raster (PNG, for slides)
and a vector PDF:

- `neocortex_grey_white.R` — was an interactive `dev.new()` + `dev.copy(png)`
  snapshot (PNG only, broke head-less). Now builds the plot in a reusable
  function and writes `neocortex_gray_white.png` + `.pdf`.
- `s4_endocranial.R` — `endocranial_region_budgets`, `endocranial_region_cost_pie_3species_volctx`,
  `volume_timeline_logage_flipped` now export PNG **and** PDF (via a
  `save_png_pdf()` helper).
- `s4_endocranial_cerebellum.R` — `cerebellum_volume_timeline_logage` now exports
  PNG + PDF.
- `figs/s4/budget_significance.png` (deck slide 69) had **no source anywhere**.
  Reconstructed as a new block in `s4_endocranial_cerebellum.R`: whole-brain
  (from `species_absolute_budgets.csv`) + cerebellum group budgets, mean +/- 95%
  CI vs MH = 1.0, with a two-sided one-sample t-test p-value. Writes
  `budget_significance.png` + `.pdf` + `budget_significance.csv`. Current-data
  p-values (~0.54 / 0.37 / 0.65 / 0.49) are close to the snapshot's
  (0.60 / 0.40 / 0.63 / 0.43); the small drift reflects data updates since the
  snapshot was taken.

## Still to resolve

- **`figs/s4/budget_residual_vs_volume.png` (slide 69) could not be
  reconstructed.** No plotting code exists in the repo or R history, and none of
  the natural definitions of "the part not explained by brain size" reproduce
  its values (±0.09 to ±0.26): in the current model the budget is almost
  perfectly determined by brain size, so any size-residual comes out ≈ 0. The
  original must have used a different decomposition or older shape/rate inputs.
  **Needs the original code or a definition of the residual to regenerate.**
- The edited scripts pass static checks (brace/paren/bracket balance, device
  open/close pairing) but could **not be executed here** (R is not installed in
  this environment). Run `s4_endocranial.R` then `s4_endocranial_cerebellum.R`
  once to confirm the new PNG/PDF pairs render.
- Tier B figures were matched by topic, not pixels, because the deck predates the
  current renders. If you re-export those slides from current scripts, the images
  will change slightly.
