# Session handoff — repo audit & reproducibility work

**Updated:** 2026-08-17

> ## 2026-08-17 — low-risk cleanup pass (this update)
>
> Done:
> - `run_all.R` rewritten to match the scripts that actually exist (the old
>   list named 6 scripts deleted in commit `7ae5165`, including the s3
>   engine + drivers; new scripts `s1c`, `s4b`, `network_residual` and both
>   s3 forks now run in a defined order). `verify_outputs.R` is no longer a
>   pipeline step — it is sourced as a **gate** after the run summary.
> - Duplicate pairs merged: `s1b_1_n`/`s1b_1_nn` →
>   `s1b_1_extract_transcriptomic_30052026.R`; `s1b_x_check_dissection_roi.R`
>   (the superset) → `s1b_2_check_dissection_roi.R`.
> - Superseded `s1b_5_n_EI_ratio_telencephalon_26052026.R` moved to
>   `scripts/archive/` and dropped from the pipeline.
> - Junk removed and git-ignored: `.__captest/`, `.Rapp.history`,
>   `.DS_Store` (tracked copies `git rm`'d).
> - README directory tree + `SCRIPTS_KEEP_LIST.md` synced to reality.
> - Orphan prune: `figs/s3/phase1/` (4 figs from the retired
>   `s3_0_missingness` diagnostic) and `checks/Stephan_primates_data_used.csv`
>   (no producer) deleted; MI sensitivity analysis confirmed retired, so the
>   Phase-1 scripts were NOT restored (they remain in git history pre-`7ae5165`
>   and in the Dropbox `archive/`). References in the s3 headers,
>   `PHASE1_missing_data_strategy.md`, keep list and `figs_inventory.csv`
>   updated to say so.
>
> Still open (larger-risk items, deferred): rebuild s3 as one engine + config
> drivers (the two forks both write to `figs/s3/all` etc. — last run wins);
> centralize the root-finder/`save_png_pdf`/`fmt_p` boilerplate in `R/`;
> remove the remaining hardcoded absolute paths (both s3 scripts'
> Dropbox `setwd`, OneDrive inputs in s1b_1/s3-VWS/qc_stephan); move
> `data_analysis/` → `tables/s1b/`; orphaned outputs
> (`checks/Stephan_primates_data_used.csv`, the 7 `figs/s3/*` config dirs
> from the deleted drivers, ~190 files).

---

**Previous update:** 2026-07-14
**Repo:** `analyses_metabol_rate_structure`
**Deliverables this repo feeds:**
`Brain energetics/Energetic Constraints on Brain Organization in Human Paleoneurology_11072026.pptx` (69 slides)
and `Brain energetics/MS Is human brain organization economical_04072026_TRACKED.docx`.
Analyses are organized as four studies: s1 (cellular composition), s2 (environmental
stress), s3 (evolutionary deviation, PGLS), s4a (fossil endocranial budgets).

---

## What has been done

**Audit of scripts vs. deck (see `SCRIPTS_KEEP_LIST.md`).**
Every embedded deck image was hashed and matched against `figs/`. 14 deck figures
map to a producing script; the manuscript embeds no script-generated figures (all
schematics/photos). Superseded/duplicate scripts were moved to
`scripts/archive/pruned_2026-07-11/` (10 files).

**Snapshot / export fixes.** Figures saved by hand now export both a raster (PNG/JPG
for slides) and a vector PDF:
- `traits_neocortex_grey_white.R` (was `dev.copy` snapshot) → PNG+PDF via a draw function.
- `s4a_endocranial.R` → PNG+PDF for `endocranial_region_budgets`,
  `endocranial_region_cost_pie_3species_volctx`, `volume_timeline_logage_flipped`
  (via a `save_png_pdf()` helper).
- `s4a_endocranial_cerebellum.R` → PNG+PDF for `cerebellum_volume_timeline_logage`,
  and a **reconstructed `budget_significance`** figure (whole-brain + cerebellum
  group budgets, mean ± 95% CI vs MH=1.0, two-sided t p-values) + CSV.

**Pipeline runner** `scripts/run_all.R` — runs all scripts in dependency order, each
in its own environment, continue-on-error, writes `logs/run_all_status_latest.csv`.
The s3 engine (`s3_predicValuesPGLS_MERGED_variant.R`) is sourced by its two `run_`
drivers, not run standalone.

**Reproducibility manifest + checker.**
- `metadata/output_manifest.csv` — 18 deliverable outputs: path, formats, producing
  script, deck slide, and a `check` mode (`file` / `script_only` / `none`).
- `scripts/verify_outputs.R` — reads the manifest and asserts outputs exist after a
  run; reports the known gap; exits non-zero on real failures.

**Repo hygiene (Part 3).**
- Deleted strays: 22 `Rplot*.png` default-device artifacts and the root
  `rcmrglc_vs_fig3a_deviation.png` / `_fit.png` (canonical copies remain in `figs/s4a/`).
- `.gitignore` extended: `Rplot*.png`, root-level `*.png/jpg/jpeg/pdf/csv`, `.RData`, `logs/`.
- Replaced the hardcoded `setwd("~/Library/CloudStorage/Dropbox/…")` in 27 scripts
  with a portable walk-up-to-`.git` root finder. Added `R/project_root.R` (reusable
  snippet) and an `.Rproj` (so `here::here()` also resolves the root).

**Full plan of record:** `metadata/REPO_AUDIT_PLAN.md`.

---

## Open items / next steps (priority order)

1. **Run the pipeline on a machine with R** (it is not installed in the Cowork
   sandbox — all edits above are verified by static syntax checks only):
   `Rscript scripts/run_all.R` then `Rscript scripts/verify_outputs.R`.
   Confirm every deck figure regenerates as PNG+PDF and the status logs are clean.
2. **Commit.** `git add -A && git commit` — stages the stray deletions and the new
   files. (The `.git` index lock encountered during the session has been cleared.)
3. **`budget_residual_vs_volume` (slide 69) still has NO producing script.** Its
   original derivation is unrecoverable; no natural size-residual reproduces its
   ±0.09–0.26 values (the model's budget is ~fully explained by size). Supply the
   original code or agree a definition, then it can be added to `s4a_endocranial.R`
   with PNG+PDF export.
4. **Clean-room reproducibility test + dead-output prune (Part 1).** Move `figs/`,
   `tables/`, `data_intermediate/` aside, run `run_all.R`, confirm manifest outputs
   reappear. Then list the ~448 files in `figs/` not in the manifest and not consumed
   downstream → archive/delete.
5. **Naming & organization (Part 2, do after the manifest is trusted).**
   - Move QC scripts to `scripts/checks/` with `check_`/`qc_` prefixes
     (`s1b_2_check_dissection_roi`, `network_residual_autocorrelation_analysis`;
     the `s3_0_missingness` / `s3_compare_stephan_vs_merged` scripts were
     retired in `7ae5165`).
   - Rebuild Study 3 as one engine in `R/` (e.g. `R/pgls_engine.R`) plus thin
     per-config drivers — the two current `s3_predicValuesPGLS_*` forks are
     near-duplicates that both write to `figs/s3/all` (last run wins).
   - Drop dated/`PATCHED` suffixes → one canonical name per analysis.
   - Shared Heiss ownership is resolved: `00_prepare_heiss_rates.R` solely writes
     `data_intermediate/heiss_2004_regions.csv`; the obsolete mixed Heiss/Stephan
     intermediate has been removed.
6. **`Stephan_primates` naming/versioning decision (pending).** Recommendation:
   keep the name (or rename once to a source-neutral name with a per-row `source`
   column) and track provenance via a column + log, not the filename.
   `build_stephan_primates_reference_sheet.R` already resolves per-datapoint
   provenance. "merged" in s3 is a separate axis (multi-source mode) — keep it.
7. **Tier B deck figures** (E:I slides 24–25, astrocyte 21–22, Study 3 37–38): the
   deck shows older renders. Re-export current script outputs if you want the deck
   to match the current analysis.

---

## Key files (created/modified this session)

- Created: `scripts/run_all.R`, `scripts/verify_outputs.R`, `R/project_root.R`,
  `analyses_metabol_rate_structure.Rproj`, `metadata/output_manifest.csv`,
  `metadata/REPO_AUDIT_PLAN.md`, `SCRIPTS_KEEP_LIST.md`,
  `scripts/archive/pruned_2026-07-11/` (10 archived scripts).
- Modified: `traits_neocortex_grey_white.R`, `s4a_endocranial.R`,
  `s4a_endocranial_cerebellum.R`, 27 scripts (`setwd` → root finder), `.gitignore`.
- Deleted: 22 `Rplot*.png`, root `rcmrglc_vs_fig3a_deviation{,_fit}.png`.

## Environment caveats

- R is not installed in the Cowork sandbox; scripts were checked statically only.
- The repo is Dropbox-synced; deleting files required elevated permission and git
  index locking was flaky — prefer running git yourself.

## 2026-10-30 — Study s1d added (neuroglial proteins vs rCMRGlc)

- New study **s1d** replicates the s1a_1 / s1a_2 design with Sjölin et al. 2022
  (Mol Brain 15:52) regional ELISA concentrations of GFAP, MBP, NFL, tau and
  UCHL1 in place of stereological cell counts/densities.
  - `scripts/s1d_1_neuroglial_proteins_metabolic_rate.R` — Spearman loop + BH,
    multivariate lm (caution: 11 regions), faceted scatter with per-panel lm;
    two measure families (µg/g wet tissue; µg/mg total protein) plus GFAP:NFL
    ratio. Sensitivity tables: direct matches only, grey matter only,
    within-donor correlations. Writes
    `data_intermediate/s1d_neuroglial_protein_regional_summary.csv`.
  - `scripts/s1d_2_neuroglial_protein_proportions.R` — per-region composition
    pies (all five / excluding MBP) and combined facet figures.
  - `metadata/anatomy/s1d_sjolin_region_map.csv` (validated by
    `00_prepare_heiss_rates.R`); registry gained `pons`, `medulla_oblongata`,
    `cervical_spinal_cord` (no Heiss rate).
  - Manifest row `sjolin_neuroglial_proteins` ->
    `data_raw/Sjolin_etal_2022_SupplementaryData.csv`; four rows added to
    `metadata/output_manifest.csv`; both scripts in `run_all.R` after s1c.
- Pipeline fix: upstream renamed `Hanson_etal_2018_Table2_derived.csv` to
  `Hanson_etal_2018_Table2.csv` (identical MD5 fa09ef00…); the manifest
  `source_relative_path` was repointed, local path unchanged. Without this the
  sync step refused to run.
- Side effect of the sync: `data_raw/Mackes_etal_2020_ERABIS_derived.csv` and
  `data_raw/volumes_wide.csv` were refreshed from upstream (s2 / traits
  inputs) — review the diff before committing.
- **Study 1 synthesis** added: `scripts/s1_synthesis_cross_study.R` reads the
  region-level tables of s1a-s1d, recomputes Spearman rho vs rCMRGlc for all
  36 predictors (BH within study x block, plus a global BH), checks the
  recomputed values against the Spearman tables s1b_6, s1b_7 and s1d_1 write
  (stops on mismatch), and writes `tables/s1_synthesis/` +
  `figs/s1_synthesis/` (rho forest, region-coverage heatmap). Only the four
  neocortical lobes are covered by all four studies.

## 2026-10-31 — Study s5 added (Coto-style analyses on Heiss 2004; humans–mice–ants project)

- New study **s5** does on the Heiss 2004 regions what Coto et al. 2026 (PNAS, *Tetramorium
  immigrans*) did on ants. Scripts (run in order; conda env `evom1-r`, from anywhere in the repo:
  `Rscript scripts/s5a_regional_neurometabolic_scaling.R` …):
  - `scripts/s5a_regional_neurometabolic_scaling.R` — regional glucose use = rCMRGlc × volume;
    SMA/OLS log–log scaling on regional volume (all / excl. white matter / grey-only / Kochiyama
    sensitivity). New map `metadata/anatomy/s5_heiss_region_volume_map.csv` (Table 1 pairings +
    s4a Kochiyama lobes; extra columns volume_set, volume_source_*, allocation_weight,
    volume_basis). Tables `tables/s5/s5a_*`, figs `figs/s5/s5a_*`.
  - `scripts/s5b_relative_brain_metabolism.R` — volume-weighted whole-brain rCMRGlc, brain W,
    relative brain mass / relative brain MR (glucose and O2 routes) with mouse and ant rows;
    reads the nonhuman repo `data_intermediate/species_brain_body_metabolism.csv`.
  - `scripts/s5c_stress_budget_preservation.R` — s2 inputs: summed budget vs summed volume change
    (normative → English / Romanian adoptees) with parametric bootstrap CIs; s2 regression in
    standardised form.
  - `scripts/s5d_cross_species_brain_mr_scaling.R` — Coto Fig. 7 replication (Coto table and
    Evo-M1 species table) with prediction-band checks for human, mouse, ant.
  - `scripts/s5e_export_comparison_inputs.R` — writes `human_comparison_rows.csv`,
    `human_regional_change_vs_rate.csv`, `human_scaling_fits.csv` to
    `../analyses_metabol_rate_structure_nonhuman/comparison/inputs/`.
  - `helpers/s5_helpers.R` — annotate_lm / save_png_pdf / make_region_palette / fit_sma_ols and
    the energy constants (mirrors the nonhuman repo helpers; ggpmisc not available).
  - Methods and all decisions: `tables/s5/METHODS_NOTE_s5.md`. Logs in `logs/s5*.log`.
- Not yet done: s5 is not in `run_all.R`, `metadata/output_manifest.csv` or `verify_outputs.R`
  (left untouched by design this session); add when the study is accepted. No existing file was
  modified.


## 2026-11-12 — Model diagnostics, s1a/s2 data flow, literature expansion, 27 new Evo-M1 items

### State right now (read first)

| What | Status |
|---|---|
| Model diagnostics (s1–s5) | In scripts; outputs from the user's 2026-11-05 10:23 `run_all.R` (116 models). |
| s1a Hanson caudate/putamen + units fix | **Done 2026-11-12.** `run_all.R` rerun 11:17–11:24 (34/34 scripts OK). `data_intermediate/s1a_stereology_comparison.csv` has 9 regions incl. Caudatum and Putamen, Hanson densities per cm3; s1a figures and diagnostics regenerated (n = 9 for neuron/glia density and glia:neuron ratio). |
| 25 new manifest rows | **Synced 2026-11-12.** All 25 files in `data_raw/`, byte-identical to their Evo-M1 sources and recorded in `data_intermediate/evo_m1_input_provenance.csv` (44 rows). Not consumed by any analysis yet. |
| 27 new Evo-M1 items | **Registered and published 2026-11-12.** 27/27 in `__ReadMe.xlsx` (one row each; Item encoded without trailing `_`); item build scripts rerun with 0 failures; 27 public TSVs in `__Public/comparative-data/`, each identical cell by cell to its paper-folder CSV; no orphan TSVs in paper folders. |

### 1. Model diagnostics (done)
- New helper `helpers/model_diagnostics.R`: `diag_stage()`, `diag_lm()`, `diag_facets()`; base graphics only.
  6 panels (resid vs fitted, Q-Q, scale-location, Cook's D with 4/n line, resid vs leverage with 2p/n and Cook contours,
  observed vs fitted), Shapiro-Wilk and Koenker Breusch-Pagan in the title. Skips saturated fits; never stops the host script.
- Outputs: `figs/<stage>/diagnostics/<label>.jpg` + `.pdf`; `tables/<stage>/diagnostics_summary.csv` (upsert by label).
  Stages: s1_synthesis, s1a, s1b, s1c, s1d, s2, s4a, s4b, s5. Calls added after every `lm()` in s1_synthesis, s1a_1, s1b_3 (x4),
  s1b_5, s1b_6, s1b_7, s1c, s1d_1, s2, s4a_endocranial, s4b, s5b, s5c, s5d and inside `helpers/s5_helpers.R::fit_sma_ols` (covers s5a/s5d OLS).
  No regression fits (not touched): s1a_2, s1b_1/2/4/5a–d, s1d_2, s4a_endocranial_cerebellum, s5e.

### 2. Script fixes (done)
- `s4a_endocranial.R`: Figure 3B lookup `f3b$Group_code` -> `f3b$Group_printed` (column renamed upstream).
- `s5c_stress_budget_preservation.R`: reads raw `Mackes_etal_2020_subcortical.csv` and `Walhovd_etal_2011_normative.csv`
  instead of non-existent `*_derived.csv`.
- `00_assemble_comparison_inputs.R` (s1a): added Caudatum (mean of Hanson TD dorsal + medial caudate) and Putamen (associative)
  rows; **all Hanson densities now x1000 (`MM3_PER_CM3`)** — Hanson prints per mm3, every other s1a density is N / cm3.
  The accumbens point was previously 1000-fold too small (explains its Cook's D 11.3 in glia density).
- s5a/s5d need `install.packages("smatr")`; s1a/s1c/s2 labels need `ggpmisc`.
- s4b ran fine in the user's runs; a reviewer flag about an s4b failure came from the scratch runner, not the script.

### 3. Diagnostics triage (decision support; numbers in session artifact `loo_robustness.csv`)
- Keep (main): s1c synaptic density vs rCMRGlc (n = 13, p = 1e-4; worst leave-one-out p 0.007; grey-only n = 12 p = 0.007);
  s5d mammal-only brain MR vs mass.
- Keep as null: s1b E:I ratio (n = 23, p ~0.9).
- Preliminary: s1d MBP/tau (white-matter driven; non-normal residuals), s2 volume change vs rCMRGlc (n = 8, p ~0.02 but
  Spearman 0.06–0.09, LOO max up to 0.18 via putamen), s1b astrocyte/OPC proportions, s5a rate vs volume (flips with white matter),
  s4b calibration (3 params on 7 taxa).
- Drop/demote: s1a p_counts (slide 58), s1a density inference (amygdala-driven; after the units fix, n = 9: glia density
  p = 0.033 but LOO max 0.33 dropping amygdala), s4a fit line (slide 43), s5d fits incl. the ant, s5b (p = 0.08), all multivariate fits,
  s5a log budget vs log volume (built-in coupling).

### 4. Literature search and Evo-M1 inventory (done)
- 104 + 32 (OpenAlex) candidates screened; all DOIs Crossref-verified; only a minority value-checked in full text.
- Already in Evo-M1 and NOT needing new papers: Andersen 2022 supplement (FDG + UCB-J, same subjects, 13 regions);
  Heiss 1991 column in `heiss_2004_regions.csv` (14 regions) and Heiss 1985 Table I (23 structures + whole brain) for x-variable
  replication; West & Gundersen 1990 (hippocampus neurons); Barger 2012/2014 (amygdala neurons); brain masses for
  Microcebus (1.77 g) and pig (112–138 g) in `__merging_brain_mass`.
- Mink 1981: 37–39/42 brain rates are model predictions (repo README: do not merge) — adds no measured species.
- Johansen Table 2 has thalamus; excluded from s1c by design, while s1d and s2 use a Heiss whole-thalamus composite (inconsistent policy).

### 5. New Evo-M1 items (built 2026-11-11; registered and TSVs written 2026-11-12)
27 items in 14 folders, each with snapshot, `.R`, CSV, README, definitions, `<Item>_registry_row.xlsx`; all builds run clean.
- s1a: Salvesen_etal_2015 (Table1, Table4, text), Karlsen_etal_2014 (Table1, Table2), Andersen_etal_1992 (Table4–7),
  Korbo__1999 (Table1, Text; PDF file was mislabelled 2006), Nykjær_etal_2017 (Table2), Wegiel_etal_2014 (Table1–3).
- s1c: Holmes_etal_2024_SupplementaryTable2, Varnäs_etal_2020_Table1, Andersen_etal_2022_TableS1/S2.
- s2: McLaughlin_etal_2014 (Table1, Table4), Sheridan_etal_2022_TableS3, Hodel_etal_2015 (InlineSupplementaryTable3/4).
- s5: Poulsen_etal_1997 (Table1, Table2), Pifferi_etal_2015_ResultsText.
- 23 publication errata recorded per the errata convention (Salvesen 3, Karlsen 2014 7, Wegiel 4, McLaughlin 1, Hodel 3,
  Poulsen/Pifferi 5). Wegiel E001 (dentate corrected volume 511 -> 460) verified against the p. 7 Results text.
- Not in manifest: Wegiel Table1 (stereology parameters), Andersen 2022 TableS2 (p-values).
- Known audit noise: Hodel and Andersen 2022 `.docx` snapshots report "snapshot MISSING" (audit regex lacks .docx).
- Human RA visual check of all snapshots still pending. Sheridan 2022 main-article PDF not in folder.
- The `build-dataset-item` skill text is stale (says header A–AN and E–M blank; live header is A–AL and HOWTO makes M manual).

### Open decisions (needed before wiring the new inputs into analyses)
1. s1a volume basis: Karlsen 2011 printed volumes are shrinkage-corrected, its densities and all new Pakkenberg-lab volumes are
   uncorrected; the current s1a divides by the corrected volume. Recommendation: recompute all densities as N / uncorrected V.
2. s1c: cross-calibrate Holmes BPND onto the Johansen scale (n 13 -> 17) or analyse it as a separate replication (recommended).
3. s2: treat BEIP (McLaughlin 2014, Sheridan 2022) as an independent replication arm rather than pooling with ERA.
4. One whole-thalamus policy across s1c, s1d and s2.

### Next steps (in order)
1. **Done 2026-11-12.** ~~`install.packages(c("smatr", "ggpmisc"))`, then `Rscript scripts/run_all.R`~~ — 25 new inputs synced; s1a outputs
   and diagnostics regenerated with the units fix.
2. **Done 2026-11-12.** ~~Paste the 27 registry rows into `__ReadMe.xlsx`, run `_tools/file_list.R`, verify K/L, rerun each item's `.R`~~ —
   27 rows registered, 27 public TSVs written and verified.
3. Visual check of snapshots; add a note to `metadata/EVO_M1_INPUTS.md` that the new rows are synced but not yet consumed.
4. Decide items 1–4 above, then extend `00_assemble_comparison_inputs.R` (s1a), the s1c/s2 scripts and the s5 species table.
5. Quick win: rerun key s1 fits with the Heiss 1991 rates as an x-variable sensitivity check.

### Session artifacts (Claude Science project, this session)
`diagnostics_all_stages.csv` (v2), `loo_robustness.csv`, `literature_candidates_all_studies.csv`, `openalex_new_candidates.csv`,
`s1a/s1c/s2/s5_*_extracted_values.csv` + notes, `evom1_inventory_reconciliation.csv` (v2), `evom1_new_items_registry_rows.xlsx`.
OpenAlex credential is now configured (Customize -> Credentials).


### 2026-11-12 continuation — implemented choices and primary-source anatomy review

User accepted: uncorrected-volume s1a densities; Holmes as a separate s1c
replication; BEIP as a separate s2 arm. User clarified that s2 asks a different
question from s1 composition/histology; do not require a common s1–s2 region set.

Implemented and run:

- KP 2011 density denominator changed to 2 × author-data hemisphere volume;
  counts unchanged. Five cortical densities increase by factors 1.71–1.90.
  Added `data_intermediate/s1a_volume_basis_audit.csv`. Current s1a remains
  nine regions; the other new stereology sources have not yet been integrated.
- `scripts/s1c_holmes_replication.R`: 12 regions without broad thalamus;
  BPND is separate from Johansen Bmax. OLS p=0.00375, Spearman p=0.00291,
  worst LOO p=0.0126. Mixed regional PVC is a material limitation:
  the PVC-only OLS sensitivity has n=7, p=0.0995.
- `scripts/s2_beip_replication.R`: 8 regions at age 8–10; 7 per adolescent
  contrast without thalamus. Childhood OLS p=0.543; adolescent CAUG p=0.0285
  and FCG p=0.0593, both weakening when the white-matter proxy is excluded
  (p=0.241 and 0.251). Reports share a cohort; no pooling across time points.
- `scripts/s1_s2_regional_sensitivity.R`: per-study thalamus/proxy sensitivities
  and same-region comparisons within research questions. Johansen–proteins
  share 8 regions; Johansen–Holmes share 9. ERA–BEIP share 7 childhood / 6
  adolescent regions. Cross-question s1–s2 intersection fits were removed.
- `scripts/s1_heiss1991_sensitivity.R`: alternate rate-source fits on identical
  complete-case regions, with/without white matter. Johansen: n=10, OLS
  p=0.000187 (2004) vs 0.000101 (1991); without white matter n=9,
  p=0.0375 vs 0.00605, but 1991 Spearman p=0.205. Avoid a blanket robustness claim.
- New scripts explicitly ordered in runner; 14 new output-manifest entries.

Anatomy review (the newest user priority):

- Added `THALAMIC_ANATOMY.md` to 12 existing Evo-M1 paper folders; verified
  written copies against staged notes. No upstream data/registry edits.
- See `metadata/anatomy/THALAMUS_SOURCE_AUDIT.md` for source pages and conclusions.
- Heiss 2004 does not measure a whole-thalamus mean. Its three-territory mean
  is a derived proxy; its medial label is not proven identical to Karlsen MDT.
- Sjölin is an unspecified thalamic tissue sample, NOT established whole-thalamus
  tissue. Checked main Methods, Fig. 1 and supplementary DOCX; nuclear location
  and dissection boundaries remain unknown. Corrected the mapping note.
- Johansen explicitly uses cerebellar cortex in Methods. Corrected its map from
  a whole-cerebellum proxy to a cerebellar-cortex match.
- Whole-thalamus source volumes remain in s2 tables and descriptive plots.
  Current metabolic fits omit the unmatched broad thalamus point and report
  proxy sensitivities. Treat this as a conservative current matching choice,
  not an irreversible final study-design decision or a common-region rule.
- Andersen 2022 uses the same Hammers thalamus VOI for FDG and SV2A: promising
  independent matched-ROI replication. FDG SUVR must not be spliced into Heiss
  absolute rates. Siletti MD+Re aggregation also needs anatomical review.

Validation: 13 affected pipeline scripts completed successfully in
`logs/nov12_analysis_expansion.log`; new replication scripts rerun with immediate
warnings in `logs/nov12_replication_validation.log` (no warnings). New figure PNGs
visually inspected. Final invariants and manifest gate recorded in
`logs/nov12_final_checks.log`.

Remaining work: 21 of the 25 new inputs still only synced; remaining stereology
and s5 expansion, full snapshot human sign-off, Andersen matched-ROI replication,
and resolution of the outstanding thalamic boundaries/rate choice. Earlier rows
in this handoff describe the state before this continuation.

## 12 Nov continuation: stricter anatomy and metabolic convergence

User explicitly requested cleaning **all studies** to use more precise anatomy,
even at the cost of dropping regions; alternate metabolic sources should supply
separate convergent analyses. User retains s2 as its own regional-volume question.
The operative policy is now `metadata/anatomy/PRIMARY_MATCH_POLICY.md`; earlier
thalamus-only decisions/counts above are superseded.

- Added inclusion decisions/reasons to nine study maps. Rejected primary mappings
  retain source measurements but receive no metabolic rate or SD. s1b exclusions
  act before pooling cells; s1a parent cortex is descriptive only alongside lobes.
- Primary region counts: s1a neuron/glia densities 8 (counts/glial subtypes 5),
  s1b 19 mapped groups, Johansen 12, Holmes 11, proteins 7, ERA 7,
  BEIP childhood 7 / adolescence 6 per contrast, s5a 8. Shared regions:
  Johansen/proteins 6; Johansen/Holmes 8; ERA/BEIP childhood 7 / adolescence 6.
  Counts are regional observations, not donors or participants.
- s4a scripts and products renamed to `*_proxy_sensitivity` / matching directories.
  Endocasts do not provide the PET grey-matter compartment; no strict primary
  metabolic estimate is claimed. s4b input paths and cross-check columns now label
  these proxy sensitivities. Primary arterial models remain intact. The network
  script depends on archived s3 anatomy and is removed from the active runner.
- s5a retains only 3.23% of brain volume. s5b therefore uses the published
  Clarke & Sokoloff whole-brain glucose rate via Karbowski (31 umol/100g/min),
  never a whole-brain extrapolation of the sparse s5a regional mean. s5c inherits
  seven strict s2 regions. s5a rate plot now uses log y to agree with its log-log fit.
- **Citation correction:** the supposed "1991" comparison column is actually
  Heiss **1984**, reference 21 of Heiss 2004 (confirmed in primary references).
  Corrected script, column and output names to 1984. Values remain the 2004
  re-tabulation, accounting for lumped constants, not newly extracted native
  1984 Table 3. Old broad thalamus/midbrain values are not exact nuclear matches.
- Added `s1c_andersen_paired_replication.R`: 13 same-atlas FDG/UCB-J regions per
  age group (11 young, 15 older adults), native GM-normalized uptake; includes
  thalamus without Heiss substitution. Young OLS p=.00749, old p=.0683;
  omitting SN gives p=.0535/.268 respectively. Rank/LOO checks retained. Shared
  normalization and regional dependence limit inference; no across-age pooling.
- Holmes strict OLS p=.00593, PVC-only n6 p=.358. BEIP strict OLS p=.580
  childhood, .241 CAUG, .251 FCG adolescent: no blanket replication claim.
  Johansen 1984-vs-2004 matched n9: OLS p=.00605/.0375, but older-source
  rank p=.205 and current-source worst LOO p=.208. Report fragility.
- Heiss1984 full native extraction is a future useful expansion (cingulate,
  lobar WM, pons, mixed midbrain, broad thalamus). Keep its own original anatomy
  and model/units. Publisher PDF fetch returned403, but primary article text was
  readable. Heiss1985 awake subjects were sleep-deprived after ~28h on duty:
  useful state sensitivity, weak choice for rested normative primary rates.
- Evo-M1: corrected Heiss2004 THALAMIC_ANATOMY.md; added ANATOMICAL_MATCHING.md
  in Heiss1984 (new notes-only folder, no local PDF), Heiss1985, Hanson2018,
  Morgan2014, Karlsen_Pakkenberg2011, Andersen2022 and Sjolin2022. Existing notes
  were checked before writing. No source measurement or registry was altered.
- Validation: 27 affected scripts reran successfully; both renamed s4a scripts,
  s4b, new mapping audit and Andersen producer also ran successfully. Both
  scientific invariant suites pass; active scripts/helpers parse; output gate
  **43 OK, 2 longstanding known gaps, no missing files/scripts**. Figure QA:
  Andersen readable; s5 annotation collision fixed. Diff whitespace check clean
  (Git LFS disabled for that read-only diff because sandbox prevents LFS temp writes).
  Logs: `logs/strict_anatomy_status.csv`, `strict_anatomy_completion_checks.log`.
- Comparison exports: three updated CSVs staged in `tables/comparison_export/`.
  Auto-review rejected copying them to sibling
  `analyses_metabol_rate_structure_nonhuman/comparison/inputs`, saying the user
  had not authorized this specific cross-project destination. An asynchronous
  approval request is pending; sibling files remain unchanged unless a later
  entry records approval and verified transfer. Do not bypass the rejection.

The 44 synced source snapshots remain original inputs. Of the 25 new synced
inputs, five are now consumed (Holmes, two McLaughlin tables, Sheridan, Andersen).
Remaining stereology and broader source ingestion work is not completed by this
anatomical cleanup. No human visual sign-off, author contact, commit or push.

### Sibling export approval and completion — 12 November 2026

The user explicitly approved updating the sibling inputs. Replaced
`human_comparison_rows.csv`, `human_regional_change_vs_rate.csv`, and
`human_scaling_fits.csv` in
`../analyses_metabol_rate_structure_nonhuman/comparison/inputs/` with the validated
local exports. All three destination files were SHA-256 verified against their
local sources. The exported thalamus rows retain volumes but have no metabolic
rate and are excluded from regression. The earlier pending transfer is resolved;
no downstream nonhuman analyses were rerun as part of this input refresh.

### Borghammer Table 2 feasibility — 12 November 2026

User asked whether Borghammer Table 2 can serve as convergence. Checked the
Borghammer2012 primary PDF Methods/Fig.2, Table1, Table2, source extraction and
TableS1 metadata, plus native Siletti dissection/donor annotations.

Verdict: useful separate normalized-rate s1b check, not a whole-thalamus solution.
Six candidate native ROI matches (Pu, CaB, GPe, GPi, VA, pure CM); VA/CM each
one donor, others three. CM-Pf is mixed, VLN boundary equivalence unresolved,
STH says nucleus and nearby (two donors), so sensitivity-only pending review.
This newly inspected STH wording also warrants revisiting its earlier Heiss
assignment; the feasibility task did not silently refit the primary analyses.
Borghammer GM = global mean; Cerebellum = central cerebellar white matter.
Controls n11; pooled absolute rates rejected by authors because of center effects.
Three normalizations are dependent; eight subcortical bilateral control ranks
agree exactly. s2/Johansen have only caudate/putamen, Holmes alsoSTN; no meaningful
broad-region replication from those tiny overlaps. Do not collapse GPe/GPi or
VA/VL/CM to manufacture broad matches, or de-normalize from TableS1 means.

Saved assessment `metadata/anatomy/BORGHAMMER_TABLE2_FEASIBILITY.md` and standalone
reproducible audit `checks/assess_borghammer_table2.R`, with two outputs in
`tables/anatomy/borghammer_*`. Audit ran successfully; source measurements,
primary fits, manifests and sibling exports were unchanged. See
`logs/borghammer_table2_feasibility.log`.

Added and byte-verified `Borghammer_etal_2012/ANATOMICAL_MATCHING.md` upstream.
Feasibility outputs checked: 30 structure-by-normalization summaries and 18 ROI-by-dataset audit rows.

### Dataset by region coverage matrix — 12 November 2026

Created a five-paper metabolic inventory (Heiss 2004/1984/1985, Borghammer 2012,
Andersen 2022), not an exhaustive Heiss bibliography. The matrix has 70 distinct
territories, 350 cells and 96 reported cells. Separate codes identify named ROIs,
spatial samples, uncertain anatomy, published-but-untranscribed means, tests-only
coverage and absent measurements. It does not declare cross-paper exact matches.

Sources contribute 26/22/24/10/14 reported territories respectively. Andersen's
fourteenth territory is HIP in S2 tests only (13 S1 means). SN is explicitly
flagged uncertain because its expansion is not defined in the publication.
1984 lentiform is kept separate from putamen, and native Table 3 coverage is
kept separate from the 2004 re-tabulation. No primary model was changed.

See `metadata/anatomy/DATASET_REGION_COVERAGE.md`, reproducible builders
`checks/build_metabolic_coverage_matrix.py` and
`checks/build_metabolic_coverage_workbook.mjs`, three CSVs under
`tables/anatomy/metabolic_dataset_*`, and the workbook under
`outputs/01a1d166-5e50-7ef2-a6aa-ee24c6b4b692/`.


## 14 November 2026 — Tüshaus, Johansen sampling, BEIP

Completed the requested three arms; detailed methods and reproduction commands
are in `metadata/anatomy/NOV14_PROTEOMICS_ATLAS.md`.

- Added `scripts/s1d_2a_extract_tushaus_proteomics.R` and
  `scripts/s1d_2b_tushaus_protein_metabolic_analysis.R`. Downloaded publisher
  Dataset EV1 and extracted 608,986 measurements from 38 tissue samples across
  13 regions in one donor. The strict shared anatomy map retains seven regions.
  Complete seven-region screens cover 6,810 LFQ and 6,797 iBAQ protein IDs.
  Published log scales are retained; regional enrichment classes are not cell
  types or cell fractions. Preserved and flagged duplicate SN replicate labels.
  No six-marker panel association survives BH correction. Protein screening is
  exploratory and cannot support donor-population inference.
- Added `scripts/s1c_sample_johansen_atlas.R` and
  `helpers/sample_johansen_atlas.py`. Sampled ten missing territories with public
  MASSP2/AAL3/Diedrichsen masks; alternative AAL3/CIT168 masks and ±1-mm shifts
  quantify sensitivity. Eight new nuclei carry small-nucleus flags. Vermis is
  descriptive only. Global cortex, medial thalamic nuclei, basal forebrain and
  centrum semiovale remain unresolved. All new estimates are excluded from
  primary analyses pending template validation.
- Critical atlas check: five overlapping atlas-derived regional means are
  30–65% below Table 2. Source smoothing/method differences and MNI template
  correspondence need review; no causal explanation is established. Final
  outputs compare 12 published regions and nine new atlas regions separately.
  No pooled 21-region regression is retained. The 21-row inventory explicitly
  forbids pooling for inference. Existing Johansen primary analysis is unchanged.
- BEIP was already implemented on Nov 12. Reran
  `scripts/s2_beip_replication.R`: separate childhood EIG–NIG (7 regions) and
  adolescent CAUG–NIG / FCG–NIG (6 each). OLS p = 0.580 / 0.241 / 0.251.
  These share a longitudinal cohort and adolescent controls; remain separate
  from ERA. No whole-thalamus proxy enters the fits.

All additions are ordered in `scripts/run_all.R` and registered in
`metadata/output_manifest.csv`. Inputs have download URLs and SHA-256 checksums
in `metadata/nov14_source_downloads.csv`; `checks/fetch_nov14_inputs.R` verifies
all 13 cached inputs or downloads missing ones. Python dependencies are isolated
in ignored `.venv-atlas`; requirements are in `helpers/requirements-atlas.txt`.
New source-folder READMEs document anatomical restrictions and atlas attribution.

Validation: extraction, protein analysis, atlas sampler, BEIP, and shared anatomy
audit executed successfully. `checks/test_nov14_proteomics_atlas.R`,
`checks/test_strict_anatomical_matching.R`, and
`checks/test_nov12_analysis_expansion.R` passed. Manifest validation: 64 OK, no
missing files/scripts; two pre-existing gaps remain (s4a budget residual and
archived s3 mosaicism). New protein, atlas, and BEIP figures were visually
inspected. Logs are under `logs/nov14_*`. No full pipeline rerun was needed.
No sibling-project exports changed in this session; earlier working-tree
changes were preserved. No commit or push performed.
