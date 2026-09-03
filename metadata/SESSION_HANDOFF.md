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
