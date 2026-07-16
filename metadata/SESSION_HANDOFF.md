# Session handoff — repo audit & reproducibility work

**Updated:** 2026-07-14
**Repo:** `analyses_metabol_rate_structure`
**Deliverables this repo feeds:**
`Brain energetics/Energetic Constraints on Brain Organization in Human Paleoneurology_11072026.pptx` (69 slides)
and `Brain energetics/MS Is human brain organization economical_04072026_TRACKED.docx`.
Analyses are organized as four studies: s1 (cellular composition), s2 (environmental
stress), s3 (evolutionary deviation, PGLS), s4 (fossil endocranial budgets).

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
- `s4_endocranial.R` → PNG+PDF for `endocranial_region_budgets`,
  `endocranial_region_cost_pie_3species_volctx`, `volume_timeline_logage_flipped`
  (via a `save_png_pdf()` helper).
- `s4_endocranial_cerebellum.R` → PNG+PDF for `cerebellum_volume_timeline_logage`,
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
  `rcmrglc_vs_fig3a_deviation.png` / `_fit.png` (canonical copies remain in `figs/s4/`).
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
   original code or agree a definition, then it can be added to `s4_endocranial.R`
   with PNG+PDF export.
4. **Clean-room reproducibility test + dead-output prune (Part 1).** Move `figs/`,
   `tables/`, `data_intermediate/` aside, run `run_all.R`, confirm manifest outputs
   reappear. Then list the ~448 files in `figs/` not in the manifest and not consumed
   downstream → archive/delete.
5. **Naming & organization (Part 2, do after the manifest is trusted).**
   - Move QC scripts to `scripts/checks/` with `check_`/`qc_` prefixes
     (`s1b_2_check_dissection_roi`, `s1b_x_check_dissection_roi`,
     `s3_0_missingness_clade_diagnostic`, `s3_compare_stephan_vs_merged`,
     `network_residual_autocorrelation_analysis`).
   - Move the s3 engine to `R/` (e.g. `R/pgls_engine.R`).
   - Drop dated/`PATCHED` suffixes → one canonical name per analysis.
   - Resolve shared ownership: `data_intermediate/Heiss_Stephan_data.csv` is written
     by `0_Heiss…` AND appended by `s1b_2`/`s1b_3` — give a single `prep_` script
     ownership.
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
- Modified: `traits_neocortex_grey_white.R`, `s4_endocranial.R`,
  `s4_endocranial_cerebellum.R`, 27 scripts (`setwd` → root finder), `.gitignore`.
- Deleted: 22 `Rplot*.png`, root `rcmrglc_vs_fig3a_deviation{,_fit}.png`.

## Environment caveats

- R is not installed in the Cowork sandbox; scripts were checked statically only.
- The repo is Dropbox-synced; deleting files required elevated permission and git
  index locking was flaky — prefer running git yourself.
