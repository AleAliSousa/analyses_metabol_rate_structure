# Repository audit & cleanup plan

**Created:** 2026-07-13
**Scope:** the `analyses_metabol_rate_structure` repo — reproducibility of
deliverable figures/tables, script naming & organization, and file hygiene.

This plan has three parts, each with a goal, a concrete method, and a
"make it a standing check" step so the audit can be re-run, not done once.

---

## Part 1 — Reproducibility: every deliverable fig/table regenerates from a script

**Goal.** For every figure or table that appears in the presentation or the
manuscript, exactly one script produces it, and running the pipeline from clean
inputs recreates it.

**Why it matters.** Right now most deck figures map to a script, but at least one
(`figs/s4/budget_residual_vs_volume.png`) has no producer, and several figures
exist in `figs/` that never reach a deliverable. Neither problem is visible
without a manifest.

**Method.**

1. **Build an output manifest** — `metadata/output_manifest.csv` with one row per
   deliverable output:
   `output_path, kind (fig/table), producing_script, consumer (deck slide / MS figure or table), status`.
   Seed it from the deck-image match already done and from the MS figure/table
   list. This is the single source of truth for "relevant" outputs.

2. **Producer check (static).** For each manifest row, confirm a script writes
   that exact path (grep for the filename stem across `scripts/`). Flag rows with
   no producer — currently `budget_residual_vs_volume.png`.

3. **Clean-room regeneration (dynamic, the real test).** Move `figs/`, `tables/`,
   and `data_intermediate/` aside into a temp folder, run `scripts/run_all.R`,
   then confirm every manifest path reappears. This proves the outputs come from
   code and raw data alone, not from stale files on disk. (This is also how the
   stray root copies in Part 3 would have been caught.)

4. **Dead-output check.** List files in `figs/` and `tables/` that are NOT in the
   manifest and are not consumed by another script — candidates for deletion or
   archiving. (There are 448 files in `figs/` today; most are per-region /
   per-config variants, so expect a long tail.)

5. **Standing check.** Add `checks/verify_outputs.R` that reads the manifest and
   asserts each output exists (and optionally is newer than its raw inputs) after
   a run. Wire it as the last step of `run_all.R` or run it in CI.

**Known gap to resolve first:** `budget_residual_vs_volume.png` — no source in
repo or history; needs the original code or an agreed definition of the residual.

---

## Part 2 — Script naming & organization

**Goal.** A reader can tell, from the name and location alone, what kind of thing
a script is, what it depends on, and whether it belongs in the deliverable
pipeline.

**Five categories** (the current repo mixes these under flat `scripts/`):

1. **Shared data-prep** — builds inputs consumed by **more than one** study
   (e.g. `Heiss_Stephan_data.csv`, the transcriptomic `.rds` extracts). These are
   NOT owned by any single study, so they must not carry a study prefix.
   *Convention:* `prep_*` (or numeric `00_*`), outputs only to
   `data_intermediate/`. Example: `prep_heiss_stephan_reference.R`,
   `prep_transcriptomic_extract_neuronal.R`.
   *Fix to flag:* today `Heiss_Stephan_data.csv` is written by `0_Heiss...` **and**
   appended by `s1b_2`/`s1b_3` — shared ownership across study scripts blurs the
   boundary. Decide on one owner (a `prep_` script) or split into clearly-named
   layers.

2. **Per-study analysis** — produces that study's figs/tables. Keep the working
   `s1a_ / s1b_ / s2_ / s3_ / s4_` prefixes plus a step number for run order.
   Outputs to `figs/<study>/` and `tables/<study>/` only.

3. **Checks / QC** — validate data or assumptions, produce **no deliverable**
   (e.g. `s1b_2_check_dissection_roi.R`, `s3_0_missingness_clade_diagnostic`,
   `s3_compare_stephan_vs_merged`, `network_residual_autocorrelation_analysis`).
   *Convention:* `check_*` / `qc_*`, kept in a `scripts/checks/` subfolder, writing
   only to `checks/`. Excluded from the deliverable pipeline; run as a separate QC
   pass.

4. **Engines / reusable code** — sourced by other scripts, never run standalone
   (e.g. `s3_predicValuesPGLS_MERGED_variant.R`, `R/plot_settings.R`).
   *Convention:* live in `R/`, named for what they provide, no side effects on
   source. Move the s3 engine to `R/` (e.g. `R/pgls_engine.R`) so it's obviously
   not a runnable step.

5. **Drivers** — orchestrate an engine across configurations
   (`s3_run_and_compare_configs.R`, `s3_run_frontal_vermis_anthro.R`).
   *Convention:* `*_run_*` prefix. These ARE pipeline steps.

**Cross-cutting rules.**

- No dates / `UPDATED` / `PATCHED` in filenames. Git holds history; superseded
  scripts go to `archive/`. (You already started this — keep it up.)
- One canonical script per output. If a method has variants, express them as a
  config passed to one engine (the s3 driver pattern), not as copied files.
- Encode run order with numeric prefixes within each study.
- Every script starts with a header comment stating: **role** (prep / analysis /
  check / engine / driver), **inputs**, **outputs**, and **who consumes** the
  outputs.

---

## Part 3 — File hygiene & output locations

**Goal.** Outputs land in predictable places; nothing is written to the repo root;
the working directory is handled the same way everywhere.

**What's wrong now (diagnosed).**

- `Rplot001.png` (root) — R's default graphics device file, produced when
  plotting code is run line-by-line without an open device. An interactive
  artifact, not a real output.
- `rcmrglc_vs_fig3a_deviation.png` and `..._fit.png` (root) — **stale** copies;
  `s4_endocranial.R` now writes both to `figs/s4/`, and the root versions differ
  (older). `..._fit.png` is committed to git.
- **Root cause:** inconsistent project-root handling — 23 script locations use
  `here()`, but ~several still hardcode
  `setwd("~/Library/CloudStorage/Dropbox/.../analyses_metabol_rate_structure")`.
  When a script is run with the wd not at root (or with old code), relative paths
  resolve to the wrong place and files land in root.

**Fixes.**

1. **Standardize the working directory.** Add an `.Rproj` at the repo root and
   have every script resolve paths with `here::here(...)` instead of `setwd(...)`.
   Remove the hardcoded Dropbox `setwd()` lines (they break on any other machine
   and are the source of misplaced files).
2. **Output-location rule.** No script writes to the repo root. Figures →
   `figs/<study>/`, tables → `tables/<study>/`, intermediates →
   `data_intermediate/`. Add this to each script header.
3. **Clean up existing strays.** Delete `Rplot001.png`; `git rm` the tracked
   `rcmrglc_vs_fig3a_deviation_fit.png` in root and delete the root
   `rcmrglc_vs_fig3a_deviation.png` (the real ones live in `figs/s4/`). Confirm
   via Part 1 clean-room rerun that nothing regenerates them in root.
4. **Extend `.gitignore`.** Add `Rplot*.png`, `/*.png`, `/*.csv` at the repo root
   (top-level only), and `.RData`, so interactive artifacts never get committed.

---

## Suggested order of work

1. Part 3 cleanup + standardize wd/`here()` — quick, removes noise, and makes the
   clean-room test in Part 1 trustworthy.
2. Part 1 manifest + clean-room regeneration — establishes the reproducibility
   baseline and surfaces every gap and dead output at once.
3. Part 2 renaming/reorganization — do last, once the manifest tells you exactly
   which scripts are deliverable-producers vs checks vs engines, so renames are
   informed and low-risk. Rename in one pass with path updates across active
   scripts (36 files reference `Stephan_primates.csv` alone — batch it).
