analyses_metabol_rate_structure/
├── README.md
├── .gitignore
├── scripts/            # analysis pipeline, run via scripts/run_all.R
│   ├── archive/        # superseded scripts kept for reference (not run)
│   └── qc_stephan/     # QC scripts for the Stephan volume compilation
├── R/                  # shared helpers (plot_settings.R, project_root.R)
├── data_raw/           # source data as published/digitised
├── data_intermediate/  # derived tables passed between scripts
├── data_analysis/      # Study 1b result tables (destination of the s1b_* scripts)
├── metadata/           # provenance, correction log, audit plans, output manifest
├── checks/             # QC / audit outputs (no deliverables)
├── tables/             # deliverable tables, one subfolder per study
├── figs/               # deliverable figures, one subfolder per study
├── network_residual_autocorrelation/  # outputs of the network QC script
└── logs/               # pipeline run logs (git-ignored)

## Running the pipeline

    Rscript scripts/run_all.R

Runs every analysis script in dependency order, then sources
`scripts/verify_outputs.R` as a gate against `metadata/output_manifest.csv`.
See `metadata/SESSION_HANDOFF.md` for current status and open items.

## Data corrections

See `metadata/CORRECTIONS_LOG.md` for the full log.

- **2026-06-03 — Stephan_primates.csv, `Telencephalon`, _Alouatta seniculus_.**
  Value was `37.39` (≈1000× too small). True value from the source table
  (Stephan et al. 1981) is **37,388**; an OCR step misread the thousands-comma
  as a decimal point (`37,388` → `37.388`), later rounded to `37.39`. Detected
  because the five fundamental brain parts did not sum to `Total_brain_net_volume`
  and the neocortex (`NeoWG` = 31,660) exceeded the stated telencephalon.
  Corrected to `37388` in `data_raw/Stephan_primates.csv` and in all 11 archived
  copies. Also removed a stray empty column from `data_raw/Stephan_primates.csv`
  (now 60 × 69). Re-run derivation scripts so downstream files pick up the fix.

## 14 November 2026 additions

- `scripts/s1d_2a_extract_tushaus_proteomics.R` extracts Tüshaus Dataset EV1;
  `scripts/s1d_2b_tushaus_protein_metabolic_analysis.R` analyzes the same seven
  matched regions across proteins. These are protein abundances from one donor.
- `scripts/s1c_sample_johansen_atlas.R` samples ten missing territories and flags
  small nuclei. Nine have named anatomical matches; vermis remains descriptive.
  New estimates remain exploratory: template correspondence is unvalidated and
  five overlap regions differ by 30–65% from Table 2. Table/map values are not pooled.
- `scripts/s2_beip_replication.R` supplies the separate BEIP arm (7/6/6 regional
  points), rerun and checked against the strict anatomical policy.

See [methods, limitations, and reproduction commands](metadata/anatomy/NOV14_PROTEOMICS_ATLAS.md).
New raw inputs have URLs and SHA-256 checksums in
`metadata/nov14_source_downloads.csv`; `checks/fetch_nov14_inputs.R` acquires or
verifies them. Atlas sampling requires the Python packages in
`helpers/requirements-atlas.txt` (local `.venv-atlas` supported).
