# =====================================================================
# run_all.R  --  run every analysis script in scripts/ in dependency order
# ---------------------------------------------------------------------
# NOTE: this runner is TEMPORARY -- it exists to build/debug the pipeline
# and will NOT be part of the published version.
#
# Usage:
#   From a terminal (project root or anywhere):
#       Rscript scripts/run_all.R
#   Or inside R / RStudio:
#       source("scripts/run_all.R")
#
# What it does
#   * Runs the scripts in the order the pipeline needs (raw-data prep ->
#     Study 1 cells -> Study 2 stress -> Study 3 evolution -> Study 4a fossils
#     -> Study 4b arterial canal -> QC checks).
#   * Each script runs in its OWN environment so leftover objects from one
#     script cannot silently feed the next.
#   * A failing script is caught, logged, and the run CONTINUES; a summary
#     table (and logs/run_all_status.csv) is printed at the end.
#   * After the summary, scripts/verify_outputs.R is sourced as a GATE: it
#     checks the deliverable manifest and exits non-zero on failure.
#
# Notes
#   * s4a_endocranial.R must run before s4a_endocranial_cerebellum.R (the
#     cerebellum script reads the whole-brain budget table the main one writes).
#   * KNOWN ISSUE (Study 3): the two s3_predicValuesPGLS_* scripts are forks
#     that BOTH write to figs/s3/all, tables/s3/all and checks/s3/all, so the
#     later one (VOLUMES_WIDE_SELECT) overwrites the earlier one's outputs.
#     A single engine + per-config drivers is the planned fix; until then the
#     order below preserves the outputs produced by the last logged run.
#   * The working directory is reset to the project root before every script,
#     because the scripts use paths relative to the project root
#     (e.g. "data_raw/...", "figs/s4a/...", "helpers/plot_settings.R").
# =====================================================================

# ---- locate the project root (parent of this scripts/ folder) -------
.this_file <- local({
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) return(normalizePath(f))
  # sourced interactively: fall back to ofile, else guess from getwd()
  of <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(of)) return(normalizePath(of))
  if (basename(getwd()) == "scripts") return(file.path(getwd(), "run_all.R"))
  file.path(getwd(), "scripts", "run_all.R")
})
SCRIPT_DIR <- dirname(.this_file)
ROOT       <- normalizePath(file.path(SCRIPT_DIR, ".."))
cat("Project root :", ROOT, "\n")
cat("Scripts dir  :", SCRIPT_DIR, "\n\n")

# ---- pipeline order --------------------------------------------------
RUN_ORDER <- c(
  # 0. raw-data preparation (Stephan volumes + Heiss/Stephan reference table)
  "0_bind_matano_1985a_to_stephan.R",
  "0_Heiss_Stephan_and_table1_30052026.R",

  # background comparative plot (slide 8)
  "traits_neocortex_grey_white.R",

  # Study 1a: stereology
  "s1a_1_stereology_cell_types_30052026.R",
  "s1a_2_stereology_proportions_30052026.R",

  # Study 1b: transcriptomic cells (extract -> map -> proportions -> analyses)
  "s1b_1_extract_transcriptomic_30052026.R",
  "s1b_2_mapping_rcmrglc_transcriptomic_cells_anatomy_21052026.R",
  "s1b_2_check_dissection_roi.R",
  "s1b_3_n_transcriptomic_neuronal_25052026.R",
  "s1b_3_n_transcriptomic_neuronal_telencephalon_25052026.R",
  "s1b_3_nn_transcriptomic_nonneuronal_25052026.R",
  "s1b_3_nn_transcriptomic_nonneuronal_telencephalon_25052026.R",
  "s1b_4_n_supercluster_rcmr_correlation_matrix_telencephalon_13062026.R",
  "s1b_5_n_EI_ratio_original_vs_jorstad_overlay_two_MSN_plots_raw_EI_only_16062026.R",
  "s1b_6_nn_type1_type2_astrocyte_compositional_rcmr_26052026.R",

  # Study 1c: synaptic density (SV2A) vs rCMRGlc
  "s1c_synaptic_density_metabolic_rate.R",

  # Study 2: environmental stress
  "s2_stress_volume_01062026.R",

  # Study 3: data prep must precede the PGLS forks -- the prep script writes
  # the species list the tree-sample script reads, and that writes the 100-tree
  # sample the PGLS scripts fit across.
  "s3_prepare_Stephan_primate_data.R",
  "s3_prepare_tree_sample.R",

  # Study 3: evolutionary deviation (PGLS). See KNOWN ISSUE in the header:
  # both forks write the same output dirs; VOLUMES_WIDE_SELECT runs last so
  # its outputs win, matching the last logged run.
  "s3_predicValuesPGLS_16062026_2deg_neocortex_cerebellum_blue_square_UPDATED_MU_PLOTS_1_PATCHED.R",
  "s3_predicValuesPGLS_16062026_2deg_neocortex_cerebellum_blue_square_VOLUMES_WIDE_SELECT.R",

  # Study 4a: fossil endocranial budgets (main before cerebellum sibling)
  "s4a_endocranial.R",
  "s4a_endocranial_cerebellum.R",

  # Study 4b: arterial-canal metabolic estimates
  "s4b_arterial_canal.R",

  # QC / robustness checks (no deliverable outputs)
  "network_residual_autocorrelation_analysis.R",

  # QC: Stephan reference sheet + raw/model comparison against
  # volumes_wide_select. The compare script's "model" mode reads
  # checks/*_data_used.csv written by the s3 scripts, so it runs after them.
  "qc_stephan/build_stephan_primates_reference_sheet.R",
  "qc_stephan/compare_stephan_vs_volumes_wide_merged_v3.R"
)

# Scripts that live in scripts/ but are NOT run here
DO_NOT_RUN <- c(
  "run_all.R",          # this file
  "verify_outputs.R"    # sourced as a gate after the run, below
)

# ---- append any script not already accounted for --------------------
# recursive so subfolder scripts (e.g. qc_stephan/) are discovered; archived
# scripts are never run.
all_scripts <- list.files(SCRIPT_DIR, pattern = "\\.R$", full.names = FALSE,
                          recursive = TRUE)
all_scripts <- all_scripts[!startsWith(all_scripts, "archive/")]
known <- c(RUN_ORDER, DO_NOT_RUN)
extra <- setdiff(all_scripts, known)
if (length(extra)) {
  cat("NOTE: these scripts are not in the known order and will run last:\n  ",
      paste(extra, collapse = "\n   "), "\n\n", sep = "")
  RUN_ORDER <- c(RUN_ORDER, extra)
}

# warn about any listed script that is missing from disk
missing <- setdiff(RUN_ORDER, all_scripts)
if (length(missing))
  cat("WARNING: listed but not found (will be skipped):\n  ",
      paste(missing, collapse = "\n   "), "\n\n", sep = "")
RUN_ORDER <- intersect(RUN_ORDER, all_scripts)

# ---- run ------------------------------------------------------------
status <- data.frame(step = integer(), script = character(),
                     result = character(), seconds = numeric(),
                     message = character(), stringsAsFactors = FALSE)

for (i in seq_along(RUN_ORDER)) {
  scr  <- RUN_ORDER[i]
  path <- file.path(SCRIPT_DIR, scr)
  cat(sprintf("\n[%2d/%2d] %s\n", i, length(RUN_ORDER), scr))
  cat(strrep("-", 70), "\n", sep = "")
  setwd(ROOT)                                  # scripts use root-relative paths
  t0  <- Sys.time()
  msg <- ""
  res <- "OK"
  tryCatch(
    sys.source(path, envir = new.env(parent = globalenv()), keep.source = FALSE),
    error = function(e) { res <<- "FAILED"; msg <<- conditionMessage(e) }
  )
  secs <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  if (res == "FAILED") cat("   >>> FAILED:", msg, "\n")
  cat(sprintf("   [%s] %.1f s\n", res, secs))
  status <- rbind(status, data.frame(step = i, script = scr, result = res,
                                      seconds = secs, message = msg,
                                      stringsAsFactors = FALSE))
}

# ---- summary --------------------------------------------------------
setwd(ROOT)
dir.create("logs", showWarnings = FALSE)
log_path <- file.path("logs",
                      format(Sys.time(), "run_all_status_%Y%m%d_%H%M%S.csv"))
write.csv(status, log_path, row.names = FALSE)
write.csv(status, file.path("logs", "run_all_status_latest.csv"), row.names = FALSE)

cat("\n", strrep("=", 70), "\n", sep = "")
cat("RUN SUMMARY  (", sum(status$result == "OK"), " OK / ",
    sum(status$result == "FAILED"), " failed of ", nrow(status), ")\n", sep = "")
cat(strrep("=", 70), "\n", sep = "")
print(status[, c("step", "script", "result", "seconds")], row.names = FALSE)
cat("\nStatus log written to:", log_path, "\n")

failed <- status$script[status$result == "FAILED"]
if (length(failed)) {
  cat("\nFAILED scripts:\n  ", paste(failed, collapse = "\n   "), "\n", sep = "")
  # non-zero exit code so shell/CI can detect failures (Rscript only)
  if (!interactive()) quit(status = 1, save = "no")
} else {
  cat("\nAll scripts completed without error.\n")
}

# ---- gate: verify deliverable outputs against the manifest ----------
# Runs after the summary so the log above is always written. verify_outputs.R
# itself exits non-zero if a deliverable file/script is missing, which makes
# the whole Rscript invocation fail -- that is the intended gate behaviour.
cat("\n", strrep("=", 70), "\n", sep = "")
cat("GATE: verify_outputs.R\n")
cat(strrep("=", 70), "\n", sep = "")
setwd(ROOT)
sys.source(file.path(SCRIPT_DIR, "verify_outputs.R"),
           envir = new.env(parent = globalenv()), keep.source = FALSE)
