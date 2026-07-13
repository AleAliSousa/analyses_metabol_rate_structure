# =====================================================================
# run_all.R  --  run every analysis script in scripts/ in dependency order
# ---------------------------------------------------------------------
# Usage:
#   From a terminal (project root or anywhere):
#       Rscript scripts/run_all.R
#   Or inside R / RStudio:
#       source("scripts/run_all.R")
#
# What it does
#   * Runs the scripts in the order the pipeline needs (raw-data prep ->
#     Study 1 cells -> Study 2 stress -> Study 3 evolution -> Study 4 fossils).
#   * Each script runs in its OWN environment so leftover objects from one
#     script cannot silently feed the next.
#   * A failing script is caught, logged, and the run CONTINUES; a summary
#     table (and logs/run_all_status.csv) is printed at the end.
#
# Notes
#   * s3_predicValuesPGLS_MERGED_variant.R is the Study-3 ENGINE. It is NOT
#     run directly -- it is source()d by s3_run_and_compare_configs.R and
#     s3_run_frontal_vermis_anthro.R, which set the per-configuration options
#     first. It therefore appears in ENGINE_ONLY below, not in RUN_ORDER.
#   * s4_endocranial.R must run before s4_endocranial_cerebellum.R (the
#     cerebellum script reads the whole-brain budget table the main one writes).
#   * The working directory is reset to the project root before every script,
#     because the scripts use paths relative to the project root
#     (e.g. "data_raw/...", "figs/s4/...", "R/plot_settings.R").
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
  "v2_0_Heiss_Stephan_and_table1_05062026.R",

  # background comparative plot (slide 8)
  "traits_neocortex_grey_white.R",

  # Study 1a: stereology
  "s1a_1_stereology_cell_types_30052026.R",
  "s1a_2_stereology_proportions_30052026.R",

  # Study 1b: transcriptomic cells (extract -> map -> proportions -> analyses)
  "s1b_1_n_extract_transcriptomic_neuronal_30052026.R",
  "s1b_1_nn_extract_transcriptomic_nonneuronal_30052026.R",
  "s1b_2_mapping_rcmrglc_transcriptomic_cells_anatomy_21052026.R",
  "s1b_2_check_dissection_roi.R",
  "s1b_x_check_dissection_roi.R",
  "s1b_3_n_transcriptomic_neuronal_25052026.R",
  "s1b_3_n_transcriptomic_neuronal_telencephalon_25052026.R",
  "s1b_3_nn_transcriptomic_nonneuronal_25052026.R",
  "s1b_3_nn_transcriptomic_nonneuronal_telencephalon_25052026.R",
  "s1b_4_n_supercluster_rcmr_correlation_matrix_telencephalon_13062026.R",
  "s1b_5_n_EI_ratio_original_vs_jorstad_overlay_two_MSN_plots_raw_EI_only_16062026.R",
  "s1b_5_n_EI_ratio_telencephalon_26052026.R",
  "s1b_6_nn_type1_type2_astrocyte_compositional_rcmr_26052026.R",

  # Study 2: environmental stress
  "s2_stress_volume_01062026.R",

  # Study 3: evolutionary deviation (diagnostics, then the config drivers that
  # source the engine)
  "s3_0_missingness_clade_diagnostic_04062026.R",
  "s3_1_phylo_multiple_imputation_04062026.R",
  "s3_compare_stephan_vs_merged.R",
  "s3_run_and_compare_configs.R",
  "s3_run_frontal_vermis_anthro.R",

  # Study 4: fossil endocranial budgets (main before cerebellum sibling)
  "s4_endocranial.R",
  "s4_endocranial_cerebellum.R"
)

# Scripts that are source()d by another script and must NOT be run on their own
ENGINE_ONLY <- c("s3_predicValuesPGLS_MERGED_variant.R")

# ---- append any top-level script not already accounted for ----------
all_scripts <- list.files(SCRIPT_DIR, pattern = "\\.R$", full.names = FALSE)
known <- c(RUN_ORDER, ENGINE_ONLY, "run_all.R")
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
cat("\nEngine sourced by drivers (not run standalone):",
    paste(ENGINE_ONLY, collapse = ", "), "\n")
cat("Status log written to:", log_path, "\n")

failed <- status$script[status$result == "FAILED"]
if (length(failed)) {
  cat("\nFAILED scripts:\n  ", paste(failed, collapse = "\n   "), "\n", sep = "")
  # non-zero exit code so shell/CI can detect failures (Rscript only)
  if (!interactive()) quit(status = 1, save = "no")
} else {
  cat("\nAll scripts completed without error.\n")
}
