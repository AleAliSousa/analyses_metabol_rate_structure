# =====================================================================
# verify_outputs.R  --  check that every deliverable output exists
# ---------------------------------------------------------------------
# Reads metadata/output_manifest.csv and confirms each deliverable figure
# (and its expected formats) is present. Intended to be run AFTER
# scripts/run_all.R, so it verifies files the pipeline just (re)generated.
#
# Usage:
#   Rscript scripts/verify_outputs.R
#   # or inside R:  source("scripts/verify_outputs.R")
#
# Manifest columns:
#   output_stem       path without extension, relative to project root
#   formats           pipe-separated extensions expected, e.g. "png|pdf"
#   kind, study       descriptive
#   producing_script  script in scripts/ that writes it
#   consumer          where it is used (deck slide / MS figure)
#   check             file        -> assert every <stem>.<fmt> exists
#                     script_only -> assert producing_script exists (used when
#                                    the deck shows an older render, so the exact
#                                    file path is not pinned)
#                     none        -> known gap; reported, not a hard failure
#   status, notes     free text
#
# Exit status (under Rscript): non-zero if any hard check fails
# (a "file" row missing a file, or a "script_only" row missing its script).
# "none" rows are reported as KNOWN GAPS but do not fail the run.
# =====================================================================

# ---- locate project root (walk up for the .git sentinel) ------------
find_repo_root <- function(start) {
  d <- normalizePath(start, mustWork = FALSE)
  repeat {
    if (file.exists(file.path(d, ".git")) ||
        length(list.files(d, pattern = "\\.Rproj$")) > 0) return(d)
    parent <- dirname(d)
    if (parent == d) stop("Could not find project root (no .git/.Rproj above ", start, ")")
    d <- parent
  }
}
.start <- local({
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else getwd()
})
ROOT <- find_repo_root(.start)
setwd(ROOT)

MANIFEST <- "metadata/output_manifest.csv"
if (!file.exists(MANIFEST)) stop("Manifest not found: ", file.path(ROOT, MANIFEST))
man <- read.csv(MANIFEST, stringsAsFactors = FALSE)
cat("Project root:", ROOT, "\n")
cat("Manifest    :", MANIFEST, " (", nrow(man), " entries)\n\n", sep = "")

# ---- evaluate each row ----------------------------------------------
res <- data.frame(consumer = character(), output = character(),
                  check = character(), result = character(),
                  detail = character(), stringsAsFactors = FALSE)

for (i in seq_len(nrow(man))) {
  r <- man[i, ]
  chk <- r$check
  if (chk == "file") {
    exts <- strsplit(r$formats, "\\|")[[1]]
    paths <- paste0(r$output_stem, ".", exts)
    miss  <- paths[!file.exists(paths)]
    result <- if (length(miss) == 0) "OK" else "MISSING_FILE"
    detail <- if (length(miss) == 0) paste(basename(paths), collapse = ", ")
              else paste("missing:", paste(miss, collapse = ", "))
    out <- r$output_stem
  } else if (chk == "script_only") {
    sp <- file.path("scripts", r$producing_script)
    result <- if (file.exists(sp)) "OK_SCRIPT" else "MISSING_SCRIPT"
    detail <- if (file.exists(sp)) paste0(r$producing_script, " (render differs from deck)")
              else paste("script not found:", r$producing_script)
    out <- r$producing_script
  } else {  # "none" -> known gap
    result <- "KNOWN_GAP"
    detail <- r$notes
    out <- r$output_stem
  }
  res <- rbind(res, data.frame(consumer = r$consumer, output = out,
                               check = chk, result = result, detail = detail,
                               stringsAsFactors = FALSE))
}

# ---- also flag manifest rows whose producing_script is missing ------
have_scripts <- list.files("scripts", pattern = "\\.R$")
bad_prod <- man$producing_script[nzchar(man$producing_script) &
                                 !(man$producing_script %in% have_scripts)]
bad_prod <- unique(bad_prod)

# ---- report ----------------------------------------------------------
dir.create("logs", showWarnings = FALSE)
ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
write.csv(res, file.path("logs", paste0("verify_outputs_", ts, ".csv")), row.names = FALSE)
write.csv(res, file.path("logs", "verify_outputs_latest.csv"), row.names = FALSE)

cat(strrep("=", 70), "\n", sep = "")
print(res, row.names = FALSE)
cat(strrep("=", 70), "\n", sep = "")

n_ok    <- sum(res$result %in% c("OK", "OK_SCRIPT"))
n_gap   <- sum(res$result == "KNOWN_GAP")
n_missf <- sum(res$result == "MISSING_FILE")
n_misss <- sum(res$result == "MISSING_SCRIPT")
cat(sprintf("OK: %d   known gaps: %d   MISSING files: %d   MISSING scripts: %d\n",
            n_ok, n_gap, n_missf, n_misss))

if (n_gap > 0)
  cat("\nKNOWN GAPS (tracked, not a failure):\n  ",
      paste(res$output[res$result == "KNOWN_GAP"], collapse = "\n   "), "\n", sep = "")
if (length(bad_prod) > 0)
  cat("\nWARNING: manifest names producing scripts that are not in scripts/:\n  ",
      paste(bad_prod, collapse = "\n   "), "\n", sep = "")

hard_fail <- n_missf + n_misss
if (hard_fail > 0) {
  cat("\nFAILURES:\n")
  print(res[res$result %in% c("MISSING_FILE", "MISSING_SCRIPT"),
            c("consumer", "output", "detail")], row.names = FALSE)
  if (!interactive()) quit(status = 1, save = "no")
} else {
  cat("\nAll deliverable outputs present (known gaps excepted).\n")
}
