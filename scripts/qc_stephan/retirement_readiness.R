# =============================================================================
# retirement_readiness.R -- can study 3 run from volumes_wide_select.csv ALONE?
#
# This is the exit test for the whole qc_stephan/ layer. That layer exists only
# to check a hand-curated file, data_raw/Stephan_primates.csv, against the merge.
# The goal is to make volumes_wide_select.csv good enough that the hand-curated
# file is no longer needed -- at which point data_raw/Stephan_primates.csv,
# scripts/qc_stephan/, checks/qc_stephan/ and metadata/qc_stephan/ are deleted
# as a unit and every variable is read straight from the merge.
#
# It reports two counts, and only one of them blocks retirement:
#
#   BLOCKERS  Stephan_primates.csv has a value where volumes_wide_select.csv has
#             none. Deleting the file would LOSE data. Every blocker must be
#             fixed upstream, in volumes_compiled_select.R -- by wiring in a
#             source, or by letting step 7 derive a bilateral partner. Fixing a
#             blocker downstream, by editing the csv, defeats the point.
#
#   differs   both files have a value and they disagree by >1%. NOT a blocker:
#             nothing is lost by deleting the csv, but somebody has decided the
#             merge is wrong, or the csv is stale, and that should be settled
#             rather than inherited. Some are expected by construction -- see
#             role = "predictor" in the crosswalk, where the two files hold
#             genuinely different constructs (net vs whole brain).
#
# READY means BLOCKERS = 0. It does not mean differs = 0.
#
# Usage: Rscript scripts/qc_stephan/retirement_readiness.R
# =============================================================================

setwd(local({ d <- normalizePath(getwd()); while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d); d }))
suppressPackageStartupMessages(library(tidyverse))
source("helpers/species_aliases.R")

wide_path    <- file.path("/Users/crossmodal/Library/CloudStorage",
                          "OneDrive-AllenInstitute/Species/Evo-M1-Trait-Data",
                          "__merging_volumes/volumes_wide_select.csv")
stephan_path <- "data_raw/Stephan_primates.csv"
xwalk_path   <- "metadata/qc_stephan/stephan_wide_crosswalk.csv"
out_dir      <- "checks/qc_stephan"

ROLES   <- c("study3", "predictor")
REL_TOL <- 0.01

nm <- function(x) { x <- trimws(as.character(x))
  suppressWarnings(as.numeric(ifelse(x %in% c("", "NA", "NaN"), NA, gsub(",", "", x)))) }
blank <- function(x) is.na(x) | trimws(as.character(x)) == ""

if (!file.exists(wide_path))
  stop("volumes_wide_select.csv not found:\n  ", wide_path, "\nIs OneDrive synced?", call. = FALSE)

stephan <- read_csv(stephan_path, show_col_types = FALSE, name_repair = "minimal") %>%
  filter(!blank(Species)) %>% mutate(tip = canon_species(Species))
wide <- read_csv(wide_path, show_col_types = FALSE, name_repair = "minimal") %>%
  filter(!blank(Species)) %>% mutate(merge_sp = gsub("\\s+", " ", trimws(gsub("_", " ", Species))))
xwalk <- read_csv(xwalk_path, show_col_types = FALSE) %>% filter(role %in% ROLES)

stephan$merge_sp <- to_merge_species(stephan$tip)

cells <- pmap_dfr(xwalk %>% select(label, stephan_var, wide_var, unit_scale, role),
                  function(label, stephan_var, wide_var, unit_scale, role) {
  if (!stephan_var %in% names(stephan)) return(NULL)
  a <- nm(stephan[[stephan_var]]) * unit_scale
  b <- if (wide_var %in% names(wide)) nm(wide[[wide_var]])[match(stephan$merge_sp, wide$merge_sp)]
       else rep(NA_real_, nrow(stephan))
  tibble(species = stephan$tip, label, stephan_var, wide_var, role,
         stephan_value = a, wide_value = b,
         verdict = case_when(
           is.na(a)              ~ "not_needed",          # csv has nothing to lose
           is.na(b)              ~ "BLOCKER",             # csv has a value the merge lacks
           abs(a - b) <= REL_TOL * pmax(abs(a), abs(b)) ~ "ok",
           TRUE                  ~ "differs"))
})

summary_tbl <- cells %>%
  group_by(variable = label, role) %>%
  summarise(wide_N   = sum(!is.na(wide_value)),
            steph_N  = sum(!is.na(stephan_value)),
            BLOCKERS = sum(verdict == "BLOCKER"),
            differs  = sum(verdict == "differs"),
            blocked_species = paste(species[verdict == "BLOCKER"], collapse = ", "),
            .groups = "drop") %>%
  arrange(desc(BLOCKERS), desc(differs))

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write_csv(summary_tbl, file.path(out_dir, "retirement_readiness_summary.csv"))
write_csv(cells %>% filter(verdict %in% c("BLOCKER", "differs")) %>% arrange(verdict, label, species),
          file.path(out_dir, "retirement_blockers.csv"))

nb <- sum(summary_tbl$BLOCKERS); nd <- sum(summary_tbl$differs)

cat("\n============================================================\n")
cat("RETIREMENT READINESS -- can study 3 run from volumes_wide_select.csv alone?\n")
cat("============================================================\n\n")
print(as.data.frame(summary_tbl %>%
        mutate(blocked_species = substr(blocked_species, 1, 46)) %>%
        select(variable, wide_N, steph_N, BLOCKERS, differs, blocked_species)),
      row.names = FALSE)

cat("\n------------------------------------------------------------\n")
if (nb == 0) {
  cat("READY. No cell exists only in Stephan_primates.csv.\n")
  cat("Deleting data_raw/Stephan_primates.csv and the qc_stephan/ layer loses nothing.\n")
  if (nd) cat(sprintf("(%d cell(s) still disagree -- settle or accept them, but they do not block.)\n", nd))
} else {
  cat(sprintf("NOT READY: %d cell(s) exist only in Stephan_primates.csv.\n", nb))
  cat("Fix each one UPSTREAM in volumes_compiled_select.R -- wire in the missing\n")
  cat("source, or let step 7 derive the bilateral partner. Do not patch the csv.\n")
  cat(sprintf("%d further cell(s) disagree by >1%%; those do not block retirement.\n", nd))
  blockers_by_sp <- cells %>% filter(verdict == "BLOCKER") %>% count(species, sort = TRUE)
  cat("\nBlocking species:\n")
  for (i in seq_len(nrow(blockers_by_sp)))
    cat(sprintf("  %-28s %d variable(s)\n", blockers_by_sp$species[i], blockers_by_sp$n[i]))
}
cat("\nWrote ", file.path(out_dir, "retirement_blockers.csv"), "\n", sep = "")
