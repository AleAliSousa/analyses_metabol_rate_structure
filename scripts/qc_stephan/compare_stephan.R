# =============================================================================
# compare_stephan.R -- is anything fishy in data_raw/Stephan_primates.csv?
#
# Compares the SOURCE files against each other:
#     data_raw/Stephan_primates.csv     (curated here, by hand)
#     volumes_wide_select.csv           (Evo-M1-Trait-Data, volumes_compiled_select.R)
#
# Source files, not model inputs, on purpose. This is the only level that can
# still see a 2x laterality slip or a 1000x unit slip; once those reach the
# fitted matrices they are baked in and look like a plain disagreement.
#
# Every disagreement is reported WITH THE SOURCE ON BOTH SIDES, so a mismatch
# reads as "Stephan took de Sousa 2010, the merge took Bauernfeind 2013" rather
# than as an anonymous gap between two numbers:
#   Stephan side  metadata/qc_stephan/Stephan_primates_references_long.csv
#                 (run build_stephan_primates_reference_sheet.R to refresh)
#   wide side     a *_Source column in the wide file if one exists, otherwise
#                 resolved from volumes_unfiltered_select.csv, the per-source
#                 ledger behind the merge
#
# Species labels come from R/species_aliases.R -- ONE table, canonical
# namespace = data_raw/species.nwk tip labels. Nothing is renamed here.
#
# Prints four blocks (FILES / SPECIES / VARIABLES / SUSPECT) and nothing else.
# Everything compared is written to checks/qc_stephan/ regardless.
#
# The variable mapping is DATA, not code: metadata/qc_stephan/stephan_wide_crosswalk.csv
# Edit that csv to add, drop or re-role a pair; this script follows it.
# =============================================================================

setwd(local({ d <- normalizePath(getwd()); while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d); d }))
suppressPackageStartupMessages(library(tidyverse))
source("R/species_aliases.R")

# ---- 0. config --------------------------------------------------------------
EVOM1_DIR    <- file.path("/Users/crossmodal/Library/CloudStorage",
                          "OneDrive-AllenInstitute/Species/Evo-M1-Trait-Data",
                          "__merging_volumes")
wide_path    <- file.path(EVOM1_DIR, "volumes_wide_select.csv")
ledger_path  <- file.path(EVOM1_DIR, "volumes_unfiltered_select.csv")
stephan_path <- "data_raw/Stephan_primates.csv"
tree_path    <- "data_raw/species.nwk"
xwalk_path   <- "metadata/qc_stephan/stephan_wide_crosswalk.csv"
refs_path    <- "metadata/qc_stephan/Stephan_primates_references_long.csv"
out_dir      <- "checks/qc_stephan"

ROLES    <- c("study3", "predictor")  # add "reference" to compare the whole file
TOL_SAME <- 1e-6                      # at or below: the same number
TOL_FLAG <- 0.01                      # above: a MATERIAL disagreement (1%)
TOP_N    <- 20                        # suspect cells printed (all are written)

# ---- 1. read ----------------------------------------------------------------
p  <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")
nm <- function(x) suppressWarnings(as.numeric(gsub(",", "", as.character(x), fixed = TRUE)))
blank <- function(x) is.na(x) | trimws(as.character(x)) == ""
stamp <- function(f) if (file.exists(f)) format(file.mtime(f), "%Y-%m-%d %H:%M") else "MISSING"
# max() over an all-NA column returns -Inf and warns; a variable with nothing
# to compare should read as blank, not as a spurious extreme.
safe_max <- function(x) if (!length(x) || all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)

for (f in c(stephan_path, xwalk_path)) if (!file.exists(f)) stop("missing: ", f, call. = FALSE)
if (!file.exists(wide_path))
  stop("volumes_wide_select.csv not found at:\n  ", wide_path,
       "\nIs OneDrive mounted / synced?", call. = FALSE)

stephan <- read_csv(stephan_path, show_col_types = FALSE, name_repair = "minimal") %>%
  filter(!blank(Species)) %>% mutate(tip = canon_species(Species))
wide <- read_csv(wide_path, show_col_types = FALSE, name_repair = "minimal") %>%
  filter(!blank(Species)) %>% mutate(merge_sp = gsub("\\s+", " ", trimws(gsub("_", " ", Species))))
xwalk <- read_csv(xwalk_path, show_col_types = FALSE) %>% filter(role %in% ROLES)

tips <- if (file.exists(tree_path) && requireNamespace("ape", quietly = TRUE))
  ape::read.tree(tree_path)$tip.label else NULL

# the per-source ledger, used only where the wide file carries no source column
ledger <- if (file.exists(ledger_path))
  read_csv(ledger_path, show_col_types = FALSE) %>%
    mutate(Value = nm(Value), Year = suppressWarnings(as.integer(Year))) %>%
    filter(!is.na(Value)) else NULL

# Stephan-side provenance, resolved per data point by the reference-sheet script
refs <- if (file.exists(refs_path))
  read_csv(refs_path, show_col_types = FALSE) %>%
    transmute(tip = canon_species(Species), stephan_var = Stephan_column,
              stephan_source = coalesce(preferred_reference, status)) %>%
    distinct(tip, stephan_var, .keep_all = TRUE) else NULL

# ---- 2. pair species --------------------------------------------------------
# Direction matters: Stephan tip -> merge label. That way the pooled merge rows
# ("Pongo sp.") are found by BOTH tips that fall under them instead of one tip
# being silently picked. A merge row claimed by >1 tip is flagged, not resolved.
pairs <- stephan %>%
  transmute(tip, merge_sp = to_merge_species(tip)) %>%
  filter(merge_sp %in% wide$merge_sp) %>%
  group_by(merge_sp) %>% mutate(pooled = n() > 1L) %>% ungroup()

only_stephan <- setdiff(stephan$tip, pairs$tip)
only_wide    <- wide$merge_sp[!wide$merge_sp %in% pairs$merge_sp]
off_tree     <- if (is.null(tips)) character(0) else setdiff(stephan$tip, tips)

# ---- 3. compare every cell --------------------------------------------------
cells <- pmap_dfr(xwalk %>% select(label, stephan_var, wide_var, unit_scale, role),
                  function(label, stephan_var, wide_var, unit_scale, role) {
  if (!stephan_var %in% names(stephan) || !wide_var %in% names(wide)) return(NULL)
  a <- nm(stephan[[stephan_var]])[match(pairs$tip, stephan$tip)] * unit_scale
  b <- nm(wide[[wide_var]])[match(pairs$merge_sp, wide$merge_sp)]
  rel <- ifelse(!is.na(a) & !is.na(b), abs(a - b) / pmax(abs(a), abs(b), 1e-300), NA_real_)
  tibble(tip = pairs$tip, merge_sp = pairs$merge_sp, pooled = pairs$pooled,
         label, stephan_var, wide_var, role,
         stephan_value = a, wide_value = b,
         ratio = ifelse(!is.na(a) & !is.na(b) & a != 0, b / a, NA_real_),
         pct   = ifelse(!is.na(a) & !is.na(b) & a != 0, 100 * (b - a) / a, NA_real_),
         rel,
         status = case_when(is.na(a) & is.na(b) ~ "both_missing",
                            is.na(b)            ~ "Stephan_only",
                            is.na(a)            ~ "wide_only",
                            rel <= TOL_SAME     ~ "same",
                            rel <= TOL_FLAG     ~ "rounding",
                            TRUE                ~ "DIFFERENT"))
}) %>% filter(status != "both_missing")

if (!nrow(cells)) stop("no comparable cells -- check ROLES and the crosswalk csv", call. = FALSE)

near <- function(x, t) !is.na(x) & abs(x / t - 1) < 0.05
cells <- cells %>% mutate(anomaly = case_when(
  near(ratio, 2)     ~ "~2x  BILATERAL vs unilateral",
  near(ratio, 0.5)   ~ "~0.5x UNILATERAL vs bilateral",
  near(ratio, 1000)  ~ "~1000x UNIT (cm3/mm3, g/mg)",
  near(ratio, 0.001) ~ "~0.001x UNIT (cm3/mm3, g/mg)",
  TRUE ~ NA_character_))

# ---- 4. attach the source on both sides -------------------------------------
cells$stephan_source <- if (is.null(refs)) NA_character_ else
  refs$stephan_source[match(paste(cells$tip, cells$stephan_var), paste(refs$tip, refs$stephan_var))]

# (a) does the wide file carry its own source columns?
src_col <- function(v) {
  cand <- c(paste0(v, "_Source"), paste0(v, "_source"), paste0("Source_", v),
            paste0(sub("_Vol\\.mm3$", "", v), "_Source"))
  hit <- cand[cand %in% names(wide)]
  if (length(hit)) hit[1] else NA_character_
}
wide_src_cols <- setNames(vapply(unique(cells$wide_var), src_col, character(1)), unique(cells$wide_var))
SRC_MODE <- if (any(!is.na(wide_src_cols))) "wide *_Source columns" else
            if (!is.null(ledger)) "volumes_unfiltered_select.csv ledger" else "unavailable"

cells$wide_source <- NA_character_
if (any(!is.na(wide_src_cols))) {
  for (v in names(wide_src_cols)) {
    if (is.na(wide_src_cols[[v]])) next
    i <- cells$wide_var == v
    cells$wide_source[i] <- as.character(wide[[wide_src_cols[[v]]]])[match(cells$merge_sp[i], wide$merge_sp)]
  }
} else if (!is.null(ledger)) {
  # winning source = the ledger row whose Value IS the value the merge published;
  # ties broken the way the merge breaks them (Selected_collection, then newest).
  L <- ledger %>% filter(Variable %in% cells$wide_var, Species %in% cells$merge_sp)
  key <- paste(L$Species, L$Variable)
  for (i in which(!is.na(cells$wide_value))) {
    j <- which(key == paste(cells$merge_sp[i], cells$wide_var[i]))
    if (!length(j)) next
    hit <- j[abs(L$Value[j] - cells$wide_value[i]) <= 0.001 * pmax(abs(L$Value[j]), abs(cells$wide_value[i]))]
    if (!length(hit)) next
    sel <- hit[L$Team[hit] == "Selected_collection"]; if (length(sel)) hit <- sel
    hit <- hit[order(-L$Year[hit])]
    cells$wide_source[i] <- if (length(hit) > 1)
      paste0(L$Source[hit[1]], " (+", length(hit) - 1, " tied)") else L$Source[hit[1]]
  }
}

# ---- 5. report --------------------------------------------------------------
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

var_summary <- cells %>%
  group_by(variable = label, role) %>%
  summarise(N = sum(status %in% c("same", "rounding", "DIFFERENT")),
            same = sum(status == "same"), round = sum(status == "rounding"),
            DIFF = sum(status == "DIFFERENT"),
            S_only = sum(status == "Stephan_only"), w_only = sum(status == "wide_only"),
            med_ratio = round(median(ratio, na.rm = TRUE), 3),
            max_pct = round(safe_max(abs(pct)), 1),
            flag = paste(unique(na.omit(anomaly)), collapse = "; "),
            .groups = "drop") %>%
  arrange(desc(DIFF), desc(max_pct))

suspect <- cells %>%
  filter(status %in% c("DIFFERENT", "Stephan_only", "wide_only") | !is.na(anomaly)) %>%
  transmute(species = tip, variable = label, Stephan = stephan_value, wide = wide_value,
            ratio = round(ratio, 3), pct = round(pct, 1), status, anomaly,
            stephan_source, wide_source, pooled) %>%
  arrange(!is.na(anomaly) == FALSE, desc(abs(replace_na(pct, 0))))

write_csv(cells,       file.path(out_dir, "compare_cells.csv"))
write_csv(var_summary, file.path(out_dir, "variable_summary.csv"))
write_csv(suspect,     file.path(out_dir, "suspect_cells.csv"))
write_csv(pairs,       file.path(out_dir, "species_pairing.csv"))

cat("\n============================================================\n")
cat("Stephan_primates.csv  vs  volumes_wide_select.csv\n")
cat("============================================================\n")

cat("\nFILES\n")
p("  stephan  %-52s %3d sp  %s", basename(stephan_path), nrow(stephan), stamp(stephan_path))
p("  wide     %-52s %3d sp  %s", basename(wide_path), nrow(wide), stamp(wide_path))
p("  sources  %-52s %s", SRC_MODE,
  if (SRC_MODE == "volumes_unfiltered_select.csv ledger") stamp(ledger_path) else "")
p("  refs     %-52s %s", basename(refs_path),
  if (is.null(refs)) "MISSING -- run build_stephan_primates_reference_sheet.R" else stamp(refs_path))

cat("\nSPECIES   (canonical namespace = species.nwk tip labels)\n")
p("  Stephan rows            %3d", nrow(stephan))
if (is.null(tips)) {
  cat("  on the tree             --   species.nwk not read (install 'ape' to enable)\n")
} else {
  p("  on the tree             %3d/%d %s", nrow(stephan) - length(off_tree), length(tips),
    if (length(off_tree)) paste0(" OFF-TREE: ", paste(off_tree, collapse = ", ")) else " OK")
}
p("  paired with wide        %3d %s", nrow(pairs),
  if (any(pairs$pooled)) sprintf("  (%d via a pooled merge row)", sum(pairs$pooled)) else "")
p("  only in Stephan         %3d %s", length(only_stephan),
  if (length(only_stephan)) paste0("  ", paste(only_stephan, collapse = ", ")) else "")
p("  only in wide            %3d %s", length(only_wide),
  if (length(only_wide)) paste0("  ", paste(only_wide, collapse = ", ")) else "")
if (any(pairs$pooled))
  for (m in unique(pairs$merge_sp[pairs$pooled]))
    p("    pooled: wide '%s' compared against %s", m,
      paste(pairs$tip[pairs$merge_sp == m], collapse = " + "))

p("\nVARIABLES  (%d compared: %s)", nrow(var_summary),
  paste(sprintf("%d %s", table(var_summary$role), names(table(var_summary$role))), collapse = ", "))
print(as.data.frame(var_summary %>% select(-role)), row.names = FALSE)

p("\nSUSPECT CELLS  (%d of %d shown; full list -> %s)",
  min(TOP_N, nrow(suspect)), nrow(suspect), file.path(out_dir, "suspect_cells.csv"))
if (nrow(suspect)) {
  print(as.data.frame(head(suspect, TOP_N) %>%
          select(species, variable, Stephan, wide, ratio, pct, stephan_source, wide_source)),
        row.names = FALSE)
} else cat("  none\n")

n_anom <- sum(!is.na(cells$anomaly))
if (n_anom)
  p("\n  %d cell(s) sit on a unit/laterality ratio -- these are the fishy ones.", n_anom)
p("\nsame <= %g | rounding <= %g | DIFFERENT above that\n", TOL_SAME, TOL_FLAG)
