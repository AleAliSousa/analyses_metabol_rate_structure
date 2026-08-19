# =============================================================================
# ARCHIVED 2026-08-18 -- superseded by scripts/qc_stephan/compare_stephan.R
#
# Replaced because it printed far more than could be read: 565 lines, two modes,
# a 60-row crosswalk inline, and its own copy of the species synonyms. In the
# rewrite the crosswalk became data (metadata/qc_stephan/stephan_wide_crosswalk.csv),
# the synonyms moved to R/species_aliases.R, and the console output is four blocks.
#
# KEPT FOR ONE REASON: mode "model", section 3 below, is still the only code that
# compares the two s3 FITTED matrices (tree restriction, effective N per
# variable, which species drop out of which analysis). The replacement is
# raw-mode only, by choice. Lift section 3 out if that question comes back.
# =============================================================================

# ============================================================
# QC: Stephan_primates.csv  vs  volumes_wide_select.csv
#
# Consolidated QC script with one shared implementation for source-file and
# model-input comparisons. The two stages remain as modes because they answer
# distinct questions:
#
#   MODE "raw"    the SOURCE files, before any analysis touches them.
#                 Species-name normalisation, duplicate-taxon coalescing, the
#                 explicit Pongo pairing, unit/laterality ratio diagnostics.
#                 This is the only level that can still see a 2x laterality or
#                 1000x unit error -- by the time the model-input matrices are
#                 written, such an error has been baked in and looks like a
#                 plain disagreement.
#
#   MODE "model"  the matrices the two s3 scripts actually fitted. Adds tree-tip
#                 restriction, effective N per variable, and which species drop
#                 out of which analysis. This is the level that explains why two
#                 fits differ even where the numbers look similar.
#
# Shared definitions (crosswalk, tolerances, classifier, summaries, and species
# normalization) are declared once and reused by both modes.
#
# Inputs
#   raw    : data_raw/Stephan_primates.csv
#            <wide_path> (volumes_wide_select.csv in the compilation repo)
#   model  : checks/Stephan_data_used.csv
#            checks/volumes_wide_select_data_used.csv
#            data_raw/species.nwk                      (optional, enables tree logic)
#
# Outputs -- every file both original scripts wrote is still written
#   checks/qc_stephan/coalesced_species_rows.csv
#   checks/qc_stephan/coalesced_row_conflicts.csv
#   checks/qc_stephan/species_only_in_volumes_wide.csv
#   checks/qc_stephan/species_only_in_stephan.csv
#   checks/qc_stephan/requested_study3_variables_missing.csv
#   checks/qc_stephan/comparison_summary.csv
#   checks/qc_stephan/species_level_differences.csv
#   checks/qc_stephan/unit_laterality_anomalies.csv
#   checks/qc_stephan/Stephan_vs_volumes_wide_side_by_side.csv
#   checks/qc_stephan/Stephan_vs_volumes_wide_cell_discrepancies.csv
#   checks/qc_stephan/Stephan_vs_volumes_wide_variable_summary.csv
#   checks/qc_stephan/Stephan_vs_volumes_wide_species_coverage.csv
#   checks/qc_stephan/compare_all_cells.csv  (both modes, one long table)
# ============================================================

setwd(local({ d <- normalizePath(getwd()); while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d); d }))  # repo root (portable; see R/project_root.R)

library(tidyverse)

# ============================================================
# 0  OPTIONS
# ============================================================

modes <- c("raw", "model")     # run either or both

wide_path    <- "/Users/crossmodal/Library/CloudStorage/OneDrive-AllenInstitute/Species/Evo-M1-Trait-Data/__merging_volumes/volumes_wide_select.csv"
stephan_path <- "data_raw/Stephan_primates.csv"
used_stephan_path <- "checks/Stephan_data_used.csv"
used_wide_path    <- "checks/volumes_wide_select_data_used.csv"
tree_file         <- "data_raw/species.nwk"

## All QC outputs live together, under the qc_stephan/ layer. That layer is
## scaffolding for checking over Stephan_primates.csv and is meant to be thrown
## away as a unit once the file is replaced with a fresh version -- keeping it
## in its own folder is what makes that a deliberate act rather than an
## accident. Paths stay repo-root-relative (see the .git walk above), so this
## script works from wherever it is invoked.
output_dir <- "checks/qc_stephan"

# --- tolerance, defined ONCE for both modes ---------------------------------
# Previously: raw used rel_tol = 0.001 as a single cut; model used REL_TOL 1e-6
# plus REL_FLAG 0.01. Three levels covers both, and the old raw behaviour is
# recoverable by setting tol_flag = 0.001.
tol_same <- 1e-6    # at or below this, the two values are the same number
tol_flag <- 0.01    # above this, call it a MATERIAL disagreement (1%)
                    # between the two: "rounding"
raw_rel_tol <- 0.001 # legacy compare_volumes_files threshold (0.1%); retained
                    # separately so raw QC still flags the same cells as before

RESTRICT_TO_TREE <- TRUE   # model mode: keep only species that can reach a fit

# Duplicate-row conflict guard: which columns make a conflict fatal.
#   "compared" = only columns this run actually compares (follows the crosswalk)
#   "all"      = every column in either file
#   character vector = an explicit set
guard_scope <- "compared"

# ============================================================
# 1  SHARED DEFINITIONS  -- the single source of truth
# ============================================================

# --- the variable crosswalk -------------------------------------------------
# role:
#   "study3"      compared, and counted in the disagreement totals
#   "report_only" compared and reported, but NOT counted as a disagreement --
#                 for pairs where the two files legitimately hold different
#                 constructs, so a difference is expected rather than a fault
#   "reference"   documentation only; never compared (kept so the mapping of the
#                 wider file is recorded in one place)
# unit_scale: multiplier applied to the STEPHAN value to reach the merge's units.
# model_stephan / model_wide: column names in the MODEL-INPUT files where they
#   differ from the raw ones; NA means "same as the raw column".
crosswalk <- tribble(
  ~label,                          ~stephan_var,                    ~wide_var,                                           ~unit_scale, ~role,        ~model_stephan,           ~model_wide,              ~note,

  # ---- Study 3 predictor -----------------------------------------------------
  # Total_brain_net_volume is Stephan's NET brain volume (ventricles and meninges
  # removed). In the merge the same term can be won on recency by Bauernfeind
  # 2013 / MacLeod 2003 whole-brain volumes, which are not the same construct --
  # hence role = "report_only": the difference is shown but is not a fault of
  # either file. In MODEL mode the column is Preferred_brain_volume, which IS the
  # fitted x-axis and IS compared, because there the two scripts' differing
  # definitions are exactly what we want to expose.
  "Preferred brain volume",        "Total_brain_net_volume",        "Total_brain_net_volume_Vol.mm3",                          1, "report_only", "Preferred_brain_volume", "Preferred_brain_volume", "raw: different constructs (net vs whole brain); model: the PGLS x-axis",

  # ---- Study 3 response regions (14) ----------------------------------------
  "Corpus geniculatum laterale",   "LGN_Sousa",                     "Corpus_geniculatum_laterale_Vol.mm3",                     1, "study3",      NA,                       NA,                       "",
  "Amygdala",                      "Amygdala",                      "Amygdala_Vol.mm3",                                        1, "study3",      NA,                       NA,                       "",
  "Pallidum",                      "Pallidum",                      "Pallidum_Vol.mm3",                                        1, "study3",      NA,                       NA,                       "",
  "Neocortex white",               "NeoW_Frahm",                    "Neocortex_white_matter_Vol.mm3",                          1, "study3",      NA,                       NA,                       "",
  "Insular cortex (grey)",         "Total_insula_volume_L",         "Insula_left_Vol.mm3",                                     1, "study3",      NA,                       NA,                       "Stephan is LEFT only; Insula_Vol.mm3 is bilateral (~2x)",
  "Nucleus subthalamicus Luysi",   "Nucleus_subthalamicus",         "Nucleus_subthalamicus_Vol.mm3",                           1, "study3",      NA,                       NA,                       "",
  "Capsula interna",               "Capsula_interna",               "Capsula_interna_Vol.mm3",                                 1, "study3",      NA,                       NA,                       "",
  "Striatum",                      "Striatum",                      "Striatum_Vol.mm3",                                        1, "study3",      NA,                       NA,                       "",
  "Area striata grey",             "ASG_Sousa",                     "Area_striata_grey_matter_Vol.mm3",                        1, "study3",      NA,                       NA,                       "de Sousa publishes 2x-left, so this is the BILATERAL merge column",
  "Neocortex grey",                "NeoG_Frahm",                    "Neocortex_grey_matter_Vol.mm3",                           1, "study3",      NA,                       NA,                       "",
  "Mesencephalon",                 "Mesencephalon",                 "Mesencephalon_Vol.mm3",                                   1, "study3",      NA,                       NA,                       "",
  "Cerebellum",                    "Cerebellum",                    "Cerebellum_Vol.mm3",                                      1, "study3",      NA,                       NA,                       "merge can take MacLeod 2003 over Stephan on recency",
  "Hippocampus",                   "Hippocampus",                   "Hippocampus_Vol.mm3",                                     1, "study3",      NA,                       NA,                       "",
  "Nucleus dentatus cerebelli",    "Lateral_cerebellar_nuclei",     "Lateral_cerebellar_nuclei_Vol.mm3",                       1, "study3",      NA,                       NA,                       "",

  # ---- documentation only: the rest of the mapping, never compared -----------
  "Body weight",                   "Body_weight",                   "Body_Mass.g",                                             1, "reference",   NA, NA, "",
  "Brain weight",                  "Brain_weight",                  "Brain_Mass.mg",                                           1, "reference",   NA, NA, "merge Bauernfeind rows are in g under an .mg term -- upstream bug",
  "Ventricles",                    "Ventricles",                    "Ventricles_Vol.mm3",                                      1, "reference",   NA, NA, "",
  "Medulla oblongata",             "Medulla_oblongata",             "Medulla_oblongata_Vol.mm3",                               1, "reference",   NA, NA, "",
  "Diencephalon",                  "Diencephalon",                  "Diencephalon_Vol.mm3",                                    1, "reference",   NA, NA, "",
  "Telencephalon",                 "Telencephalon",                 "Telencephalon_Vol.mm3",                                   1, "reference",   NA, NA, "",
  "Pons",                          "Pons",                          "Pons_Vol.mm3",                                            1, "reference",   NA, NA, "",
  "Bulbus olfactorius",            "Bulbus_olfactorius",            "Bulbus_olfactorius_Vol.mm3",                              1, "reference",   NA, NA, "",
  "Bulbus olfactorius accessorius","Bulbus_olfactorius_accessorius","Bulbus_olfactorius_accessorius_Vol.mm3",                  1, "reference",   NA, NA, "",
  "Nucleus tractus olfactorius",   "Nucleus_tractus_olfactorius",   "Nucleus_tractus_olfactorius_Vol.mm3",                     1, "reference",   NA, NA, "",
  "Lobus piriformis",              "Lobus_piriformis",              "Lobus_piriformis_Vol.mm3",                                1, "reference",   NA, NA, "",
  "Septum",                        "Septum",                        "Septum_Vol.mm3",                                          1, "reference",   NA, NA, "",
  "Schizocortex",                  "Schizocortex",                  "Schizo_cortex_Vol.mm3",                                   1, "reference",   NA, NA, "",
  "Palaeocortex",                  "Palaeocortex",                  "Palaeocortex_Vol.mm3",                                    1, "reference",   NA, NA, "",
  "Neocortex (grey+white)",        "NeoWG",                         "Neocortex_Vol.mm3",                                       1, "reference",   NA, NA, "",
  "Complexus centromedialis",      "Complexus_centromedialis",      "Complexus_centromedialis_Vol.mm3",                        1, "reference",   NA, NA, "",
  "Complexus corticobasolateralis","Complexus_cortico_basolateralis","Complexus_corticobasolateralis_Vol.mm3",                 1, "reference",   NA, NA, "",
  "Nucleus amygdalae basalis",     "Nucleus_amygdalae_basalis_parsmagnocellularis", "Nucleus_amygdalae_basalis_pars_magnocellularis_Vol.mm3", 1, "reference", NA, NA, "",
  "Epithalamus",                   "Epithalamus",                   "Epithalamus_Vol.mm3",                                     1, "reference",   NA, NA, "",
  "Thalamus",                      "Thalamus",                      "Thalamus_Vol.mm3",                                        1, "reference",   NA, NA, "",
  "Hypothalamus",                  "Hypothalamus",                  "Hypothalamus_Vol.mm3",                                    1, "reference",   NA, NA, "",
  "Subthalamus",                   "Subthalamus",                   "Subthalamus_Vol.mm3",                                     1, "reference",   NA, NA, "",
  "Tractus opticus",               "Tractus_opticus",               "Tractus_opticus_Vol.mm3",                                 1, "reference",   NA, NA, "",
  "Cerebellar nuclei total",       "Cerebellar_nuclei_total",       "Cerebellar_nuclei_total_Vol.mm3",                         1, "reference",   NA, NA, "",
  "Interpositus cerebellar nuclei","Interpositus_cerebellar_nuclei","Interpositus_cerebellar_nuclei_Vol.mm3",                  1, "reference",   NA, NA, "",
  "Medial cerebellar nuclei",      "Medial_cerebellar_nuclei",      "Medial_cerebellar_nuclei_Vol.mm3",                        1, "reference",   NA, NA, "",
  "Complexus vestibularis",        "Complexus.vestibularis",        "Complexus_vestibularis_unilateral_Vol.mm3",               1, "reference",   NA, NA, "Stephan is UNILATERAL",
  "N. vestibularis superior",      "Nucleus.vestibularis.superior", "Nucleus_vestibularis_superior_unilateral_Vol.mm3",        1, "reference",   NA, NA, "Stephan is UNILATERAL",
  "N. vestibularis lateralis",     "Nucleus.vestibularis.lateralis","Nucleus_vestibularis_lateralis_unilateral_Vol.mm3",       1, "reference",   NA, NA, "Stephan is UNILATERAL",
  "N. vestibularis medialis",      "Nucleus.vestibularis.medialis", "Nucleus_vestibularis_medialis_unilateral_Vol.mm3",        1, "reference",   NA, NA, "Stephan is UNILATERAL",
  "N. vestibularis descendens",    "Nucleus.vestibularis.descendens","Nucleus_vestibularis_descendens_unilateral_Vol.mm3",     1, "reference",   NA, NA, "Stephan is UNILATERAL",
  "N. septalis triangularis",      "Nucleus.septalis.triangularis", "Nucleus_septalis_triangularis_Vol.mm3",                   1, "reference",   NA, NA, "",
  "Corpus subfornicale",           "Corpus.subfornicale",           "Corpus_subfornicale_Vol.mm3",                             1, "reference",   NA, NA, "",
  "N. habenularis medialis",       "Nucleus.habenularis.medialis",  "Nucleus_habenularis_medialis_Vol.mm3",                    1, "reference",   NA, NA, "",
  "Corpus pineale",                "Corpus.pineale",                "Corpus_pineal_Vol.mm3",                                   1, "reference",   NA, NA, "",
  "Corpus subcommissurale",        "Corpus.subcommissurale",        "Corpus_subcommissurale_Vol.mm3",                          1, "reference",   NA, NA, "",
  "Granular insula (L)",           "Granular_volume_L",             "Granular_insular_cortex_left_Vol.mm3",                    1, "reference",   NA, NA, "Stephan is LEFT only",
  "Dysgranular insula (L)",        "Dysgranular_volume_L",          "Dysgranular_insular_cortex_left_Vol.mm3",                 1, "reference",   NA, NA, "Stephan is LEFT only",
  "Agranular insula (L)",          "Agranular_volume_L",            "Agranular_insular_cortex_left_Vol.mm3",                   1, "reference",   NA, NA, "Stephan is LEFT only",
  "Fronto-insular (L)",            "FI_volume_L",                   "fronto_insular_cortex_left_Vol.mm3",                      1, "reference",   NA, NA, "Stephan is LEFT only",
  "Frontal motor grey",            "Frontal.motor.Gray",            "FrontalMotor_grey_matter_Vol.mm3",                     1000, "reference",   NA, NA, "Stephan in cm3, merge in mm3",
  "Frontal motor white",           "Frontal.motor",                 "FrontalMotor_white_matter_Vol.mm3",                    1000, "reference",   NA, NA, "Stephan in cm3, merge in mm3"
)
stopifnot(!anyDuplicated(crosswalk$label), !anyDuplicated(crosswalk$stephan_var))

# The set this run compares. Comment a line out of `crosswalk` above, or change
# its role, and everything downstream follows -- there is no second list.
compared <- crosswalk %>% filter(role %in% c("study3", "report_only"))

# --- species normalisation --------------------------------------------------
# Underscores -> spaces, then superseded labels rewritten to the accepted name.
# Pongo is deliberately NOT handled here; see pongo_pairs.
taxon_synonyms <- c(
  "Lagothrix lagothricha" = "Lagothrix lagotricha",
  "Callicebus moloch"     = "Plecturocebus moloch",
  "Gorilla sp."           = "Gorilla gorilla",
  "Tarsius sp."           = "Tarsius syrichta"
)
norm_species <- function(x) {
  x <- str_squish(str_replace_all(as.character(x), "_", " "))
  hit <- x %in% names(taxon_synonyms)
  x[hit] <- unname(taxon_synonyms[x[hit]])
  x
}
is_blank <- function(x) is.na(x) | trimws(as.character(x)) == ""
num <- function(x) suppressWarnings(as.numeric(gsub(",", "", as.character(x), fixed = TRUE)))

# The merge holds one "Pongo sp." row -- a sensu lato mean of Bornean, Sumatran
# and hybrid material -- while Stephan splits abelii from pygmaeus. Compare
# against BOTH, flagged, rather than silently pinning it to one.
pongo_pairs <- tribble(
  ~wide_species, ~stephan_species,  ~pair_note,
  "Pongo sp.",   "Pongo abelii",    "merge 'Pongo sp.' is a sensu lato mean (Bornean + Sumatran + hybrid); Stephan row is Sumatran",
  "Pongo sp.",   "Pongo pygmaeus",  "merge 'Pongo sp.' is a sensu lato mean (Bornean + Sumatran + hybrid); Stephan row is Bornean"
)

# --- ONE cell classifier, used by both modes --------------------------------
classify <- function(a, b) {
  denom <- pmax(abs(a), abs(b))
  rel <- ifelse(!is.na(a) & !is.na(b), ifelse(denom == 0, 0, abs(a - b) / denom), NA_real_)
  status <- case_when(
    is.na(a) &  is.na(b) ~ "both_missing",
    !is.na(a) & is.na(b) ~ "Stephan_only",
    is.na(a) & !is.na(b) ~ "volumes_wide_only",
    rel <= tol_same      ~ "same",
    rel <= tol_flag      ~ "rounding",
    TRUE                 ~ "DIFFERENT")
  tibble(status = status, rel = rel,
         pct_diff = ifelse(!is.na(a) & !is.na(b) & a != 0, 100 * (b - a) / a, NA_real_),
         ratio    = ifelse(!is.na(a) & !is.na(b) & a != 0, b / a, NA_real_))
}

# --- ONE ratio-anomaly rule -------------------------------------------------
near_ratio <- function(x, target, tol = 0.05) !is.na(x) & abs(x / target - 1) < tol
anomaly_of <- function(ratio) case_when(
  near_ratio(ratio, 1000)  ~ "merge is ~1000x Stephan (cm3/mm3 or g/mg)",
  near_ratio(ratio, 0.001) ~ "merge is ~1/1000x Stephan (cm3/mm3 or g/mg)",
  near_ratio(ratio, 2)     ~ "merge is ~2x Stephan (bilateral vs unilateral)",
  near_ratio(ratio, 0.5)   ~ "merge is ~0.5x Stephan (unilateral vs bilateral)",
  TRUE ~ NA_character_)

# --- ONE per-variable summariser --------------------------------------------
safe_max <- function(x) if (length(x) == 0 || all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
safe_cor <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  if (sum(ok) < 3 || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  cor(x[ok], y[ok])
}
summarise_cells <- function(cells, keep = rep(TRUE, nrow(cells))) {
  cells[keep, , drop = FALSE] %>%
    group_by(label, stephan_var, wide_var, unit_scale, role) %>%
    summarise(
      n_both         = sum(status %in% c("same", "rounding", "DIFFERENT")),
      N_Stephan      = sum(!is.na(stephan_value)),
      N_volumes_wide = sum(!is.na(wide_value)),
      n_same         = sum(status == "same"),
      n_rounding     = sum(status == "rounding"),
      n_DIFFERENT    = sum(status == "DIFFERENT"),
      n_agree_0_1pct  = sum(!is.na(rel) & rel <= raw_rel_tol),
      n_disagree_0_1pct = sum(!is.na(rel) & rel > raw_rel_tol),
      n_Stephan_only = sum(status == "Stephan_only"),
      n_wide_only    = sum(status == "volumes_wide_only"),
      cor            = safe_cor(stephan_value, wide_value),
      median_ratio   = median(ratio, na.rm = TRUE),
      median_abs_pct = median(abs(pct_diff), na.rm = TRUE),
      max_abs_pct    = safe_max(abs(pct_diff)),
      max_abs_diff   = safe_max(abs(wide_value - stephan_value)),
      .groups = "drop") %>%
    arrange(desc(max_abs_pct))
}

# --- duplicate-taxon coalescing (raw mode) ----------------------------------
coalesce_duplicate_rows <- function(dat, label) {
  dups <- unique(dat$Species_std[duplicated(dat$Species_std)])
  if (!length(dups)) return(list(data = dat, report = tibble(), conflicts = tibble()))
  report <- tibble(); conflicts <- tibble(file = character(), Species_std = character(),
                                          column = character(), values = character())
  keep_rows <- rep(TRUE, nrow(dat))
  for (sp in dups) {
    idx <- which(dat$Species_std == sp); block <- dat[idx, , drop = FALSE]
    report <- bind_rows(report, tibble(
      file = label, Species_std = sp, n_rows_merged = length(idx),
      labels_merged = paste(block$Species, collapse = " | "),
      populated_fields_per_row = paste(
        vapply(seq_len(nrow(block)), function(i) sum(!is_blank(unlist(block[i, ]))), integer(1)),
        collapse = " | ")))
    for (col in setdiff(names(dat), c("Species", "Species_std"))) {
      vals <- block[[col]]; good <- unique(vals[!is_blank(vals)])
      if (length(good) > 1L)
        conflicts <- bind_rows(conflicts, tibble(file = label, Species_std = sp, column = col,
                                                 values = paste(format(good), collapse = " | ")))
      if (length(good) >= 1L) dat[[col]][idx[1]] <- good[1]
    }
    dat[["Species"]][idx[1]] <- sp; keep_rows[idx[-1]] <- FALSE
  }
  list(data = dat[keep_rows, , drop = FALSE], report = report, conflicts = conflicts)
}

# --- tree tips (model mode) -------------------------------------------------
read_tree_tips <- function() {
  if (!file.exists(tree_file)) {
    warning("Tree file not found at ", tree_file, "; species-coverage section skipped."); return(NULL)
  }
  if (!requireNamespace("ape", quietly = TRUE)) {
    warning("Package 'ape' not available; species-coverage section skipped."); return(NULL)
  }
  ape::read.tree(tree_file)$tip.label
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
all_cells_both <- list()

# ============================================================
# 2  MODE "raw" -- the source files
# ============================================================
if ("raw" %in% modes) {
  cat("\n============================================================\n")
  cat("MODE raw : source files, before any analysis\n")
  cat("============================================================\n")
  if (!file.exists(wide_path))
    stop("volumes_wide_select.csv not found at:\n  ", wide_path, "\nIs OneDrive mounted / synced?")
  if (!file.exists(stephan_path))
    stop("Stephan_primates.csv not found at:\n  ", file.path(getwd(), stephan_path))
  cat("wide    : ", wide_path, "  (", format(file.mtime(wide_path), "%Y-%m-%d %H:%M"), ")\n", sep = "")
  cat("stephan : ", stephan_path, "  (", format(file.mtime(stephan_path), "%Y-%m-%d %H:%M"), ")\n", sep = "")

  wide <- read_csv(wide_path, show_col_types = FALSE, name_repair = "minimal") %>%
    filter(!is_blank(Species)) %>% mutate(Species_std = norm_species(Species))
  stephan <- read_csv(stephan_path, show_col_types = FALSE, name_repair = "minimal") %>%
    filter(!is_blank(Species)) %>% mutate(Species_std = norm_species(Species))
  stopifnot("Species" %in% names(wide), "Species" %in% names(stephan))

  wide_co <- coalesce_duplicate_rows(wide, "volumes_wide_select")
  step_co <- coalesce_duplicate_rows(stephan, "Stephan_primates")
  coalesce_report    <- bind_rows(wide_co$report, step_co$report)
  coalesce_conflicts <- bind_rows(wide_co$conflicts, step_co$conflicts)
  write_csv(coalesce_report, file.path(output_dir, "coalesced_species_rows.csv"))
  if (nrow(coalesce_report))
    for (i in seq_len(nrow(coalesce_report)))
      message("Coalesced ", coalesce_report$file[i], ": ", coalesce_report$labels_merged[i],
              " -> ", coalesce_report$Species_std[i])
  wide <- wide_co$data; stephan <- step_co$data
  wide_names <- names(wide); stephan_names <- names(stephan)

  # which compared pairs actually exist in both files
  cw <- compared %>%
    mutate(missing_from = case_when(
      !(wide_var %in% wide_names) & !(stephan_var %in% stephan_names) ~ "both",
      !(wide_var %in% wide_names) ~ "volumes_wide_select",
      !(stephan_var %in% stephan_names) ~ "Stephan",
      TRUE ~ NA_character_))
  write_csv(cw %>% filter(!is.na(missing_from)) %>% select(label, stephan_var, wide_var, missing_from),
            file.path(output_dir, "requested_study3_variables_missing.csv"))
  if (any(!is.na(cw$missing_from)))
    warning("Dropping ", sum(!is.na(cw$missing_from)), " pair(s) whose columns are absent:\n",
            paste0("  ", cw$stephan_var[!is.na(cw$missing_from)], " <-> ",
                   cw$wide_var[!is.na(cw$missing_from)], collapse = "\n"))
  cw <- cw %>% filter(is.na(missing_from)) %>% select(-missing_from)
  stopifnot(nrow(cw) > 0)

  # duplicate-row conflict guard, resolved against the FINAL compared set
  guarded <- function(all_cols, compared_cols) {
    if (identical(guard_scope, "all")) all_cols
    else if (identical(guard_scope, "compared")) compared_cols
    else guard_scope
  }
  if (nrow(coalesce_conflicts))
    coalesce_conflicts <- coalesce_conflicts %>%
      mutate(fatal = (file == "volumes_wide_select" & column %in% guarded(wide_names, cw$wide_var)) |
                     (file == "Stephan_primates"    & column %in% guarded(stephan_names, cw$stephan_var))) %>%
      arrange(desc(fatal), file, Species_std, column)
  write_csv(coalesce_conflicts, file.path(output_dir, "coalesced_row_conflicts.csv"))
  n_fatal <- if (nrow(coalesce_conflicts)) sum(coalesce_conflicts$fatal) else 0L
  if (n_fatal > 0) {
    print(coalesce_conflicts %>% filter(fatal), n = Inf)
    stop("Rows that collapse onto one accepted name hold conflicting values in ", n_fatal,
         " compared column(s) -- see ", file.path(output_dir, "coalesced_row_conflicts.csv"),
         ". Resolve in the source file before comparing.")
  }

  # species pairing: identity, plus the explicit Pongo one-to-many
  species_pairs <- bind_rows(
    tibble(wide_species = intersect(wide$Species_std, stephan$Species_std)) %>%
      mutate(stephan_species = wide_species, pair_note = NA_character_),
    pongo_pairs %>% filter(wide_species %in% wide$Species_std,
                           stephan_species %in% stephan$Species_std)) %>%
    arrange(wide_species, stephan_species)
  write_csv(tibble(Species = sort(setdiff(wide$Species_std, species_pairs$wide_species))),
            file.path(output_dir, "species_only_in_volumes_wide.csv"))
  write_csv(tibble(Species = sort(setdiff(stephan$Species_std, species_pairs$stephan_species))),
            file.path(output_dir, "species_only_in_stephan.csv"))
  cat("\nSpecies pairs compared: ", nrow(species_pairs), " (",
      sum(!is.na(species_pairs$pair_note)), " flagged)\n", sep = "")

  raw_cells <- pmap_dfr(cw %>% select(label, stephan_var, wide_var, unit_scale, role),
                        function(label, stephan_var, wide_var, unit_scale, role) {
    s <- tibble(stephan_species = stephan$Species_std,
                stephan_value = num(stephan[[stephan_var]]) * unit_scale)
    w <- tibble(wide_species = wide$Species_std, wide_value = num(wide[[wide_var]]))
    d <- species_pairs %>% left_join(s, by = "stephan_species") %>% left_join(w, by = "wide_species")
    bind_cols(
      tibble(mode = "raw", Species = d$wide_species, stephan_species = d$stephan_species,
             pair_note = d$pair_note, label = label, stephan_var = stephan_var,
             wide_var = wide_var, unit_scale = unit_scale, role = role,
             stephan_value = d$stephan_value, wide_value = d$wide_value, on_tree = NA_character_),
      classify(d$stephan_value, d$wide_value))
  }) %>% filter(status != "both_missing")

  raw_cells$anomaly <- anomaly_of(raw_cells$ratio)
  all_cells_both[["raw"]] <- raw_cells

  both_present <- raw_cells$status %in% c("same", "rounding", "DIFFERENT")
  write_csv(summarise_cells(raw_cells, both_present), file.path(output_dir, "comparison_summary.csv"))
  # Preserve compare_volumes_files.R's original 0.1% reporting threshold.
  # The status classifier above still distinguishes numerical equality, small
  # differences (<=1%), and material differences (>1%).
  write_csv(raw_cells %>%
              filter(!is.na(rel), rel > raw_rel_tol) %>%
              mutate(match_cat = cut(abs(pct_diff), c(-Inf, .5, 2, 5, 100, Inf),
                                     labels = c("exact(<0.5%)","close(0.5-2%)","minor(2-5%)",
                                                "notable(5-100%)","gross(>100%)"))) %>%
              arrange(desc(abs(pct_diff))),
            file.path(output_dir, "species_level_differences.csv"))
  write_csv(raw_cells %>% filter(!is.na(anomaly)) %>%
              select(Species, stephan_species, label, stephan_var, wide_var,
                     stephan_value, wide_value, ratio, anomaly) %>%
              arrange(anomaly, stephan_var, Species),
            file.path(output_dir, "unit_laterality_anomalies.csv"))

  cat("Cells compared      :", sum(both_present), "\n")
  cat("  same              :", sum(raw_cells$status == "same"), "\n")
  cat("  rounding          :", sum(raw_cells$status == "rounding"), "\n")
  cat("  DIFFERENT         :", sum(raw_cells$status == "DIFFERENT"),
      " (of which", sum(raw_cells$status == "DIFFERENT" & raw_cells$role == "report_only"),
      "are report_only pairs)\n")
  cat("  >0.1% (legacy raw):", sum(!is.na(raw_cells$rel) & raw_cells$rel > raw_rel_tol), "\n")
  cat("Unit/laterality     :", sum(!is.na(raw_cells$anomaly)), "\n")
  print(summarise_cells(raw_cells, both_present) %>%
          select(label, n_both, n_same, n_DIFFERENT, median_ratio, max_abs_pct, role), n = Inf)
}

# ============================================================
# 3  MODE "model" -- the matrices the s3 scripts fitted
# ============================================================
if ("model" %in% modes) {
  cat("\n============================================================\n")
  cat("MODE model : the matrices the two s3 scripts actually fitted\n")
  cat("============================================================\n")
  for (f in c(used_stephan_path, used_wide_path))
    if (!file.exists(f)) stop("Model-input file not found: ", f)

  ustephan <- read.csv(used_stephan_path, check.names = FALSE, stringsAsFactors = FALSE)
  uwide    <- read.csv(used_wide_path,    check.names = FALSE, stringsAsFactors = FALSE)
  tree_tips <- read_tree_tips()

  # model-input column names: the crosswalk's model_* override where given
  cwm <- compared %>%
    mutate(s_col = coalesce(model_stephan, stephan_var),
           w_col = coalesce(model_wide, wide_var))
  miss_s <- setdiff(c("Species", cwm$s_col), names(ustephan))
  miss_w <- setdiff(c("Species", cwm$w_col), names(uwide))
  if (length(miss_s)) stop("Missing column(s) from ", used_stephan_path, ": ", paste(miss_s, collapse = ", "))
  if (length(miss_w)) stop("Missing column(s) from ", used_wide_path, ": ", paste(miss_w, collapse = ", "))
  for (nm in c("ustephan", "uwide")) {
    d <- get(nm)
    if (anyDuplicated(d$Species))
      warning("Duplicate Species rows in ", nm, ": ",
              paste(unique(d$Species[duplicated(d$Species)]), collapse = ", "))
  }

  all_species_raw <- unique(c(ustephan$Species, uwide$Species))
  all_species <- if (RESTRICT_TO_TREE && !is.null(tree_tips))
    all_species_raw[all_species_raw %in% tree_tips] else all_species_raw
  idx_s <- match(all_species, ustephan$Species); idx_w <- match(all_species, uwide$Species)
  on_tree <- if (is.null(tree_tips)) rep(NA_character_, length(all_species))
             else ifelse(all_species %in% tree_tips, "yes", "PRUNED")

  model_cells <- pmap_dfr(cwm %>% select(label, stephan_var, wide_var, unit_scale, role, s_col, w_col),
                          function(label, stephan_var, wide_var, unit_scale, role, s_col, w_col) {
    a <- num(ustephan[[s_col]][idx_s]) * unit_scale
    b <- num(uwide[[w_col]][idx_w])
    bind_cols(
      tibble(mode = "model", Species = all_species, stephan_species = all_species,
             pair_note = NA_character_, label = label, stephan_var = s_col, wide_var = w_col,
             unit_scale = unit_scale, role = role, stephan_value = a, wide_value = b,
             on_tree = on_tree),
      classify(a, b))
  }) %>% filter(status != "both_missing")
  model_cells$anomaly <- anomaly_of(model_cells$ratio)
  all_cells_both[["model"]] <- model_cells

  # side-by-side wide view, one block of columns per variable (as before)
  comparison <- tibble(Species = all_species)
  if (!is.null(tree_tips)) comparison$on_tree <- on_tree
  comparison$in_Stephan      <- ifelse(all_species %in% ustephan$Species, "yes", "-")
  comparison$in_volumes_wide <- ifelse(all_species %in% uwide$Species, "yes", "-")
  for (i in seq_len(nrow(cwm))) {
    base <- gsub("^_|_$", "", gsub("[^A-Za-z0-9]+", "_", cwm$label[i]))
    b_i <- model_cells %>% filter(label == cwm$label[i])
    j <- match(all_species, b_i$Species)
    comparison[[paste0(base, "_Stephan")]]      <- b_i$stephan_value[j]
    comparison[[paste0(base, "_volumes_wide")]] <- b_i$wide_value[j]
    comparison[[paste0(base, "_pct_diff")]]     <- round(b_i$pct_diff[j], 2)
    comparison[[paste0(base, "_status")]]       <- b_i$status[j]
  }
  write.csv(comparison, file.path(output_dir, "Stephan_vs_volumes_wide_side_by_side.csv"),
            row.names = FALSE, na = "")

  cell_disc <- model_cells %>%
    filter(status %in% c("Stephan_only", "volumes_wide_only", "DIFFERENT")) %>%
    transmute(Species, variable = label, col_Stephan = stephan_var, col_wide = wide_var,
              Stephan = stephan_value, volumes_wide = wide_value,
              pct_diff = round(pct_diff, 2), status, on_tree, role) %>%
    arrange(on_tree != "yes", status != "DIFFERENT", desc(abs(replace_na(pct_diff, 0))))
  write.csv(cell_disc, file.path(output_dir, "Stephan_vs_volumes_wide_cell_discrepancies.csv"),
            row.names = FALSE, na = "")

  in_model <- if (is.null(tree_tips)) rep(TRUE, nrow(model_cells)) else model_cells$on_tree == "yes"
  var_summary <- summarise_cells(model_cells, in_model) %>%
    select(variable = label, N_Stephan, N_volumes_wide, n_shared = n_both,
           n_Stephan_only, n_wide_only, n_DIFFERENT, max_abs_pct_diff = max_abs_pct, role)
  write.csv(var_summary, file.path(output_dir, "Stephan_vs_volumes_wide_variable_summary.csv"),
            row.names = FALSE, na = "")

  species_coverage <- tibble(
    Species = all_species_raw,
    in_Stephan = all_species_raw %in% ustephan$Species,
    in_volumes_wide = all_species_raw %in% uwide$Species,
    on_tree = if (is.null(tree_tips)) NA else all_species_raw %in% tree_tips) %>%
    mutate(verdict = case_when(
      !is.na(on_tree) & !on_tree ~ "off-tree - out of scope, resolve upstream",
      in_Stephan & !in_volumes_wide ~ "on tree but absent from volumes_wide - LOST",
      !in_Stephan & in_volumes_wide ~ "extra in volumes_wide",
      TRUE ~ "in both")) %>%
    arrange(match(verdict, c("on tree but absent from volumes_wide - LOST",
                             "extra in volumes_wide", "in both",
                             "off-tree - out of scope, resolve upstream")), Species)
  write.csv(species_coverage, file.path(output_dir, "Stephan_vs_volumes_wide_species_coverage.csv"),
            row.names = FALSE, na = "")

  cat("Rows read      : Stephan ", nrow(ustephan), " | volumes_wide ", nrow(uwide), "\n", sep = "")
  if (!is.null(tree_tips)) {
    cat("Tree tips      : ", length(tree_tips), "\n", sep = "")
    cat("Modelled       : Stephan ", sum(ustephan$Species %in% tree_tips),
        " | volumes_wide ", sum(uwide$Species %in% tree_tips), "\n", sep = "")
  }
  lost <- species_coverage$Species[species_coverage$verdict == "on tree but absent from volumes_wide - LOST"]
  if (length(lost)) {
    cat("\n--- Tree species in Stephan but ABSENT from volumes_wide (dropped from every wide fit) ---\n")
    for (sp in lost) cat("      ", sp, "\n")
  }
  offtree <- species_coverage$Species[species_coverage$verdict == "off-tree - out of scope, resolve upstream"]
  if (length(offtree))
    cat("\n--- Off-tree rows:", length(offtree), "--- (never modelled; resolve names upstream)\n")
  cat("\n--- Per-variable effective N (tree tips only) ---\n")
  print(as.data.frame(var_summary), row.names = FALSE)
  onmodel <- cell_disc %>% filter(is.na(on_tree) | on_tree == "yes", status == "DIFFERENT")
  cat("\n--- Cell discrepancies affecting the models:", nrow(onmodel), "---\n")
  if (nrow(onmodel)) print(as.data.frame(head(onmodel, 20)), row.names = FALSE)
}

# ============================================================
# 4  Combined long table
# ============================================================
if (length(all_cells_both)) {
  write_csv(bind_rows(all_cells_both), file.path(output_dir, "compare_all_cells.csv"))
  cat("\nWrote ", file.path(output_dir, "compare_all_cells.csv"), " (",
      nrow(bind_rows(all_cells_both)), " rows, modes: ",
      paste(names(all_cells_both), collapse = " + "), ")\n", sep = "")
}
cat("\nTolerances: same <= ", tol_same, " | rounding <= ", tol_flag,
    " | DIFFERENT above that\n", sep = "")
