# ============================================================
# QC: compare the manually curated data_raw/Stephan_primates.csv against the
# script-built merge, __merging_volumes/volumes_wide_select.csv.
#
# Cell by cell (species x structure), reports both values, absolute and percent
# difference, and the wide/Stephan ratio.  The ratio column is the point of this
# script: a ratio near 2 means a laterality mismatch, near 1000 a unit mismatch,
# near 1 genuine agreement.
#
# Outputs (checks/compare_volumes_files/):
#   coalesced_species_rows.csv          rows merged because two labels = one taxon
#   species_only_in_volumes_wide.csv    species present only in the merge
#   species_only_in_stephan.csv         species present only in Stephan
#   variables_only_in_volumes_wide.csv  unmapped merge columns
#   variables_only_in_stephan.csv       unmapped Stephan columns
#   comparison_summary.csv              one row per mapped variable
#   species_level_differences.csv       one row per disagreeing cell
#   unit_laterality_anomalies.csv       cells whose ratio implies a unit/2x error
#
# Run from anywhere inside the repo.
# ============================================================

setwd(local({ d <- normalizePath(getwd()); while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d); d }))  # repo root (portable; see R/project_root.R)

library(tidyverse)

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

wide_path <- "/Users/crossmodal/Library/CloudStorage/OneDrive-AllenInstitute/Species/Evo-M1-Trait-Data/__merging_volumes/volumes_wide_select.csv"
stephan_path <- "data_raw/Stephan_primates.csv"

out_dir <- "checks/compare_volumes_files"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Relative tolerance. Cells differing by more than this fraction are reported.
rel_tol <- 0.001   # 0.1%

if (!file.exists(wide_path)) {
  stop("volumes_wide_select.csv not found at:\n  ", wide_path,
       "\nIs OneDrive mounted / the file synced?")
}
if (!file.exists(stephan_path)) {
  stop("Stephan_primates.csv not found at:\n  ", file.path(getwd(), stephan_path))
}

# ------------------------------------------------------------
# Read files
# ------------------------------------------------------------

wide <- read_csv(wide_path, show_col_types = FALSE, name_repair = "minimal")
stephan <- read_csv(stephan_path, show_col_types = FALSE, name_repair = "minimal")

stopifnot("Species" %in% names(wide), "Species" %in% names(stephan))

is_blank <- function(x) is.na(x) | trimws(as.character(x)) == ""

# ------------------------------------------------------------
# Standardise species names
#
# Underscores -> spaces, then rewrite superseded labels to the accepted name.
# Pongo is deliberately NOT handled here -- see extra_species_pairs below.
# ------------------------------------------------------------

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

wide <- wide %>%
  filter(!is_blank(Species)) %>%
  mutate(Species_std = norm_species(Species))

stephan <- stephan %>%
  filter(!is_blank(Species)) %>%
  mutate(Species_std = norm_species(Species))

# ------------------------------------------------------------
# Coalesce rows that collapse onto one accepted name
#
# volumes_wide_select.csv carries the same taxon twice in at least one case
# (Callicebus moloch, populated; Plecturocebus moloch, near-empty).  Merge them
# field by field, stop if the two rows actually disagree, and record what was
# merged so the split in the source file stays visible.
# ------------------------------------------------------------

coalesce_duplicate_rows <- function(dat, label) {
  dups <- unique(dat$Species_std[duplicated(dat$Species_std)])
  if (length(dups) == 0L) {
    return(list(data = dat, report = tibble(), conflicts = tibble()))
  }

  report <- tibble()
  conflicts <- tibble()
  keep_rows <- rep(TRUE, nrow(dat))

  for (sp in dups) {
    idx <- which(dat$Species_std == sp)
    block <- dat[idx, , drop = FALSE]

    report <- bind_rows(report, tibble(
      file = label,
      Species_std = sp,
      n_rows_merged = length(idx),
      labels_merged = paste(block$Species, collapse = " | "),
      populated_fields_per_row = paste(
        vapply(seq_len(nrow(block)),
               function(i) sum(!is_blank(unlist(block[i, ]))), integer(1)),
        collapse = " | ")
    ))

    # The label columns are expected to differ -- that is the whole point of the
    # merge -- so they are recorded in labels_merged, not treated as conflicts.
    for (col in setdiff(names(dat), c("Species", "Species_std"))) {
      vals <- block[[col]]
      good <- unique(vals[!is_blank(vals)])
      if (length(good) > 1L) {
        conflicts <- bind_rows(conflicts, tibble(
          file = label, Species_std = sp, column = col,
          values = paste(format(good), collapse = " | ")))
      }
      if (length(good) >= 1L) dat[[col]][idx[1]] <- good[1]
    }
    dat[["Species"]][idx[1]] <- sp
    keep_rows[idx[-1]] <- FALSE
  }

  list(data = dat[keep_rows, , drop = FALSE], report = report, conflicts = conflicts)
}

wide_co <- coalesce_duplicate_rows(wide, "volumes_wide_select")
stephan_co <- coalesce_duplicate_rows(stephan, "Stephan_primates")

coalesce_report <- bind_rows(wide_co$report, stephan_co$report)
coalesce_conflicts <- bind_rows(wide_co$conflicts, stephan_co$conflicts)

write_csv(coalesce_report, file.path(out_dir, "coalesced_species_rows.csv"))

if (nrow(coalesce_conflicts) > 0) {
  print(coalesce_conflicts, n = Inf)
  stop("Rows that collapse onto one accepted name hold conflicting values ",
       "(see printout above). Resolve in the source file before comparing.")
}

if (nrow(coalesce_report) > 0) {
  message("Coalesced ", nrow(coalesce_report), " duplicated taxon label(s):")
  for (i in seq_len(nrow(coalesce_report))) {
    message("  ", coalesce_report$file[i], ": ", coalesce_report$labels_merged[i],
            " -> ", coalesce_report$Species_std[i],
            "  (populated fields per row: ",
            coalesce_report$populated_fields_per_row[i], ")")
  }
}

wide <- wide_co$data
stephan <- stephan_co$data

# ------------------------------------------------------------
# Species pairing
#
# Identity pairs for names that match after normalisation, plus explicit pairs
# where the two files disagree about how finely a taxon is split.
#
# Pongo: the merge holds a single "Pongo sp." row, which is a sensu lato mean
# (Bornean + Sumatran + hybrid material).  Stephan splits abelii from pygmaeus.
# Compare against BOTH, flagged, rather than silently pinning it to one.
# ------------------------------------------------------------

extra_species_pairs <- tribble(
  ~wide_species, ~stephan_species,  ~pair_note,
  "Pongo sp.",   "Pongo abelii",    "merge 'Pongo sp.' is a sensu lato mean (Bornean + Sumatran + hybrid); Stephan row is Sumatran",
  "Pongo sp.",   "Pongo pygmaeus",  "merge 'Pongo sp.' is a sensu lato mean (Bornean + Sumatran + hybrid); Stephan row is Bornean"
) %>%
  filter(wide_species %in% wide$Species_std,
         stephan_species %in% stephan$Species_std)

species_pairs <- bind_rows(
  tibble(wide_species = intersect(wide$Species_std, stephan$Species_std)) %>%
    mutate(stephan_species = wide_species, pair_note = NA_character_),
  extra_species_pairs
) %>%
  arrange(wide_species, stephan_species)

# ------------------------------------------------------------
# Species coverage
# ------------------------------------------------------------

species_only_wide <- setdiff(wide$Species_std, species_pairs$wide_species)
species_only_stephan <- setdiff(stephan$Species_std, species_pairs$stephan_species)

cat("\nSpecies pairs compared:", nrow(species_pairs),
    " (", sum(!is.na(species_pairs$pair_note)), " flagged )\n", sep = "")
cat("Only in volumes_wide_select:", length(species_only_wide), "\n")
cat("Only in Stephan_primates:", length(species_only_stephan), "\n")

write_csv(tibble(Species = sort(species_only_wide)),
          file.path(out_dir, "species_only_in_volumes_wide.csv"))
write_csv(tibble(Species = sort(species_only_stephan)),
          file.path(out_dir, "species_only_in_stephan.csv"))

# ------------------------------------------------------------
# Variable crosswalk
#
# unit_scale = multiplier applied to the STEPHAN value to put it in the merge's
# units.  1 everywhere except the Smaers cortical columns, which Stephan holds
# in cm^3 while the merge holds mm^3.
#
# Laterality: Stephan's vestibular columns are unilateral, so they are compared
# against the merge's *_unilateral_* columns (the bilateral merge columns are
# exactly 2x).  Likewise Stephan's insula columns are LEFT only, so they are
# compared against *_left_* -- Insula_Vol.mm3 is bilateral (~2x).
# ------------------------------------------------------------

crosswalk <- tribble(
  ~stephan_var,                                    ~wide_var,                                                ~unit_scale,

  # --- body and brain size ---
  "Body_weight",                                   "Body_Mass.g",                                                 1,
  "Brain_weight",                                  "Brain_Mass.mg",                                               1,
  "Total_brain_net_volume",                        "Total_brain_net_volume_Vol.mm3",                              1,
  "Ventricles",                                    "Ventricles_Vol.mm3",                                          1,

  # --- major subdivisions ---
  "Medulla_oblongata",                             "Medulla_oblongata_Vol.mm3",                                   1,
  "Cerebellum",                                    "Cerebellum_Vol.mm3",                                          1,
  "Mesencephalon",                                 "Mesencephalon_Vol.mm3",                                       1,
  "Diencephalon",                                  "Diencephalon_Vol.mm3",                                        1,
  "Telencephalon",                                 "Telencephalon_Vol.mm3",                                       1,
  "Pons",                                          "Pons_Vol.mm3",                                                1,

  # --- olfactory ---
  "Bulbus_olfactorius",                            "Bulbus_olfactorius_Vol.mm3",                                  1,
  "Bulbus_olfactorius_accessorius",                "Bulbus_olfactorius_accessorius_Vol.mm3",                      1,
  "Nucleus_tractus_olfactorius",                   "Nucleus_tractus_olfactorius_Vol.mm3",                         1,

  # --- telencephalic ---
  "Lobus_piriformis",                              "Lobus_piriformis_Vol.mm3",                                    1,
  "Septum",                                        "Septum_Vol.mm3",                                              1,
  "Striatum",                                      "Striatum_Vol.mm3",                                            1,
  "Schizocortex",                                  "Schizo_cortex_Vol.mm3",                                       1,
  "Hippocampus",                                   "Hippocampus_Vol.mm3",                                         1,
  "Palaeocortex",                                  "Palaeocortex_Vol.mm3",                                        1,
  "NeoWG",                                         "Neocortex_Vol.mm3",                                           1,
  "NeoG_Frahm",                                    "Neocortex_grey_matter_Vol.mm3",                               1,
  "NeoW_Frahm",                                    "Neocortex_white_matter_Vol.mm3",                              1,

  # --- amygdala ---
  "Amygdala",                                      "Amygdala_Vol.mm3",                                            1,
  "Complexus_centromedialis",                      "Complexus_centromedialis_Vol.mm3",                            1,
  "Complexus_cortico_basolateralis",               "Complexus_corticobasolateralis_Vol.mm3",                      1,
  "Nucleus_amygdalae_basalis_parsmagnocellularis", "Nucleus_amygdalae_basalis_pars_magnocellularis_Vol.mm3",      1,

  # --- diencephalic ---
  "Epithalamus",                                   "Epithalamus_Vol.mm3",                                         1,
  "Thalamus",                                      "Thalamus_Vol.mm3",                                            1,
  "Hypothalamus",                                  "Hypothalamus_Vol.mm3",                                        1,
  "Subthalamus",                                   "Subthalamus_Vol.mm3",                                         1,
  "Pallidum",                                      "Pallidum_Vol.mm3",                                            1,
  "Nucleus_subthalamicus",                         "Nucleus_subthalamicus_Vol.mm3",                               1,
  "Capsula_interna",                               "Capsula_interna_Vol.mm3",                                     1,
  "Tractus_opticus",                               "Tractus_opticus_Vol.mm3",                                     1,

  # --- visual (de Sousa) ---
  "ASG_Sousa",                                     "Area_striata_grey_matter_Vol.mm3",                            1,
  "LGN_Sousa",                                     "Corpus_geniculatum_laterale_Vol.mm3",                         1,

  # --- cerebellar nuclei (Matano 1985a) ---
  "Cerebellar_nuclei_total",                       "Cerebellar_nuclei_total_Vol.mm3",                             1,
  "Interpositus_cerebellar_nuclei",                "Interpositus_cerebellar_nuclei_Vol.mm3",                      1,
  "Lateral_cerebellar_nuclei",                     "Lateral_cerebellar_nuclei_Vol.mm3",                           1,
  "Medial_cerebellar_nuclei",                      "Medial_cerebellar_nuclei_Vol.mm3",                            1,

  # --- vestibular complex: Stephan is UNILATERAL ---
  "Complexus.vestibularis",                        "Complexus_vestibularis_unilateral_Vol.mm3",                   1,
  "Nucleus.vestibularis.superior",                 "Nucleus_vestibularis_superior_unilateral_Vol.mm3",            1,
  "Nucleus.vestibularis.lateralis",                "Nucleus_vestibularis_lateralis_unilateral_Vol.mm3",           1,
  "Nucleus.vestibularis.medialis",                 "Nucleus_vestibularis_medialis_unilateral_Vol.mm3",            1,
  "Nucleus.vestibularis.descendens",               "Nucleus_vestibularis_descendens_unilateral_Vol.mm3",          1,

  # --- periventricular organs ---
  "Nucleus.septalis.triangularis",                 "Nucleus_septalis_triangularis_Vol.mm3",                       1,
  "Corpus.subfornicale",                           "Corpus_subfornicale_Vol.mm3",                                 1,
  "Nucleus.habenularis.medialis",                  "Nucleus_habenularis_medialis_Vol.mm3",                        1,
  "Corpus.pineale",                                "Corpus_pineal_Vol.mm3",                                       1,
  "Corpus.subcommissurale",                        "Corpus_subcommissurale_Vol.mm3",                              1,

  # --- insula: Stephan is LEFT hemisphere only ---
  "Total_insula_volume_L",                         "Insula_left_Vol.mm3",                                         1,
  "Granular_volume_L",                             "Granular_insular_cortex_left_Vol.mm3",                        1,
  "Dysgranular_volume_L",                          "Dysgranular_insular_cortex_left_Vol.mm3",                     1,
  "Agranular_volume_L",                            "Agranular_insular_cortex_left_Vol.mm3",                       1,
  "FI_volume_L",                                   "fronto_insular_cortex_left_Vol.mm3",                          1,

  # --- Smaers cortical parcels: Stephan in cm^3, merge in mm^3 ---
  "Frontal.motor.Gray",                            "FrontalMotor_grey_matter_Vol.mm3",                         1000,
  "Frontal.motor",                                 "FrontalMotor_white_matter_Vol.mm3",                        1000
)

# Keep only pairs where BOTH columns actually exist.
# NOTE: this must not be written as `filter(wide %in% names(wide))` -- inside
# filter(), `wide` resolves to the crosswalk's own column, not to the data
# frame, so names() returns NULL and every row is silently dropped.
wide_names <- names(wide)
stephan_names <- names(stephan)

crosswalk_missing <- crosswalk %>%
  filter(!(wide_var %in% wide_names) | !(stephan_var %in% stephan_names)) %>%
  mutate(missing_from = case_when(
    !(wide_var %in% wide_names) & !(stephan_var %in% stephan_names) ~ "both",
    !(wide_var %in% wide_names) ~ "volumes_wide_select",
    TRUE ~ "Stephan"
  ))

if (nrow(crosswalk_missing) > 0) {
  warning("Dropping ", nrow(crosswalk_missing), " crosswalk pair(s) whose columns are absent:\n",
          paste0("  ", crosswalk_missing$stephan_var, " <-> ",
                 crosswalk_missing$wide_var, "  (missing from ",
                 crosswalk_missing$missing_from, ")", collapse = "\n"))
}

crosswalk <- crosswalk %>%
  filter(wide_var %in% wide_names, stephan_var %in% stephan_names)

stopifnot(nrow(crosswalk) > 0)

# ------------------------------------------------------------
# Variable coverage
# ------------------------------------------------------------

non_data_cols <- c("Species", "Species_std", "Species_Stephan_file",
                   "Code_number_Stephan", "Species_Matano1985a",
                   "code_Matano1985a", "order")

write_csv(
  tibble(variable = sort(setdiff(wide_names, c(crosswalk$wide_var, non_data_cols)))),
  file.path(out_dir, "variables_only_in_volumes_wide.csv")
)
write_csv(
  tibble(variable = sort(setdiff(stephan_names, c(crosswalk$stephan_var, non_data_cols)))),
  file.path(out_dir, "variables_only_in_stephan.csv")
)

# ------------------------------------------------------------
# Cell-by-cell comparison
#
# Built long, one variable at a time, through species_pairs.  Do NOT
# inner_join() the two wide frames and then reach for paste0(var, ".wide"):
# dplyr's `suffix` is applied only to columns whose names collide, and almost
# none of these do, so those lookups return NULL.
# ------------------------------------------------------------

num <- function(x) suppressWarnings(as.numeric(gsub(",", "", as.character(x), fixed = TRUE)))

all_cells <- pmap_dfr(crosswalk, function(stephan_var, wide_var, unit_scale) {
  s_dat <- tibble(stephan_species = stephan$Species_std,
                  stephan_value = num(stephan[[stephan_var]]) * unit_scale)
  w_dat <- tibble(wide_species = wide$Species_std,
                  wide_value = num(wide[[wide_var]]))

  species_pairs %>%
    left_join(s_dat, by = "stephan_species") %>%
    left_join(w_dat, by = "wide_species") %>%
    filter(!is.na(stephan_value), !is.na(wide_value)) %>%
    transmute(
      Species = wide_species,
      stephan_species,
      pair_note,
      stephan_var = .env$stephan_var,
      wide_var = .env$wide_var,
      unit_scale = .env$unit_scale,
      stephan_value,
      wide_value,
      difference = wide_value - stephan_value,
      abs_difference = abs(difference),
      pct_diff = if_else(stephan_value == 0, NA_real_,
                         100 * (wide_value - stephan_value) / stephan_value),
      ratio = if_else(stephan_value == 0, NA_real_, wide_value / stephan_value)
    )
})

# ------------------------------------------------------------
# Per-variable summary
# ------------------------------------------------------------

safe_max <- function(x) if (length(x) == 0 || all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
safe_cor <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  if (sum(ok) < 3 || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  cor(x[ok], y[ok])
}

comparison_results <- all_cells %>%
  group_by(stephan_var, wide_var, unit_scale) %>%
  summarise(
    n_both          = n(),
    n_agree         = sum(abs(pct_diff) <= rel_tol * 100, na.rm = TRUE),
    n_disagree      = sum(abs(pct_diff) > rel_tol * 100, na.rm = TRUE),
    cor             = safe_cor(stephan_value, wide_value),
    median_ratio    = median(ratio, na.rm = TRUE),
    median_abs_pct  = median(abs(pct_diff), na.rm = TRUE),
    max_abs_pct     = safe_max(abs(pct_diff)),
    mean_abs_diff   = mean(abs_difference, na.rm = TRUE),
    max_abs_diff    = safe_max(abs_difference),
    .groups = "drop"
  ) %>%
  arrange(desc(max_abs_pct))

write_csv(comparison_results, file.path(out_dir, "comparison_summary.csv"))

# ------------------------------------------------------------
# Species-level disagreements (relative tolerance)
# ------------------------------------------------------------

all_diffs <- all_cells %>%
  filter(is.na(pct_diff) | abs(pct_diff) > rel_tol * 100) %>%
  mutate(match_cat = cut(abs(pct_diff), c(-Inf, 0.5, 2, 5, 100, Inf),
                         labels = c("exact(<0.5%)", "close(0.5-2%)", "minor(2-5%)",
                                    "notable(5-100%)", "gross(>100%)"))) %>%
  arrange(desc(abs(pct_diff)))

write_csv(all_diffs, file.path(out_dir, "species_level_differences.csv"))

# ------------------------------------------------------------
# Unit / laterality anomalies
#
# A ratio clustering near 1000, 0.001, 2 or 0.5 is a systematic error, not a
# measurement disagreement.  Surfaced separately so it does not hide in the tail
# of species_level_differences.csv.
# ------------------------------------------------------------

near_ratio <- function(x, target, tol = 0.05) !is.na(x) & abs(x / target - 1) < tol

anomalies <- all_cells %>%
  mutate(anomaly = case_when(
    near_ratio(ratio, 1000)  ~ "merge is ~1000x Stephan (cm3/mm3 or g/mg)",
    near_ratio(ratio, 0.001) ~ "merge is ~1/1000x Stephan (cm3/mm3 or g/mg)",
    near_ratio(ratio, 2)     ~ "merge is ~2x Stephan (bilateral vs unilateral)",
    near_ratio(ratio, 0.5)   ~ "merge is ~0.5x Stephan (unilateral vs bilateral)",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(anomaly)) %>%
  select(Species, stephan_species, stephan_var, wide_var,
         stephan_value, wide_value, ratio, anomaly) %>%
  arrange(anomaly, stephan_var, Species)

write_csv(anomalies, file.path(out_dir, "unit_laterality_anomalies.csv"))

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

cat("\nVariables compared:", nrow(comparison_results), "\n")
cat("Cells compared:", nrow(all_cells), "\n")
cat("Cells agreeing within ", rel_tol * 100, "%: ",
    sum(abs(all_cells$pct_diff) <= rel_tol * 100, na.rm = TRUE), "\n", sep = "")
cat("Cells disagreeing:", nrow(all_diffs), "\n")
cat("Unit / laterality anomalies:", nrow(anomalies), "\n")
cat("\nOutputs written to:", normalizePath(out_dir), "\n\n")

print(comparison_results, n = Inf)

if (nrow(anomalies) > 0) {
  cat("\n--- Unit / laterality anomalies ---\n")
  print(anomalies %>% count(stephan_var, wide_var, anomaly), n = Inf)
}
