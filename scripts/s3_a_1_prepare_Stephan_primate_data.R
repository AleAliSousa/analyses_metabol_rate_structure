# ============================================================
# Study 3 - Prepare Stephan primate brain-region data
#
# Purpose:
#   1. Clean the Stephan primate dataset.
#   2. Create a unified preferred total-brain-volume variable.
#   3. Standardise anatomical terms across the Stephan and Heiss datasets.
#   4. Create analysis-ready wide and long data files.
#   5. Export the species list for phylogenetic-tree preparation/checking.
#
# This script performs data preparation only.
# It does not fit models, run analyses, or create plots.
# ============================================================

library(here)
library(tidyverse)

options(scipen = 999)

# ============================================================
# FILE PATHS
# ============================================================
input_stephan <- here("data_raw", "Stephan_primates.csv")
input_heiss   <- here("data_intermediate", "Heiss_Stephan_data.csv")
output_dir    <- here("data_intermediate")

output_wide_csv    <- file.path(output_dir, "s3_Stephan_primate_data_prepared.csv")
output_wide_rds    <- file.path(output_dir, "s3_Stephan_primate_data_prepared.rds")
output_long_csv    <- file.path(output_dir, "s3_Stephan_primate_regions_long.csv")
output_long_rds    <- file.path(output_dir, "s3_Stephan_primate_regions_long.rds")
output_species_csv <- file.path(output_dir, "s3_Stephan_species_list.csv")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# ANATOMICAL TERM CROSSWALK
# ============================================================
region_crosswalk <- tribble(
  ~raw_region,                    ~Structure,
  "LGN_Sousa",                   "Corpus geniculatum laterale",
  "Amygdala",                    "Amygdala",
  "Pallidum",                    "Pallidum",
  "NeoW_Frahm",                  "Neocortex white",
  "Total_insula_volume_L",       "Insular cortex (grey)",
  "Nucleus_subthalamicus",       "Nucleus subthalamicus Luysi",
  "Capsula_interna",             "Capsula interna",
  "Striatum",                    "Striatum",
  "ASG_Sousa",                   "Area striata grey",
  "NeoG_Frahm",                  "Neocortex grey",
  "Mesencephalon",               "Mesencephalon",
  "Cerebellum",                  "Cerebellum",
  "Hippocampus",                 "Hippocampus",
  "Lateral_cerebellar_nuclei",  "Nucleus dentatus cerebelli"
)

target_cols <- region_crosswalk$raw_region

# ============================================================
# LOAD AND CHECK INPUTS
# ============================================================
if (!file.exists(input_stephan)) stop("File not found: ", input_stephan)
if (!file.exists(input_heiss)) stop("File not found: ", input_heiss)

Stephan_primates <- read_csv(input_stephan, show_col_types = FALSE)
heiss_stephan_tbl <- read_csv(input_heiss, show_col_types = FALSE)

required_stephan <- c(
  "Species", "Brain_volume", "Brainvol", "Total_brain_net_volume",
  target_cols
)
required_heiss <- c("volume_term", "rCMRGlc_mean_both_hemispheres")

missing_stephan <- setdiff(required_stephan, names(Stephan_primates))
missing_heiss <- setdiff(required_heiss, names(heiss_stephan_tbl))

if (length(missing_stephan) > 0) {
  stop("Missing columns in Stephan_primates.csv: ",
       paste(missing_stephan, collapse = ", "))
}
if (length(missing_heiss) > 0) {
  stop("Missing columns in Heiss_Stephan_data.csv: ",
       paste(missing_heiss, collapse = ", "))
}

# ============================================================
# CLEAN AND PREPARE WIDE DATA
# ============================================================
s3_Stephan_prepared <- Stephan_primates %>%
  select(-any_of(c("X", "order"))) %>%
  filter(!is.na(Species), trimws(Species) != "") %>%
  mutate(
    Species = trimws(as.character(Species)),
    Preferred_brain_volume = coalesce(
      as.numeric(Brain_volume),
      as.numeric(Brainvol),
      as.numeric(Total_brain_net_volume)
    )
  ) %>%
  select(Species, Preferred_brain_volume, all_of(target_cols), everything())

missing_predictor_species <- s3_Stephan_prepared %>%
  filter(is.na(Preferred_brain_volume)) %>%
  pull(Species)

if (length(missing_predictor_species) > 0) {
  warning(
    length(missing_predictor_species),
    " species lack Brain_volume, Brainvol, and Total_brain_net_volume: ",
    paste(missing_predictor_species, collapse = ", ")
  )
}

# ============================================================
# MATCH ANATOMICAL TERMS AND CREATE LONG DATA
# ============================================================
heiss_region_data <- heiss_stephan_tbl %>%
  transmute(
    Structure = trimws(as.character(volume_term)),
    rCMRGlc = as.numeric(rCMRGlc_mean_both_hemispheres)
  ) %>%
  distinct(Structure, .keep_all = TRUE)

s3_Stephan_regions_long <- s3_Stephan_prepared %>%
  select(Species, Preferred_brain_volume, all_of(target_cols)) %>%
  pivot_longer(
    cols = all_of(target_cols),
    names_to = "raw_region",
    values_to = "Region_volume"
  ) %>%
  left_join(region_crosswalk, by = "raw_region") %>%
  left_join(heiss_region_data, by = "Structure") %>%
  arrange(Species, match(raw_region, target_cols))

# Report anatomical terms that did not match the Heiss table.
unmatched_structures <- s3_Stephan_regions_long %>%
  distinct(raw_region, Structure, rCMRGlc) %>%
  filter(is.na(rCMRGlc)) %>%
  select(raw_region, Structure)

if (nrow(unmatched_structures) > 0) {
  warning(
    "No rCMRGlc match was found for: ",
    paste(unmatched_structures$Structure, collapse = ", ")
  )
}

# ============================================================
# SPECIES LIST FOR TREE PREPARATION/CHECKING
# ============================================================
s3_Stephan_species <- s3_Stephan_prepared %>%
  distinct(Species) %>%
  arrange(Species)

# ============================================================
# SAVE OUTPUTS
# ============================================================
write_csv(s3_Stephan_prepared, output_wide_csv, na = "")
saveRDS(s3_Stephan_prepared, output_wide_rds)

write_csv(s3_Stephan_regions_long, output_long_csv, na = "")
saveRDS(s3_Stephan_regions_long, output_long_rds)

write_csv(s3_Stephan_species, output_species_csv, na = "")

message("Saved prepared wide data: ", output_wide_csv)
message("Saved prepared long data: ", output_long_csv)
message("Saved species list: ", output_species_csv)
