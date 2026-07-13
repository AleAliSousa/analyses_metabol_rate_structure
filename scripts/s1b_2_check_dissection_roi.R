setwd(local({ d <- normalizePath(getwd()); while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d); d }))  # repo root (portable; replaces hardcoded path -- see R/project_root.R)

############################################################
# Create CSVs to inspect obs terms that map to rCMRGlc terms
############################################################

## Install and Load up packages
library(ggplot2)
library(ggpmisc)
library(readxl)
library(tidyverse)

data_dir <- "data_intermediate"
out_dir <- "checks/anatomy_rule_audit"

## cortical ROIs / Lobes are here: https://github.com/linnarsson-lab/adult-human-brain/blob/main/notebooks/Revision/RevisionFig2.ipynb
## Also See Fig 2 in the original paper: https://www.science.org/doi/10.1126/science.add7046

############################
## Load saved obs metadata
############################
obs1 <- readRDS("data_intermediate/linnarsson_adult_human_brain_obs_metadata_nonneuronal.rds") %>%
  as_tibble() %>%
  mutate(obs_dataset = "nonneuronal")

obs2 <- readRDS("data_intermediate/linnarsson_adult_human_brain_obs_metadata_neuronal.rds") %>%
  as_tibble() %>%
  mutate(obs_dataset = "neuronal")

# Check if the columns are identical
setdiff(colnames(obs1), colnames(obs2))
setdiff(colnames(obs2), colnames(obs1))

obs <- bind_rows(obs1, obs2)

# inspect the anndata object
colnames(obs)
#####################################
# 1. Count the number of observations for each observed roi-dissection pairing.
# This produces a compact long-format table with one row per observed pairing.
roi_dissection_counts <- obs %>%                 
  count(roi, dissection, name = "n") %>%         # Count rows for each unique roi–dissection combination; store the count in a column named n
  arrange(roi, desc(n))                          # Sort by roi, and within each roi put the most common dissection assignments first
showrows<-nrow(roi_dissection_counts)            # Get the exact number of rows that will be printed

# Print exactly all rows in the summarized count table
roi_dissection_counts %>%
  print(n = showrows)

# 2. Count the number of unique rois associated with each dissection, and vice versa. 
# This produces a compact summary table with one row per unique roi and one row per unique dissection.
# Add explicit missingness flags for roi and dissection
obs_check <- obs %>%
  mutate(
    roi_missing = is.na(roi) | str_trim(as.character(roi)) == "",
    dissection_missing = is.na(dissection) | str_trim(as.character(dissection)) == ""
  )

# Count rows with roi present but dissection missing, and vice versa
missing_pair_summary <- obs_check %>%
  summarise(
    n_total = n(),
    n_roi_present_dissection_missing = sum(!roi_missing & dissection_missing),
    n_roi_missing_dissection_present = sum(roi_missing & !dissection_missing),
    n_both_missing = sum(roi_missing & dissection_missing),
    n_both_present = sum(!roi_missing & !dissection_missing)
  )

missing_pair_summary
#n_total n_roi_present_dissection_missing n_roi_missing_dissection_present n_both_missing n_both_present
#<int>                            <int>                            <int>          <int>          <int>
#  1 3369219                             8275                                0              0        3360944
### Conclusion: 
### - Some rows for roi lack dissection
### - All of the dissecytion rows have roi 

# To print the exact affected roi values where dissection is missing:
obs_check %>%
  filter(!roi_missing & dissection_missing) %>%
  count(roi, name = "n_rows") %>%
  arrange(desc(n_rows)) %>%
  print(n = Inf)

# And the reverse, if any:
obs_check %>%
  filter(roi_missing & !dissection_missing) %>%
  count(dissection, name = "n_rows") %>%
  arrange(desc(n_rows)) %>%
  print(n = Inf)

# Now test many-to-one in both directions, using only rows where both fields are present:
pair_map <- obs_check %>%
  filter(!roi_missing, !dissection_missing) %>%
  distinct(roi, dissection)

# One roi mapped to multiple dissections
roi_to_many_dissections <- pair_map %>%
  count(roi, name = "n_dissections") %>%
  filter(n_dissections > 1) %>%
  arrange(desc(n_dissections), roi)

# Print the actual (roi, dissection) pairs for those ROIs
pair_map %>%
  filter(roi %in% roi_to_many_dissections) %>%
  arrange(roi, dissection) %>%
  print(n = Inf, width = Inf)

# One dissection mapped to multiple rois
dissection_to_many_rois <- pair_map %>%
  count(dissection, name = "n_rois") %>%
  filter(n_rois > 1) %>%
  arrange(desc(n_rois), dissection)

roi_to_many_dissections
dissection_to_many_rois

# Filter the original obs data to get the full rows for those cases. This allows us to inspect all metadata fields for those specific observations.
rows_roi_to_many_dissections <- obs_check %>%
  filter(!roi_missing, !dissection_missing, roi %in% roi_to_many_dissections$roi)

# ROI → many dissections (compact list)
rows_roi_to_many_dissections %>%
  distinct(roi, dissection) %>%
  arrange(roi, dissection) %>%
  print(n = Inf, width = Inf)

# Now do the same for dissection → many ROIs
rows_dissection_to_many_rois <- obs_check %>%
  filter(!roi_missing, !dissection_missing, dissection %in% dissection_to_many_rois$dissection)

# Dissection → many ROIs (compact list)
rows_dissection_to_many_rois %>%
  distinct(roi, dissection) %>%
  arrange(roi, dissection) %>%
  print(n = Inf)

# Coverage check: For each ROI that appears in rows with dissection blank, does that ROI also appear somewhere else with a non-blank dissection?
 # Uniqueness check: For those ROIs, when dissection is present, does the ROI map to exactly one dissection (unique), or to multiple?

############################################################
# CHECK: For blank dissection rows, does ROI appear elsewhere,
# and is ROI->dissection mapping unique?
############################################################

# 1) ROIs that have blank dissection rows
rois_with_blank_dissection <- obs_check %>%
  filter(!roi_missing, dissection_missing) %>%
  distinct(roi)

# 2) For those ROIs, check whether they appear anywhere with a non-blank dissection
roi_blank_dissection_coverage <- rois_with_blank_dissection %>%
  left_join(
    obs_check %>%
      filter(!roi_missing, !dissection_missing) %>%
      distinct(roi) %>%
      mutate(appears_with_dissection_elsewhere = TRUE),
    by = "roi"
  ) %>%
  mutate(
    appears_with_dissection_elsewhere = tidyr::replace_na(appears_with_dissection_elsewhere, FALSE)
  ) %>%
  arrange(appears_with_dissection_elsewhere, roi)

cat("\n=== ROIs in blank-dissection rows: do they appear elsewhere with dissection present? ===\n")
roi_blank_dissection_coverage %>% print(n = Inf)

# 3) For the subset that DO appear elsewhere, check if roi->dissection is unique or many-to-one
roi_blank_dissection_mapping <- pair_map %>%                         # already only both-present
  semi_join(rois_with_blank_dissection, by = "roi") %>%             # only ROIs that appear in blank-dissection rows
  count(roi, name = "n_dissections_when_present") %>%
  arrange(desc(n_dissections_when_present), roi)

cat("\n=== For those ROIs (when dissection is present), how many distinct dissections do they map to? ===\n")
roi_blank_dissection_mapping %>% print(n = Inf)

# 4) Explicit lists (the answer you likely want)
rois_only_in_blank_rows <- roi_blank_dissection_coverage %>%
  filter(!appears_with_dissection_elsewhere) %>%
  pull(roi)

rois_blank_rows_but_unique_mapping <- roi_blank_dissection_mapping %>%
  filter(n_dissections_when_present == 1) %>%
  pull(roi)

rois_blank_rows_and_ambiguous_mapping <- roi_blank_dissection_mapping %>%
  filter(n_dissections_when_present > 1) %>%
  pull(roi)

cat("\n=== SUMMARY ===\n")
cat("ROIs that occur ONLY when dissection is blank (no nonblank dissection rows exist):\n")
print(rois_only_in_blank_rows)

cat("\nROIs that have blank-dissection rows BUT map uniquely to 1 dissection when present:\n")
print(rois_blank_rows_but_unique_mapping)

cat("\nROIs that have blank-dissection rows AND map to MULTIPLE dissections when present (ambiguous):\n")
print(rois_blank_rows_and_ambiguous_mapping)

# 5) (Optional) Show the actual dissection(s) for ambiguous ROIs
if (length(rois_blank_rows_and_ambiguous_mapping) > 0) {
  cat("\n=== Ambiguous ROI -> dissection pairs (when dissection is present) ===\n")
  pair_map %>%
    filter(roi %in% rois_blank_rows_and_ambiguous_mapping) %>%
    arrange(roi, dissection) %>%
    print(n = Inf, width = Inf)
}

# RESULTS SUMMARY:
# All ROIs that lack a dissection assignment also occur elsewhere with a defined dissection except two cases (Human A35-36, Human Gpe).
# For all ROIs with both missing and present dissection values, the mapping is uniquely defined (each ROI maps to exactly one dissection when observed).
# No ROI shows one-to-many ambiguity in the subset relevant to missing dissection rows.