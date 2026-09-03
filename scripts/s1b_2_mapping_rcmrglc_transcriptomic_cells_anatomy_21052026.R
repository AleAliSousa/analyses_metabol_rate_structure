# =====================================================================
# Validate the single authoritative s1b Linnarsson-ROI -> Heiss map.
#
# Mapping decisions live only in:
#   metadata/anatomy/s1b_linnarsson_roi_map.csv
#
# This script performs coverage and consistency checks against the extracted
# neuronal and non-neuronal observation metadata. It does not discover or
# silently rewrite anatomical mappings.
# =====================================================================

setwd(local({
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d)
  d
}))

suppressPackageStartupMessages(library(tidyverse))
source("helpers/read_heiss_rates.R")

map_path <- "metadata/anatomy/s1b_linnarsson_roi_map.csv"
coverage_path <- "data_intermediate/s1b_roi_mapping_coverage.csv"

obs_neuronal <- readRDS(
  "data_intermediate/linnarsson_adult_human_brain_obs_metadata_neuronal.rds"
) %>%
  as_tibble() %>%
  mutate(obs_dataset = "neuronal")
obs_nonneuronal <- readRDS(
  "data_intermediate/linnarsson_adult_human_brain_obs_metadata_nonneuronal.rds"
) %>%
  as_tibble() %>%
  mutate(obs_dataset = "nonneuronal")

obs <- bind_rows(obs_neuronal, obs_nonneuronal) %>%
  mutate(roi = str_squish(as.character(roi)))

roi_map <- read_s1b_roi_map(map_path)

if (anyDuplicated(roi_map$roi)) {
  stop("The authoritative s1b map assigns an ROI more than once.", call. = FALSE)
}
if (anyNA(roi_map$rcmr_value)) {
  stop("Every s1b mapping row must resolve to one Heiss rate.", call. = FALSE)
}

scope_check <- obs %>%
  distinct(roi, ROIGroupCoarse) %>%
  inner_join(
    as_tibble(roi_map) %>%
      select(
        roi, anatomy_id, is_cerebral_cortex, is_neocortex,
        is_telencephalon
      ),
    by = "roi"
  ) %>%
  mutate(
    expected_is_neocortex = ROIGroupCoarse == "Cerebral cortex",
    expected_is_cerebral_cortex =
      expected_is_neocortex | anatomy_id %in% c("amygdala", "hippocampus"),
    expected_is_telencephalon = ROIGroupCoarse %in%
      c("Cerebral cortex", "Cerebral nuclei", "Hippocampus")
  )

scope_mismatches <- scope_check %>%
  filter(
    is_neocortex != expected_is_neocortex |
      is_cerebral_cortex != expected_is_cerebral_cortex |
      is_telencephalon != expected_is_telencephalon
  )
if (nrow(scope_mismatches)) {
  stop(
    "The s1b cortical-scope flags disagree with the source anatomical ",
    "annotations for ROI(s): ",
    paste(unique(scope_mismatches$roi), collapse = ", "),
    call. = FALSE
  )
}

acc_check <- roi_map %>%
  filter(roi %in% c("Human A25", "Human A32"))
if (nrow(acc_check) != 2L || any(!acc_check$is_neocortex)) {
  stop(
    "Human A25 and Human A32 must both be retained in the s1b neocortex scope.",
    call. = FALSE
  )
}

orphan_map_rows <- setdiff(roi_map$roi, unique(obs$roi))
if (length(orphan_map_rows)) {
  stop(
    "Authoritative s1b ROI(s) absent from the extracted observations: ",
    paste(orphan_map_rows, collapse = ", "),
    call. = FALSE
  )
}

coverage <- obs %>%
  count(obs_dataset, roi, name = "n_obs_rows") %>%
  left_join(
    as_tibble(roi_map) %>%
      select(
        roi, anatomy_id, anatomy_group,
        is_cerebral_cortex, is_neocortex, is_telencephalon
      ),
    by = "roi"
  ) %>%
  mutate(mapping_status = if_else(
    is.na(anatomy_id), "unmapped", "mapped"
  )) %>%
  arrange(obs_dataset, mapping_status, roi)

write.csv(coverage, coverage_path, row.names = FALSE, na = "")

coverage_summary <- coverage %>%
  group_by(obs_dataset, mapping_status) %>%
  summarise(
    n_rois = n(),
    n_obs_rows = sum(n_obs_rows),
    .groups = "drop"
  )

print(coverage_summary, n = Inf)
message(
  "Validated ", nrow(roi_map), " authoritative s1b ROI mappings; wrote ",
  coverage_path
)
