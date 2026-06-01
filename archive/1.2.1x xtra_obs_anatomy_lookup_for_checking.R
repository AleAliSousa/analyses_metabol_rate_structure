setwd("~/Library/CloudStorage/Dropbox/COLLABORATIVE/Do expensive brain regions increase less in humans/analyses_metabol_rate_structure")

############################################################
# Create CSVs to inspect obs terms that map to rCMRGlc terms
############################################################

## Install and Load up packages
library(ggplot2)
library(ggpmisc)
library(readxl)
library(tidyverse)

data_dir <- "data"
out_dir <- "output/anatomy_rule_audit"

## cortical ROIs / Lobes are here: https://github.com/linnarsson-lab/adult-human-brain/blob/main/notebooks/Revision/RevisionFig2.ipynb
## Also See Fig 2 in the original paper: https://www.science.org/doi/10.1126/science.add7046

############################
## Load saved obs metadata
############################
obs1 <- readRDS("data/linnarsson_adult_human_brain_obs_metadata_nonneuronal.rds") %>%
  as_tibble() %>%
  mutate(obs_dataset = "nonneuronal")

obs2 <- readRDS("data/linnarsson_adult_human_brain_obs_metadata_neuronal.rds") %>%
  as_tibble() %>%
  mutate(obs_dataset = "neuronal")

# Check if the columns are identical
setdiff(colnames(obs1), colnames(obs2))
setdiff(colnames(obs2), colnames(obs1))

obs <- bind_rows(obs1, obs2)

cols_keep <- c(
  "dissection",
  "ROIGroup",
  "ROIGroupCoarse",
  "ROIGroupFine",
  "roi"
)

obs_anatomy_lookup <- obs %>%
  mutate(across(all_of(cols_keep), ~ str_squish(as.character(.x)))) %>%
  select(all_of(cols_keep)) %>%
  filter(if_any(all_of(cols_keep), ~ !is.na(.x) & .x != "")) %>%
  distinct() %>%
  arrange(dissection, roi, ROIGroup, ROIGroupCoarse, ROIGroupFine)

cat("Number of unique anatomy rows:", nrow(obs_anatomy_lookup), "\n")


write_csv(
  obs_anatomy_lookup,
  "output/obs_unique_anatomy_lookup.csv"
)

filtered_obs <- obs %>%
  filter(grepl("cingulate", dissection, ignore.case = TRUE))

unique <- data.frame(unique(filtered_obs$roi))
unique_dis <- data.frame(unique(filtered_obs$dissection))