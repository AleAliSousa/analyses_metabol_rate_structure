library(dplyr)
library(stringr)
library(readr)

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