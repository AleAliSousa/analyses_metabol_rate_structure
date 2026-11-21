############################################################
# mouse_s1c_0 -- Turn the SV2A Mouse Synaptome Atlas download (Grant lab,
# University of Edinburgh) into the regional table read by
# mouse_s1c_synaptic_density_vs_rcmr.R.
#
# Input : data_raw/SV2A_Mouse_Synaptome_Atlas_Data.xlsx
#         one row per mouse x atlas region (33 mice, both sexes, ages 6W/6M/18M,
#         140 regions). Columns: mouseID, sex, age, region_renamed,
#         count_for_stats, mean, count_sum, main_region.
# Map   : metadata/mouse_sv2a_region_to_cmrgl_region_map.csv
# Output: data_raw/mouse_synaptic_density_by_region.csv   (analysis sets)
#         tables/mouse_s1c/S1_layer_profile.csv
#
# Aggregation: layers/lobules -> area (equal weight) -> CMRgl region (equal
# area weight), within each mouse; then mean across mice of the age group.
# Measures are taken as provided; their meaning is NOT yet confirmed:
#   count_for_stats -> treated as the density-like measure (primary)
#   mean            -> second measure (probably puncta intensity)
# count_sum (total puncta) is not used: it scales with the sampled area.
############################################################
.a <- commandArgs(FALSE); .f <- sub("^--file=", "", .a[grepl("^--file=", .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else "analyses_metabol_rate_structure_nonhuman/scripts", "common.R"))
source(nh("scripts/mouse_s1b_helpers.R"))
suppressPackageStartupMessages(library(readxl))
dir.create(nh("tables/mouse_s1c"), recursive = TRUE, showWarnings = FALSE)

x <- read_excel(nh("data_raw/SV2A_Mouse_Synaptome_Atlas_Data.xlsx"))
req <- c("mouseID", "sex", "age", "region_renamed", "count_for_stats", "mean")
if (length(m <- setdiff(req, names(x)))) stop("Atlas file lacks column(s): ", paste(m, collapse = ", "), call. = FALSE)
if (anyDuplicated(x[c("mouseID", "region_renamed")])) stop("Duplicate mouse x region rows.", call. = FALSE)
if (!all(x$age %in% c("6W", "6M", "18M"))) stop("Unexpected age labels.", call. = FALSE)
if (anyNA(x$count_for_stats) || anyNA(x$mean)) stop("Missing measurements in atlas file.", call. = FALSE)

map <- read_csv(nh("metadata/mouse_sv2a_region_to_cmrgl_region_map.csv"), show_col_types = FALSE)
if (length(m <- setdiff(map$sv2a_region, x$region_renamed))) stop("Mapped region(s) absent from atlas: ", paste(m, collapse = ", "), call. = FALSE)
rates <- read_mouse_rates()
if (length(m <- setdiff(map$cmrgl_region, rates$cmrgl_region))) stop("Mapped region(s) without a rate: ", paste(m, collapse = ", "), call. = FALSE)

long <- x %>% pivot_longer(c(count_for_stats, mean), names_to = "measure", values_to = "value") %>%
  inner_join(map, by = c("region_renamed" = "sv2a_region"), relationship = "many-to-many")

per_mouse <- long %>%
  group_by(mouseID, sex, age, measure, cmrgl_region, primary_included, nested_in_other_region, area) %>%
  summarise(value = mean(value), .groups = "drop") %>%
  group_by(mouseID, sex, age, measure, cmrgl_region, primary_included, nested_in_other_region) %>%
  summarise(value = mean(value), n_areas = n(), .groups = "drop")

age_sets <- c("6M", "6W", "18M")
by_age <- per_mouse %>% mutate(analysis_set = paste0("age_", age)) 
all_ages <- per_mouse %>% mutate(analysis_set = "all_ages")
sets <- bind_rows(by_age, all_ages) %>%
  group_by(analysis_set, measure, cmrgl_region, primary_included, nested_in_other_region) %>%
  summarise(synaptic_density = mean(value), synaptic_density_sd = sd(value), n_animals = n_distinct(mouseID), .groups = "drop") %>%
  mutate(measure = recode(measure, count_for_stats = "SV2A_count_for_stats", mean = "SV2A_mean"),
         analysis_set = paste(analysis_set, measure, sep = "|"),
         unit = "as provided (unit unconfirmed)",
         source = "SV2A Mouse Synaptome Atlas (Grant lab, Univ. Edinburgh)",
         comparability_note = "Mouse age group pooled over sex; layers/lobules and areas equal-weighted")
write_csv(sets, nh("data_raw/mouse_synaptic_density_by_region.csv"))
message("Wrote ", nrow(sets), " rows in ", n_distinct(sets$analysis_set), " analysis sets")

# S1 layer profile: SS layers vs the other isocortical layers (per age, count_for_stats)
iso_layers <- x %>% filter(grepl("^(SS|MO|VIS|PTLp|RSPd|RSPv|FRP|ORB)", region_renamed), !grepl("^MOB", region_renamed)) %>%
  mutate(area = sub("(1|23|4|5|6a|6b)$", "", region_renamed), layer = sub("^.*?(1|23|4|5|6a|6b)$", "\\1", region_renamed, perl = TRUE),
         group = if_else(area == "SS", "S1 (SS)", "other isocortex"))
prof <- iso_layers %>% group_by(age, group, layer) %>%
  summarise(count_for_stats = mean(count_for_stats), mean = mean(mean), n_mice = n_distinct(mouseID), .groups = "drop")
write_csv(prof, nh("tables/mouse_s1c/S1_layer_profile.csv"))
