setwd("~/Library/CloudStorage/Dropbox/COLLABORATIVE/Do expensive brain regions increase less in humans/analyses_metabol_rate_structure")

source("R/0.01_plot_settings.R")

library(readxl)
library(tidyverse)

############################
## Load saved obs metadata
############################
obs <- readRDS("data/linnarsson_adult_human_brain_obs_metadata_neuronal.rds")

# neurons only
obs <- obs %>%
  filter(cell_type == "neuron") %>%
  mutate(
    roi = stringr::str_squish(roi),
    supercluster_term = stringr::str_squish(as.character(supercluster_term))
  ) %>%
  filter(!is.na(supercluster_term), supercluster_term != "")

######################################################
# Read rCMRGlc values from Heiss et al. 2004
######################################################
heiss_stephan_tbl <- read.csv("data/Heiss_Stephan_data.csv")

rcmr <- heiss_stephan_tbl %>%
  transmute(
    rcmr_term = str_trim(as.character(term_3)),
    rcmr_value = round(as.numeric(rCMRGlc_mean_both_hemispheres), 1)
  ) %>%
  filter(
    !is.na(rcmr_term),
    rcmr_term != "",
    !str_detect(str_to_lower(rcmr_term), "average"),
    !is.na(rcmr_value)
  )

#####################################
# Anatomical grouping of ROIs
#####################################
anatomy_rules <- readr::read_csv(
  "data/rcmr_roi_relationship.csv",
  show_col_types = FALSE
) %>%
  transmute(
    anatomy_group = rcmr_term,
    roi = rois
  ) %>%
  separate_rows(roi, sep = "\\s*\\|\\|\\s*") %>%
  mutate(roi = stringr::str_squish(roi)) %>%
  distinct(roi, anatomy_group)

obs <- obs %>%
  select(-any_of("anatomy_group")) %>%
  left_join(anatomy_rules, by = "roi") %>%
  mutate(anatomy_group = coalesce(anatomy_group, "Unmapped"))

###################################################################
# Telencephalon classification
###################################################################
telencephalon_coarse_groups <- c("Cerebral cortex", "Cerebral nuclei", "Hippocampus")

telencephalon_table <- obs %>%
  filter(anatomy_group != "", anatomy_group != "Unmapped") %>%
  count(anatomy_group, ROIGroupCoarse, name = "n_cells") %>%
  group_by(anatomy_group) %>%
  mutate(
    n_total = sum(n_cells),
    frac = n_cells / n_total
  ) %>%
  slice_max(n_cells, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    anatomy_group,
    dominant_ROIGroupCoarse = ROIGroupCoarse,
    n_cells_dominant = n_cells,
    n_cells_total = n_total,
    dominant_fraction = round(frac, 3),
    is_telencephalon = dominant_ROIGroupCoarse %in% telencephalon_coarse_groups,
    division = ifelse(is_telencephalon, "Telencephalon", "Non-telencephalon")
  ) %>%
  arrange(desc(is_telencephalon), anatomy_group)

write.csv(
  telencephalon_table,
  "data/telencephalon_classification_neuronal.csv",
  row.names = FALSE
)

###################################################################
# Mean participant-level proportions by anatomy_group over ALL
# neuronal supercluster_term labels
###################################################################
obs_supercluster <- obs %>%
  filter(anatomy_group != "", anatomy_group != "Unmapped")

superclusters <- obs_supercluster %>%
  distinct(supercluster_term)

donor_region <- obs_supercluster %>%
  distinct(anatomy_group, donor_id)

supercluster_counts <- obs_supercluster %>%
  count(anatomy_group, donor_id, supercluster_term, name = "n_cells")

supercluster_table_long <- donor_region %>%
  crossing(superclusters) %>%
  left_join(
    supercluster_counts,
    by = c("anatomy_group", "donor_id", "supercluster_term")
  ) %>%
  mutate(n_cells = replace_na(n_cells, 0)) %>%
  group_by(anatomy_group, donor_id) %>%
  mutate(
    donor_total_cells = sum(n_cells),
    donor_prop = n_cells / donor_total_cells
  ) %>%
  ungroup() %>%
  group_by(anatomy_group, supercluster_term) %>%
  summarise(
    p_cells = mean(donor_prop, na.rm = TRUE),
    n_donors = n_distinct(donor_id),
    .groups = "drop"
  )

# Wide table with one p_cells_* column per original neuronal supercluster_term.
supercluster_proportion_table <- supercluster_table_long %>%
  select(anatomy_group, supercluster_term, p_cells) %>%
  pivot_wider(
    names_from = supercluster_term,
    values_from = p_cells,
    values_fill = 0,
    names_prefix = "p_cells_",
    names_repair = "unique"
  )

analysis_df <- supercluster_proportion_table %>%
  inner_join(rcmr, by = c("anatomy_group" = "rcmr_term")) %>%
  left_join(
    telencephalon_table %>% select(anatomy_group, is_telencephalon, division),
    by = "anatomy_group"
  ) %>%
  filter(!is.na(division))

cat("\n================ Regions in analysis by division ================\n")
print(analysis_df %>% count(division, name = "n_regions") %>% arrange(desc(n_regions)))
cat("\nRegion list by division:\n")
print(
  analysis_df %>% select(division, anatomy_group, rcmr_value) %>% arrange(division, anatomy_group),
  n = Inf
)

write.csv(
  analysis_df,
  "data/supercluster_proportions_by_region_with_rcmr_telencephalon_split.csv",
  row.names = FALSE
)

###################################################################
# Correlations: rCMRGlc against every neuronal supercluster_term,
# separately in telencephalon and non-telencephalon.
###################################################################
all_supercluster_predictors <- names(analysis_df)[grepl("^p_cells_", names(analysis_df))]

safe_cor_one <- function(df, variable, label, method = "spearman") {
  x <- df$rcmr_value
  y <- df[[variable]]
  ok <- complete.cases(x, y)
  n <- sum(ok)
  supercluster <- sub("^p_cells_", "", variable)

  if (n < 4 || dplyr::n_distinct(x[ok]) < 2 || dplyr::n_distinct(y[ok]) < 2) {
    return(tibble(
      division = label,
      supercluster_term = supercluster,
      variable = variable,
      rho = NA_real_,
      p_value = NA_real_,
      n_regions = n,
      n_nonzero_regions = sum(y[ok] > 0),
      mean_region_prop = mean(y[ok], na.rm = TRUE),
      sd_region_prop = sd(y[ok], na.rm = TRUE),
      status = "not_tested_constant_or_low_n"
    ))
  }

  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = method, exact = FALSE))
  tibble(
    division = label,
    supercluster_term = supercluster,
    variable = variable,
    rho = unname(ct$estimate),
    p_value = ct$p.value,
    n_regions = n,
    n_nonzero_regions = sum(y[ok] > 0),
    mean_region_prop = mean(y[ok], na.rm = TRUE),
    sd_region_prop = sd(y[ok], na.rm = TRUE),
    status = "tested"
  )
}

run_supercluster_correlations <- function(df, label, predictors) {
  bind_rows(lapply(predictors, function(v) safe_cor_one(df, v, label))) %>%
    mutate(
      p_adj_BH = ifelse(status == "tested", p.adjust(p_value, method = "BH"), NA_real_)
    ) %>%
    arrange(p_value, supercluster_term)
}

cor_tel <- run_supercluster_correlations(
  analysis_df %>% filter(is_telencephalon),
  "Telencephalon",
  all_supercluster_predictors
)

cor_nontel <- run_supercluster_correlations(
  analysis_df %>% filter(!is_telencephalon),
  "Non-telencephalon",
  all_supercluster_predictors
)

cor_long <- bind_rows(cor_tel, cor_nontel) %>%
  arrange(supercluster_term, division)

# Primary matrix requested: rho values, one row per supercluster_term,
# one column per brain-region division.
rho_matrix <- cor_long %>%
  select(supercluster_term, division, rho) %>%
  pivot_wider(names_from = division, values_from = rho) %>%
  arrange(supercluster_term)

# Companion matrices/tables for interpretation.
p_matrix <- cor_long %>%
  select(supercluster_term, division, p_value) %>%
  pivot_wider(names_from = division, values_from = p_value) %>%
  arrange(supercluster_term)

padj_matrix <- cor_long %>%
  select(supercluster_term, division, p_adj_BH) %>%
  pivot_wider(names_from = division, values_from = p_adj_BH) %>%
  arrange(supercluster_term)

n_regions_matrix <- cor_long %>%
  select(supercluster_term, division, n_regions) %>%
  pivot_wider(names_from = division, values_from = n_regions) %>%
  arrange(supercluster_term)

n_nonzero_regions_matrix <- cor_long %>%
  select(supercluster_term, division, n_nonzero_regions) %>%
  pivot_wider(names_from = division, values_from = n_nonzero_regions) %>%
  arrange(supercluster_term)

# One combined matrix with rho, nominal p, BH-adjusted p, and nonzero-region counts.
cor_matrix_full <- cor_long %>%
  select(
    supercluster_term, division, rho, p_value, p_adj_BH,
    n_regions, n_nonzero_regions, mean_region_prop, sd_region_prop, status
  ) %>%
  pivot_wider(
    names_from = division,
    values_from = c(rho, p_value, p_adj_BH, n_regions, n_nonzero_regions,
                    mean_region_prop, sd_region_prop, status),
    names_glue = "{division}_{.value}"
  ) %>%
  arrange(supercluster_term)

write.csv(cor_long, "data/spearman_correlations_rcmr_by_supercluster_term_long.csv", row.names = FALSE)
write.csv(rho_matrix, "data/spearman_rho_matrix_rcmr_by_supercluster_term_telencephalon_split.csv", row.names = FALSE)
write.csv(p_matrix, "data/spearman_pvalue_matrix_rcmr_by_supercluster_term_telencephalon_split.csv", row.names = FALSE)
write.csv(padj_matrix, "data/spearman_BH_padj_matrix_rcmr_by_supercluster_term_telencephalon_split.csv", row.names = FALSE)
write.csv(n_regions_matrix, "data/n_regions_matrix_rcmr_by_supercluster_term_telencephalon_split.csv", row.names = FALSE)
write.csv(n_nonzero_regions_matrix, "data/n_nonzero_regions_matrix_rcmr_by_supercluster_term_telencephalon_split.csv", row.names = FALSE)
write.csv(cor_matrix_full, "data/spearman_correlation_full_matrix_rcmr_by_supercluster_term_telencephalon_split.csv", row.names = FALSE)

cat("\n================ Spearman rho matrix: rCMRGlc ~ supercluster proportions ================\n")
print(rho_matrix, n = Inf)

cat("\n================ Full correlation matrix ================\n")
print(cor_matrix_full, n = Inf)

cat("\nDone. Outputs:\n")
cat("  data/supercluster_proportions_by_region_with_rcmr_telencephalon_split.csv\n")
cat("  data/spearman_correlations_rcmr_by_supercluster_term_long.csv\n")
cat("  data/spearman_rho_matrix_rcmr_by_supercluster_term_telencephalon_split.csv\n")
cat("  data/spearman_pvalue_matrix_rcmr_by_supercluster_term_telencephalon_split.csv\n")
cat("  data/spearman_BH_padj_matrix_rcmr_by_supercluster_term_telencephalon_split.csv\n")
cat("  data/n_regions_matrix_rcmr_by_supercluster_term_telencephalon_split.csv\n")
cat("  data/n_nonzero_regions_matrix_rcmr_by_supercluster_term_telencephalon_split.csv\n")
cat("  data/spearman_correlation_full_matrix_rcmr_by_supercluster_term_telencephalon_split.csv\n")
