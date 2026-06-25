setwd("~/Library/CloudStorage/Dropbox/COLLABORATIVE/Do expensive brain regions increase less in humans/analyses_metabol_rate_structure")

# ------------------------------------------------------------
# Re-run the correlation pipeline. Sourcing this also runs its own
# source("R/plot_settings.R"), so we inherit region_palette, canonical_region(),
# scale_color_regions(), theme_facet_compact(), facet_dims(), AND the objects
# analysis_df, all_supercluster_predictors, run_supercluster_correlations().
# Edit the path below if your correlation script lives elsewhere / is renamed.
# ------------------------------------------------------------

out_dir <- "figs/s1b"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

source("R/plot_settings.R")

library(readxl)
library(tidyverse)

############################
## Load saved obs metadata
############################
obs <- readRDS("data_intermediate/linnarsson_adult_human_brain_obs_metadata_neuronal.rds")

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
heiss_stephan_tbl <- read.csv("data_intermediate/Heiss_Stephan_data.csv")

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
  "data_intermediate/rcmr_roi_relationship.csv",
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
  "data_analysis/telencephalon_classification_neuronal.csv",
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
  "data_analysis/supercluster_proportions_by_region_with_rcmr_telencephalon_split.csv",
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

write.csv(cor_long, "data_analysis/spearman_correlations_rcmr_by_supercluster_term_long.csv", row.names = FALSE)
write.csv(rho_matrix, "data_analysis/spearman_rho_matrix_rcmr_by_supercluster_term_telencephalon_split.csv", row.names = FALSE)
write.csv(p_matrix, "data_analysis/spearman_pvalue_matrix_rcmr_by_supercluster_term_telencephalon_split.csv", row.names = FALSE)
write.csv(padj_matrix, "data_analysis/spearman_BH_padj_matrix_rcmr_by_supercluster_term_telencephalon_split.csv", row.names = FALSE)
write.csv(n_regions_matrix, "data_analysis/n_regions_matrix_rcmr_by_supercluster_term_telencephalon_split.csv", row.names = FALSE)
write.csv(n_nonzero_regions_matrix, "data_analysis/n_nonzero_regions_matrix_rcmr_by_supercluster_term_telencephalon_split.csv", row.names = FALSE)
write.csv(cor_matrix_full, "data_analysis/spearman_correlation_full_matrix_rcmr_by_supercluster_term_telencephalon_split.csv", row.names = FALSE)

cat("\n================ Spearman rho matrix: rCMRGlc ~ supercluster proportions ================\n")
print(rho_matrix, n = Inf)

cat("\n================ Full correlation matrix ================\n")
print(cor_matrix_full, n = Inf)

cat("\nDone. Outputs:\n")
cat("  data_analysis/supercluster_proportions_by_region_with_rcmr_telencephalon_split.csv\n")
cat("  data_analysis/spearman_correlations_rcmr_by_supercluster_term_long.csv\n")
cat("  data_analysis/spearman_rho_matrix_rcmr_by_supercluster_term_telencephalon_split.csv\n")
cat("  data_analysis/spearman_pvalue_matrix_rcmr_by_supercluster_term_telencephalon_split.csv\n")
cat("  data_analysis/spearman_BH_padj_matrix_rcmr_by_supercluster_term_telencephalon_split.csv\n")
cat("  data_analysis/n_regions_matrix_rcmr_by_supercluster_term_telencephalon_split.csv\n")
cat("  data_analysis/n_nonzero_regions_matrix_rcmr_by_supercluster_term_telencephalon_split.csv\n")
cat("  data_analysis/spearman_correlation_full_matrix_rcmr_by_supercluster_term_telencephalon_split.csv\n")


# ------------------------------------------------------------
# Long form: one row per (region, supercluster_term).
# analysis_df is the in-memory wide table, so supercluster names are intact
# (no read.csv name-mangling).
# ------------------------------------------------------------
panel_long <- analysis_df %>%
  pivot_longer(
    cols = all_of(all_supercluster_predictors),
    names_to = "variable",
    values_to = "proportion"
  ) %>%
  mutate(
    supercluster_term = sub("^p_cells_", "", variable),
    region_color = canonical_region(anatomy_group)
  )

# Surface (don't silently grey out) any region lacking a palette colour.
missing_palette <- setdiff(unique(panel_long$region_color), names(region_palette))
if (length(missing_palette)) {
  warning("Regions without a palette colour (drawn in grey): ",
          paste(missing_palette, collapse = ", "))
}

# Fixed facet order, shared across all three figures so panels line up.
supercluster_levels <- sort(unique(panel_long$supercluster_term))
panel_long$supercluster_term <- factor(panel_long$supercluster_term,
                                       levels = supercluster_levels)

# ------------------------------------------------------------
# Spearman rho + nominal p per supercluster, computed per kind.
# Reuses run_supercluster_correlations() from the sourced script so the values
# match cor_matrix_full exactly for the telencephalon / non-telencephalon kinds.
# ------------------------------------------------------------
stats_for <- function(df_wide, label) {
  run_supercluster_correlations(df_wide, label, all_supercluster_predictors) %>%
    select(supercluster_term, rho, p_value)
}

stats_tel <- stats_for(analysis_df %>% filter(is_telencephalon),  "Telencephalon")
stats_non <- stats_for(analysis_df %>% filter(!is_telencephalon), "Non-telencephalon")
stats_all <- stats_for(analysis_df,                                "All")

fmt_p <- function(p) {
  ifelse(is.na(p), "P = NA",
         ifelse(p < 0.001, "P < 0.001",
                paste0("P = ", formatC(p, digits = 3, format = "f"))))
}

make_labels <- function(stats_df) {
  stats_df %>%
    transmute(
      supercluster_term = factor(supercluster_term, levels = supercluster_levels),
      label = ifelse(is.na(rho),
                     "\u03C1 = NA",
                     paste0("\u03C1 = ", formatC(rho, digits = 2, format = "f"),
                            ", ", fmt_p(p_value)))
    )
}

# ------------------------------------------------------------
# Plot builder (matches slides 16-18).
#   - single blue lm line per panel (NOT coloured by region)
#   - points coloured by region
#   - rho + P annotation in the top-right of each panel
# ------------------------------------------------------------
line_col   <- "#377EB8"   # blue lm line (slides 16-18)
ribbon_col <- "grey70"

make_panel_fig <- function(df_long, stats_df, title) {
  labels_df <- make_labels(stats_df)
  ggplot(df_long, aes(x = proportion, y = rcmr_value)) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                color = line_col, fill = ribbon_col, linewidth = 0.8, na.rm = TRUE) +
    geom_point(aes(color = region_color), size = 2.3, na.rm = TRUE) +
    geom_text(data = labels_df, aes(x = Inf, y = Inf, label = label),
              inherit.aes = FALSE, hjust = 1.05, vjust = 1.4, size = 3) +
    facet_wrap(~ supercluster_term, scales = "free_x") +
    scale_color_regions(name = "Region", na.value = "grey60") +
    labs(
      title = title,
      x = "Cell-type proportion",
      y = "rCMRGlc (\u00B5mol/100 g/min.)"
    ) +
    theme_facet_compact()
}

# ------------------------------------------------------------
# Build and save one dense figure per kind.
# ------------------------------------------------------------
n_panels <- length(supercluster_levels)
dims <- facet_dims(n_panels, ncol = 6, panel_w = 3.0, panel_h = 2.6)

save_fig <- function(p, name) {
  ggsave(file.path(out_dir, paste0(name, ".png")), p,
         width = dims$width, height = dims$height, dpi = 300,
         limitsize = FALSE, bg = "white")
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p,
         width = dims$width, height = dims$height, limitsize = FALSE)
}

p_tel <- make_panel_fig(
  panel_long %>% filter(is_telencephalon), stats_tel,
  "Telencephalon regions: rCMRGlc vs neuronal supercluster proportion"
)
p_non <- make_panel_fig(
  panel_long %>% filter(!is_telencephalon), stats_non,
  "Non-telencephalon regions: rCMRGlc vs neuronal supercluster proportion"
)
p_all <- make_panel_fig(
  panel_long, stats_all,
  "All regions: rCMRGlc vs neuronal supercluster proportion"
)

save_fig(p_tel, "scatter_rcmr_supercluster_telencephalon")
save_fig(p_non, "scatter_rcmr_supercluster_non_telencephalon")
save_fig(p_all, "scatter_rcmr_supercluster_all_regions")

message(sprintf("Saved 3 figures (png+pdf) to %s | %d panels per figure (%d cols x %d rows)",
                out_dir, n_panels, dims$ncol, dims$nrow))