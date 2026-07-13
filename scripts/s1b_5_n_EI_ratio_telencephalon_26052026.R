setwd(local({ d <- normalizePath(getwd()); while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d); d }))  # repo root (portable; replaces hardcoded path -- see R/project_root.R)

source("R/plot_settings.R")

## Install and Load up packages
library(ggplot2)
library(ggpmisc)
library(readxl)
library(tidyverse)

# Make sure output directories exist
if (!dir.exists("data")) dir.create("data", recursive = TRUE)
if (!dir.exists("figs")) dir.create("figs", recursive = TRUE)

############################
## Load saved obs metadata
############################
obs <- readRDS("data_intermediate/linnarsson_adult_human_brain_obs_metadata_neuronal.rds")

# filter to neurons only
obs <- obs %>%
  filter(cell_type == "neuron") %>%
  mutate(
    roi = stringr::str_squish(roi),
    supercluster_term = stringr::str_squish(as.character(supercluster_term))
  )

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
# Anatomical grouping of rois
#####################################
anatomy_rules <- readr::read_csv(
  "data_intermediate/rcmr_roi_relationship.csv",
  show_col_types = FALSE
) %>%
  dplyr::transmute(
    anatomy_group = rcmr_term,
    roi = rois
  ) %>%
  tidyr::separate_rows(roi, sep = "\\s*\\|\\|\\s*") %>%
  dplyr::mutate(roi = stringr::str_squish(roi)) %>%
  dplyr::distinct(roi, anatomy_group)

obs <- obs %>%
  dplyr::select(-dplyr::any_of("anatomy_group")) %>%
  dplyr::left_join(anatomy_rules, by = "roi") %>%
  dplyr::mutate(anatomy_group = dplyr::coalesce(anatomy_group, "Unmapped"))

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
    frac    = n_cells / n_total
  ) %>%
  slice_max(n_cells, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    anatomy_group,
    dominant_ROIGroupCoarse = ROIGroupCoarse,
    n_cells_dominant        = n_cells,
    n_cells_total           = n_total,
    dominant_fraction       = round(frac, 3),
    is_telencephalon        = dominant_ROIGroupCoarse %in% telencephalon_coarse_groups,
    division                = ifelse(is_telencephalon, "Telencephalon", "Non-telencephalon")
  ) %>%
  arrange(desc(is_telencephalon), anatomy_group)

cat("\n================ Telencephalon classification ================\n")
print(telencephalon_table, n = Inf)

write.csv(
  telencephalon_table,
  "data_analysis/telencephalon_classification_neuronal.csv",
  row.names = FALSE
)

###################################################################
# Cell type proportion calculations for E:I ratio
###################################################################
# E:I ratio is calculated as:
#   p_cells_Excitatory_projection / p_cells_Inhibitory_interneuron
# This matches the previous E and I proportion endpoints. MSN-type
# inhibitory projection neurons are retained in the denominator diagnostics
# below but are NOT included in the primary ratio, because the earlier plots
# treated the pan-region inhibitory endpoint as inhibitory interneuron.

obs <- obs %>%
  mutate(
    cell_category = case_when(
      supercluster_term %in% c(
        "Amygdala excitatory",
        "Deep-layer corticothalamic and 6b",
        "Deep-layer intratelencephalic",
        "Deep-layer near-projecting",
        "Upper-layer intratelencephalic",
        "Thalamic excitatory",
        "Hippocampal CA1-3",
        "Hippocampal CA4",
        "Hippocampal dentate gyrus",
        "Mammillary body",
        "Lower rhombic lip",
        "Upper rhombic lip"
      ) ~ "Excitatory_projection",
      supercluster_term %in% c(
        "CGE interneuron",
        "MGE interneuron",
        "LAMP5-LHX6 and Chandelier",
        "Cerebellar inhibitory",
        "Midbrain-derived inhibitory"
      ) ~ "Inhibitory_interneuron",
      supercluster_term %in% c(
        "Medium spiny neuron",
        "Eccentric medium spiny neuron"
      ) ~ "Inhibitory_projection_MSN",
      supercluster_term == "Splatter" ~ "Splatter",
      supercluster_term == "Miscellaneous" ~ "Miscellaneous",
      TRUE ~ "Other_neuron"
    )
  )

cat("\n================ Cell category counts ================\n")
print(table(obs$cell_category))

obs_celltype <- obs %>%
  filter(anatomy_group != "", anatomy_group != "Unmapped")

celltypes <- obs_celltype %>% distinct(cell_category)
donor_region <- obs_celltype %>% distinct(anatomy_group, donor_id)

celltype_counts <- obs_celltype %>%
  count(anatomy_group, donor_id, cell_category, name = "n_cells")

celltype_table_long <- donor_region %>%
  crossing(celltypes) %>%
  left_join(
    celltype_counts,
    by = c("anatomy_group", "donor_id", "cell_category")
  ) %>%
  mutate(n_cells = replace_na(n_cells, 0)) %>%
  group_by(anatomy_group, donor_id) %>%
  mutate(
    donor_total_cells = sum(n_cells),
    donor_prop = n_cells / donor_total_cells
  ) %>%
  ungroup() %>%
  group_by(anatomy_group, cell_category) %>%
  summarise(
    p_cells = mean(donor_prop, na.rm = TRUE),
    n_donors = n_distinct(donor_id),
    .groups = "drop"
  )

celltype_proportion_table <- celltype_table_long %>%
  select(anatomy_group, cell_category, p_cells) %>%
  pivot_wider(
    names_from  = cell_category,
    values_from = p_cells,
    values_fill = 0,
    names_prefix = "p_cells_"
  )

analysis_df <- celltype_proportion_table %>%
  inner_join(rcmr, by = c("anatomy_group" = "rcmr_term")) %>%
  left_join(
    telencephalon_table %>% select(anatomy_group, is_telencephalon, division),
    by = "anatomy_group"
  ) %>%
  filter(!is.na(division)) %>%
  mutate(
    p_cells_Excitatory_projection = coalesce(p_cells_Excitatory_projection, 0),
    p_cells_Inhibitory_interneuron = coalesce(p_cells_Inhibitory_interneuron, 0),
    p_cells_Inhibitory_projection_MSN = coalesce(p_cells_Inhibitory_projection_MSN, 0),
    p_cells_Inhibitory_total_with_MSN = p_cells_Inhibitory_interneuron + p_cells_Inhibitory_projection_MSN,
    EI_ratio = if_else(
      !is.na(p_cells_Inhibitory_interneuron) & p_cells_Inhibitory_interneuron > 0,
      p_cells_Excitatory_projection / p_cells_Inhibitory_interneuron,
      NA_real_
    ),
    log2_EI_ratio = log2(EI_ratio)
  )

cat("\n================ Regions in analysis by division ================\n")
print(analysis_df %>% count(division, name = "n_regions") %>% arrange(desc(n_regions)))
cat("\nRegion list by division with E:I ratio:\n")
print(
  analysis_df %>%
    select(division, anatomy_group, rcmr_value,
           p_cells_Excitatory_projection, p_cells_Inhibitory_interneuron,
           EI_ratio, log2_EI_ratio) %>%
    arrange(division, anatomy_group),
  n = Inf
)

write.csv(
  analysis_df,
  "data_analysis/neuronal_EI_ratio_by_region_with_rcmr_telencephalon_split.csv",
  row.names = FALSE
)

###########
# Analyses
###########
predictor <- "EI_ratio"

safe_spearman <- function(df, label, variable = "EI_ratio") {
  x <- df$rcmr_value
  y <- df[[variable]]
  ok <- complete.cases(x, y)
  n <- sum(ok)

  if (n < 4 || dplyr::n_distinct(x[ok]) < 2 || dplyr::n_distinct(y[ok]) < 2) {
    return(tibble(
      subset = label,
      variable = variable,
      rho = NA_real_,
      p_value = NA_real_,
      p_adj_BH = NA_real_,
      n_regions = n,
      n_regions_with_finite_ratio = sum(is.finite(y[ok])),
      mean_EI_ratio = mean(y[ok], na.rm = TRUE),
      sd_EI_ratio = sd(y[ok], na.rm = TRUE),
      status = "not_tested_constant_or_low_n"
    ))
  }

  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  tibble(
    subset = label,
    variable = variable,
    rho = unname(ct$estimate),
    p_value = ct$p.value,
    p_adj_BH = NA_real_,
    n_regions = n,
    n_regions_with_finite_ratio = sum(is.finite(y[ok])),
    mean_EI_ratio = mean(y[ok], na.rm = TRUE),
    sd_EI_ratio = sd(y[ok], na.rm = TRUE),
    status = "tested"
  )
}

run_lm_subset <- function(df, label, variable = "EI_ratio") {
  df2 <- df %>% filter(is.finite(.data[[variable]]), !is.na(rcmr_value))
  cat("\n========== LINEAR MODEL:", label, "(n =", nrow(df2), "regions) ==========" , "\n")
  if (nrow(df2) < 4 || dplyr::n_distinct(df2[[variable]]) < 2) {
    cat("  Skipping LM: insufficient finite observations or constant predictor.\n")
    return(invisible(NULL))
  }
  fit <- lm(reformulate(variable, response = "rcmr_value"), data = df2)
  print(summary(fit))
  invisible(fit)
}

make_plot_subset <- function(df, label, file_slug, variable = "EI_ratio") {
  plot_df <- df %>%
    filter(is.finite(.data[[variable]]), !is.na(rcmr_value)) %>%
    mutate(EI_ratio_label = "E:I ratio")

  if (nrow(plot_df) < 3 || dplyr::n_distinct(plot_df[[variable]]) < 2) {
    cat("  Skipping plot for ", label, ": n=", nrow(plot_df), " finite regions.\n", sep = "")
    return(invisible(NULL))
  }

  region_scale <- tryCatch({
    check_region_palette(plot_df, region_col = "anatomy_group")
    present <- intersect(region_order, unique(as.character(plot_df$anatomy_group)))
    plot_df$anatomy_group <- factor(plot_df$anatomy_group, levels = present)
    ggplot2::scale_color_manual(values = region_palette, drop = TRUE)
  }, error = function(e) {
    message("Region palette unavailable for ", label, " (", e$message,
            "); using default ggplot palette.")
    NULL
  })

  p <- ggplot(plot_df, aes(x = .data[[variable]], y = rcmr_value, color = anatomy_group)) +
    geom_point(size = 2.8, alpha = 0.85) +
    geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "steelblue") +
    stat_poly_eq(
      aes(label = paste(after_stat(rr.label), after_stat(p.value.label), sep = "*\", \"*")),
      formula = y ~ x,
      parse = TRUE,
      label.x = "right",
      label.y = "top",
      size = 4,
      color = "black"
    ) +
    labs(
      title = paste0("Neuronal E:I ratio vs rCMRGlc — ", label),
      subtitle = paste0("n = ", nrow(plot_df), " regions; E:I = excitatory projection proportion / inhibitory interneuron proportion"),
      x = "Neuronal E:I ratio",
      y = "rCMRGlc (µmol/100 g/min.)",
      color = "Region"
    ) +
    theme_classic(base_size = 14)
  if (!is.null(region_scale)) p <- p + region_scale

  print(p)

  ggsave(
    filename = paste0("figs/s1b/p_neuronal_EI_ratio_", file_slug, ".pdf"),
    plot = p, width = 8, height = 7, units = "in"
  )
  ggsave(
    filename = paste0("figs/s1b/p_neuronal_EI_ratio_", file_slug, ".jpg"),
    plot = p, width = 8, height = 7, units = "in", dpi = 300
  )
  invisible(p)
}

###################################################################
# Run analyses for each subset
###################################################################
tel_df    <- analysis_df %>% filter(is_telencephalon)
nontel_df <- analysis_df %>% filter(!is_telencephalon)

fit_tel    <- run_lm_subset(tel_df,    "Telencephalon")
fit_nontel <- run_lm_subset(nontel_df, "Non-telencephalon")
fit_all    <- run_lm_subset(analysis_df, "All mapped regions")

cor_combined <- bind_rows(
  safe_spearman(tel_df, "Telencephalon"),
  safe_spearman(nontel_df, "Non-telencephalon"),
  safe_spearman(analysis_df, "All mapped regions")
) %>%
  mutate(p_adj_BH = ifelse(status == "tested", p.adjust(p_value, method = "BH"), NA_real_))

cat("\n========== SPEARMAN CORRELATIONS: rCMRGlc vs neuronal E:I ratio ==========" , "\n")
print(cor_combined)

write.csv(
  cor_combined,
  "data_analysis/spearman_correlations_rcmr_neuronal_EI_ratio_telencephalon_split.csv",
  row.names = FALSE
)

p_tel    <- make_plot_subset(tel_df,    "Telencephalon",     "telencephalon")
p_nontel <- make_plot_subset(nontel_df, "Non-telencephalon", "nontelencephalon")
p_all    <- make_plot_subset(analysis_df, "All mapped regions", "all_regions")

###################################################################
# Combined plots: telencephalon vs non-telencephalon
###################################################################
plot_df_all <- analysis_df %>%
  filter(is.finite(EI_ratio), !is.na(rcmr_value), !is.na(division))

if (nrow(plot_df_all) < 4 || dplyr::n_distinct(plot_df_all$EI_ratio) < 2) {
  cat("\nSkipping combined plots: insufficient finite E:I ratio values.\n")
} else {

  region_scale_all <- tryCatch({
    check_region_palette(plot_df_all, region_col = "anatomy_group")
    present_all <- intersect(region_order, unique(as.character(plot_df_all$anatomy_group)))
    plot_df_all$anatomy_group <- factor(plot_df_all$anatomy_group, levels = present_all)
    ggplot2::scale_color_manual(values = region_palette, drop = TRUE)
  }, error = function(e) {
    message("Region palette unavailable for combined plot (", e$message,
            "); using default palette.")
    NULL
  })

  # (a) Grid plot: one row per division, each region colored separately
  p_combined <- ggplot(plot_df_all, aes(x = EI_ratio, y = rcmr_value, color = anatomy_group)) +
    geom_point(size = 2.8, alpha = 0.85) +
    geom_smooth(aes(group = division), method = "lm", se = TRUE, color = "steelblue") +
    stat_poly_eq(
      aes(label = paste(after_stat(rr.label), after_stat(p.value.label), sep = "*\", \"*")),
      formula = y ~ x,
      parse = TRUE,
      label.x = "right",
      label.y = "top",
      size = 3.5,
      color = "black"
    ) +
    facet_grid(division ~ ., scales = "free_x") +
    labs(
      title = "Neuronal E:I ratio vs rCMRGlc by telencephalon vs non-telencephalon",
      subtitle = "E:I = excitatory projection proportion / inhibitory interneuron proportion",
      x = "Neuronal E:I ratio",
      y = "rCMRGlc (µmol/100 g/min.)",
      color = "Region"
    ) +
    theme_classic(base_size = 12)
  if (!is.null(region_scale_all)) p_combined <- p_combined + region_scale_all

  print(p_combined)

  ggsave(
    filename = "figs/s1b/p_neuronal_EI_ratio_by_telencephalon.pdf",
    plot = p_combined, width = 9, height = 9, units = "in"
  )
  ggsave(
    filename = "figs/s1b/p_neuronal_EI_ratio_by_telencephalon.jpg",
    plot = p_combined, width = 9, height = 9, units = "in", dpi = 300
  )

  # (b) Overlay plot: both divisions on same axes, separate regression lines
  p_overlay <- ggplot(plot_df_all,
                      aes(x = EI_ratio, y = rcmr_value, color = division, fill = division)) +
    geom_point(size = 2.6, alpha = 0.85) +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.15) +
    stat_poly_eq(
      aes(label = paste(after_stat(rr.label), after_stat(p.value.label), sep = "*\", \"*")),
      formula = y ~ x,
      parse = TRUE,
      label.x = "right",
      size = 3.2,
      show.legend = FALSE
    ) +
    scale_color_manual(values = c("Telencephalon" = "#D7263D",
                                  "Non-telencephalon" = "#1B998B")) +
    scale_fill_manual(values  = c("Telencephalon" = "#D7263D",
                                  "Non-telencephalon" = "#1B998B")) +
    labs(
      title = "Neuronal E:I ratio vs rCMRGlc -- Telencephalon vs Non-telencephalon",
      subtitle = "E:I = excitatory projection proportion / inhibitory interneuron proportion",
      x = "Neuronal E:I ratio",
      y = "rCMRGlc (µmol/100 g/min.)",
      color = "Division",
      fill  = "Division"
    ) +
    theme_classic(base_size = 13) +
    theme(legend.position = "top")

  print(p_overlay)

  ggsave(
    filename = "figs/s1b/p_neuronal_EI_ratio_overlay_telencephalon.pdf",
    plot = p_overlay, width = 8, height = 7, units = "in"
  )
  ggsave(
    filename = "figs/s1b/p_neuronal_EI_ratio_overlay_telencephalon.jpg",
    plot = p_overlay, width = 8, height = 7, units = "in", dpi = 300
  )
}

cat("\nDone. Outputs:\n")
cat("  data_analysis/telencephalon_classification_neuronal.csv\n")
cat("  data_analysis/neuronal_EI_ratio_by_region_with_rcmr_telencephalon_split.csv\n")
cat("  data_analysis/spearman_correlations_rcmr_neuronal_EI_ratio_telencephalon_split.csv\n")
cat("  figs/s1b/p_neuronal_EI_ratio_telencephalon.{pdf,jpg}\n")
cat("  figs/s1b/p_neuronal_EI_ratio_nontelencephalon.{pdf,jpg}\n")
cat("  figs/s1b/p_neuronal_EI_ratio_all_regions.{pdf,jpg}\n")
cat("  figs/s1b/p_neuronal_EI_ratio_by_telencephalon.{pdf,jpg}\n")
cat("  figs/s1b/p_neuronal_EI_ratio_overlay_telencephalon.{pdf,jpg}\n")
