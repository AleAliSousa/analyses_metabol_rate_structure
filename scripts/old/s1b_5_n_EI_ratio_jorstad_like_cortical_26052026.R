setwd("~/Library/CloudStorage/Dropbox/COLLABORATIVE/Do expensive brain regions increase less in humans/analyses_metabol_rate_structure")

source("R/plot_settings.R")

## Install and Load up packages
library(ggplot2)
library(ggpmisc)
library(readxl)
library(tidyverse)
library(tibble)

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
    supercluster_term = stringr::str_squish(as.character(supercluster_term)),
    ROIGroup = stringr::str_squish(as.character(ROIGroup)),
    ROIGroupCoarse = stringr::str_squish(as.character(ROIGroupCoarse)),
    ROIGroupFine = stringr::str_squish(as.character(ROIGroupFine)),
    dissection = stringr::str_squish(as.character(dissection))
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
# Jorstad-like cortical E:I ratio definition
###################################################################
# This is the closest supercluster-level match to the Jorstad et al.
# neocortical E:I definition, using cortical excitatory neuron families
# and cortical interneuron families only.
#
# Siletti supercluster_term                     Jorstad-like family
# -----------------------------------------------------------------------------
# Upper-layer intratelencephalic               L2/3 IT + L4 IT
# Deep-layer intratelencephalic                L5 IT + L6 IT + L6 IT Car3
# Deep-layer corticothalamic and 6b            L6 CT + L6b
# Deep-layer near-projecting                   L5/6 NP
# CGE interneuron                              LAMP5 + SNCG + VIP + PAX6
# MGE interneuron                              PVALB + SST + SST CHODL
# LAMP5-LHX6 and Chandelier                    LAMP5 LHX6 + Chandelier
#
# Excluded from this cortical Jorstad-like E:I calculation:
# hippocampal, amygdala, thalamic, mammillary body, rhombic lip,
# cerebellar inhibitory, midbrain-derived inhibitory, medium spiny neurons,
# Splatter, and Miscellaneous.

jorstad_like_E <- c(
  "Upper-layer intratelencephalic",
  "Deep-layer intratelencephalic",
  "Deep-layer corticothalamic and 6b",
  "Deep-layer near-projecting"
)

jorstad_like_I <- c(
  "CGE interneuron",
  "MGE interneuron",
  "LAMP5-LHX6 and Chandelier"
)

siletti_jorstad_crosswalk <- tibble::tribble(
  ~supercluster_term,                         ~EI_class_jorstad_like, ~jorstad_equivalent,
  "Upper-layer intratelencephalic",           "E",                    "L2/3 IT + L4 IT",
  "Deep-layer intratelencephalic",            "E",                    "L5 IT + L6 IT + L6 IT Car3",
  "Deep-layer corticothalamic and 6b",        "E",                    "L6 CT + L6b",
  "Deep-layer near-projecting",               "E",                    "L5/6 NP",
  "CGE interneuron",                          "I",                    "LAMP5 + SNCG + VIP + PAX6",
  "MGE interneuron",                          "I",                    "PVALB + SST + SST CHODL",
  "LAMP5-LHX6 and Chandelier",                "I",                    "LAMP5 LHX6 + Chandelier"
)

write.csv(
  siletti_jorstad_crosswalk,
  "data_analysis/siletti_to_jorstad_like_cortical_EI_crosswalk.csv",
  row.names = FALSE
)

cat("\n================ Jorstad-like cortical E:I crosswalk ================\n")
print(siletti_jorstad_crosswalk, n = Inf)

###################################################################
# Identify cortical rCMR/anatomy groups
###################################################################
# The main criterion is dominant Siletti ROIGroupCoarse == "Cerebral cortex".
# A permissive text backup is included because local metadata versions can differ.

cortex_annotation_regex <- "cortex|cortical|neocortex|gyrus|precentral|postcentral|frontal|temporal|parietal|occipital|cingulate|insula"

cortex_table <- obs %>%
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
    is_cortex_by_dominant_group = dominant_ROIGroupCoarse == "Cerebral cortex"
  ) %>%
  left_join(
    obs %>%
      filter(anatomy_group != "", anatomy_group != "Unmapped") %>%
      group_by(anatomy_group) %>%
      summarise(
        any_cortical_text = any(
          str_detect(str_to_lower(paste(ROIGroup, ROIGroupCoarse, ROIGroupFine, dissection, sep = " ")),
                     cortex_annotation_regex),
          na.rm = TRUE
        ),
        .groups = "drop"
      ),
    by = "anatomy_group"
  ) %>%
  mutate(
    is_cortex = is_cortex_by_dominant_group | any_cortical_text,
    cortex_call = ifelse(is_cortex, "Cortex", "Non-cortex")
  ) %>%
  arrange(desc(is_cortex), anatomy_group)

cat("\n================ Cortex classification ================\n")
print(cortex_table, n = Inf)

write.csv(
  cortex_table,
  "data_analysis/cortex_classification_neuronal.csv",
  row.names = FALSE
)

###################################################################
# Audit which Siletti numeric clusters enter the Jorstad-like classes
###################################################################
cluster_membership_check <- obs %>%
  inner_join(siletti_jorstad_crosswalk, by = "supercluster_term") %>%
  count(EI_class_jorstad_like, jorstad_equivalent, supercluster_term,
        cluster_id, subcluster_id, sort = TRUE) %>%
  arrange(EI_class_jorstad_like, jorstad_equivalent, supercluster_term,
          cluster_id, subcluster_id)

cat("\n================ Numeric cluster/subcluster membership in Jorstad-like E:I classes ================\n")
print(tibble::as_tibble(cluster_membership_check), n = Inf)

write.csv(
  cluster_membership_check,
  "data_analysis/siletti_cluster_subcluster_membership_jorstad_like_cortical_EI.csv",
  row.names = FALSE
)

###################################################################
# Calculate cortex-only Jorstad-like E:I ratio per rCMR/anatomy group
###################################################################
obs_cortex_ei <- obs %>%
  inner_join(cortex_table %>% filter(is_cortex) %>% select(anatomy_group), by = "anatomy_group") %>%
  inner_join(siletti_jorstad_crosswalk, by = "supercluster_term")

cat("\n================ Jorstad-like cortical E/I cell counts ================\n")
x <- obs_cortex_ei %>%
  count(EI_class_jorstad_like, supercluster_term, sort = TRUE)

View(x)
# Donor-balanced proportions: calculate E and I fractions within the included
# cortical E/I neuron set for each donor and region, then average donors.
ei_classes <- tibble(EI_class_jorstad_like = c("E", "I"))
donor_region <- obs_cortex_ei %>% distinct(anatomy_group, donor_id)

ei_counts <- obs_cortex_ei %>%
  count(anatomy_group, donor_id, EI_class_jorstad_like, name = "n_cells")

ei_table_long <- donor_region %>%
  crossing(ei_classes) %>%
  left_join(ei_counts, by = c("anatomy_group", "donor_id", "EI_class_jorstad_like")) %>%
  mutate(n_cells = replace_na(n_cells, 0)) %>%
  group_by(anatomy_group, donor_id) %>%
  mutate(
    donor_total_EI_cells = sum(n_cells),
    donor_prop = if_else(donor_total_EI_cells > 0, n_cells / donor_total_EI_cells, NA_real_)
  ) %>%
  ungroup() %>%
  group_by(anatomy_group, EI_class_jorstad_like) %>%
  summarise(
    p_cells = mean(donor_prop, na.rm = TRUE),
    n_donors = n_distinct(donor_id[donor_total_EI_cells > 0]),
    n_cells_total = sum(n_cells, na.rm = TRUE),
    .groups = "drop"
  )

ei_region <- ei_table_long %>%
  select(anatomy_group, EI_class_jorstad_like, p_cells) %>%
  pivot_wider(
    names_from = EI_class_jorstad_like,
    values_from = p_cells,
    values_fill = NA_real_,
    names_prefix = "p_cells_jorstad_like_"
  ) %>%
  left_join(
    ei_table_long %>%
      group_by(anatomy_group) %>%
      summarise(
        n_donors_EI = max(n_donors, na.rm = TRUE),
        n_cells_EI_total = sum(n_cells_total, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "anatomy_group"
  ) %>%
  inner_join(rcmr, by = c("anatomy_group" = "rcmr_term")) %>%
  left_join(cortex_table, by = "anatomy_group") %>%
  mutate(
    EI_ratio_jorstad_like = if_else(
      !is.na(p_cells_jorstad_like_I) & p_cells_jorstad_like_I > 0,
      p_cells_jorstad_like_E / p_cells_jorstad_like_I,
      NA_real_
    ),
    log2_EI_ratio_jorstad_like = log2(EI_ratio_jorstad_like)
  ) %>%
  filter(is_cortex)

cat("\n================ Cortex-only Jorstad-like E:I table ================\n")
print(
  ei_region %>%
    select(anatomy_group, rcmr_value,
           p_cells_jorstad_like_E, p_cells_jorstad_like_I,
           EI_ratio_jorstad_like, log2_EI_ratio_jorstad_like,
           n_donors_EI, n_cells_EI_total,
           dominant_ROIGroupCoarse, dominant_fraction) %>%
    arrange(anatomy_group),
  n = Inf
)

write.csv(
  ei_region,
  "data_analysis/jorstad_like_cortical_EI_ratio_by_region_with_rcmr.csv",
  row.names = FALSE
)

###########
# Analyses
###########
safe_spearman <- function(df, label, variable = "EI_ratio_jorstad_like") {
  x <- df$rcmr_value
  y <- df[[variable]]
  ok <- complete.cases(x, y) & is.finite(x) & is.finite(y)
  n <- sum(ok)

  if (n < 4 || dplyr::n_distinct(x[ok]) < 2 || dplyr::n_distinct(y[ok]) < 2) {
    return(tibble(
      subset = label,
      variable = variable,
      rho = NA_real_,
      p_value = NA_real_,
      n_regions = n,
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
    n_regions = n,
    mean_EI_ratio = mean(y[ok], na.rm = TRUE),
    sd_EI_ratio = sd(y[ok], na.rm = TRUE),
    status = "tested"
  )
}

run_lm_subset <- function(df, label, variable = "EI_ratio_jorstad_like") {
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

make_cortical_plot <- function(df, variable = "EI_ratio_jorstad_like", file_slug = "cortex_jorstad_like") {
  plot_df <- df %>%
    filter(is.finite(.data[[variable]]), !is.na(rcmr_value))

  if (nrow(plot_df) < 3 || dplyr::n_distinct(plot_df[[variable]]) < 2) {
    cat("  Skipping cortex-only plot: n=", nrow(plot_df), " finite regions.\n", sep = "")
    return(invisible(NULL))
  }

  region_scale <- tryCatch({
    check_region_palette(plot_df, region_col = "anatomy_group")
    present <- intersect(region_order, unique(as.character(plot_df$anatomy_group)))
    plot_df$anatomy_group <- factor(plot_df$anatomy_group, levels = present)
    ggplot2::scale_color_manual(values = region_palette, drop = TRUE)
  }, error = function(e) {
    message("Region palette unavailable for cortex-only plot (", e$message,
            "); using default ggplot palette.")
    NULL
  })

  p <- ggplot(plot_df, aes(x = .data[[variable]], y = rcmr_value, color = anatomy_group)) +
    geom_point(size = 3.1, alpha = 0.9) +
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
      title = "Cortex-only Jorstad-like E:I ratio vs rCMRGlc",
      subtitle = paste0(
        "n = ", nrow(plot_df),
        " cortical rCMR regions; E:I uses only cortical IT/CT/6b/NP excitatory and CGE/MGE/LAMP5-LHX6/Chandelier interneuron superclusters"
      ),
      x = "Jorstad-like cortical E:I ratio",
      y = "rCMRGlc (µmol/100 g/min.)",
      color = "Cortical region"
    ) +
    theme_classic(base_size = 13)
  if (!is.null(region_scale)) p <- p + region_scale

  print(p)

  ggsave(
    filename = paste0("figs/s1b/p_", file_slug, ".pdf"),
    plot = p, width = 8.5, height = 7, units = "in"
  )
  ggsave(
    filename = paste0("figs/s1b/p_", file_slug, ".jpg"),
    plot = p, width = 8.5, height = 7, units = "in", dpi = 300
  )
  invisible(p)
}

fit_cortex <- run_lm_subset(ei_region, "Cortex-only Jorstad-like E:I")

cor_cortex <- safe_spearman(ei_region, "Cortex-only Jorstad-like E:I")
cat("\n========== SPEARMAN CORRELATION: rCMRGlc vs cortex-only Jorstad-like E:I ratio ==========" , "\n")
print(cor_cortex)

write.csv(
  cor_cortex,
  "data_analysis/spearman_correlation_rcmr_jorstad_like_cortical_EI_ratio.csv",
  row.names = FALSE
)

p_cortex <- make_cortical_plot(ei_region)

###################################################################
# Optional diagnostic: log2 transformed ratio
###################################################################
if (nrow(ei_region %>% filter(is.finite(log2_EI_ratio_jorstad_like))) >= 3) {
  fit_log <- run_lm_subset(ei_region, "Cortex-only log2 Jorstad-like E:I", "log2_EI_ratio_jorstad_like")
  cor_log <- safe_spearman(ei_region, "Cortex-only log2 Jorstad-like E:I", "log2_EI_ratio_jorstad_like")
  write.csv(
    cor_log,
    "data_analysis/spearman_correlation_rcmr_log2_jorstad_like_cortical_EI_ratio.csv",
    row.names = FALSE
  )
}

cat("\nDone. Outputs:\n")
cat("  data_analysis/siletti_to_jorstad_like_cortical_EI_crosswalk.csv\n")
cat("  data_analysis/cortex_classification_neuronal.csv\n")
cat("  data_analysis/siletti_cluster_subcluster_membership_jorstad_like_cortical_EI.csv\n")
cat("  data_analysis/jorstad_like_cortical_EI_ratio_by_region_with_rcmr.csv\n")
cat("  data_analysis/spearman_correlation_rcmr_jorstad_like_cortical_EI_ratio.csv\n")
cat("  figs/s1b/p_cortex_jorstad_like.{pdf,jpg}\n")
