setwd("~/Library/CloudStorage/Dropbox/COLLABORATIVE/Do expensive brain regions increase less in humans/analyses_metabol_rate_structure")

# Compare broad all-region E:I ratios with a Jorstad-like cortical E:I ratio.
#
# This rewritten version makes two primary overlay plots using the raw, unlogged E:I ratio on the x-axis:
#   1. Broad all-region E:I with medium spiny neurons excluded from the I denominator.
#   2. Broad all-region E:I with medium spiny neurons included in the I denominator.
#
# The blue Jorstad-like cortical points/line are unchanged in both plots.
# Plot captions were removed to keep the figure visually cleaner; key definitions are in titles/subtitles.
# Main plot outputs:
#   figs/s1b/p_EI_original_all_regions_with_jorstad_cortex_overlay_raw_EI.{pdf,jpg}
#   figs/s1b/p_EI_original_with_MSN_all_regions_with_jorstad_cortex_overlay_raw_EI.{pdf,jpg}
#   data_analysis/EI_original_no_MSN_with_MSN_vs_jorstad_comparison_table.csv
#   data_analysis/spearman_EI_original_no_MSN_with_MSN_vs_jorstad_comparison.csv
#   data_analysis/lm_EI_original_no_MSN_with_MSN_vs_jorstad_comparison.csv
#   data_analysis/broad_EI_MSN_sensitivity_table.csv
#   data_analysis/siletti_to_jorstad_like_cortical_EI_crosswalk.csv
#   data_analysis/siletti_cluster_subcluster_membership_jorstad_like_cortical_EI.csv

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggpmisc)
})

has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)

# Optional local plot settings. The script still runs if unavailable.
if (file.exists("R/plot_settings.R")) {
  source("R/plot_settings.R")
}

if (!dir.exists("data_analysis")) dir.create("data_analysis", recursive = TRUE)
if (!dir.exists("figs/s1b")) dir.create("figs/s1b", recursive = TRUE)

DEF_BROAD_NO_MSN <- "Broad all-region E:I (MSN excluded)"
DEF_BROAD_WITH_MSN <- "Broad all-region E:I (MSN included)"
DEF_JORSTAD <- "Jorstad-like cortical E:I"

############################
## Load Siletti obs metadata
############################
obs <- readRDS("data_intermediate/linnarsson_adult_human_brain_obs_metadata_neuronal.rds")

obs <- obs %>%
  filter(cell_type == "neuron") %>%
  mutate(
    across(
      any_of(c("roi", "supercluster_term", "ROIGroup", "ROIGroupCoarse", "ROIGroupFine", "dissection", "donor_id")),
      ~ stringr::str_squish(as.character(.x))
    )
  )

######################################################
# Read rCMRGlc values from Heiss/Stephan table
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
# Map Siletti ROI to rCMR anatomy terms
#####################################
anatomy_rules <- readr::read_csv(
  "data_intermediate/rcmr_roi_relationship.csv",
  show_col_types = FALSE
) %>%
  transmute(
    anatomy_group = rcmr_term,
    roi = rois
  ) %>%
  tidyr::separate_rows(roi, sep = "\\s*\\|\\|\\s*") %>%
  mutate(roi = stringr::str_squish(roi)) %>%
  distinct(roi, anatomy_group)

obs <- obs %>%
  select(-any_of("anatomy_group")) %>%
  left_join(anatomy_rules, by = "roi") %>%
  mutate(anatomy_group = coalesce(anatomy_group, "Unmapped"))

###################################################################
# Broad all-region E:I category definition
###################################################################
# Broad E numerator:
#   all Excitatory_projection superclusters listed below.
# Broad I denominator, no-MSN version:
#   all Inhibitory_interneuron superclusters listed below.
# Broad I denominator, MSN-included version:
#   Inhibitory_interneuron + Inhibitory_projection_MSN.
#
# The two broad definitions are kept side by side so that the effect of
# including MSNs can be inspected directly.

broad_E_terms <- c(
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
)

broad_I_interneuron_terms <- c(
  "CGE interneuron",
  "MGE interneuron",
  "LAMP5-LHX6 and Chandelier",
  "Cerebellar inhibitory",
  "Midbrain-derived inhibitory"
)

broad_I_MSN_terms <- c(
  "Medium spiny neuron",
  "Eccentric medium spiny neuron"
)

obs <- obs %>%
  mutate(
    cell_category_recomputed = case_when(
      supercluster_term %in% broad_E_terms ~ "Excitatory_projection",
      supercluster_term %in% broad_I_interneuron_terms ~ "Inhibitory_interneuron",
      supercluster_term %in% broad_I_MSN_terms ~ "Inhibitory_projection_MSN",
      supercluster_term == "Splatter" ~ "Splatter",
      supercluster_term == "Miscellaneous" ~ "Miscellaneous",
      TRUE ~ "Other_neuron"
    )
  )

cat("\n================ Recomputed broad cell category counts ================\n")
print(as_tibble(obs %>% count(cell_category_recomputed, sort = TRUE)), n = Inf)

###################################################################
# Cortex classification
###################################################################
# Prefer dominant ROIGroupCoarse == "Cerebral cortex" when available.
# Backup text regex protects against local metadata naming differences.

cortex_annotation_regex <- "cortex|cortical|neocortex|gyrus|precentral|postcentral|frontal|temporal|parietal|occipital|cingulate|insula"

telencephalon_coarse_groups <- c("Cerebral cortex", "Cerebral nuclei", "Hippocampus")

region_classification <- obs %>%
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
    is_cortex_by_dominant_group = dominant_ROIGroupCoarse == "Cerebral cortex",
    is_telencephalon = dominant_ROIGroupCoarse %in% telencephalon_coarse_groups
  ) %>%
  left_join(
    obs %>%
      filter(anatomy_group != "", anatomy_group != "Unmapped") %>%
      group_by(anatomy_group) %>%
      summarise(
        any_cortical_text = any(
          str_detect(
            str_to_lower(paste(ROIGroup, ROIGroupCoarse, ROIGroupFine, dissection, sep = " ")),
            cortex_annotation_regex
          ),
          na.rm = TRUE
        ),
        .groups = "drop"
      ),
    by = "anatomy_group"
  ) %>%
  mutate(
    is_cortex = is_cortex_by_dominant_group | any_cortical_text,
    cortex_call = if_else(is_cortex, "Cortex", "Non-cortex"),
    division = if_else(is_telencephalon, "Telencephalon", "Non-telencephalon")
  ) %>%
  arrange(desc(is_cortex), anatomy_group)

cat("\n================ Region classification ================\n")
print(as_tibble(region_classification), n = Inf)

write.csv(region_classification, "data_analysis/region_classification_for_EI_comparison.csv", row.names = FALSE)

###################################################################
# Jorstad-like cortical E:I definition
###################################################################
# Closest Siletti supercluster-level match to Jorstad et al. neocortical E:I.
# Jorstad excitatory subclasses: L2/3 IT, L4 IT, L5 IT, L6 IT, L6 IT Car3,
# L5 ET, L5/6 NP, L6b, L6 CT.
# Jorstad inhibitory subclasses: LAMP5 LHX6, LAMP5, SNCG, VIP, PAX6,
# Chandelier, PVALB, SST, SST CHODL.

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

cluster_membership_check <- obs %>%
  inner_join(siletti_jorstad_crosswalk, by = "supercluster_term") %>%
  count(EI_class_jorstad_like, jorstad_equivalent, supercluster_term,
        cluster_id, subcluster_id, sort = TRUE) %>%
  arrange(EI_class_jorstad_like, jorstad_equivalent, supercluster_term,
          cluster_id, subcluster_id)

write.csv(
  cluster_membership_check,
  "data_analysis/siletti_cluster_subcluster_membership_jorstad_like_cortical_EI.csv",
  row.names = FALSE
)

###################################################################
# Helper: donor-balanced regional proportions
###################################################################
# We compute within-donor regional proportions first, then average donors.
# This avoids letting donors with more nuclei dominate regional estimates.

calculate_region_props <- function(obs_in, class_col, included_classes, value_prefix) {
  stopifnot(class_col %in% colnames(obs_in))

  class_tbl <- tibble(class_value = included_classes)

  obs_use <- obs_in %>%
    filter(
      anatomy_group != "",
      anatomy_group != "Unmapped",
      .data[[class_col]] %in% included_classes
    ) %>%
    mutate(class_value = .data[[class_col]])

  donor_region <- obs_use %>% distinct(anatomy_group, donor_id)

  counts <- obs_use %>%
    count(anatomy_group, donor_id, class_value, name = "n_cells")

  long <- donor_region %>%
    tidyr::crossing(class_tbl) %>%
    left_join(counts, by = c("anatomy_group", "donor_id", "class_value")) %>%
    mutate(n_cells = replace_na(n_cells, 0)) %>%
    group_by(anatomy_group, donor_id) %>%
    mutate(
      donor_total_cells = sum(n_cells),
      donor_prop = if_else(donor_total_cells > 0, n_cells / donor_total_cells, NA_real_)
    ) %>%
    ungroup() %>%
    group_by(anatomy_group, class_value) %>%
    summarise(
      p_cells = mean(donor_prop, na.rm = TRUE),
      n_donors = n_distinct(donor_id[donor_total_cells > 0]),
      n_cells_total = sum(n_cells, na.rm = TRUE),
      .groups = "drop"
    )

  wide <- long %>%
    select(anatomy_group, class_value, p_cells) %>%
    pivot_wider(
      names_from = class_value,
      values_from = p_cells,
      values_fill = 0,
      names_prefix = value_prefix
    ) %>%
    left_join(
      long %>%
        group_by(anatomy_group) %>%
        summarise(
          n_donors_EI = max(n_donors, na.rm = TRUE),
          n_cells_EI_total = sum(n_cells_total, na.rm = TRUE),
          .groups = "drop"
        ),
      by = "anatomy_group"
    )

  wide
}

###################################################################
# Broad all-region E:I, calculated both without and with MSNs
###################################################################
original_wide <- calculate_region_props(
  obs_in = obs,
  class_col = "cell_category_recomputed",
  included_classes = c("Excitatory_projection", "Inhibitory_interneuron", "Inhibitory_projection_MSN"),
  value_prefix = "p_original_"
) %>%
  inner_join(rcmr, by = c("anatomy_group" = "rcmr_term")) %>%
  left_join(region_classification, by = "anatomy_group") %>%
  mutate(
    p_original_Excitatory_projection = coalesce(p_original_Excitatory_projection, 0),
    p_original_Inhibitory_interneuron = coalesce(p_original_Inhibitory_interneuron, 0),
    p_original_Inhibitory_projection_MSN = coalesce(p_original_Inhibitory_projection_MSN, 0),
    p_original_Inhibitory_all_with_MSN =
      p_original_Inhibitory_interneuron + p_original_Inhibitory_projection_MSN,
    EI_ratio_no_MSN = if_else(
      p_original_Inhibitory_interneuron > 0,
      p_original_Excitatory_projection / p_original_Inhibitory_interneuron,
      NA_real_
    ),
    log2_EI_ratio_no_MSN = if_else(EI_ratio_no_MSN > 0, log2(EI_ratio_no_MSN), NA_real_),
    # Legacy aliases: these retain the previous MSN-excluded definition.
    EI_ratio = EI_ratio_no_MSN,
    log2_EI_ratio = log2_EI_ratio_no_MSN,
    EI_ratio_with_MSN = if_else(
      p_original_Inhibitory_all_with_MSN > 0,
      p_original_Excitatory_projection / p_original_Inhibitory_all_with_MSN,
      NA_real_
    ),
    log2_EI_ratio_with_MSN = if_else(EI_ratio_with_MSN > 0, log2(EI_ratio_with_MSN), NA_real_),
    MSN_fraction_of_broad_I = if_else(
      p_original_Inhibitory_all_with_MSN > 0,
      p_original_Inhibitory_projection_MSN / p_original_Inhibitory_all_with_MSN,
      NA_real_
    ),
    delta_log2_EI_with_MSN_minus_no_MSN =
      log2_EI_ratio_with_MSN - log2_EI_ratio_no_MSN
  )

write.csv(original_wide, "data_analysis/original_broad_EI_ratio_all_regions_with_rcmr.csv", row.names = FALSE)

original_no_msn_long <- original_wide %>%
  transmute(
    anatomy_group,
    rcmr_value,
    EI_definition = DEF_BROAD_NO_MSN,
    EI_ratio = EI_ratio_no_MSN,
    log2_EI_ratio = log2_EI_ratio_no_MSN,
    n_donors_EI,
    n_cells_EI_total,
    cortex_call,
    is_cortex,
    division,
    MSN_included_in_denominator = FALSE,
    p_E = p_original_Excitatory_projection,
    p_I = p_original_Inhibitory_interneuron,
    p_I_interneuron = p_original_Inhibitory_interneuron,
    p_I_MSN = p_original_Inhibitory_projection_MSN,
    p_I_all_with_MSN = p_original_Inhibitory_all_with_MSN,
    EI_ratio_no_MSN,
    log2_EI_ratio_no_MSN,
    EI_ratio_with_MSN,
    log2_EI_ratio_with_MSN,
    MSN_fraction_of_broad_I,
    delta_log2_EI_with_MSN_minus_no_MSN
  )

original_with_msn_long <- original_wide %>%
  transmute(
    anatomy_group,
    rcmr_value,
    EI_definition = DEF_BROAD_WITH_MSN,
    EI_ratio = EI_ratio_with_MSN,
    log2_EI_ratio = log2_EI_ratio_with_MSN,
    n_donors_EI,
    n_cells_EI_total,
    cortex_call,
    is_cortex,
    division,
    MSN_included_in_denominator = TRUE,
    p_E = p_original_Excitatory_projection,
    p_I = p_original_Inhibitory_all_with_MSN,
    p_I_interneuron = p_original_Inhibitory_interneuron,
    p_I_MSN = p_original_Inhibitory_projection_MSN,
    p_I_all_with_MSN = p_original_Inhibitory_all_with_MSN,
    EI_ratio_no_MSN,
    log2_EI_ratio_no_MSN,
    EI_ratio_with_MSN,
    log2_EI_ratio_with_MSN,
    MSN_fraction_of_broad_I,
    delta_log2_EI_with_MSN_minus_no_MSN
  )

msn_sensitivity_tbl <- original_wide %>%
  transmute(
    anatomy_group,
    rcmr_value,
    cortex_call,
    is_cortex,
    division,
    p_E = p_original_Excitatory_projection,
    p_I_interneuron = p_original_Inhibitory_interneuron,
    p_I_MSN = p_original_Inhibitory_projection_MSN,
    p_I_all_with_MSN = p_original_Inhibitory_all_with_MSN,
    MSN_fraction_of_broad_I,
    EI_ratio_no_MSN,
    log2_EI_ratio_no_MSN,
    EI_ratio_with_MSN,
    log2_EI_ratio_with_MSN,
    delta_log2_EI_with_MSN_minus_no_MSN
  ) %>%
  arrange(desc(MSN_fraction_of_broad_I), anatomy_group)

write.csv(msn_sensitivity_tbl, "data_analysis/broad_EI_MSN_sensitivity_table.csv", row.names = FALSE)

###################################################################
# Jorstad-like cortex-only E:I
###################################################################
obs_jorstad_cortex <- obs %>%
  inner_join(region_classification %>% filter(is_cortex) %>% select(anatomy_group), by = "anatomy_group") %>%
  inner_join(siletti_jorstad_crosswalk, by = "supercluster_term")

cat("\n================ Jorstad-like cortical E/I cell counts ================\n")
print(
  as_tibble(obs_jorstad_cortex %>% count(EI_class_jorstad_like, supercluster_term, sort = TRUE)),
  n = Inf
)

jorstad_wide <- calculate_region_props(
  obs_in = obs_jorstad_cortex,
  class_col = "EI_class_jorstad_like",
  included_classes = c("E", "I"),
  value_prefix = "p_jorstad_like_"
) %>%
  inner_join(rcmr, by = c("anatomy_group" = "rcmr_term")) %>%
  left_join(region_classification, by = "anatomy_group") %>%
  filter(is_cortex) %>%
  mutate(
    p_jorstad_like_E = coalesce(p_jorstad_like_E, 0),
    p_jorstad_like_I = coalesce(p_jorstad_like_I, 0),
    EI_ratio = if_else(
      p_jorstad_like_I > 0,
      p_jorstad_like_E / p_jorstad_like_I,
      NA_real_
    ),
    log2_EI_ratio = if_else(EI_ratio > 0, log2(EI_ratio), NA_real_),
    EI_definition = DEF_JORSTAD
  )

write.csv(jorstad_wide, "data_analysis/jorstad_like_cortical_EI_ratio_by_region_with_rcmr.csv", row.names = FALSE)

jorstad_long <- jorstad_wide %>%
  transmute(
    anatomy_group,
    rcmr_value,
    EI_definition,
    EI_ratio,
    log2_EI_ratio,
    n_donors_EI,
    n_cells_EI_total,
    cortex_call,
    is_cortex,
    division,
    MSN_included_in_denominator = NA,
    p_E = p_jorstad_like_E,
    p_I = p_jorstad_like_I,
    p_I_interneuron = p_jorstad_like_I,
    p_I_MSN = NA_real_,
    p_I_all_with_MSN = NA_real_,
    EI_ratio_no_MSN = NA_real_,
    log2_EI_ratio_no_MSN = NA_real_,
    EI_ratio_with_MSN = NA_real_,
    log2_EI_ratio_with_MSN = NA_real_,
    MSN_fraction_of_broad_I = NA_real_,
    delta_log2_EI_with_MSN_minus_no_MSN = NA_real_
  )

###################################################################
# Combined comparison table
###################################################################
comparison_df <- bind_rows(
  original_no_msn_long,
  original_with_msn_long,
  jorstad_long
) %>%
  filter(is.finite(EI_ratio), is.finite(rcmr_value), !is.na(rcmr_value))

write.csv(
  comparison_df,
  "data_analysis/EI_original_no_MSN_with_MSN_vs_jorstad_comparison_table.csv",
  row.names = FALSE
)
# Compatibility copy using the previous summary-table filename.
write.csv(
  comparison_df,
  "data_analysis/EI_original_vs_jorstad_comparison_table.csv",
  row.names = FALSE
)

###################################################################
# Correlation and LM summaries
###################################################################
safe_spearman <- function(df, variable = "log2_EI_ratio") {
  x <- df[[variable]]
  y <- df$rcmr_value
  ok <- complete.cases(x, y) & is.finite(x) & is.finite(y)

  if (sum(ok) < 4 || n_distinct(x[ok]) < 2 || n_distinct(y[ok]) < 2) {
    return(tibble(
      variable = variable,
      n_regions = sum(ok),
      rho = NA_real_,
      p_value = NA_real_,
      status = "not_tested_constant_or_low_n"
    ))
  }

  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  tibble(
    variable = variable,
    n_regions = sum(ok),
    rho = unname(ct$estimate),
    p_value = ct$p.value,
    status = "tested"
  )
}

spearman_tbl <- comparison_df %>%
  group_by(EI_definition) %>%
  group_modify(~ bind_rows(
    safe_spearman(.x, "EI_ratio"),
    safe_spearman(.x, "log2_EI_ratio")
  )) %>%
  ungroup() %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH"))

cat("\n================ Spearman comparison ================\n")
print(as_tibble(spearman_tbl), n = Inf)

write.csv(
  spearman_tbl,
  "data_analysis/spearman_EI_original_no_MSN_with_MSN_vs_jorstad_comparison.csv",
  row.names = FALSE
)
# Compatibility copy using the previous Spearman-summary filename.
write.csv(
  spearman_tbl,
  "data_analysis/spearman_EI_original_vs_jorstad_comparison.csv",
  row.names = FALSE
)

lm_summary_one <- function(df, variable = "log2_EI_ratio") {
  df2 <- df %>% filter(is.finite(.data[[variable]]), !is.na(rcmr_value))
  if (nrow(df2) < 4 || n_distinct(df2[[variable]]) < 2) {
    return(tibble(
      variable = variable,
      n_regions = nrow(df2),
      beta = NA_real_,
      se = NA_real_,
      t_value = NA_real_,
      p_value = NA_real_,
      r_squared = NA_real_,
      adj_r_squared = NA_real_,
      status = "not_tested_constant_or_low_n"
    ))
  }
  fit <- lm(reformulate(variable, response = "rcmr_value"), data = df2)
  sm <- summary(fit)
  coef_tbl <- coef(sm)
  tibble(
    variable = variable,
    n_regions = nrow(df2),
    beta = coef_tbl[variable, "Estimate"],
    se = coef_tbl[variable, "Std. Error"],
    t_value = coef_tbl[variable, "t value"],
    p_value = coef_tbl[variable, "Pr(>|t|)"],
    r_squared = sm$r.squared,
    adj_r_squared = sm$adj.r.squared,
    status = "tested"
  )
}

lm_tbl <- comparison_df %>%
  group_by(EI_definition) %>%
  group_modify(~ bind_rows(
    lm_summary_one(.x, "EI_ratio"),
    lm_summary_one(.x, "log2_EI_ratio")
  )) %>%
  ungroup() %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH"))

cat("\n================ LM comparison ================\n")
print(as_tibble(lm_tbl), n = Inf)

write.csv(
  lm_tbl,
  "data_analysis/lm_EI_original_no_MSN_with_MSN_vs_jorstad_comparison.csv",
  row.names = FALSE
)
# Compatibility copy using the previous LM-summary filename.
write.csv(
  lm_tbl,
  "data_analysis/lm_EI_original_vs_jorstad_comparison.csv",
  row.names = FALSE
)

###################################################################
# Plot helpers
###################################################################
plot_df_original_no_msn <- comparison_df %>% filter(EI_definition == DEF_BROAD_NO_MSN)
plot_df_original_with_msn <- comparison_df %>% filter(EI_definition == DEF_BROAD_WITH_MSN)
plot_df_jorstad <- comparison_df %>% filter(EI_definition == DEF_JORSTAD)

can_add_fit_label <- function(df) {
  nrow(df) >= 4 &&
    n_distinct(df$EI_ratio[is.finite(df$EI_ratio)]) >= 2 &&
    n_distinct(df$rcmr_value[is.finite(df$rcmr_value)]) >= 2
}

make_overlay_plot <- function(broad_df,
                              plot_title,
                              plot_subtitle,
                              output_stem) {
  p <- ggplot() +
    geom_point(
      data = broad_df,
      aes(x = EI_ratio, y = rcmr_value),
      color = "grey70",
      alpha = 0.65,
      size = 2.4
    ) +
    geom_smooth(
      data = broad_df,
      aes(x = EI_ratio, y = rcmr_value),
      method = "lm",
      se = TRUE,
      color = "grey45",
      fill = "grey85",
      linewidth = 0.8
    ) +
    geom_point(
      data = plot_df_jorstad,
      aes(x = EI_ratio, y = rcmr_value),
      color = "#0072B2",
      alpha = 0.95,
      size = 3.1
    ) +
    geom_smooth(
      data = plot_df_jorstad,
      aes(x = EI_ratio, y = rcmr_value),
      method = "lm",
      se = TRUE,
      color = "#0072B2",
      fill = "#D8ECF7",
      linewidth = 1.0
    ) +
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      x = "E:I ratio",
      y = "rCMRGlc (umol/100 g/min.)"
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 13)
    )

  # Per-line regression equations + R^2 + p, colour-matched to each fit.
  if (can_add_fit_label(broad_df)) {
    p <- p +
      stat_poly_eq(
        data = broad_df,
        aes(
          x = EI_ratio,
          y = rcmr_value,
          label = paste(after_stat(eq.label), after_stat(rr.label), after_stat(p.value.label),
                        sep = "*\", \"*")
        ),
        formula = y ~ x,
        parse = TRUE,
        color = "grey45",
        label.x = "left",
        label.y = 0.99,
        size = 3
      )
  }

  if (can_add_fit_label(plot_df_jorstad)) {
    p <- p +
      stat_poly_eq(
        data = plot_df_jorstad,
        aes(
          x = EI_ratio,
          y = rcmr_value,
          label = paste(after_stat(eq.label), after_stat(rr.label), after_stat(p.value.label),
                        sep = "*\", \"*")
        ),
        formula = y ~ x,
        parse = TRUE,
        color = "#0072B2",
        label.x = "left",
        label.y = 0.91,
        size = 3
      )
  }

  if (has_ggrepel) {
    p <- p +
      ggrepel::geom_text_repel(
        data = plot_df_jorstad,
        aes(x = EI_ratio, y = rcmr_value, label = anatomy_group),
        color = "#0072B2",
        size = 3.0,
        max.overlaps = Inf,
        min.segment.length = 0,
        box.padding = 0.3
      )
  } else {
    message("Package ggrepel is not installed; skipping cortical region labels on overlay plot.")
  }

  print(p)
  ggsave(paste0(output_stem, ".pdf"), p, width = 8.5, height = 7, units = "in")
  ggsave(paste0(output_stem, ".jpg"), p, width = 8.5, height = 7, units = "in", dpi = 300)

  p
}

###################################################################
# Plot 1: broad all-region E:I with MSNs excluded
###################################################################
p_overlay_no_msn <- make_overlay_plot(
  broad_df = plot_df_original_no_msn,
  plot_title = "Regional glucose metabolism vs raw neuronal E:I ratio:\nbroad all-region, MSN excluded, with Jorstad-like cortical overlay",
  plot_subtitle = "X-axis is raw E:I, not log2-transformed. Each cortical region appears twice: grey = broad all-region E:I; blue = Jorstad-like cortical E:I.",
  output_stem = "figs/s1b/p_EI_original_all_regions_with_jorstad_cortex_overlay_raw_EI"
)

###################################################################
# Plot 2: broad all-region E:I with MSNs included
###################################################################
p_overlay_with_msn <- make_overlay_plot(
  broad_df = plot_df_original_with_msn,
  plot_title = "Regional glucose metabolism vs raw neuronal E:I ratio:\nbroad all-region, MSN included, with Jorstad-like cortical overlay",
  plot_subtitle = "X-axis is raw E:I, not log2-transformed. Grey = broad all-region E:I with MSN included; blue = Jorstad-like cortical E:I.",
  output_stem = "figs/s1b/p_EI_original_with_MSN_all_regions_with_jorstad_cortex_overlay_raw_EI"
)

cat("\nDone. Key outputs:\n")
cat("  data_analysis/EI_original_no_MSN_with_MSN_vs_jorstad_comparison_table.csv\n")
cat("  data_analysis/broad_EI_MSN_sensitivity_table.csv\n")
cat("  data_analysis/spearman_EI_original_no_MSN_with_MSN_vs_jorstad_comparison.csv\n")
cat("  data_analysis/lm_EI_original_no_MSN_with_MSN_vs_jorstad_comparison.csv\n")
cat("  figs/s1b/p_EI_original_all_regions_with_jorstad_cortex_overlay_raw_EI.{pdf,jpg}\n")
cat("  figs/s1b/p_EI_original_with_MSN_all_regions_with_jorstad_cortex_overlay_raw_EI.{pdf,jpg}\n")
