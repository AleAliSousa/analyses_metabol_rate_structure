setwd("~/Library/CloudStorage/Dropbox/COLLABORATIVE/Do expensive brain regions increase less in humans/analyses_metabol_rate_structure")

source("R/plot_settings.R")

## Install and Load up packages
library(ggplot2)
library(ggpmisc)
library(readxl)
library(tidyverse)

############################
## Load saved obs metadata
############################
obs <- readRDS("data_intermediate/linnarsson_adult_human_brain_obs_metadata_nonneuronal.rds")

######################################################
# Read rCMRGlc values from Heiss et al. 2004
######################################################

# Read Table with rCMRGlc values
heiss_stephan_tbl <- read.csv("data_intermediate/Heiss_Stephan_data.csv")

# Keep only the relevant columns and no average columns; ensure rcmr_value is numeric and rounded to 1 decimal place
rcmr <- heiss_stephan_tbl %>%
  transmute(
    rcmr_term = str_trim(as.character(term_3)),     # term_3 is the original list from Heiss
    rcmr_value = round(as.numeric(rCMRGlc_mean_both_hemispheres), 1)
  ) %>%
  filter(
    !is.na(rcmr_term),
    rcmr_term != "",
    !str_detect(str_to_lower(rcmr_term), "average"),     # Drop rows where the term contains "average"
    !is.na(rcmr_value)
  )

head(rcmr)

#####################################
# Anatomical grouping of rois
#####################################

# Read the roi-to-rCMR mapping table directly.
# Expected columns:
#   rcmr_term
#   rois: one or more exact obs$roi strings separated by "||"

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
  dplyr::mutate(roi = stringr::str_squish(roi)) %>%
  dplyr::left_join(anatomy_rules, by = "roi") %>%
  dplyr::mutate(anatomy_group = dplyr::coalesce(anatomy_group, "Unmapped"))

###################################################################
# Telencephalon classification
###################################################################
# Telencephalon (forebrain - rostral): cerebral cortex, cerebral nuclei
# (basal ganglia + amygdala + claustrum + basal forebrain), and hippocampus.
# Non-telencephalon: diencephalon (thalamus, hypothalamus, epithalamus),
# midbrain, pons, medulla (myelencephalon), cerebellum, spinal cord.
#
# We classify each anatomy_group by the dominant ROIGroupCoarse of its
# constituent cells. This is data-driven and robust to label variants.

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

# Persist the classification table
write.csv(
  telencephalon_table,
  "data_analysis/telencephalon_classification.csv",
  row.names = FALSE
)


###################################################################
# Astrocyte Type 1 / Type 2 proportion calculations
# Mean participant-level proportions by anatomy_group
###################################################################
# The denominator is the full nonneuronal compartment in each
# donor × anatomy_group sample. The only numerators retained are
# Type 1 and Type 2 astrocytes.
#
# Cluster mapping inferred by overlaying Fig. 5A onto Fig. S14A:
#   Type 1 astrocytes: clusters 52, 53, 54, 55, 56, 57, 58, 61
#   Type 2 astrocytes: clusters 59, 60, 62, 63, 64
#
# If you later obtain an official machine-readable Type 1 / Type 2
# annotation, replace these two vectors and rerun the script.
###################################################################

type1_astro_cluster_ids <- c(52, 53, 54, 55, 56, 57, 58, 61)
type2_astro_cluster_ids <- c(59, 60, 62, 63, 64)
all_astro_type_cluster_ids <- c(type1_astro_cluster_ids, type2_astro_cluster_ids)

astro_type_levels <- c("Astrocyte_Type1", "Astrocyte_Type2")

# Keep mapped regions only; obs is already nonneuronal metadata.
obs_celltype <- obs %>%
  filter(anatomy_group != "", anatomy_group != "Unmapped") %>%
  mutate(
    cluster_id_int = as.integer(stringr::str_extract(as.character(cluster_id), "\\d+")),
    is_astrocyte = stringr::str_detect(
      as.character(supercluster_term),
      stringr::regex("astro", ignore_case = TRUE)
    ),
    astro_type = dplyr::case_when(
      is_astrocyte & cluster_id_int %in% type1_astro_cluster_ids ~ "Astrocyte_Type1",
      is_astrocyte & cluster_id_int %in% type2_astro_cluster_ids ~ "Astrocyte_Type2",
      TRUE ~ NA_character_
    ),
    astro_type = factor(astro_type, levels = astro_type_levels)
  )

# QC: inspect which astrocyte clusters are assigned to each type and region.
astro_cluster_qc <- obs_celltype %>%
  filter(is_astrocyte) %>%
  count(
    astro_type,
    cluster_id_int,
    supercluster_term,
    ROIGroupCoarse,
    anatomy_group,
    name = "n_cells",
    sort = TRUE
  )

cat("\n================ Astrocyte cluster Type 1 / Type 2 QC ================\n")
print(
  astro_cluster_qc %>%
    group_by(astro_type, cluster_id_int) %>%
    summarise(n_cells = sum(n_cells), .groups = "drop") %>%
    arrange(astro_type, cluster_id_int),
  n = Inf
)

unassigned_astro_clusters <- obs_celltype %>%
  filter(is_astrocyte, is.na(astro_type)) %>%
  count(cluster_id_int, supercluster_term, name = "n_cells", sort = TRUE)

if (nrow(unassigned_astro_clusters) > 0) {
  warning(
    "Some astrocyte clusters were not assigned to Type 1 or Type 2. ",
    "Inspect data_analysis/astrocyte_type_unassigned_clusters.csv."
  )
}

write.csv(
  astro_cluster_qc,
  "data_analysis/astrocyte_type_cluster_qc_by_region.csv",
  row.names = FALSE
)
write.csv(
  unassigned_astro_clusters,
  "data_analysis/astrocyte_type_unassigned_clusters.csv",
  row.names = FALSE
)

# Donor × region denominators: all nonneuronal cells in the mapped sample.
donor_region <- obs_celltype %>%
  distinct(anatomy_group, donor_id)

donor_region_denominators <- obs_celltype %>%
  count(anatomy_group, donor_id, name = "donor_total_nonneuronal_cells")

# Numerators: Type 1 and Type 2 astrocytes only.
astrotype_counts <- obs_celltype %>%
  filter(!is.na(astro_type)) %>%
  count(anatomy_group, donor_id, astro_type, name = "n_cells")

# Add zero counts for astrocyte types absent from a donor-region sample,
# then calculate donor-level proportions relative to all nonneuronal cells.
celltype_table_long <- donor_region %>%
  tidyr::crossing(astro_type = factor(astro_type_levels, levels = astro_type_levels)) %>%
  left_join(
    astrotype_counts,
    by = c("anatomy_group", "donor_id", "astro_type")
  ) %>%
  left_join(
    donor_region_denominators,
    by = c("anatomy_group", "donor_id")
  ) %>%
  mutate(
    n_cells = tidyr::replace_na(n_cells, 0L),
    donor_prop = n_cells / donor_total_nonneuronal_cells
  ) %>%
  group_by(anatomy_group, astro_type) %>%
  summarise(
    p_cells = mean(donor_prop, na.rm = TRUE),
    n_donors = dplyr::n_distinct(donor_id),
    .groups = "drop"
  ) %>%
  rename(cell_type = astro_type)

# Wide table:
# p_* columns are mean participant-level proportions among all nonneuronal cells.
celltype_proportion_table <- celltype_table_long %>%
  select(anatomy_group, cell_type, p_cells) %>%
  pivot_wider(
    names_from  = cell_type,
    values_from = p_cells,
    values_fill = 0,
    names_prefix = "p_"
  )

write.csv(
  celltype_proportion_table,
  "data_analysis/type1_type2_astrocyte_mean_proportion_by_region.csv",
  row.names = FALSE
)

###################################################################
# Bind astrocyte Type 1 / Type 2 proportions to rCMRGlc table
# and telencephalon flag
###################################################################

analysis_df <- celltype_proportion_table %>%
  inner_join(rcmr, by = c("anatomy_group" = "rcmr_term")) %>%
  left_join(
    telencephalon_table %>% select(anatomy_group, is_telencephalon, division),
    by = "anatomy_group"
  )

cat("\n================ Regions in analysis by division ================\n")
print(
  analysis_df %>%
    count(division, name = "n_regions") %>%
    arrange(desc(n_regions))
)
cat("\nRegion list by division:\n")
print(
  analysis_df %>%
    select(division, anatomy_group, rcmr_value, starts_with("p_")) %>%
    arrange(division, anatomy_group),
  n = Inf
)

###########
# Analyses
###########

#############################
## PREDICTOR SETUP
#############################

predictors <- names(analysis_df)[grepl("^p_Astrocyte_Type", names(analysis_df))]

cat("\nPredictors used for Type 1 / Type 2 astrocyte analyses:\n")
print(predictors)

# Keep predictors that have sufficient non-missing values and nonzero variance
# within a subset. Zeros are allowed; constant predictors are not analyzable.
predictors_for_subset <- function(df, base_predictors) {
  keep <- vapply(base_predictors, function(v) {
    x <- df[[v]]
    x <- x[is.finite(x)]
    length(x) >= 4 && dplyr::n_distinct(x) > 1
  }, logical(1))
  base_predictors[keep]
}

###################################################################
# Helper functions: per-subset analysis (lm + correlations + plot)
###################################################################

run_lm_subset <- function(df, predictors, label) {
  n <- nrow(df)
  k <- length(predictors)
  df_resid <- n - k - 1
  cat("\n========== MULTIVARIATE LM:", label,
      "(n =", n, "regions, ", k, "predictors, df_residual =", df_resid, ") ==========\n")
  if (k == 0) {
    cat("  Skipping LM: no predictors with nonzero variance in this subset.\n")
    return(invisible(NULL))
  }
  if (df_resid < 2) {
    cat("  WARNING: df_residual =", df_resid,
        "-- model would be saturated. Refusing to fit.\n")
    return(invisible(NULL))
  }
  fit <- lm(
    reformulate(predictors, response = "rcmr_value"),
    data = df
  )
  print(summary(fit))
  invisible(fit)
}

run_cor_subset <- function(df, predictors, label) {
  cat("\n========== SPEARMAN CORRELATIONS:", label, "(n =", nrow(df), "regions) ==========\n")
  if (nrow(df) < 4 || length(predictors) == 0) {
    cat("  Skipping correlations: n=", nrow(df), ", predictors=", length(predictors), ".\n", sep = "")
    return(invisible(NULL))
  }
  cor_results <- do.call(
    rbind,
    lapply(predictors, function(v) {
      x <- df[[v]]
      ok <- complete.cases(df$rcmr_value, x)
      if (sum(ok) < 4 || dplyr::n_distinct(x[ok]) < 2) {
        return(data.frame(
          subset = label,
          variable = v,
          r = NA_real_,
          p_value = NA_real_,
          n = sum(ok)
        ))
      }
      ct <- suppressWarnings(cor.test(
        df$rcmr_value[ok], x[ok],
        method = "spearman", exact = FALSE
      ))
      data.frame(
        subset   = label,
        variable = v,
        r        = unname(ct$estimate),
        p_value  = ct$p.value,
        n        = sum(ok)
      )
    })
  ) %>%
    dplyr::arrange(p_value) %>%
    dplyr::mutate(p_adj_BH = p.adjust(p_value, method = "BH"))
  print(cor_results)
  invisible(cor_results)
}

make_plot_subset <- function(df, predictors, label, file_slug) {
  if (nrow(df) < 3 || length(predictors) == 0) {
    cat("  Skipping plot for ", label, ": n=", nrow(df), " regions, ",
        length(predictors), " predictors.\n", sep = "")
    return(invisible(NULL))
  }
  plot_df <- df %>%
    dplyr::select(anatomy_group, rcmr_value, all_of(predictors)) %>%
    tidyr::pivot_longer(
      cols = all_of(predictors),
      names_to = "predictor",
      values_to = "prop"
    ) %>%
    mutate(
      predictor = dplyr::recode(
        predictor,
        "p_Astrocyte_Type1" = "Type 1 astrocytes",
        "p_Astrocyte_Type2" = "Type 2 astrocytes"
      )
    )

  region_scale <- tryCatch({
    check_region_palette(df, region_col = "anatomy_group")
    present <- intersect(region_order, unique(as.character(df$anatomy_group)))
    plot_df$anatomy_group <- factor(plot_df$anatomy_group, levels = present)
    ggplot2::scale_color_manual(values = region_palette, drop = TRUE)
  }, error = function(e) {
    message("Region palette unavailable for ", label, " (", e$message,
            "); using default ggplot palette.")
    NULL
  })

  p <- ggplot(plot_df, aes(x = prop, y = rcmr_value, color = anatomy_group)) +
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
    facet_wrap(~ predictor, scales = "free_x") +
    labs(
      title = paste0("Astrocyte Type 1 / Type 2 proportions vs rCMRGlc — ", label),
      subtitle = paste0("n = ", nrow(df), " regions; ", length(predictors), " predictors"),
      x = "Mean cell type proportion among nonneuronal cells",
      y = "rCMRGlc (µmol/100 g/min.)",
      color = "Region"
    ) +
    theme_classic(base_size = 14)
  if (!is.null(region_scale)) p <- p + region_scale

  print(p)

  ggsave(
    filename = paste0("figs/s1b/p_type1_type2_astrocytes_", file_slug, ".pdf"),
    plot = p, width = 10, height = 7, units = "in"
  )
  ggsave(
    filename = paste0("figs/s1b/p_type1_type2_astrocytes_", file_slug, ".jpg"),
    plot = p, width = 10, height = 7, units = "in", dpi = 300
  )
  invisible(p)
}

###################################################################
# Run analyses for each subset
###################################################################

tel_df    <- analysis_df %>% filter(is_telencephalon)
nontel_df <- analysis_df %>% filter(!is_telencephalon)

predictors_tel    <- predictors_for_subset(tel_df,    predictors)
predictors_nontel <- predictors_for_subset(nontel_df, predictors)
predictors_union  <- union(predictors_tel, predictors_nontel)

cat("\n================ Predictors after per-subset variance check ================\n")
cat("Telencephalon predictors (", length(predictors_tel), "):\n", sep = "")
print(predictors_tel)
cat("Non-telencephalon predictors (", length(predictors_nontel), "):\n", sep = "")
print(predictors_nontel)
cat("Union used for combined plots (", length(predictors_union), "):\n", sep = "")
print(predictors_union)
cat("Dropped in Tel because constant or too sparse:    ", paste(setdiff(predictors, predictors_tel),    collapse = ", "), "\n")
cat("Dropped in Non-tel because constant or too sparse:", paste(setdiff(predictors, predictors_nontel), collapse = ", "), "\n")

write.csv(
  tel_df,
  "data_analysis/type1_type2_astrocyte_proportions_by_region_telencephalon.csv",
  row.names = FALSE
)
write.csv(
  nontel_df,
  "data_analysis/type1_type2_astrocyte_proportions_by_region_nontelencephalon.csv",
  row.names = FALSE
)

fit_tel    <- run_lm_subset(tel_df,    predictors_tel,    "Telencephalon")
fit_nontel <- run_lm_subset(nontel_df, predictors_nontel, "Non-telencephalon")

cor_tel    <- run_cor_subset(tel_df,    predictors_tel,    "Telencephalon")
cor_nontel <- run_cor_subset(nontel_df, predictors_nontel, "Non-telencephalon")

cor_combined <- dplyr::bind_rows(cor_tel, cor_nontel)
write.csv(
  cor_combined,
  "data_analysis/spearman_correlations_type1_type2_astrocytes_by_telencephalon.csv",
  row.names = FALSE
)

p_tel    <- make_plot_subset(tel_df,    predictors_tel,    "Telencephalon",     "telencephalon")
p_nontel <- make_plot_subset(nontel_df, predictors_nontel, "Non-telencephalon", "nontelencephalon")

###################################################################
# Combined plots
###################################################################

if (length(predictors_union) > 0) {
  plot_df_all <- analysis_df %>%
    dplyr::select(anatomy_group, rcmr_value, division, all_of(predictors_union)) %>%
    tidyr::pivot_longer(
      cols = all_of(predictors_union),
      names_to = "predictor",
      values_to = "prop"
    ) %>%
    mutate(
      predictor = dplyr::recode(
        predictor,
        "p_Astrocyte_Type1" = "Type 1 astrocytes",
        "p_Astrocyte_Type2" = "Type 2 astrocytes"
      )
    )

  region_scale_all <- tryCatch({
    check_region_palette(analysis_df, region_col = "anatomy_group")
    present_all <- intersect(region_order, unique(as.character(analysis_df$anatomy_group)))
    plot_df_all$anatomy_group <- factor(plot_df_all$anatomy_group, levels = present_all)
    ggplot2::scale_color_manual(values = region_palette, drop = TRUE)
  }, error = function(e) {
    message("Region palette unavailable for combined plot (", e$message,
            "); using default palette.")
    NULL
  })

  p_combined <- ggplot(plot_df_all, aes(x = prop, y = rcmr_value, color = anatomy_group)) +
    geom_point(size = 2.4, alpha = 0.85) +
    geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "steelblue") +
    stat_poly_eq(
      aes(label = paste(after_stat(rr.label), after_stat(p.value.label), sep = "*\", \"*")),
      formula = y ~ x,
      parse = TRUE,
      label.x = "right",
      label.y = "top",
      size = 3,
      color = "black"
    ) +
    facet_grid(division ~ predictor, scales = "free_x") +
    labs(
      title = "Telencephalic (Type 1) vs non-telencephalic (Type 2) astrocyte proportions\nvs regional glucose metabolism (rCMRGlc)",
      subtitle = "Columns: astrocyte type. Rows: telencephalon vs non-telencephalon. Points are brain regions.",
      x = "Mean proportion among nonneuronal cells",
      y = "rCMRGlc (µmol/100 g/min.)",
      color = "Region",
      caption = paste(
        "Siletti et al. (2023) grouped 13 astrocyte clusters into telencephalic (Type 1) and non-telencephalic (Type 2)",
        "types, each with GFAP-low and -high populations. Type 1 includes human cortical populations: WIF1+ gray-matter,",
        "TNC+ white-matter and LMO2+ interlaminar astrocytes. Proportions are each type as a fraction of all nonneuronal",
        "cells. Lines are OLS fits with 95% CI (R^2 and p shown). In our data, telencephalic regions show a positive",
        "relationship between rCMRGlc and the Type 1 / all-nonneuronal proportion; the Type 2 relationship is negative and",
        "borderline (not quite significant). Non-telencephalic regions show no such relationship.",
        sep = "\n"
      )
    ) +
    theme_facet_compact(12)
  if (!is.null(region_scale_all)) p_combined <- p_combined + region_scale_all

  print(p_combined)

  grid_w <- length(unique(plot_df_all$predictor)) * 3.3 + 2.6
  ggsave(
    filename = "figs/s1b/p_type1_type2_astrocytes_by_telencephalon.pdf",
    plot = p_combined, width = grid_w, height = 7.8, units = "in"
  )
  ggsave(
    filename = "figs/s1b/p_type1_type2_astrocytes_by_telencephalon.jpg",
    plot = p_combined, width = grid_w, height = 7.8, units = "in", dpi = 300
  )

  p_overlay <- ggplot(plot_df_all,
                      aes(x = prop, y = rcmr_value, color = division, fill = division)) +
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
    facet_wrap(~ predictor, scales = "free_x") +
    scale_color_manual(values = c("Telencephalon" = "#D7263D",
                                  "Non-telencephalon" = "#1B998B")) +
    scale_fill_manual(values  = c("Telencephalon" = "#D7263D",
                                  "Non-telencephalon" = "#1B998B")) +
    labs(
      title = "Type 1 / Type 2 astrocyte proportions vs rCMRGlc",
      subtitle = "Separate regression lines for telencephalon and non-telencephalon",
      x = "Mean cell type proportion among nonneuronal cells",
      y = "rCMRGlc (µmol/100 g/min.)",
      color = "Division",
      fill  = "Division"
    ) +
    theme_classic(base_size = 13) +
    theme(legend.position = "top")

  print(p_overlay)

  ggsave(
    filename = "figs/s1b/p_type1_type2_astrocytes_overlay_telencephalon.pdf",
    plot = p_overlay, width = 10, height = 7, units = "in"
  )
  ggsave(
    filename = "figs/s1b/p_type1_type2_astrocytes_overlay_telencephalon.jpg",
    plot = p_overlay, width = 10, height = 7, units = "in", dpi = 300
  )
} else {
  cat("\nSkipping combined plots: no predictors with nonzero variance in either subset.\n")
}

###################################################################
# ADDITIONAL ANALYSIS: Type 1 / Type 2 astrocyte composition
# Denominator is Type 1 + Type 2 astrocytes, so proportions sum to 1
# within each anatomy_group. This is analogous to the Type 1 / Type 2
# oligodendrocyte and OPC composition panel in the reference figure.
###################################################################

# Recalculate donor-level Type1/Type2 composition with denominator restricted
# to assigned Type 1 + Type 2 astrocytes, then average across donors.
astrotype_composition_long <- donor_region %>%
  tidyr::crossing(astro_type = factor(astro_type_levels, levels = astro_type_levels)) %>%
  left_join(
    astrotype_counts,
    by = c("anatomy_group", "donor_id", "astro_type")
  ) %>%
  mutate(n_cells = tidyr::replace_na(n_cells, 0L)) %>%
  group_by(anatomy_group, donor_id) %>%
  mutate(
    donor_total_type1_type2_astrocytes = sum(n_cells),
    donor_comp = dplyr::if_else(
      donor_total_type1_type2_astrocytes > 0,
      n_cells / donor_total_type1_type2_astrocytes,
      NA_real_
    )
  ) %>%
  ungroup() %>%
  group_by(anatomy_group, astro_type) %>%
  summarise(
    comp_mean = mean(donor_comp, na.rm = TRUE),
    comp_sd = sd(donor_comp, na.rm = TRUE),
    n_donors = sum(is.finite(donor_comp)),
    mean_assigned_astro_n = mean(donor_total_type1_type2_astrocytes, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    astro_type_label = dplyr::recode(
      as.character(astro_type),
      "Astrocyte_Type1" = "Type 1 astrocytes",
      "Astrocyte_Type2" = "Type 2 astrocytes"
    )
  )

astrotype_composition_wide <- astrotype_composition_long %>%
  select(anatomy_group, astro_type, comp_mean) %>%
  tidyr::pivot_wider(
    names_from = astro_type,
    values_from = comp_mean,
    values_fill = NA_real_,
    names_prefix = "comp_"
  ) %>%
  mutate(
    comp_Type1_minus_Type2 = comp_Astrocyte_Type1 - comp_Astrocyte_Type2,
    comp_Type1_over_Type2 = comp_Astrocyte_Type1 / comp_Astrocyte_Type2
  )

astrotype_composition_analysis <- astrotype_composition_wide %>%
  inner_join(rcmr, by = c("anatomy_group" = "rcmr_term")) %>%
  left_join(
    telencephalon_table %>% select(anatomy_group, is_telencephalon, division),
    by = "anatomy_group"
  ) %>%
  arrange(desc(comp_Astrocyte_Type1))

write.csv(
  astrotype_composition_long,
  "data_analysis/type1_type2_astrocyte_composition_long_by_region.csv",
  row.names = FALSE
)
write.csv(
  astrotype_composition_analysis,
  "data_analysis/type1_type2_astrocyte_composition_with_rcmr_by_region.csv",
  row.names = FALSE
)

cat("\n================ Type 1 / Type 2 astrocyte composition ================\n")
print(
  astrotype_composition_analysis %>%
    select(anatomy_group, division, rcmr_value,
           comp_Astrocyte_Type1, comp_Astrocyte_Type2,
           comp_Type1_minus_Type2, comp_Type1_over_Type2) %>%
    tibble::as_tibble(),
  n = Inf
)

# Correlations between rCMRGlc and Type1/Type2 composition.
composition_cor <- astrotype_composition_analysis %>%
  summarise(
    n = sum(complete.cases(rcmr_value, comp_Astrocyte_Type1)),
    spearman_rho_type1 = suppressWarnings(cor(rcmr_value, comp_Astrocyte_Type1,
                                              method = "spearman", use = "complete.obs")),
    pearson_r_type1 = suppressWarnings(cor(rcmr_value, comp_Astrocyte_Type1,
                                           method = "pearson", use = "complete.obs")),
    spearman_p_type1 = suppressWarnings(cor.test(rcmr_value, comp_Astrocyte_Type1,
                                                 method = "spearman", exact = FALSE)$p.value),
    pearson_p_type1 = suppressWarnings(cor.test(rcmr_value, comp_Astrocyte_Type1,
                                                method = "pearson")$p.value)
  )
write.csv(
  composition_cor,
  "data_analysis/type1_type2_astrocyte_composition_rcmr_correlations.csv",
  row.names = FALSE
)
cat("\nComposition correlations with rCMRGlc:\n")
print(composition_cor)

###################################################################
# Figure A: signed paired bar plot like the oligodendrocyte panel
# Type 1 is plotted upward; Type 2 is plotted downward. Since the
# denominator is Type 1 + Type 2 astrocytes, the absolute values sum
# to 1 for every region.
###################################################################

composition_bar_df <- astrotype_composition_long %>%
  inner_join(
    astrotype_composition_analysis %>% select(anatomy_group, rcmr_value, division),
    by = "anatomy_group"
  ) %>%
  mutate(
    signed_comp = dplyr::if_else(astro_type == "Astrocyte_Type1", comp_mean, -comp_mean),
    anatomy_group = forcats::fct_reorder(anatomy_group, rcmr_value),
    astro_type_label = factor(astro_type_label,
                              levels = c("Type 1 astrocytes", "Type 2 astrocytes"))
  )

p_type1_type2_signed_bar <- ggplot(
  composition_bar_df,
  aes(x = anatomy_group, y = signed_comp, fill = rcmr_value)
) +
  geom_col(width = 0.85, color = NA) +
  geom_hline(yintercept = 0, linewidth = 0.35) +
  facet_grid(astro_type_label ~ ., scales = "free_y") +
  scale_y_continuous(
    labels = function(x) abs(x),
    breaks = c(-1, -0.5, 0, 0.5, 1)
  ) +
  scale_fill_viridis_c(option = "magma", direction = -1) +
  labs(
    title = "Relative composition of Type 1 and Type 2 astrocytes by region",
    subtitle = "Within each region, Type 1 + Type 2 astrocyte fractions sum to 1; regions ordered by rCMRGlc",
    x = NULL,
    y = "Fraction of assigned Type 1 + Type 2 astrocytes",
    fill = "rCMRGlc\n(µmol/100 g/min.)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    strip.background = element_blank(),
    strip.text.y = element_text(angle = 0)
  )

print(p_type1_type2_signed_bar)
ggsave(
  filename = "figs/s1b/p_type1_type2_astrocyte_composition_signed_bar_rcmr_ordered.pdf",
  plot = p_type1_type2_signed_bar, width = 11, height = 5.5, units = "in"
)
ggsave(
  filename = "figs/s1b/p_type1_type2_astrocyte_composition_signed_bar_rcmr_ordered.jpg",
  plot = p_type1_type2_signed_bar, width = 11, height = 5.5, units = "in", dpi = 300
)

###################################################################
# Figure B: stacked composition bars plus rCMRGlc points/line.
# This shows Type1/Type2 composition and metabolism in the same panel
# without using a dual y-axis for the bars themselves.
###################################################################

composition_stacked_df <- astrotype_composition_long %>%
  inner_join(
    astrotype_composition_analysis %>% select(anatomy_group, rcmr_value, division),
    by = "anatomy_group"
  ) %>%
  mutate(
    anatomy_group = forcats::fct_reorder(anatomy_group, rcmr_value),
    astro_type_label = factor(astro_type_label,
                              levels = c("Type 2 astrocytes", "Type 1 astrocytes"))
  )

rcmr_overlay_df <- astrotype_composition_analysis %>%
  mutate(
    anatomy_group = forcats::fct_reorder(anatomy_group, rcmr_value),
    rcmr_scaled_to_fraction = scales::rescale(rcmr_value, to = c(0, 1))
  )

p_type1_type2_stacked_with_rcmr <- ggplot(
  composition_stacked_df,
  aes(x = anatomy_group, y = comp_mean, fill = astro_type_label)
) +
  geom_col(width = 0.85, color = "white", linewidth = 0.2) +
  geom_line(
    data = rcmr_overlay_df,
    aes(x = anatomy_group, y = rcmr_scaled_to_fraction, group = 1),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.6
  ) +
  geom_point(
    data = rcmr_overlay_df,
    aes(x = anatomy_group, y = rcmr_scaled_to_fraction),
    inherit.aes = FALSE,
    color = "black",
    size = 2
  ) +
  scale_fill_manual(
    values = c("Type 1 astrocytes" = "#9E2F7F", "Type 2 astrocytes" = "#F0A202")
  ) +
  labs(
    title = "Type 1 / Type 2 astrocyte composition with rCMRGlc overlay",
    subtitle = "Bars show astrocyte composition; black points/line show rCMRGlc rescaled to 0–1 and regions ordered by rCMRGlc",
    x = NULL,
    y = "Fraction of assigned Type 1 + Type 2 astrocytes",
    fill = "Astrocyte type"
  ) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

print(p_type1_type2_stacked_with_rcmr)
ggsave(
  filename = "figs/s1b/p_type1_type2_astrocyte_composition_stacked_with_rcmr_overlay.pdf",
  plot = p_type1_type2_stacked_with_rcmr, width = 11, height = 5.5, units = "in"
)
ggsave(
  filename = "figs/s1b/p_type1_type2_astrocyte_composition_stacked_with_rcmr_overlay.jpg",
  plot = p_type1_type2_stacked_with_rcmr, width = 11, height = 5.5, units = "in", dpi = 300
)

###################################################################
# Figure C: direct rCMRGlc association plot.
# Because Type 1 and Type 2 fractions sum to 1, plotting Type 1 is
# sufficient; Type 2 is exactly 1 - Type 1.
###################################################################

region_scale_comp <- tryCatch({
  check_region_palette(astrotype_composition_analysis, region_col = "anatomy_group")
  present <- intersect(region_order, unique(as.character(astrotype_composition_analysis$anatomy_group)))
  astrotype_composition_analysis$anatomy_group <- factor(astrotype_composition_analysis$anatomy_group,
                                                         levels = present)
  ggplot2::scale_color_manual(values = region_palette, drop = TRUE)
}, error = function(e) {
  message("Region palette unavailable for composition rCMR plot (", e$message,
          "); using default ggplot palette.")
  NULL
})

p_type1_fraction_vs_rcmr <- ggplot(
  astrotype_composition_analysis,
  aes(x = comp_Astrocyte_Type1, y = rcmr_value, color = anatomy_group)
) +
  geom_point(size = 3, alpha = 0.9) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black") +
  ggpmisc::stat_poly_eq(
    aes(label = paste(after_stat(rr.label), after_stat(p.value.label), sep = "*\", \"*")),
    formula = y ~ x,
    parse = TRUE,
    label.x = "right",
    label.y = "top",
    size = 4,
    color = "black"
  ) +
  labs(
    title = "rCMRGlc versus Type 1 astrocyte fraction",
    subtitle = "Type 1 fraction is computed within assigned Type 1 + Type 2 astrocytes; Type 2 fraction = 1 - Type 1",
    x = "Type 1 astrocyte fraction among Type 1 + Type 2 astrocytes",
    y = "rCMRGlc (µmol/100 g/min.)",
    color = "Region"
  ) +
  theme_classic(base_size = 13)
if (!is.null(region_scale_comp)) p_type1_fraction_vs_rcmr <- p_type1_fraction_vs_rcmr + region_scale_comp

print(p_type1_fraction_vs_rcmr)
ggsave(
  filename = "figs/s1b/p_type1_astrocyte_fraction_vs_rcmr.pdf",
  plot = p_type1_fraction_vs_rcmr, width = 9, height = 6, units = "in"
)
ggsave(
  filename = "figs/s1b/p_type1_astrocyte_fraction_vs_rcmr.jpg",
  plot = p_type1_fraction_vs_rcmr, width = 9, height = 6, units = "in", dpi = 300
)

###################################################################
# Figure D: telencephalon split for the compositional Type 1 fraction.
###################################################################

p_type1_fraction_vs_rcmr_division <- ggplot(
  astrotype_composition_analysis,
  aes(x = comp_Astrocyte_Type1, y = rcmr_value, color = division, fill = division)
) +
  geom_point(size = 3, alpha = 0.9) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15) +
  ggpmisc::stat_poly_eq(
    aes(label = paste(after_stat(rr.label), after_stat(p.value.label), sep = "*\", \"*")),
    formula = y ~ x,
    parse = TRUE,
    label.x = "right",
    size = 3.5,
    show.legend = FALSE
  ) +
  facet_wrap(~ division) +
  scale_color_manual(values = c("Telencephalon" = "#D7263D",
                                "Non-telencephalon" = "#1B998B")) +
  scale_fill_manual(values  = c("Telencephalon" = "#D7263D",
                                "Non-telencephalon" = "#1B998B")) +
  labs(
    title = "rCMRGlc versus Type 1 astrocyte fraction by telencephalon status",
    subtitle = "Type 2 fraction is the complement of Type 1 within assigned astrocyte types",
    x = "Type 1 astrocyte fraction among Type 1 + Type 2 astrocytes",
    y = "rCMRGlc (µmol/100 g/min.)",
    color = "Division",
    fill = "Division"
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top")

print(p_type1_fraction_vs_rcmr_division)
ggsave(
  filename = "figs/s1b/p_type1_astrocyte_fraction_vs_rcmr_by_telencephalon.pdf",
  plot = p_type1_fraction_vs_rcmr_division, width = 10, height = 5.8, units = "in"
)
ggsave(
  filename = "figs/s1b/p_type1_astrocyte_fraction_vs_rcmr_by_telencephalon.jpg",
  plot = p_type1_fraction_vs_rcmr_division, width = 10, height = 5.8, units = "in", dpi = 300
)


####################
####################

###################################################################
# Figure: Type 1 / Type 2 astrocyte proportions with rCMRGlc overlay
# Split into telencephalon and non-telencephalon panels.
# Regions are ordered by increasing rCMRGlc within each division.
###################################################################

if (!requireNamespace("tidytext", quietly = TRUE)) {
  install.packages("tidytext")
}
library(tidytext)

###################################################################
# Build shared region order
###################################################################

region_order_df <- astrotype_composition_analysis %>%
  distinct(anatomy_group, division, rcmr_value) %>%
  arrange(division, rcmr_value) %>%
  group_by(division) %>%
  mutate(region_order = row_number()) %>%
  ungroup()

###################################################################
# Composition data for bars
###################################################################

composition_stacked_df <- astrotype_composition_long %>%
  inner_join(
    astrotype_composition_analysis %>%
      select(anatomy_group, rcmr_value, division),
    by = "anatomy_group"
  ) %>%
  left_join(
    region_order_df,
    by = c("anatomy_group", "division", "rcmr_value")
  ) %>%
  mutate(
    anatomy_group_ordered = tidytext::reorder_within(
      anatomy_group,
      region_order,
      division
    ),
    astro_type_label = factor(
      astro_type_label,
      levels = c("Type 2 astrocytes", "Type 1 astrocytes")
    )
  )

###################################################################
# rCMRGlc data for black overlay line
###################################################################

rcmr_overlay_df <- astrotype_composition_analysis %>%
  distinct(anatomy_group, division, rcmr_value) %>%
  left_join(
    region_order_df,
    by = c("anatomy_group", "division", "rcmr_value")
  ) %>%
  group_by(division) %>%
  mutate(
    rcmr_scaled_to_fraction = scales::rescale(rcmr_value, to = c(0, 1)),
    anatomy_group_ordered = tidytext::reorder_within(
      anatomy_group,
      region_order,
      division
    )
  ) %>%
  arrange(division, region_order) %>%
  ungroup()

###################################################################
# Plot
###################################################################

p_type1_type2_dodged_with_rcmr_split <- ggplot(
  composition_stacked_df,
  aes(
    x = anatomy_group_ordered,
    y = comp_mean,
    fill = astro_type_label
  )
) +
  geom_col(
    position = position_dodge(width = 0.65),
    width = 0.45,
    color = "white",
    linewidth = 0.25
  ) +
  geom_line(
    data = rcmr_overlay_df,
    aes(
      x = anatomy_group_ordered,
      y = rcmr_scaled_to_fraction,
      group = division
    ),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.75
  ) +
  geom_point(
    data = rcmr_overlay_df,
    aes(
      x = anatomy_group_ordered,
      y = rcmr_scaled_to_fraction
    ),
    inherit.aes = FALSE,
    color = "black",
    size = 2
  ) +
  facet_wrap(
    ~ division,
    scales = "free_x",
    ncol = 1
  ) +
  tidytext::scale_x_reordered() +
  scale_fill_manual(
    values = c(
      "Type 1 astrocytes" = "#9E2F7F",
      "Type 2 astrocytes" = "#F0A202"
    )
  ) +
  labs(
    title = "Type 1 / Type 2 astrocyte proportions with rCMRGlc overlay",
    subtitle = "Regions ordered by increasing rCMRGlc within each division; black line shows rCMRGlc rescaled within division",
    x = NULL,
    y = "Fraction of Type 1 + Type 2 astrocytes",
    fill = "Astrocyte type"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

print(p_type1_type2_dodged_with_rcmr_split)

ggsave(
  filename = "figs/s1b/p_type1_type2_astrocyte_dodged_with_rcmr_split.pdf",
  plot = p_type1_type2_dodged_with_rcmr_split,
  width = 12,
  height = 8,
  units = "in"
)

ggsave(
  filename = "figs/s1b/p_type1_type2_astrocyte_dodged_with_rcmr_split.jpg",
  plot = p_type1_type2_dodged_with_rcmr_split,
  width = 12,
  height = 8,
  units = "in",
  dpi = 300
)
