setwd(local({ d <- normalizePath(getwd()); while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d); d }))  # repo root (portable; replaces hardcoded path -- see R/project_root.R)

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
# Cell type proportion calculations
# Mean participant-level proportions by anatomy_group
###################################################################

# Keep mapped regions only
obs_celltype <- obs %>%
  filter(anatomy_group != "", anatomy_group != "Unmapped")

# All cell types observed in the analysis universe
celltypes <- obs_celltype %>%
  distinct(supercluster_term)

# All observed donor × region combinations
donor_region <- obs_celltype %>%
  distinct(anatomy_group, donor_id)

# Raw counts per donor × region × cell type
celltype_counts <- obs_celltype %>%
  count(anatomy_group, donor_id, supercluster_term, name = "n_cells")

# Add zero counts for cell types absent from a donor-region sample,
# then calculate donor-level proportions
celltype_table_long <- donor_region %>%
  crossing(celltypes) %>%
  left_join(
    celltype_counts,
    by = c("anatomy_group", "donor_id", "supercluster_term")
  ) %>%
  mutate(n_cells = replace_na(n_cells, 0)) %>%
  group_by(anatomy_group, donor_id) %>%
  mutate(
    donor_total_cells = sum(n_cells),
    donor_prop = n_cells / donor_total_cells
  ) %>%
  ungroup() %>%

  # Region-level mean of donor-level proportions
  group_by(anatomy_group, supercluster_term) %>%
  summarise(
    p_cells = mean(donor_prop, na.rm = TRUE),
    .groups = "drop"
  )

# Wide table:
# p_cells_* columns are mean participant-level cell type proportions
celltype_proportion_table <- celltype_table_long %>%
  pivot_wider(
    names_from  = supercluster_term,
    values_from = p_cells,
    values_fill = 0,
    names_prefix = "p_cells_"
  )

###################################################################
# Bind cell type proportions to rCMRGlc table AND telencephalon flag
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
    select(division, anatomy_group, rcmr_value) %>%
    arrange(division, anatomy_group),
  n = Inf
)

###########
# Analyses
###########

#############################
## AUTOMATIC PREDICTOR SETUP
#############################

# Explicit exclusions: known anatomically-restricted markers (not pan-regional axes).
# NOTE: actual column names use the "p_cells_" prefix (from the pivot_wider with
# names_sep="_"). The earlier "p_*" entries silently never matched -- now fixed.
exclude <- c("p_cells_Bergmann glia",   # cerebellum-specific
             "p_cells_Choroid plexus",  # ventricle-associated; sampling artefact
             "p_cells_Ependymal"        # ventricle-lining; anatomically restricted
)

all_p_cells <- names(analysis_df)[grepl("^p_cells_", names(analysis_df))]
predictors  <- setdiff(all_p_cells, exclude)

unmatched <- setdiff(exclude, all_p_cells)
if (length(unmatched) > 0) {
  warning("These exclusions did not match any p_cells_ column: ",
          paste(unmatched, collapse = ", "))
}

cat("\nPredictors after EXPLICIT exclusion (Bergmann glia, Choroid plexus, Ependymal):\n")
print(predictors)

# Wrap any predictor containing spaces in backticks (formula-safe)
predictors_backticked <- ifelse(
  grepl("\\s", predictors),
  paste0("`", predictors, "`"),
  predictors
)

# ----------------------------------------------------------------
# DYNAMIC per-subset exclusion (region-specific markers):
# A predictor with any zero proportion within a subset is treated as
# region-specific within that subset and excluded from the analyses/plots
# for that subset. Rationale: zeros indicate the cell type does not exist
# in some regions, which violates the pan-region predictor assumption.
# ----------------------------------------------------------------
predictors_for_subset <- function(df, base_predictors) {
  keep <- vapply(base_predictors, function(v) all(df[[v]] > 0), logical(1))
  base_predictors[keep]
}

# ----------------------------------------------------------------
# CORE predictor set for the multivariate LM
# ----------------------------------------------------------------
# With ~11-12 regions per subset, a 9-10 predictor model is saturated
# (df_residual = 0-2). Restrict to the four canonical glial classes,
# which together account for the vast majority of parenchymal nonneuronal
# cells and are interpretable. This gives df_residual = n - 5, i.e. ~6-7
# for each subset -- enough for meaningful inference.
lm_predictors <- c(
  "p_cells_Astrocyte",
  "p_cells_Oligodendrocyte",
  "p_cells_Oligodendrocyte precursor",
  "p_cells_Microglia"
)
lm_predictors_backticked <- ifelse(
  grepl("\\s", lm_predictors),
  paste0("`", lm_predictors, "`"),
  lm_predictors
)
cat("\nCORE predictors used for the multivariate LM (per subset):\n")
print(lm_predictors)

###################################################################
# Helper functions: per-subset analysis (lm + correlations + plot)
###################################################################

run_lm_subset <- function(df, predictors_bt, label) {
  n   <- nrow(df)
  k   <- length(predictors_bt)
  df_resid <- n - k - 1
  cat("\n========== MULTIVARIATE LM:", label,
      "(n =", n, "regions, ", k, "predictors, df_residual =", df_resid, ") ==========\n")
  if (df_resid < 2) {
    cat("  WARNING: df_residual =", df_resid,
        "-- model would be saturated. Refusing to fit.\n")
    return(invisible(NULL))
  }
  fit <- lm(
    reformulate(predictors_bt, response = "rcmr_value"),
    data = df
  )
  print(summary(fit))
  invisible(fit)
}

run_cor_subset <- function(df, predictors, label) {
  cat("\n========== SPEARMAN CORRELATIONS:", label, "(n =", nrow(df), "regions) ==========\n")
  if (nrow(df) < 4) {
    cat("  Skipping correlations: only", nrow(df), "regions.\n")
    return(invisible(NULL))
  }
  cor_results <- do.call(
    rbind,
    lapply(predictors, function(v) {
      ct <- suppressWarnings(cor.test(
        df$rcmr_value, df[[v]],
        method = "spearman", exact = FALSE
      ))
      data.frame(
        subset   = label,
        variable = v,
        r        = unname(ct$estimate),
        p_value  = ct$p.value,
        n        = sum(complete.cases(df$rcmr_value, df[[v]]))
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
    )

  # Region palette: use the project palette but drop unused levels so the
  # legend only shows regions actually plotted in this subset.
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
    facet_wrap(~ predictor, scales = "free_x", ncol = facet_ncol(length(predictors))) +
    labs(
      title = paste0("Nonneuronal cell-type proportions vs rCMRGlc \u2014 ", label),
      subtitle = paste0("n = ", nrow(df), " regions; ", length(predictors), " predictors"),
      x = "Nonneuronal Cell-type proportion",
      y = "rCMRGlc (\u00b5mol/100 g/min.)",
      color = "Region"
    ) +
    theme_facet_compact(12)
  if (!is.null(region_scale)) p <- p + region_scale

  print(p)

  fd <- facet_dims(length(predictors))
  ggsave(
    filename = paste0("figs/s1b/p_nonneuronal_", file_slug, ".pdf"),
    plot = p, width = fd$width, height = fd$height, units = "in"
  )
  ggsave(
    filename = paste0("figs/s1b/p_nonneuronal_", file_slug, ".jpg"),
    plot = p, width = fd$width, height = fd$height, units = "in", dpi = 300
  )
  invisible(p)
}

###################################################################
# Run analyses for each subset
###################################################################

tel_df    <- analysis_df %>% filter(is_telencephalon)
nontel_df <- analysis_df %>% filter(!is_telencephalon)

# Compute per-subset and shared predictor sets
predictors_tel    <- predictors_for_subset(tel_df,    predictors)
predictors_nontel <- predictors_for_subset(nontel_df, predictors)
predictors_common <- intersect(predictors_tel, predictors_nontel)

cat("\n================ Predictors after dynamic per-subset exclusion ================\n")
cat("Telencephalon predictors (", length(predictors_tel), "):\n", sep = "")
print(predictors_tel)
cat("Non-telencephalon predictors (", length(predictors_nontel), "):\n", sep = "")
print(predictors_nontel)
cat("Common (used for overlay/grid plots) (", length(predictors_common), "):\n", sep = "")
print(predictors_common)
cat("Dropped in Tel due to zeros:    ", paste(setdiff(predictors, predictors_tel),    collapse = ", "), "\n")
cat("Dropped in Non-tel due to zeros:", paste(setdiff(predictors, predictors_nontel), collapse = ", "), "\n")

# Write per-subset proportion tables for the record
write.csv(
  tel_df,
  "data_analysis/nonneuronal_celltype_proportions_by_region_telencephalon.csv",
  row.names = FALSE
)
write.csv(
  nontel_df,
  "data_analysis/nonneuronal_celltype_proportions_by_region_nontelencephalon.csv",
  row.names = FALSE
)

fit_tel    <- run_lm_subset(tel_df,    lm_predictors_backticked, "Telencephalon")
fit_nontel <- run_lm_subset(nontel_df, lm_predictors_backticked, "Non-telencephalon")

cor_tel    <- run_cor_subset(tel_df,    predictors_tel,    "Telencephalon")
cor_nontel <- run_cor_subset(nontel_df, predictors_nontel, "Non-telencephalon")

# Combined correlation table (long)
cor_combined <- dplyr::bind_rows(cor_tel, cor_nontel)
write.csv(
  cor_combined,
  "data_analysis/spearman_correlations_by_telencephalon.csv",
  row.names = FALSE
)

p_tel    <- make_plot_subset(tel_df,    predictors_tel,    "Telencephalon",     "telencephalon")
p_nontel <- make_plot_subset(nontel_df, predictors_nontel, "Non-telencephalon", "nontelencephalon")

###################################################################
# Combined plots (use predictors_common: present and non-zero in both subsets)
###################################################################

# (a) Grid plot: predictor (columns) x division (rows) -- separate panels
plot_df_all <- analysis_df %>%
  dplyr::select(anatomy_group, rcmr_value, division, all_of(predictors_common)) %>%
  tidyr::pivot_longer(
    cols = all_of(predictors_common),
    names_to = "predictor",
    values_to = "prop"
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
  geom_smooth(aes(group = division), method = "lm", se = TRUE, color = "steelblue") +
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
    title = "Nonneuronal cell-type proportions vs rCMRGlc by telencephalon vs non-telencephalon",
    subtitle = paste0(length(predictors_common), " predictors shared across both divisions"),
    x = "Nonneuronal Cell-type proportion",
    y = "rCMRGlc (\u00b5mol/100 g/min.)",
    color = "Region"
  ) +
  theme_facet_compact(11)
if (!is.null(region_scale_all)) p_combined <- p_combined + region_scale_all

print(p_combined)

# 2 division rows x predictors_common columns; size width by columns so panels
# are not squashed, with a fixed taller height for the two rows.
grid_w <- length(predictors_common) * 3.3 + 2.6
ggsave(
  filename = "figs/s1b/p_nonneuronal_by_telencephalon.pdf",
  plot = p_combined, width = grid_w, height = 9, units = "in"
)
ggsave(
  filename = "figs/s1b/p_nonneuronal_by_telencephalon.jpg",
  plot = p_combined, width = grid_w, height = 9, units = "in", dpi = 300
)

# (b) Overlay plot: one panel per predictor, both divisions on same axes,
#     separate regression lines + per-line R^2 / p annotations.
p_overlay <- ggplot(plot_df_all,
                    aes(x = prop, y = rcmr_value, color = division, fill = division)) +
  geom_point(size = 2.4, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15) +
  stat_poly_eq(
    aes(label = paste(after_stat(rr.label), after_stat(p.value.label), sep = "*\", \"*")),
    formula = y ~ x,
    parse = TRUE,
    label.x = "right",
    size = 3,
    show.legend = FALSE
  ) +
  facet_wrap(~ predictor, scales = "free_x", ncol = facet_ncol(length(predictors_common))) +
  scale_color_manual(values = c("Telencephalon" = "#D7263D",
                                "Non-telencephalon" = "#1B998B")) +
  scale_fill_manual(values  = c("Telencephalon" = "#D7263D",
                                "Non-telencephalon" = "#1B998B")) +
  labs(
    title = "Nonneuronal cell-type proportion vs rCMRGlc -- Telencephalon vs Non-telencephalon",
    subtitle = paste0(length(predictors_common), " predictors shared across both divisions"),
    x = "Nonneuronal Cell-type proportion",
    y = "rCMRGlc (\u00b5mol/100 g/min.)",
    color = "Division",
    fill  = "Division"
  ) +
  theme_facet_compact(12) +
  theme(legend.position = "top")

print(p_overlay)

fd_overlay <- facet_dims(length(predictors_common))
ggsave(
  filename = "figs/s1b/p_nonneuronal_overlay_telencephalon.pdf",
  plot = p_overlay, width = fd_overlay$width, height = fd_overlay$height, units = "in"
)
ggsave(
  filename = "figs/s1b/p_nonneuronal_overlay_telencephalon.jpg",
  plot = p_overlay, width = fd_overlay$width, height = fd_overlay$height, units = "in", dpi = 300
)

cat("\nDone. Outputs:\n")
cat("  data_analysis/telencephalon_classification.csv\n")
cat("  data_analysis/nonneuronal_celltype_proportions_by_region_telencephalon.csv\n")
cat("  data_analysis/nonneuronal_celltype_proportions_by_region_nontelencephalon.csv\n")
cat("  data_analysis/spearman_correlations_by_telencephalon.csv\n")
cat("  figs/s1b/p_nonneuronal_telencephalon.{pdf,jpg}\n")
cat("  figs/s1b/p_nonneuronal_nontelencephalon.{pdf,jpg}\n")
cat("  figs/s1b/p_nonneuronal_by_telencephalon.{pdf,jpg}\n")
cat("  figs/s1b/p_nonneuronal_overlay_telencephalon.{pdf,jpg}\n")
