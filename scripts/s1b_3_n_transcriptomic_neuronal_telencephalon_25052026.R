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
obs <- readRDS("data_intermediate/linnarsson_adult_human_brain_obs_metadata_neuronal.rds")

# inspect the anndata object
colnames(obs)
#####################################

# filter to neurons only
obs <- obs %>%
  filter(cell_type == "neuron")

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
  "data_analysis/telencephalon_classification_neuronal.csv",
  row.names = FALSE
)

###################################################################
# Cell type proportion calculations (ALL major neuronal categories)
###################################################################
# ---- Define major neuronal categories ----
# These are mutually exclusive and exhaustive for neurons
obs <- obs %>%
  mutate(
    cell_category = case_when(
      # ---- Excitatory projection neurons ----
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
      # ---- Inhibitory interneurons ----
      supercluster_term %in% c(
        "CGE interneuron",
        "MGE interneuron",
        "LAMP5-LHX6 and Chandelier",
        "Cerebellar inhibitory",
        "Midbrain-derived inhibitory"
      ) ~ "Inhibitory_interneuron",
      # ---- Inhibitory projection neurons (basal ganglia principal cells) ----
      supercluster_term %in% c(
        "Medium spiny neuron",
        "Eccentric medium spiny neuron"
      ) ~ "Inhibitory_projection_MSN",
      # ---- Splatter ----
      supercluster_term == "Splatter" ~ "Splatter",
      # ---- Miscellaneous ----
      supercluster_term == "Miscellaneous" ~ "Miscellaneous",
      # ---- Fallback (should be rare) ----
      TRUE ~ "Other_neuron"
    )
  )
# ---- Sanity check: every neuron assigned ----
table(obs$cell_category)

###################################################################
# Mean participant-level proportions by anatomy_group (over cell_category)
###################################################################

# Keep mapped regions only
obs_celltype <- obs %>%
  filter(anatomy_group != "", anatomy_group != "Unmapped")

# All cell categories observed in the analysis universe
celltypes <- obs_celltype %>%
  distinct(cell_category)

# All observed donor × region combinations
donor_region <- obs_celltype %>%
  distinct(anatomy_group, donor_id)

# Raw counts per donor × region × cell category
celltype_counts <- obs_celltype %>%
  count(anatomy_group, donor_id, cell_category, name = "n_cells")

# Add zero counts for categories absent from a donor-region sample,
# then calculate donor-level proportions
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

  # Region-level mean of donor-level proportions
  group_by(anatomy_group, cell_category) %>%
  summarise(
    p_cells = mean(donor_prop, na.rm = TRUE),
    .groups = "drop"
  )

# Wide table:
# p_cells_* columns are mean participant-level proportions by cell_category
celltype_proportion_table <- celltype_table_long %>%
  pivot_wider(
    names_from  = cell_category,
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

# Explicit exclusions for neuronal cell_category predictors:
# - Splatter / Miscellaneous: QC/annotation residuals, not real cell classes.
# - Inhibitory_projection_MSN: basal-ganglia (striatal) projection neurons;
#   essentially a "basal ganglia vs not" indicator rather than a smooth
#   pan-region cell-type axis.
exclude <- c("p_cells_Splatter",
             "p_cells_Miscellaneous",
             "p_cells_Inhibitory_projection_MSN"
)

all_p_cells <- names(analysis_df)[grepl("^p_cells_", names(analysis_df))]
predictors  <- setdiff(all_p_cells, exclude)

unmatched <- setdiff(exclude, all_p_cells)
if (length(unmatched) > 0) {
  warning("These exclusions did not match any p_cells_ column: ",
          paste(unmatched, collapse = ", "))
}

cat("\nPredictors after EXPLICIT exclusion (Splatter, Miscellaneous, MSN):\n")
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
# Using the cell_category collapse, the two canonical pan-region neuronal
# axes are excitatory projection neurons and inhibitory interneurons.
# With ~11-12 regions per subset and 2 predictors, df_residual = n - 3,
# i.e. ~8-9 -- comfortable for inference.
lm_predictors <- c(
  "p_cells_Excitatory_projection",
  "p_cells_Inhibitory_interneuron"
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
  if (length(predictors) == 0) {
    cat("  Skipping correlations: no predictors survived dynamic exclusion.\n")
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
      title = paste0("Neuronal cell-type proportions vs rCMRGlc \u2014 ", label),
      subtitle = paste0("n = ", nrow(df), " regions; ", length(predictors), " predictors"),
      x = "Neuronal Cell-type proportion",
      y = "rCMRGlc (\u00b5mol/100 g/min.)",
      color = "Region"
    ) +
    theme_facet_compact(12)
  if (!is.null(region_scale)) p <- p + region_scale

  print(p)

  fd <- facet_dims(length(predictors))
  ggsave(
    filename = paste0("figs/s1b/p_neuronal_", file_slug, ".pdf"),
    plot = p, width = fd$width, height = fd$height, units = "in"
  )
  ggsave(
    filename = paste0("figs/s1b/p_neuronal_", file_slug, ".jpg"),
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
  "data_analysis/neuronal_celltype_proportions_by_region_telencephalon.csv",
  row.names = FALSE
)
write.csv(
  nontel_df,
  "data_analysis/neuronal_celltype_proportions_by_region_nontelencephalon.csv",
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
  "data_analysis/spearman_correlations_by_telencephalon_neuronal.csv",
  row.names = FALSE
)

p_tel    <- make_plot_subset(tel_df,    predictors_tel,    "Telencephalon",     "telencephalon")
p_nontel <- make_plot_subset(nontel_df, predictors_nontel, "Non-telencephalon", "nontelencephalon")

###################################################################
# Combined plots (use predictors_common: present and non-zero in both subsets)
###################################################################

if (length(predictors_common) == 0) {
  cat("\nSkipping combined plots: no predictors are non-zero in BOTH divisions.\n")
} else {

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
    title = "Neuronal cell-type proportions vs rCMRGlc by telencephalon vs non-telencephalon",
    subtitle = paste0(length(predictors_common), " predictors shared across both divisions"),
    x = "Neuronal Cell-type proportion",
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
  filename = "figs/s1b/p_neuronal_by_telencephalon.pdf",
  plot = p_combined, width = grid_w, height = 9, units = "in"
)
ggsave(
  filename = "figs/s1b/p_neuronal_by_telencephalon.jpg",
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
    title = "Neuronal cell-type proportion vs rCMRGlc -- Telencephalon vs Non-telencephalon",
    subtitle = paste0(length(predictors_common), " predictors shared across both divisions"),
    x = "Neuronal Cell-type proportion",
    y = "rCMRGlc (\u00b5mol/100 g/min.)",
    color = "Division",
    fill  = "Division"
  ) +
  theme_facet_compact(12) +
  theme(legend.position = "top")

print(p_overlay)

fd_overlay <- facet_dims(length(predictors_common))
ggsave(
  filename = "figs/s1b/p_neuronal_overlay_telencephalon.pdf",
  plot = p_overlay, width = fd_overlay$width, height = fd_overlay$height, units = "in"
)
ggsave(
  filename = "figs/s1b/p_neuronal_overlay_telencephalon.jpg",
  plot = p_overlay, width = fd_overlay$width, height = fd_overlay$height, units = "in", dpi = 300
)

}  # end of `if (length(predictors_common) > 0)` guard

cat("\nDone. Outputs:\n")
cat("  data_analysis/telencephalon_classification_neuronal.csv\n")
cat("  data_analysis/neuronal_celltype_proportions_by_region_telencephalon.csv\n")
cat("  data_analysis/neuronal_celltype_proportions_by_region_nontelencephalon.csv\n")
cat("  data_analysis/spearman_correlations_by_telencephalon_neuronal.csv\n")
cat("  figs/s1b/p_neuronal_telencephalon.{pdf,jpg}\n")
cat("  figs/s1b/p_neuronal_nontelencephalon.{pdf,jpg}\n")
cat("  figs/s1b/p_neuronal_by_telencephalon.{pdf,jpg}\n")
cat("  figs/s1b/p_neuronal_overlay_telencephalon.{pdf,jpg}\n")
