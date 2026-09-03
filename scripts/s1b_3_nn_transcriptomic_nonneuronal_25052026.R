setwd(local({ d <- normalizePath(getwd()); while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d); d }))  # repo root (portable; replaces hardcoded path -- see helpers/project_root.R)

source("helpers/plot_settings.R")

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

# Load shared anatomy readers and the source-independent Heiss preparation.
source("helpers/read_heiss_rates.R")

#####################################
# Anatomical grouping of rois
#####################################

# The one authoritative s1b ROI crosswalk is
# metadata/anatomy/s1b_linnarsson_roi_map.csv.
obs <- attach_s1b_anatomy(obs)

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
# p_* columns are mean participant-level cell type proportions
celltype_proportion_table <- celltype_table_long %>%
  pivot_wider(
    names_from  = supercluster_term,
    values_from = p_cells,
    values_fill = 0,
    names_prefix = "p_"
  )

write.csv(
  celltype_proportion_table,
  "data_analysis/nonneuronal_celltype_p_by_region.csv",
  row.names = FALSE
)

###################################################################
# Bind cell type proportions to rCMRGlc table
###################################################################

analysis_df <- celltype_proportion_table %>%
  join_heiss_rates()

###########
# Analyses
###########

#############################
## AUTOMATIC PREDICTOR SETUP
#############################

exclude <- c("p_Bergmann glia",	#Not a good pan-region predictor. It is essentially cerebellum-specific: nonzero in only 7/23 regions and dominant in cerebellar cortex/vermis.
             "p_Choroid plexus",	#I would not use as a main parenchymal predictor. Likely reflects sampling/proximity/contamination or ventricular-associated tissue rather than regional neural composition.
             "p_Ependymal"	#Same issue: anatomically restricted and ventricle-associated, not a general regional cell-type axis.
)

# Detect all p_ columns and optionally exclude some
all_p <- names(analysis_df)[grepl("^p_", names(analysis_df))]
predictors <- setdiff(all_p, exclude)

# Print exclusions that were not actually among the candidate p_ columns.
predictors <- setdiff(all_p, exclude)
setdiff(exclude, all_p)

cat("Using predictors:\n")
print(predictors)

# Check that predictors match all p_ columns except those explicitly excluded

cat("\nCheck: predictors match all p_ columns except excluded:\n")

if (setequal(predictors, setdiff(all_p, exclude))) {
  cat(" ✓ OK\n")
} else {
  cat(" ⚠ Mismatch detected\n")
}

# Wrap any predictor containing spaces in backticks
predictors_backticked <- ifelse(
  grepl("\\s", predictors),
  paste0("`", predictors, "`"),
  predictors
)

cat("\nPredictors (formula-safe):\n")
print(predictors_backticked)

###########################
## MULTIVARIATE REGRESSION
###########################

fit <- lm(
  reformulate(predictors_backticked, response = "rcmr_value"),
  data = analysis_df
)
summary(fit)

###########################
## CORRELATIONS (LOOP)
###########################

cor_results <- do.call(
  rbind,
  lapply(predictors, function(v) {
    ct <- cor.test(
      analysis_df$rcmr_value,
      analysis_df[[v]],
      method = "spearman",
      exact = FALSE
    )
    data.frame(
      variable = v,
      r        = unname(ct$estimate),
      p_value  = ct$p.value,
      n        = sum(complete.cases(analysis_df$rcmr_value, analysis_df[[v]]))
    )
  })
) %>%
  dplyr::arrange(p_value) %>%
  dplyr::mutate(p_adj_BH = p.adjust(p_value, method = "BH"))

cor_results

###########################
## FACETED SCATTERPLOT
###########################

plot_df <- analysis_df %>%
  dplyr::select(anatomy_group, rcmr_value, all_of(predictors)) %>%
  tidyr::pivot_longer(
    cols = all_of(predictors),
    names_to = "predictor",
    values_to = "proportion"
  )

check_region_palette(analysis_df, region_col = "anatomy_group")
analysis_df <- set_region_order(analysis_df, region_col = "anatomy_group")

fd_nn <- facet_dims(length(predictors))

p_nonneuronal<-ggplot(plot_df, aes(x = proportion, y = rcmr_value, color = anatomy_group)) +
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
  facet_wrap(~ predictor, scales = "free_x", ncol = fd_nn$ncol) +
  #scale_color_manual(values = pal) +
  scale_color_regions() +
  labs(
    x = "Nonneuronal Cell-type mean proportion across donors",
    y = "rCMRGlc (µmol/100 g/min.)",
    color = "Region"
  ) +
  theme_facet_compact(12)

p_nonneuronal

ggsave(
  filename = "figs/s1b/p_nonneuronal.pdf",
  plot = p_nonneuronal,
  width = fd_nn$width,
  height = fd_nn$height,
  units = "in"
)

ggsave(
  filename = "figs/s1b/p_nonneuronal.jpg",
  plot = p_nonneuronal,
  width = fd_nn$width,
  height = fd_nn$height,
  units = "in",
  dpi = 300
)
