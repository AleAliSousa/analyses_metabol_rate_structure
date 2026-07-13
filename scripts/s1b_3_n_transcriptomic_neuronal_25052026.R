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

# filter to neurons only
obs_neuron <- obs %>%
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

obs_neuron <- obs_neuron %>%
  dplyr::select(-dplyr::any_of("anatomy_group")) %>%
  dplyr::mutate(roi = stringr::str_squish(roi)) %>%
  dplyr::left_join(anatomy_rules, by = "roi") %>%
  dplyr::mutate(anatomy_group = dplyr::coalesce(anatomy_group, "Unmapped"))

###################################################################
# Cell type proportion calculations (ALL major neuronal categories)
###################################################################
# ---- Define major neuronal categories ----
# These are mutually exclusive and exhaustive for neurons
obs_neuron <- obs_neuron %>%
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
table(obs_neuron$cell_category)
# ---- Region × cell-category counts and proportions ----

###################################################################
# Cell type proportion calculations
# Mean participant-level proportions by anatomy_group
# (UPDATED: using cell_category)
###################################################################

# Keep mapped regions only
obs_celltype <- obs_neuron %>%
  filter(anatomy_group != "", anatomy_group != "Unmapped")

# All categories observed
celltypes <- obs_celltype %>%
  distinct(cell_category)

# All observed donor × region combinations
donor_region <- obs_celltype %>%
  distinct(anatomy_group, donor_id)

# Raw counts per donor × region × category
celltype_counts <- obs_celltype %>%
  count(anatomy_group, donor_id, cell_category, name = "n_cells")

# Add zero counts, then compute donor-level proportions
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
  
  # Region-level mean
  group_by(anatomy_group, cell_category) %>%
  summarise(
    p_cells = mean(donor_prop, na.rm = TRUE),
    .groups = "drop"
  )

# Wide table
celltype_proportion_table <- celltype_table_long %>%
  pivot_wider(
    names_from  = cell_category,
    values_from = p_cells,
    values_fill = 0,
    names_prefix = "p_"
  )

write.csv(
  celltype_proportion_table,
  "data_analysis/neuronal_category_p_by_region.csv",
  row.names = FALSE
)

# ---- QC: do neuron proportions sum to ~1 by region? ----

qc_prop_sum <- celltype_proportion_table %>%
  mutate(
    p_sum = rowSums(
      dplyr::select(., dplyr::starts_with("p_")),
      na.rm = TRUE
    )
  ) %>%
  select(anatomy_group, p_sum) %>%
  arrange(desc(abs(p_sum - 1)))

qc_prop_sum

###################################################################
# Bind cell type proportions to rCMRGlc table
###################################################################

# join to rCMRGlc table
analysis_df <- celltype_proportion_table %>%
  inner_join(rcmr, by = c("anatomy_group" = "rcmr_term"))

###########
# Analyses
###########

#############################
## AUTOMATIC PREDICTOR SETUP
#############################

exclude <- c("p_Splatter" #Exclude from biology. This is QC/noise-like annotation, not a meaningful cell class.
             ,"p_Miscellaneous"	#Exclude from biological interpretation. It is an annotation residual.
             ,"p_Inhibitory_projection_MSN"	#Use cautiously. It is mostly a basal ganglia/striatal signal. Across all regions, the regression mostly tests “basal ganglia-like vs not,” not a smooth cell-type gradient.
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

###########################
## MULTIVARIATE REGRESSION
###########################

fit <- lm(
  reformulate(predictors, response = "rcmr_value"),
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
    values_to = "prop"
  )

check_region_palette(analysis_df, region_col = "anatomy_group")
analysis_df <- set_region_order(analysis_df, region_col = "anatomy_group")

fd_n <- facet_dims(length(predictors))

p_neuronal<-ggplot(plot_df, aes(x = prop, y = rcmr_value, color = anatomy_group)) +
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
  facet_wrap(~ predictor, scales = "free_x", ncol = fd_n$ncol) +
  #scale_color_manual(values = pal) +
  scale_color_regions() +
  labs(
    x = "Neuronal Cell-type proportion across donors",
    y = "rCMRGlc (µmol/100 g/min.)",
    color = "Region"
  ) +
  theme_facet_compact(12)

p_neuronal

ggsave(
  filename = "figs/s1b/p_neuronal.pdf",
  plot = p_neuronal,
  width = fd_n$width,
  height = fd_n$height,
  units = "in"
)

ggsave(
  filename = "figs/s1b/p_neuronal.jpg",
  plot = p_neuronal,
  width = fd_n$width,
  height = fd_n$height,
  units = "in",
  dpi = 300
)
