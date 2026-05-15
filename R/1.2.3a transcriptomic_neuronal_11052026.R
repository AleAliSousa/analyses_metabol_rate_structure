setwd("~/Library/CloudStorage/Dropbox/COLLABORATIVE/Do expensive brain regions increase less in humans/analyses_metabol_rate_structure")

## Install and Load up packages
library(ggplot2)
library(ggpmisc)
library(readxl)
library(tidyverse)

############################
## Load saved obs metadata
############################
obs <- readRDS("data/linnarsson_adult_human_brain_obs_metadata_neuronal.rds")

# inspect the anndata object
colnames(obs)
#####################################

# filter to neurons only
obs_neuron <- obs %>%
  filter(cell_type == "neuron")

######################################################
# Read rCMRGlc values from Heiss et al. 2004
######################################################

# Read Table with rCMRGlc values
heiss_stephan_tbl <- read.csv("data/Heiss_Stephan_data.csv")

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
# Anatomical grouping of dissections
#####################################

# Read the dissection-to-rCMR mapping table directly.
# Expected columns:
#   rcmr_term
#   dissection_terms: one or more exact obs$dissection strings separated by "||"

anatomy_rules <- readr::read_csv(
  "data/rcmr_dissection_terms_by_rcmr_term.csv",
  show_col_types = FALSE
) %>%
  dplyr::transmute(
    anatomy_group = rcmr_term,
    dissection = dissection_terms
  ) %>%
  tidyr::separate_rows(dissection, sep = "\\s*\\|\\|\\s*") %>%
  dplyr::mutate(dissection = stringr::str_squish(dissection)) %>%
  dplyr::distinct(dissection, anatomy_group)

obs_neuron <- obs_neuron %>%
  dplyr::select(-dplyr::any_of("anatomy_group")) %>%
  dplyr::mutate(dissection = stringr::str_squish(dissection)) %>%
  dplyr::left_join(anatomy_rules, by = "dissection") %>%
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

celltype_table_long <- obs_neuron %>%
  filter(anatomy_group != "", anatomy_group != "Unmapped") %>%
  count(anatomy_group, cell_category, name = "n_cells") %>%
  group_by(anatomy_group) %>%
  mutate(
    n_all_neurons = sum(n_cells),
    p_cells = n_cells / n_all_neurons
  ) %>%
  ungroup()

celltype_proportion_table <- celltype_table_long %>%
  # Pivot counts
  select(anatomy_group, cell_category, n_cells) %>%
  tidyr::pivot_wider(
    names_from  = cell_category,
    values_from = n_cells,
    values_fill = 0,
    names_prefix = "n_"
  ) %>%
  # Join proportions
  left_join(
    celltype_table_long %>%
      select(anatomy_group, cell_category, p_cells) %>%
      tidyr::pivot_wider(
        names_from  = cell_category,
        values_from = p_cells,
        values_fill = 0,
        names_prefix = "p_"
      ),
    by = "anatomy_group"
  ) %>%
  # Add total neuron count only (no forced proportions)
  left_join(
    celltype_table_long %>%
      distinct(anatomy_group, n_all_neurons),
    by = "anatomy_group"
  )

write.csv(
  celltype_proportion_table,
  "data/neuronal_celltype_proportions_by_region.csv",
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

exclude <- c("p_all_neurons", "p_all_nonneurons", "p_sum")

# Detect all p_ columns and optionally exclude some
all_p <- names(analysis_df)[grepl("^p_", names(analysis_df))]
predictors <- setdiff(all_p, exclude)

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

# Build high-contrast palette
regions <- sort(unique(plot_df$anatomy_group))
base_cols <- c(
  "#000000", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
)
pal <- grDevices::colorRampPalette(base_cols)(length(regions))
names(pal) <- regions

p_neuronal<-ggplot(plot_df, aes(x = prop, y = rcmr_value, color = anatomy_group)) +
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
  scale_color_manual(values = pal) +
  labs(
    x = "Neuronal Cell-type proportion",
    y = "rCMRGlc (µmol/100 g/min.)",
    color = "Region"
  ) +
  theme_classic(base_size = 14)

p_neuronal

ggsave(
  filename = "figs/p_neuronal.pdf",
  plot = p_neuronal,
  width = 10,
  height = 7,
  units = "in"
)

ggsave(
  filename = "figs/p_neuronal.jpg",
  plot = p_neuronal,
  width = 10,
  height = 7,
  units = "in",
  dpi = 300
)
########################################################
# Are there differences between cortical and noncortical?
########################################################

library(broom)

############################################################
# 1) Add cortical vs non-cortical class at region level
############################################################

cortical_groups <- anatomy_rules %>%
  dplyr::filter(stringr::str_detect(dissection, stringr::fixed("Cerebral cortex (Cx)"))) %>%
  dplyr::distinct(anatomy_group) %>%
  dplyr::pull(anatomy_group)

analysis_df <- analysis_df %>%
  mutate(
    region_class = if_else(
      anatomy_group %in% cortical_groups,
      "Cortical",
      "Non-cortical"
    )
  )

table(analysis_df$region_class)

############################################################
# 2) Multivariate regression separately by group
############################################################

fit_multi_by_group <- function(df, predictors) {
  # guard against overfitting/singular model in small n
  if (nrow(df) <= length(predictors) + 1) {
    return(tibble(
      term = NA_character_,
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      n = nrow(df),
      note = "Too few regions for multivariate model"
    ))
  }
  
  lm(reformulate(predictors, response = "rcmr_value"), data = df) %>%
    broom::tidy() %>%
    mutate(n = nrow(df), note = NA_character_)
}

multi_cortical <- analysis_df %>%
  filter(region_class == "Cortical") %>%
  fit_multi_by_group(predictors) %>%
  mutate(region_class = "Cortical")

multi_noncortical <- analysis_df %>%
  filter(region_class == "Non-cortical") %>%
  fit_multi_by_group(predictors) %>%
  mutate(region_class = "Non-cortical")

multi_results <- bind_rows(multi_cortical, multi_noncortical)

multi_results

############################################################
# 3) Per-predictor regressions separately by group
############################################################

plot_df <- analysis_df %>%
  select(anatomy_group, region_class, rcmr_value, all_of(predictors)) %>%
  pivot_longer(
    cols = all_of(predictors),
    names_to = "predictor",
    values_to = "prop"
  )

uni_results <- plot_df %>%
  group_by(region_class, predictor) %>%
  group_modify(~{
    dat <- .x %>% filter(is.finite(prop), is.finite(rcmr_value))
    
    if (nrow(dat) < 3 || sd(dat$prop) == 0) {
      return(tibble(
        n = nrow(dat),
        intercept = NA_real_,
        slope = NA_real_,
        r2 = NA_real_,
        p_value = NA_real_
      ))
    }
    
    fit <- lm(rcmr_value ~ prop, data = dat)
    sm <- summary(fit)
    
    tibble(
      n = nrow(dat),
      intercept = coef(fit)[["(Intercept)"]],
      slope = coef(fit)[["prop"]],
      r2 = sm$r.squared,
      p_value = sm$coefficients["prop", "Pr(>|t|)"]
    )
  }) %>%
  ungroup() %>%
  group_by(region_class) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  arrange(predictor, region_class)

uni_results
############################################################
# 4) Faceted plot with TWO lines per facet
#    (solid = cortical, dashed = non-cortical)
############################################################

p_2_neuronal<-ggplot(plot_df, aes(x = prop, y = rcmr_value, color = region_class)) +
  geom_point(size = 2.6, alpha = 0.85) +
  geom_smooth(aes(linetype = region_class), method = "lm", se = TRUE) +
  ggpmisc::stat_poly_eq(
    aes(
      group = region_class,
      label = paste(after_stat(rr.label), after_stat(p.value.label), sep = "*\", \"*")
    ),
    formula = y ~ x,
    parse = TRUE,
    label.x = 0.02,   # or "left"
    label.y = 0.98,   # or "top"
    vstep = 0.10,
    show.legend = FALSE
  ) +
  facet_wrap(~ predictor, scales = "free_x") +
  scale_linetype_manual(values = c("Cortical" = "solid", "Non-cortical" = "dashed")) +
  labs(
    x = "Cell-type proportion",
    y = "rCMRGlc (µmol/100 g/min.)",
    color = "Region class",
    linetype = "Region class"
  ) +
  theme_classic(base_size = 14)

p_2_neuronal

ggsave(
  filename = "figs/p_2_neuronal.pdf",
  plot = p_2_neuronal,
  width = 10,
  height = 7,
  units = "in"
)

ggsave(
  filename = "figs/p_2_neuronal.jpg",
  plot = p_2_neuronal,
  width = 10,
  height = 7,
  units = "in",
  dpi = 300
)