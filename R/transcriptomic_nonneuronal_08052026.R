setwd("~/Library/CloudStorage/Dropbox/COLLABORATIVE/Do expensive brain regions increase less in humans/analyses_metabol_rate_structure")

## Install and Load up packages
library(ggplot2)
library(ggpmisc)
library(readxl)
library(tidyverse)

############################
## Load saved obs metadata
############################
obs <- readRDS("data/linnarsson_adult_human_brain_obs_metadata_nonneuronal.rds")

# inspect the anndata object
colnames(obs)
#####################################

#####################################
# Anatomical grouping of dissections
#####################################

# --- NEW: derive cortical lobe for cerebral cortex dissections ---
anatomy_rules <- tribble(
  ~pattern, ~anatomy_group,
  
  # ===================== CORTEX =====================
  
  "^Cerebral cortex \\(Cx\\).*cuneus|lingual|occipital|peristriate|prostriata|visual|V1|V2|A19|MT",
  "Occipital lobe",
  
  "^Cerebral cortex \\(Cx\\).*temporal|parahippocampal|entorhinal|perirhinal|fusiform|occipitotemporal|hippocamp",
  "Temporal lobe",
  
  "^Cerebral cortex \\(Cx\\).*frontal|precentral|orbit|gyrus rectus|subcallosal|A44|A45|A46|A13|A14|A24|A25|A32",
  "Frontal lobe",
  
  "^Cerebral cortex \\(Cx\\).*parietal|postcentral|supramarginal|supraparietal|precuneus|operculum|A40|A43|A5|A7|S1",
  "Parietal lobe",
  
  "^Cerebral cortex \\(Cx\\).*insula|insular",
  "Insular lobe",
  
  "^Cerebral cortex \\(Cx\\).*cingulate|retrosplenial|A23|A29|A30",
  "Cingulate / limbic",
  
  # ===================== LIMBIC =====================
  
  "Amygdaloid|Amygdala|Basolateral nuclear group|Central nuclear group|Corticomedial nuclear group",
  "Corpus amygdaloideum",
  
  "Hippocampal CA|Dentate gyrus",
  "Hippocampus",
  
  # ===================== BASAL GANGLIA =====================
  
  "Caudate",
  "Caudatum",
  
  "Putamen",
  "Putamen",
  
  "Accumbens",
  "Nucleus accumbens",
  
  "Pallidum|Globus pallidus",
  "Pallidum",
  
  # ===================== THALAMUS =====================
  
  "Medial thalam",
  "Nucleus medial thalami",
  
  "Lateral geniculate",
  "Corpus geniculatum laterale",
  
  "Medial geniculate",
  "Corpus geniculatum mediale",
  
  "Subthalamic",
  "Nucleus subthalamicus",
  
  # ===================== MIDBRAIN =====================
  
  "Inferior colliculus",
  "Colliculus inferior",
  
  "Superior colliculus",
  "Colliculus superior",
  
  "Substantia nigra",
  "Substantia nigra",
  
  "Red nucleus",
  "Nucleus ruber",
  
  # ===================== CEREBELLUM =====================
  
  "Cerebellar cortex",
  "Cerebellar cortex",
  
  "Vermis",
  "Vermis",
  
  "Dentate nucleus",
  "Nucleus dentatus cerebelli",
  
  # ===================== WHITE MATTER =====================
  
  "Centrum semiovale",
  "Centrum semiovale",
  
  "Internal capsule",
  "Capsula interna"
)

# ---- Apply rules (first match wins) ----
obs <- obs %>%
  mutate(anatomy_group = NA_character_)

for (i in seq_len(nrow(anatomy_rules))) {
  obs$anatomy_group[
    is.na(obs$anatomy_group) &
      str_detect(obs$dissection,
                 regex(anatomy_rules$pattern[i], ignore_case = TRUE))
  ] <- anatomy_rules$anatomy_group[i]
}

obs$anatomy_group[is.na(obs$anatomy_group)] <- "Unmapped"

# inspect the distribution of anatomy groups
unmapped_df <- obs %>%
  filter(anatomy_group == "Unmapped") %>%
  select(
    dissection,
    ROIGroupFine,
    ROIGroupCoarse,
    roi
  ) %>%
  distinct()

unmapped_df

###################################################################
# Cell type proportion calculations (ALL nonneuronal categories)
###################################################################

# ---- Region × cell-category proportions ----
celltype_table_long <- obs %>%
  filter(anatomy_group != "", anatomy_group != "Unmapped") %>%
  count(anatomy_group, supercluster_term, name = "n_cells") %>%
  group_by(anatomy_group) %>%
  mutate(
    n_all_nonneurons = sum(n_cells),
    p_cells = n_cells / n_all_nonneurons
  ) %>%
  ungroup()

celltype_proportion_table <- celltype_table_long %>%
  # Pivot counts
  select(anatomy_group, supercluster_term, n_cells) %>%
  tidyr::pivot_wider(
    names_from  = supercluster_term,
    values_from = n_cells,
    values_fill = 0,
    names_prefix = "n_"
  ) %>%
  # Join proportions (pivoted separately)
  left_join(
    celltype_table_long %>%
      select(anatomy_group, supercluster_term, p_cells) %>%
      tidyr::pivot_wider(
        names_from  = supercluster_term,
        values_from = p_cells,
        values_fill = 0,
        names_prefix = "p_"
      ),
    by = "anatomy_group"
  ) %>%
  # Add totals and force “all nonneurons proportion = 1”
  left_join(
    celltype_table_long %>%
      distinct(anatomy_group, n_all_nonneurons),
    by = "anatomy_group"
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

######################################################
# Read and Bind rCMRGlc values from Heiss et al. 2004
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

# join to rCMRGlc table
analysis_df <- celltype_proportion_table %>%
  inner_join(rcmr, by = c("anatomy_group" = "rcmr_term"))

###########
# Analyses
###########

#############################
## AUTOMATIC PREDICTOR SETUP
#############################

exclude <- c("p_Bergmann glia", "p_all_nonneurons", "p_sum")

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

p_nonneuronal<-ggplot(plot_df, aes(x = prop, y = rcmr_value, color = anatomy_group)) +
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
    x = "Nonneuronal Cell-type proportion",
    y = "rCMRGlc (µmol/100 g/min.)",
    color = "Region"
  ) +
  theme_classic(base_size = 14)


p_nonneuronal

ggsave(
  filename = "figs/p_nonneuronal.pdf",
  plot = p_nonneuronal,
  width = 10,
  height = 7,
  units = "in"
)

ggsave(
  filename = "figs/p_nonneuronal.jpg",
  plot = p_nonneuronal,
  width = 10,
  height = 7,
  units = "in",
  dpi = 300
)