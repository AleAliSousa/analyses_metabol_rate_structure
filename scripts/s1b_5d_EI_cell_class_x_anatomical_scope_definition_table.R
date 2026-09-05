setwd(local({
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d)
  d
}))

# Produce a compact definition-and-results table for the three biologically
# distinct s1b_5 comparisons plus the MSN-included sensitivity definition.

suppressPackageStartupMessages(library(tidyverse))
source("helpers/read_heiss_rates.R")
source("helpers/s1b_EI.R")

lm_path <-
  "data_analysis/lm_EI_original_no_MSN_with_MSN_vs_jorstad_comparison.csv"
spearman_path <-
  "data_analysis/spearman_EI_original_no_MSN_with_MSN_vs_jorstad_comparison.csv"

if (!file.exists(lm_path) || !file.exists(spearman_path)) {
  stop(
    "Run scripts/s1b_5_n_EI_ratio_original_vs_jorstad_overlay_two_MSN_plots_raw_EI_only_16062026.R first.",
    call. = FALSE
  )
}

definitions <- tribble(
  ~comparison_id, ~EI_definition, ~cell_class_scope, ~cell_class_definition, ~anatomical_region_scope, ~regions_included, ~MSN_treatment, ~comparison_type,
  "s1b_5a", "Jorstad-like neocortical E:I", "Neocortical", "All nine Jorstad excitatory subclasses (including L5 ET) divided by all nine Jorstad inhibitory subclasses", "Neocortex", "Frontal, insular, occipital, parietal, and temporal lobes", "Not part of this crosswalk", "Jorstad-aligned primary comparison",
  "s1b_5b", "Siletti cerebral-cortex scope with neocortical-class crosswalk", "Neocortical", "Same Jorstad-like neocortical crosswalk used in s1b_5a", "Cerebral cortex", "Five neocortical lobes plus Corpus amygdaloideum and Hippocampus", "Not part of this crosswalk", "Anatomical-scope extrapolation",
  "s1b_5c_primary", "Siletti cerebral-cortex-class E:I (off-scope classes and MSN excluded)", "Cerebral cortex", "Siletti neocortical, amygdalar, hippocampal, and cluster-specific L5 ET excitatory classes divided by cortical interneuron classes; thalamic, hypothalamic, midbrain, hindbrain, and cerebellar classes excluded", "Cerebral cortex", "Five neocortical lobes plus Corpus amygdaloideum and Hippocampus", "Excluded from inhibitory denominator", "Matched cell-class and anatomical scope",
  "s1b_5c_sensitivity", "Siletti cerebral-cortex-class E:I (off-scope classes excluded; MSN included)", "Cerebral cortex", "Same cerebral-cortex Siletti classes as s1b_5c_primary", "Cerebral cortex", "Five neocortical lobes plus Corpus amygdaloideum and Hippocampus", "Included in inhibitory denominator", "MSN sensitivity analysis"
)

lm_results <- read.csv(lm_path, check.names = FALSE) %>%
  transmute(
    EI_definition,
    n_regions,
    lm_beta = beta,
    lm_r_squared = r_squared,
    lm_p_unadjusted = p_value,
    lm_p_BH_across_six_definitions = p_adj_BH
  )
spearman_results <- read.csv(spearman_path, check.names = FALSE) %>%
  transmute(
    EI_definition,
    spearman_rho = rho,
    spearman_p_unadjusted = p_value,
    spearman_p_BH_across_six_definitions = p_adj_BH
  )

comparison_table <- definitions %>%
  left_join(lm_results, by = "EI_definition") %>%
  left_join(spearman_results, by = "EI_definition") %>%
  mutate(
    lm_significant_unadjusted_0_05 = lm_p_unadjusted < 0.05,
    lm_significant_after_BH_0_05 = lm_p_BH_across_six_definitions < 0.05,
    spearman_significant_unadjusted_0_05 = spearman_p_unadjusted < 0.05,
    spearman_significant_after_BH_0_05 =
      spearman_p_BH_across_six_definitions < 0.05,
    conclusion = case_when(
      comparison_id == "s1b_5b" ~
        "Significant only for the unadjusted linear model; not significant after BH correction or by Spearman correlation",
      TRUE ~ "No statistically significant relationship"
    )
  )

if (anyNA(comparison_table$n_regions)) {
  stop("At least one definition was not found in the s1b_5 statistics.", call. = FALSE)
}

write.csv(
  comparison_table,
  "data_analysis/s1b_5_EI_cell_class_x_anatomical_scope_definition_table.csv",
  row.names = FALSE
)

# Make every inclusion and exclusion visible. The Jorstad-like crosswalk is
# used by s1b_5b; s1b_5c uses the restricted cerebral-cortex terms above.
obs <- readRDS(
  "data_intermediate/linnarsson_adult_human_brain_obs_metadata_neuronal.rds"
) %>%
  filter(cell_type == "neuron") %>%
  mutate(
    roi = stringr::str_squish(as.character(roi)),
    supercluster_term = stringr::str_squish(as.character(supercluster_term)),
    donor_id = stringr::str_squish(as.character(donor_id))
  )
obs <- attach_s1b_anatomy(obs)

jorstad_supercluster_map <- s1b_neocortical_ei_crosswalk() %>%
  filter(mapping_level == "supercluster") %>%
  select(
    supercluster_term,
    s1b_5b_neocortical_crosswalk_EI = EI_class,
    jorstad_equivalent = class_equivalent
  )

supercluster_inclusion <- obs %>%
  count(supercluster_term, name = "all_region_cells") %>%
  left_join(
    obs %>%
      filter(is_cerebral_cortex %in% TRUE) %>%
      count(supercluster_term, name = "cerebral_cortex_region_cells"),
    by = "supercluster_term"
  ) %>%
  mutate(cerebral_cortex_region_cells = replace_na(cerebral_cortex_region_cells, 0L)) %>%
  left_join(jorstad_supercluster_map, by = "supercluster_term") %>%
  mutate(
    s1b_5b_neocortical_crosswalk_EI = case_when(
      !is.na(s1b_5b_neocortical_crosswalk_EI) ~ s1b_5b_neocortical_crosswalk_EI,
      supercluster_term == "Miscellaneous" ~ "Partial E",
      TRUE ~ "Neither"
    ),
    jorstad_equivalent = case_when(
      !is.na(jorstad_equivalent) ~ jorstad_equivalent,
      supercluster_term == "Miscellaneous" ~ "L5 ET clusters 113/114/117/118 only",
      TRUE ~ NA_character_
    ),
    s1b_5c_cerebral_cortex_primary_EI = case_when(
      supercluster_term == "Miscellaneous" ~ "Partial E",
      supercluster_term %in% s1b_siletti_cerebral_cortex_e_terms() ~ "E",
      supercluster_term %in% s1b_siletti_cerebral_cortex_i_terms() ~ "I",
      TRUE ~ "Neither"
    ),
    s1b_5c_MSN_sensitivity_EI = case_when(
      supercluster_term %in% s1b_siletti_broad_i_msn_terms() ~ "I",
      TRUE ~ s1b_5c_cerebral_cortex_primary_EI
    ),
    classification_note = case_when(
      supercluster_term == "Miscellaneous" ~
        "Only L5ET clusters 113/114/117/118 are E; all remaining Miscellaneous cells are excluded",
      supercluster_term %in% c(
        "Thalamic excitatory", "Mammillary body", "Lower rhombic lip",
        "Upper rhombic lip", "Cerebellar inhibitory",
        "Midbrain-derived inhibitory"
      ) ~ "Excluded as outside the cerebral-cortex cell-class scope",
      supercluster_term %in% s1b_siletti_broad_i_msn_terms() ~
        "Excluded from the primary definition; included as I only in the MSN sensitivity",
      supercluster_term == "Splatter" ~ "Excluded from both E and I",
      TRUE ~ NA_character_
    )
  ) %>%
  arrange(supercluster_term) %>%
  as_tibble()

write.csv(
  supercluster_inclusion,
  "data_analysis/s1b_5_siletti_supercluster_EI_inclusion_comparison.csv",
  row.names = FALSE
)

print(comparison_table, width = Inf)
print(supercluster_inclusion, n = Inf, width = Inf)
