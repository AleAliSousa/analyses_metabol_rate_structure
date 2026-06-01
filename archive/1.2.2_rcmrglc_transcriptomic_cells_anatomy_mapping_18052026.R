setwd("~/Library/CloudStorage/Dropbox/COLLABORATIVE/Do expensive brain regions increase less in humans/analyses_metabol_rate_structure")

############################################################
# Create CSVs to inspect obs terms that map to rCMRGlc terms
############################################################

## Install and Load up packages
library(ggplot2)
library(ggpmisc)
library(readxl)
library(tidyverse)

data_dir <- "data"
out_dir <- "output/anatomy_rule_audit"

## cortical ROIs / Lobes are here: https://github.com/linnarsson-lab/adult-human-brain/blob/main/notebooks/Revision/RevisionFig2.ipynb
## Also See Fig 2 in the original paper: https://www.science.org/doi/10.1126/science.add7046

# Note from script that preceded this one:
# Likely interpretations:
#   Human Gpe vs Human GPe → case inconsistency
#   Human A35-36 vs Human A35-A36 → naming inconsistency

############################
## Load saved obs metadata
############################
obs1 <- readRDS("data/linnarsson_adult_human_brain_obs_metadata_nonneuronal.rds") %>%
  as_tibble() %>%
  mutate(obs_dataset = "nonneuronal")

obs2 <- readRDS("data/linnarsson_adult_human_brain_obs_metadata_neuronal.rds") %>%
  as_tibble() %>%
  mutate(obs_dataset = "neuronal")

# Check if the columns are identical
setdiff(colnames(obs1), colnames(obs2))
setdiff(colnames(obs2), colnames(obs1))

obs <- bind_rows(obs1, obs2)

# inspect the anndata object
colnames(obs)
#####################################

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

############################################################
# Generate review CSVs for mapping obs terms to rCMRGlc terms
############################################################

match_cols <- c("ROIGroup", "ROIGroupCoarse", "ROIGroupFine", "roi", "dissection")

# -----------------------------
# 1. rCMRGlc terms
# -----------------------------

rcmr_terms <- rcmr %>%
  distinct(rcmr_term) %>%
  mutate(rcmr_term = str_squish(as.character(rcmr_term))) %>%
  filter(!is.na(rcmr_term), rcmr_term != "") %>%
  pull(rcmr_term)

# -----------------------------
# 2. Normalize obs
# -----------------------------
# Standardize anatomy strings for reliable matching
obs_norm <- obs %>%
  mutate(
    across(all_of(match_cols), ~ str_squish(as.character(.x))),
    across(all_of(match_cols), ~ na_if(.x, ""))
  )
# Create a unique key for each observed anatomy context, and count how many obs rows fall into each context.
obs_key <- obs_norm %>%
  count(
    ROIGroup,
    ROIGroupCoarse,
    ROIGroupFine,
    roi,
    dissection,
    name = "n_obs_rows"
  ) %>%
  mutate(.obs_key_id = row_number())

# -----------------------------
# Anatomy universes being audited
# -----------------------------
# obs: original row-level metadata.
# obs_key: unique anatomy-context keys:
#   ROIGroup × ROIGroupCoarse × ROIGroupFine × dissection.
# dissection_universe: unique dissection strings with all observed ROI contexts.

dissection_universe <- obs_key %>%
  filter(!is.na(dissection), dissection != "") %>%
  group_by(dissection) %>%
  summarise(
    n_roi_contexts = n(),
    n_obs_rows_total = sum(n_obs_rows),
    roi_contexts = paste(
      sort(unique(paste(ROIGroup, ROIGroupCoarse, ROIGroupFine, roi, sep = " / "))),
      collapse = " || "
    ),
    .groups = "drop"
  ) %>%
  arrange(dissection)

write_csv(
  dissection_universe,
  file.path(out_dir, "all_obs_unique_dissection_terms.csv")
)

# -----------------------------
# Explicit ambiguous anatomy exclusions
# These are valid obs dissections, but should not receive unattended rCMR assignment.
# -----------------------------

ambiguous_dissection_terms <- c(
  "Midbrain (M) - Substantia Nigra, Red Nucleus, and nearby nuclei - SN-RN"
)

ambiguous_anatomy <- obs_key %>%
  filter(dissection %in% ambiguous_dissection_terms) %>%
  select(.obs_key_id) %>%
  mutate(ambiguous_anatomy_flag = TRUE)

obs_terms_long <- obs_key %>%
  select(.obs_key_id, all_of(match_cols)) %>%
  pivot_longer(
    cols = all_of(match_cols),
    names_to = "match_column",
    values_to = "match_value"
  ) %>%
  filter(!is.na(match_value), match_value != "")

column_priority <- tibble(
  match_column = match_cols,
  column_rank = c(1L, 1L, 2L, 3L),
  column_order = c(1L, 2L, 3L, 4L)
)

# -----------------------------
# 3. Query rules for discovery
#    Check the output
#    If they are not inclusive/exclusive enough update them later
# -----------------------------

mk_rule <- function(rcmr_term,
                    match_column,
                    pattern,
                    match_mode = "regex_ci",
                    rule_rank = 10L,
                    rule_label = NA_character_) {
  if (is.na(rule_label)) rule_label <- pattern
  
  tibble(
    rcmr_term = rcmr_term,
    match_column = match_column,
    pattern = pattern,
    match_mode = match_mode,
    rule_rank = rule_rank,
    rule_label = rule_label
  )
}

# Exact rCMR term matches in any queried column
exact_rcmr_rules <- tidyr::crossing(
  rcmr_term = rcmr_terms,
  match_column = match_cols
) %>%
  transmute(
    rcmr_term,
    match_column,
    pattern = rcmr_term,
    match_mode = "exact_ci",
    rule_rank = 0L,
    rule_label = "exact rCMR term"
  )

# Broader discovery rules.
# These are for generating candidates, not for final unattended mapping.
# NOTE: match_mode="regex_ci" already makes matching case-insensitive.

manual_query_rules <- tibble::tribble(
  ~rcmr_term,                         ~pattern,                                                                                                                                      ~rule_rank, ~rule_label,
  
  # ===================== CORTEX (Cx) =====================
  "Insular lobe",
  "^Cerebral cortex \\(Cx\\).*(insula|insular|short insular|long insular|granular insular|dysgranular insular|frontal agranular insular)",
  10L, "Cx insular keywords",
  
  "Occipital lobe",
  "^Cerebral cortex \\(Cx\\).*(cuneus|lingual|occip|peristriate|prostriata|area\\s*19)",
  20L, "Cx occipital keywords",
  
  "Temporal lobe",
  "^Cerebral cortex \\(Cx\\).*(temporal|fusiform|parahippocampal|perirhinal|entorhinal)",
  30L, "Cx temporal keywords",
  
  "Frontal lobe",
  "^Cerebral cortex \\(Cx\\).*(frontal|precentral|orbital|orbitofrontal|rostral gyrus|subcallosal|gyrus rectus|motor)",
  40L, "Cx frontal keywords",
  
  "Parietal lobe",
  "^Cerebral cortex \\(Cx\\).*(parietal|postcentral|somatosensory|supraparietal|supramarginal|parietal operculum)",
  50L, "Cx parietal keywords",
  
  
  # ===================== LIMBIC =====================
  "Hippocampus",
  "(head of hippocampus|body of hippocampus|tail of hippocampus|hippocampus|dentate gyrus|subicular cortex)",
  10L, "hippocampus words only",
  
  "Corpus amygdaloideum",
  "(amygdaloid complex|basolateral nuclear group|central nuclear group|corticomedial nuclear group|extended amygdala|bed nucleus of stria terminalis)",
  10L, "amygdala/extended amygdala words",
  
  "Basal forebrain",
  "(basal forebrain|septal nuclei|substantia innominata)",
  10L, "BF/SEP/SI words",
  
  
  # ===================== BASAL GANGLIA =====================
  "Caudatum",
  "(caudate|body of the caudate)",
  10L, "caudate words",
  
  "Putamen",
  "(putamen)",
  10L, "putamen words",
  
  "Nucleus accumbens",
  "(nucleus accumbens|accumbens)",
  10L, "accumbens words",
  
  "Pallidum",
  "(globus pallidus)",
  10L, "pallidum words",
  
  
  # ===================== THALAMUS / SUBTHALAMUS =====================
  "Corpus geniculatum laterale",
  "(lateral geniculate nucleus|lateral geniculate)",
  10L, "LG words",
  
  "Corpus geniculatum mediale",
  "(medial geniculate nuclei|medial geniculate)",
  10L, "MG words",
  
  "Nucleus subthalamicus",
  "(subthalamic nucleus|subthalamic)",
  10L, "subthalamic words",
  
  "Nucleus medial thalami",
  "(medial nuclear complex|mediodorsal|reuniens|medioventral)",
  20L, "MNC/MD/Re words",
  
  
  # ===================== MIDBRAIN =====================
  "Substantia nigra",
  "(substantia\\s+nigra)",
  10L, "SN words only",
  
  "Nucleus ruber",
  "(red\\s+nucleus|nucleus\\s+ruber)",
  20L, "RN words only",
  
  "Colliculus inferior",
  "(inferior\\s+colliculus)",
  10L, "IC words only",
  
  "Colliculus superior",
  "(superior\\s+colliculus)",
  10L, "SC words only",
  
  
  # ===================== CEREBELLUM =====================
  "Vermis",
  "(vermis|cerebellar vermis)",
  10L, "vermis words",
  
  "Nucleus dentatus cerebelli",
  "(dentate|dentate nucleus|nucleus dentatus|cerebellar deep nuclei)",
  20L, "dentate words",
  
  "Cerebellar cortex",
  "(cerebellar cortex|lateral hemisphere of cerebellum)",
  30L, "cerebellar cortex words",
  
  
  # ===================== WHITE MATTER =====================
  "Centrum semiovale",
  "(centrum semiovale)",
  10L, "centrum semiovale",
  
  "Capsula interna",
  "(internal capsule|capsula interna)",
  10L, "internal capsule",
) %>%
  dplyr::mutate(
    match_column = "dissection",
    match_mode   = "regex_ci"
  ) %>%
  dplyr::select(rcmr_term, match_column, pattern, match_mode, rule_rank, rule_label) %>%
  dplyr::filter(rcmr_term %in% rcmr_terms)

rcmr_query_rules <- bind_rows(exact_rcmr_rules, manual_query_rules) %>%
  left_join(column_priority, by = "match_column") %>%
  arrange(column_rank, column_order, rule_rank, rcmr_term)


# 2) How to see which obs_key rows link to which rules
# Keep your current obs_key structure, but match only on dissection
# A) Make a “rule hits” table (all matches, not just the winner)
rule_hits <- obs_terms_long %>%
  inner_join(rcmr_query_rules, by = "match_column") %>%
  rowwise() %>%
  mutate(
    rule_match = case_when(
      match_mode == "exact_ci" ~ str_to_lower(match_value) == str_to_lower(pattern),
      match_mode == "regex_ci" ~ str_detect(match_value, regex(pattern, ignore_case = TRUE)),
      TRUE ~ NA
    )
  ) %>%
  ungroup() %>%
  filter(rule_match) %>%
  select(.obs_key_id, rcmr_term, match_column, match_value, rule_label, column_rank, column_order, rule_rank)
