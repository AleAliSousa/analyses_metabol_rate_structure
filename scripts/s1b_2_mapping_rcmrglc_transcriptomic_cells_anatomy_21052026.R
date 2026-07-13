setwd(local({ d <- normalizePath(getwd()); while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d); d }))  # repo root (portable; replaces hardcoded path -- see R/project_root.R)

############################################################
# Create CSVs to inspect obs terms that map to rCMRGlc terms
############################################################

## Install and Load up packages
library(ggplot2)
library(ggpmisc)
library(readxl)
library(tidyverse)

## cortical ROIs / Lobes are here: https://github.com/linnarsson-lab/adult-human-brain/blob/main/notebooks/Revision/RevisionFig2.ipynb
## Also See Fig 2 in the original paper: https://www.science.org/doi/10.1126/science.add7046

# NOTE from script that preceded this one:
# Likely interpretations:
#   Human Gpe vs Human GPe → case inconsistency
#   Human A35-36 vs Human A35-A36 → naming inconsistency

# I. Load and normalize anatomy/rCMR input data
############################
## Load saved obs metadata
############################
obs1 <- readRDS("data_intermediate/linnarsson_adult_human_brain_obs_metadata_nonneuronal.rds") %>%
  as_tibble() %>%
  mutate(obs_dataset = "nonneuronal")

obs2 <- readRDS("data_intermediate/linnarsson_adult_human_brain_obs_metadata_neuronal.rds") %>%
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

############################################################
# Generate review CSVs for mapping obs terms to rCMRGlc terms
############################################################

# -----------------------------
# Match-column specification
# -----------------------------
# Define anatomy columns and their matching priority in one place.
# This avoids length mismatches if columns are added or removed.

match_spec <- tribble(
  ~match_column,      ~column_rank, ~column_order,
  "ROIGroup",                1L,            1L,
  "ROIGroupCoarse",          2L,            2L,
  "ROIGroupFine",            3L,            3L,
  "roi",                     4L,            4L,
  "tissue",                  5L,            5L,
  "dissection",              6L,            6L
)

match_cols <- match_spec$match_column
column_priority <- match_spec

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
    across(all_of(match_cols)),
    name = "n_obs_rows"
  ) %>%
  mutate(.obs_key_id = row_number())

# -----------------------------
# Audit cortical contexts missed by "Cerebral cortex (Cx)"-anchored rules
# -----------------------------
# These rows have cortical ROI/tissue annotations, but dissection uses another
# prefix, e.g. Perirhinal cortex or Paleocortex. They may need explicit rules.
cortical_context_not_cx_prefixed <- obs_key %>%
  filter(
    if_any(
      c(ROIGroup, ROIGroupCoarse, ROIGroupFine),
      ~ .x == "Cerebral cortex"
    ) |
      str_detect(tissue, regex("Brodmann|cerebral cortex", ignore_case = TRUE)),
    !str_detect(dissection, regex("^Cerebral cortex \\(Cx\\)", ignore_case = TRUE))
  ) %>%
  arrange(dissection, roi, tissue) %>%
  select(
    .obs_key_id,
    ROIGroup,
    ROIGroupCoarse,
    ROIGroupFine,
    roi,
    tissue,
    dissection,
    n_obs_rows
  )

print(cortical_context_not_cx_prefixed, n = Inf)

# NOTE: Human A35r is BA35/perirhinal cortex that lacks "Cerebral cortex (Cx)" dissection prefix and a Temporal lobe candidate.

# -----------------------------
# Anatomy universes being audited
# -----------------------------
# obs: original row-level metadata.
# obs_key: unique anatomy-context keys:
#   ROIGroup × ROIGroupCoarse × ROIGroupFine × roi × dissection × tissue.
# dissection_universe: unique dissection strings with all observed ROI contexts.

dissection_universe <- obs_key %>%
  filter(!is.na(dissection), dissection != "") %>%
  group_by(dissection) %>%
  summarise(
    n_roi_contexts = n(),
    n_obs_rows_total = sum(n_obs_rows),
    roi_contexts = paste(
      sort(unique(paste(ROIGroup, ROIGroupCoarse, ROIGroupFine, roi, tissue, sep = " / "))),
      collapse = " || "
    ),
    .groups = "drop"
  ) %>%
  arrange(dissection)

# Print dissection terms that occur in more than one ROI context.
# These are the cases where a single dissection string is not uniquely associated
# with one ROIGroup / ROIGroupCoarse / ROIGroupFine / roi / tissue context.
dissection_universe %>%
  filter(n_roi_contexts > 1)

# II Define rules and audit dissection-based mappings

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

# -----------------------------
# Reshape the anatomy fields into long format so rules can match against any selected anatomy column:
# -----------------------------
obs_terms_long <- obs_key %>%
  select(.obs_key_id, all_of(match_cols)) %>%
  pivot_longer(
    cols = all_of(match_cols),
    names_to = "match_column",
    values_to = "match_value"
  ) %>%
  filter(!is.na(match_value), match_value != "")

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

# Exact rCMR term matches
exact_rcmr_rules <- tibble(
  rcmr_term = rcmr_terms,
  match_column = "dissection",
  pattern = rcmr_terms,
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
  "^Cerebral cortex \\(Cx\\).*(cuneus|lingual|occipital gyrus|peristriate|prostriata|area\\s*19)",
  20L, "Cx occipital keywords",
  
  "Temporal lobe",
  "^Cerebral cortex \\(Cx\\).*(temporal|fusiform|parahippocampal|perirhinal|entorhinal)",
  30L, "Cx temporal keywords",
  
  # Non-Cx-prefixed temporal cortex case.
  # Human A35r is BA35/perirhinal cortex, but its dissection label starts with
  # "Perirhinal cortex" rather than "Cerebral cortex (Cx)".
  "Temporal lobe",
  "^Perirhinal cortex \\(area 35\\).*A35r$",
  31L, "A35r perirhinal BA35",
  
  "Frontal lobe",
  "^Cerebral cortex \\(Cx\\).*(frontal gyrus|precentral|orbital|orbitofrontal|rostral gyrus|subcallosal|gyrus rectus|motor)", 
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
  10L, "internal capsule"
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

# -----------------------------
# Dissection-only rule hits
# -----------------------------
# For this audit, only the dissection field determines whether an obs_key row
# matched an rCMR term. Other anatomy fields are retained only as context.

dissection_rule_hits <- obs_terms_long %>%
  filter(
    match_column == "dissection",
    !match_value %in% ambiguous_dissection_terms
  ) %>%
  inner_join(
    rcmr_query_rules %>%
      filter(match_column == "dissection"),
    by = "match_column",
    relationship = "many-to-many"
  ) %>%
  rowwise() %>%
  mutate(
    rule_match = case_when(
      match_mode == "exact_ci" ~ str_to_lower(match_value) == str_to_lower(pattern),
      match_mode == "regex_ci" ~ str_detect(match_value, regex(pattern, ignore_case = TRUE)),
      TRUE ~ FALSE
    )
  ) %>%
  ungroup() %>%
  filter(rule_match) %>%
  select(
    .obs_key_id,
    rcmr_term,
    match_value,
    rule_label,
    rule_rank
  ) %>%
  left_join(
    obs_key %>%
      select(.obs_key_id, all_of(match_cols), n_obs_rows),
    by = ".obs_key_id"
  ) %>%
  arrange(rcmr_term, dissection, roi, tissue, rule_rank)

# Matched dissection cases, one row per obs_key_id × rCMR term.
dissection_hit_cases <- dissection_rule_hits %>%
  group_by(
    rcmr_term,
    .obs_key_id,
    ROIGroup,
    ROIGroupCoarse,
    ROIGroupFine,
    roi,
    tissue,
    dissection,
    n_obs_rows
  ) %>%
  summarise(
    rule_labels = paste(sort(unique(rule_label)), collapse = " || "),
    matched_values = paste(sort(unique(match_value)), collapse = " || "),
    n_rules_hit = n_distinct(rule_label),
    .groups = "drop"
  ) %>%
  arrange(rcmr_term, dissection, roi, tissue)

View(dissection_hit_cases)

# Rows whose dissection matched more than one rCMR term.
multi_rcmr_dissection_hits <- dissection_hit_cases %>%
  group_by(.obs_key_id, ROIGroup, ROIGroupCoarse, ROIGroupFine, roi, tissue, dissection, n_obs_rows) %>%
  filter(n_distinct(rcmr_term) > 1) %>%
  arrange(dissection, roi, tissue, rcmr_term) %>%
  ungroup()

print(multi_rcmr_dissection_hits, n = Inf)
View(multi_rcmr_dissection_hits)

# Master dissection-only matched/unmatched table.
dissection_match_status <- obs_key %>%
  left_join(
    dissection_rule_hits %>%
      group_by(.obs_key_id) %>%
      summarise(
        matched_rcmr_terms = paste(sort(unique(rcmr_term)), collapse = " || "),
        matched_rule_labels = paste(sort(unique(rule_label)), collapse = " || "),
        n_dissection_rules_hit = n_distinct(rule_label),
        .groups = "drop"
      ),
    by = ".obs_key_id"
  ) %>%
  mutate(
    dissection_match_status = case_when(
      is.na(dissection) | dissection == "" ~ "missing_dissection",
      is.na(matched_rcmr_terms) ~ "dissection_unmatched",
      TRUE ~ "dissection_matched"
    )
  )

unmatched_dissection_cases <- dissection_match_status %>%
  filter(dissection_match_status %in% c("dissection_unmatched", "missing_dissection")) %>%
  select(
    .obs_key_id,
    dissection_match_status,
    all_of(match_cols),
    n_obs_rows
  ) %>%
  arrange(dissection_match_status, dissection, roi, tissue)

View(unmatched_dissection_cases)

dissection_match_status %>%
  count(dissection_match_status) %>%
  print(n = Inf)

stopifnot(
  n_distinct(dissection_match_status$.obs_key_id) == nrow(obs_key)
)

# III. Build final rcmr_term → roi relationship, add manual patches, compare to lobe dictionary, save

# Extract the relationship between rcmr_term and roi from the dissection_rule_hits table
rcmr_roi_relationship_raw <- dissection_rule_hits %>%
  select(rcmr_term, roi) %>%
  distinct() %>%
  group_by(rcmr_term) %>%
  summarise(
    rois = paste(sort(unique(roi)), collapse = " || "),
    #n_rois = n_distinct(roi),
    .groups = "drop"
  ) %>%
  arrange(rcmr_term)
View(rcmr_roi_relationship_raw)

######### ------ ------ ######### 
# Check missing-dissection ROIs that occur only in missing_dissection cases
missing_dissection_only_rois_strict <- dissection_match_status %>%
  group_by(roi) %>%
  summarise(
    n_missing_dissection = sum(dissection_match_status == "missing_dissection"),
    n_dissection_unmatched = sum(dissection_match_status == "dissection_unmatched"),
    n_dissection_matched = sum(dissection_match_status == "dissection_matched"),
    n_obs_rows_total = sum(n_obs_rows),
    .groups = "drop"
  ) %>%
  filter(
    n_missing_dissection > 0,
    n_dissection_unmatched == 0,
    n_dissection_matched == 0
  ) %>%
  arrange(roi)

print(missing_dissection_only_rois_strict, n = Inf)
View(missing_dissection_only_rois_strict)

# Manually append ROI spelling/name variants before saving
rcmr_roi_relationship <- rcmr_roi_relationship_raw %>%
  mutate(
    rois = case_when(
      rcmr_term == "Pallidum" ~ paste(rois, "Human Gpe", sep = " || "),
      rcmr_term == "Temporal lobe" ~ paste(rois, "Human A35-36", sep = " || "),
      TRUE ~ rois
    )
  )
######### ------ ------ ######### 

# =========================
# Compare Linnarsson lobe_dict vs RCMR roi relationships
# =========================
# Source:
# Linnarsson Lab GitHub
# https://github.com/linnarsson-lab/adult-human-brain/blob/main/notebooks/Revision/RevisionFig2.ipynb
# Accessed: 2026-05-21
lobe_dict <- list(
  Frontal   = c("A13", "A14", "A32", "A44-A45", "A46", "M1C", "A25", "FI"),
  Parietal  = c("A5-A7", "A40", "A43", "S1C"),
  Limbic    = c("ACC", "A23", "A29-A30", "TH-TL"),
  Occipital = c("A19", "Pro", "V1C", "V2"),
  Temporal  = c("A1C", "A38", "ITG", "Idg", "Ig", "MTG", "STG", "TF", "A35-A36", "LEC", "MEC", "A35r"),
  Paleo     = c("AON", "Pir")
)

# 1) Make the lobe_dict long
lobe_long <- data.frame(
  lobe = rep(names(lobe_dict), lengths(lobe_dict)),
  roi_code = unlist(lobe_dict, use.names = FALSE),
  stringsAsFactors = FALSE
)

# 2) Make rcmr_roi_relationship long
rcmr_long <- do.call(rbind, lapply(seq_len(nrow(rcmr_roi_relationship)), function(i) {
  rois <- trimws(unlist(strsplit(rcmr_roi_relationship$rois[i], "\\|\\|")))
  rois <- trimws(gsub("^\\s*Human\\s+", "", rois, ignore.case = TRUE))
  
  data.frame(
    lobe = sub("\\s+lobe$", "", rcmr_roi_relationship$rcmr_term[i], ignore.case = TRUE),
    roi_code = rois,
    stringsAsFactors = FALSE
  )
}))

# 0) Define which lobes you are comparing (scope)
keep_lobes <- c("Frontal", "Parietal", "Occipital", "Temporal", "Insular")
rcmr_long <- rcmr_long[rcmr_long$lobe %in% keep_lobes, , drop = FALSE]

# 1) Pair-level differences (not "missing/extra")
key_lobe  <- paste(lobe_long$lobe, lobe_long$roi_code, sep = " :: ")
key_rcmr  <- paste(rcmr_long$lobe, rcmr_long$roi_code, sep = " :: ")

# in lobe_dict but not in rcmr
missing_keys <- setdiff(key_lobe, key_rcmr)

# in rcmr but not in lobe_dict
extra_keys   <- setdiff(key_rcmr, key_lobe)

# Turn back into data frames
missing_in_rcmr <- data.frame(
  lobe = sub(" :: .*", "", missing_keys),
  roi_code = sub(".* :: ", "", missing_keys),
  stringsAsFactors = FALSE
)

extra_in_rcmr <- data.frame(
  lobe = sub(" :: .*", "", extra_keys),
  roi_code = sub(".* :: ", "", extra_keys),
  stringsAsFactors = FALSE
)

missing_in_rcmr
extra_in_rcmr

# For the cases that are in both, check if the lobe assignment matches
dict <- setNames(unique(lobe_long[,c("roi_code","lobe")]), c("roi_code","lobe_dict"))
rcmr <- setNames(unique(rcmr_long[,c("roi_code","lobe")]), c("roi_code","lobe_rcmr"))
conflicts <- merge(dict, rcmr, by="roi_code")
conflicts <- conflicts[conflicts$lobe_dict != conflicts$lobe_rcmr, ]
conflicts

# NOTE: Expected lobe-assignment differences between sources:
# FI: Frontal -> Insular;  denotes frontal agranular insular cortex
# Idg: Temporal -> Insular;  denotes dysgranular insular cortex
# Ig: Temporal -> Insular; denotes granular insular cortex
# TH-TL: Limbic -> Temporal; denotes posterior parahippocampal cortex.
#   Heiss used FSL cortical assignment which would follow gyral anatomy, so this is
#   treated as temporal cortex here rather than as a limbic-system functional category.

# Get truly missing codes
missing_true <- missing_in_rcmr[!missing_in_rcmr$roi_code %in% extra_in_rcmr$roi_code, ]
missing_true

# Get the obs rows that correspond to the missing codes
hits <- obs[obs$roi %in% paste("Human", missing_true$roi_code), c("roi","dissection")]
unique(hits)

## NOTE:: these regions that were not included are cingulate cortex and paleocortex 

## Conclusion: we are happy with the mapping of rCMRGlc terms to obs dissection terms, and the coverage of the mapping across the observed anatomy contexts. The few unmatched contexts are either ambiguous (e.g. "Midbrain (M)
# Save rcmr_roi_relationship to on output so it can be used in the final mapping script.
write_csv(rcmr_roi_relationship, "data_intermediate/rcmr_roi_relationship.csv")
