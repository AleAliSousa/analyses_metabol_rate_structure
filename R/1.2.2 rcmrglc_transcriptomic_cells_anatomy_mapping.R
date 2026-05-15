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

match_cols <- c("ROIGroup", "ROIGroupCoarse", "ROIGroupFine", "dissection")

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

obs_norm <- obs %>%
  mutate(
    across(all_of(match_cols), ~ str_squish(as.character(.x))),
    across(all_of(match_cols), ~ na_if(.x, ""))
  )

obs_key <- obs_norm %>%
  count(
    ROIGroup,
    ROIGroupCoarse,
    ROIGroupFine,
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
      sort(unique(paste(ROIGroup, ROIGroupCoarse, ROIGroupFine, sep = " / "))),
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
#    These are intentionally broader than final anatomy rules.
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
manual_query_rules <- bind_rows(
  
  # ===================== CORTEX =====================
  # Only query dissection for lobe calls. ROIGroup == "Cerebral cortex"
  # is not specific enough.
  
  mk_rule(
    "Insular lobe",
    "dissection",
    "^Cerebral cortex \\(Cx\\).*?(insula|insular|Short insular|Long insular|Granular insular|Dysgranular insular|Frontal agranular insular|\\bIg\\b|\\bIdg\\b|\\bFI\\b)",
    rule_rank = 10L,
    rule_label = "cortical insular pattern"
  ),
  
  mk_rule(
    "Occipital lobe",
    "dissection",
    "^Cerebral cortex \\(Cx\\).*?(cuneus|lingual|occipital|peristriate|prostriata|visual|Primary Visual|Lingual gyrus \\(LiG\\)|\\bV1C?\\b|\\bV2\\b|\\bA19\\b|\\bMT\\b|\\bSOG\\b|\\bPro\\b)",
    rule_rank = 20L,
    rule_label = "cortical occipital pattern" ## Avoid bare \\bLiG\\b because case-insensitive matching also hits LIG = Long insular gyri.
  ),
  
  mk_rule(
    "Temporal lobe",
    "dissection",
    "^Cerebral cortex \\(Cx\\).*?(temporal|parahippocampal|entorhinal|perirhinal|fusiform|occipitotemporal|auditory|\\bTTG\\b|\\bSTG\\b|\\bMTG\\b|\\bITG\\b|\\bTP\\b|\\bA1C\\b|\\bA35\\b|\\bA36\\b|\\bA38\\b|\\bLEC\\b|\\bMEC\\b|\\bFuGt\\b|\\bPPH\\b|\\bAPH\\b|\\bPRG\\b)",
    rule_rank = 30L,
    rule_label = "cortical temporal pattern"
  ),
  
  mk_rule(
    "Frontal lobe",
    "dissection",
    "^Cerebral cortex \\(Cx\\)(?!.*\\binsular\\b).*?(frontal|precentral|orbit|gyrus rectus|subcallosal|subgenual|rostral gyrus|motor cortex|\\bMFC\\b|\\bOFC\\b|\\bPrCG\\b|\\bPOrG\\b|\\bSCG\\b|\\bMFG\\b|\\bIFG\\b|\\bRoG\\b|\\bReG\\b|\\bM1C\\b|\\bA44\\b|\\bA45\\b|\\bA46\\b|\\bA13\\b|\\bA14\\b|\\bA24\\b|\\bA25\\b|\\bA32\\b)",
    rule_rank = 40L,
    rule_label = "cortical frontal pattern" #excluding insular cortex
  ), 
    
  mk_rule(
    "Parietal lobe",
    "dissection",
    "^Cerebral cortex \\(Cx\\).*?(parietal|postcentral|supramarginal|supraparietal|somatosensory|operculum|\\bPoCG\\b|\\bSMG\\b|\\bSPL\\b|\\bPaO\\b|\\bS1C\\b|\\bA40\\b|\\bA43\\b|\\bA5\\b|\\bA7\\b)",
    rule_rank = 50L,
    rule_label = "cortical parietal pattern"
  ),
  
  # ===================== LIMBIC =====================
  
  mk_rule(
    "Hippocampus",
    c("ROIGroup", "ROIGroupCoarse", "ROIGroupFine", "dissection"),
    "\\bhippocampus\\b|Head of hippocampus|Body of hippocampus|Tail of Hippocampus|\\bHippocampal\\b|Dentate gyrus|\\bDG\\b|\\bCA1\\b|\\bCA2\\b|\\bCA3\\b|\\bCA4\\b|Subicular cortex",
    rule_rank = 10L,
    rule_label = "hippocampus group or dissection pattern"
  ),
  
  mk_rule(
    "Corpus amygdaloideum",
    c("ROIGroupFine", "dissection"),
    "Amygdaloid complex|Amygdala|Basolateral nuclear group|Central nuclear group|Corticomedial nuclear group|\\bAMY\\b|\\bBLN\\b|\\bCEN\\b|\\bCMN\\b",
    rule_rank = 10L,
    rule_label = "amygdaloid group or dissection pattern"
  ),
  
  mk_rule(
    "Basal forebrain",
    c("ROIGroupFine", "dissection"),
    "^Basal forebrain$|Basal forebrain \\(BF\\)|septal nuclei|substantia innominata|\\bSEP\\b|\\bSI\\b",
    rule_rank = 10L,
    rule_label = "basal forebrain group or dissection pattern"
  ),
  
  # ===================== BASAL GANGLIA =====================
  
  mk_rule(
    "Caudatum",
    "dissection",
    "Caudate|Body of the Caudate|\\bCaB\\b",
    rule_rank = 10L,
    rule_label = "caudate pattern"
  ),
  
  mk_rule(
    "Putamen",
    "dissection",
    "Putamen|\\bPu\\b",
    rule_rank = 10L,
    rule_label = "putamen pattern"
  ),
  
  mk_rule(
    "Nucleus accumbens",
    "dissection",
    "Accumbens|Nucleus Accumbens|\\bNAC\\b",
    rule_rank = 10L,
    rule_label = "nucleus accumbens pattern"
  ),
  
  mk_rule(
    "Pallidum",
    "dissection",
    "Pallidum|Globus pallidus|\\bGP\\b|\\bGPe\\b|\\bGPi\\b",
    rule_rank = 10L,
    rule_label = "pallidum pattern"
  ),
  
  # ===================== THALAMUS =====================
  
  mk_rule(
    "Corpus geniculatum laterale",
    "dissection",
    "lateral geniculate|\\bLG\\b",
    rule_rank = 10L,
    rule_label = "lateral geniculate pattern"
  ),
  
  mk_rule(
    "Corpus geniculatum mediale",
    "dissection",
    "medial geniculate|\\bMG\\b",
    rule_rank = 10L,
    rule_label = "medial geniculate pattern"
  ),
  
  mk_rule(
    "Nucleus subthalamicus",
    "dissection",
    "Subthalamic|\\bSTH\\b",
    rule_rank = 10L,
    rule_label = "subthalamic pattern"
  ),
  
  mk_rule(
    "Nucleus medial thalami",
    "dissection",
    "medial nuclear complex|mediodorsal|reuniens|medioventral|\\bMNC\\b|\\bMD\\b|\\bRe\\b",
    rule_rank = 20L,
    rule_label = "medial thalamic pattern"
  ),
  
  # ===================== MIDBRAIN =====================
  
  mk_rule(
    "Substantia nigra",
    "dissection",
    "Substantia Nigra|\\bSN\\b",
    rule_rank = 10L,
    rule_label = "substantia nigra pattern"
  ),
  
  mk_rule(
    "Nucleus ruber",
    "dissection",
    "^Midbrain \\(RN\\) - Red Nucleus - RN$",
    rule_rank = 20L,
    rule_label = "red nucleus exact dissection"
  ),
  
  mk_rule(
    "Colliculus inferior",
    "dissection",
    "Inferior colliculus|\\bIC\\b",
    rule_rank = 10L,
    rule_label = "inferior colliculus pattern"
  ),
  
  mk_rule(
    "Colliculus superior",
    "dissection",
    "Superior colliculus|\\bSC\\b",
    rule_rank = 10L,
    rule_label = "superior colliculus pattern"
  ),
  
  # ===================== CEREBELLUM =====================
  
  mk_rule(
    "Vermis",
    "dissection",
    "Vermis|Cerebellar Vermis|\\bCBV\\b",
    rule_rank = 10L,
    rule_label = "vermis pattern"
  ),
  
  mk_rule(
    "Nucleus dentatus cerebelli",
    "dissection",
    "Dentate nucleus|Nucleus dentatus|Cerebellar deep nuclei|\\bCbDN\\b",
    rule_rank = 20L,
    rule_label = "dentate or cerebellar deep nuclei pattern"
  ),
  
  mk_rule(
    "Cerebellar cortex",
    "dissection",
    "Cerebellar cortex|Lateral hemisphere of cerebellum|\\bCBL\\b",
    rule_rank = 30L,
    rule_label = "cerebellar cortex or lateral hemisphere pattern"
  ),
  
  # ===================== WHITE MATTER =====================
  
  mk_rule(
    "Centrum semiovale",
    "dissection",
    "Centrum semiovale",
    rule_rank = 10L,
    rule_label = "centrum semiovale pattern"
  ),
  
  mk_rule(
    "Capsula interna",
    "dissection",
    "Internal capsule|Capsula interna",
    rule_rank = 10L,
    rule_label = "internal capsule pattern"
  )
) %>%
  filter(rcmr_term %in% rcmr_terms)

rcmr_query_rules <- bind_rows(exact_rcmr_rules, manual_query_rules) %>%
  left_join(column_priority, by = "match_column") %>%
  arrange(column_rank, column_order, rule_rank, rcmr_term)

write_csv(
  rcmr_query_rules,
  file.path(out_dir, "rcmr_query_rules_used_for_candidate_generation.csv")
)

# -----------------------------
# 4. Apply discovery rules
# -----------------------------

match_one_rule <- function(rcmr_term,
                           match_column,
                           pattern,
                           match_mode,
                           rule_rank,
                           rule_label,
                           column_rank,
                           column_order) {
  vals <- obs_terms_long %>%
    filter(.data$match_column == !!match_column)
  
  if (nrow(vals) == 0) {
    return(tibble())
  }
  
  hit <- switch(
    match_mode,
    exact_ci = str_to_lower(vals$match_value) == str_to_lower(pattern),
    regex_ci = str_detect(vals$match_value, regex(pattern, ignore_case = TRUE)),
    stop("Unknown match_mode: ", match_mode)
  )
  
  vals %>%
    filter(hit) %>%
    mutate(
      rcmr_term = rcmr_term,
      matched_pattern = pattern,
      match_mode = match_mode,
      rule_rank = rule_rank,
      rule_label = rule_label,
      column_rank = column_rank,
      column_order = column_order
    )
}

candidate_matches <- pmap_dfr(
  rcmr_query_rules,
  match_one_rule
) %>%
  left_join(obs_key, by = ".obs_key_id") %>%
  arrange(
    .obs_key_id,
    column_rank,
    column_order,
    rule_rank,
    rcmr_term,
    match_column,
    match_value
  )

# All candidate matches, including multiple candidates per obs row
write_csv(
  candidate_matches,
  file.path(data_dir, "rcmr_obs_candidate_matches_long.csv")
)

# A compact term-level table: all matched values found in the four queried columns
matched_terms_by_column <- candidate_matches %>%
  distinct(
    rcmr_term,
    match_column,
    match_value,
    match_mode,
    rule_label,
    matched_pattern
  ) %>%
  arrange(rcmr_term, match_column, match_value)

write_csv(
  matched_terms_by_column,
  file.path(out_dir, "rcmr_matched_terms_by_column.csv")
)

# -----------------------------
# 5. First-pass one-term-per-row assignment
# -----------------------------

first_assignment <- candidate_matches %>%
  group_by(.obs_key_id) %>%
  arrange(column_rank, column_order, rule_rank, rcmr_term, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    .obs_key_id,
    rcmr_term,
    assigned_column = match_column,
    assigned_value = match_value,
    assigned_match_mode = match_mode,
    assigned_rule_label = rule_label,
    assigned_pattern = matched_pattern,
    assigned_column_rank = column_rank,
    assigned_column_order = column_order,
    assigned_rule_rank = rule_rank
  )

candidate_summary <- candidate_matches %>%
  group_by(.obs_key_id) %>%
  summarise(
    n_candidate_rcmr_terms = n_distinct(rcmr_term),
    candidate_rcmr_terms = paste(sort(unique(rcmr_term)), collapse = " | "),
    candidate_evidence = paste(
      sort(unique(paste0(match_column, "=", match_value, " -> ", rcmr_term))),
      collapse = " || "
    ),
    conflict_flag = n_candidate_rcmr_terms > 1,
    .groups = "drop"
  )

obs_match_audit <- obs_key %>%
  left_join(first_assignment, by = ".obs_key_id") %>%
  left_join(candidate_summary, by = ".obs_key_id") %>%
  left_join(ambiguous_anatomy, by = ".obs_key_id") %>%
  mutate(
    ambiguous_anatomy_flag = replace_na(ambiguous_anatomy_flag, FALSE),
    n_candidate_rcmr_terms = replace_na(n_candidate_rcmr_terms, 0L),
    candidate_rcmr_terms = replace_na(candidate_rcmr_terms, ""),
    candidate_evidence = replace_na(candidate_evidence, ""),
    conflict_flag = replace_na(conflict_flag, FALSE),
    rcmr_term = if_else(ambiguous_anatomy_flag, NA_character_, rcmr_term),
    match_status = case_when(
      ambiguous_anatomy_flag ~ "excluded_ambiguous_anatomy",
      is.na(rcmr_term) ~ "unmatched",
      conflict_flag ~ "assigned_with_multiple_candidates",
      TRUE ~ "assigned_single_candidate"
    )
  ) %>%
  arrange(
    match_status,
    rcmr_term,
    assigned_column_rank,
    assigned_rule_rank,
    ROIGroup,
    ROIGroupCoarse,
    ROIGroupFine,
    dissection
  ) %>%
  select(
    rcmr_term,
    match_status,
    conflict_flag,
    n_candidate_rcmr_terms,
    candidate_rcmr_terms,
    assigned_column,
    assigned_value,
    assigned_match_mode,
    assigned_rule_label,
    assigned_pattern,
    candidate_evidence,
    n_obs_rows,
    ROIGroup,
    ROIGroupCoarse,
    ROIGroupFine,
    dissection
  )

write_csv(
  obs_match_audit,
  file.path(out_dir, "rcmr_obs_match_audit_first_pass.csv")
)

# -----------------------------
# 6. Main review table for exact dissection-based rules
# -----------------------------

dissection_lookup_review <- obs_match_audit %>%
  filter(!is.na(dissection), dissection != "") %>%
  mutate(
    review_keep = !is.na(rcmr_term),
    review_corrected_rcmr_term = rcmr_term,
    review_note = ""
  ) %>%
  select(
    review_keep,
    review_corrected_rcmr_term,
    review_note,
    rcmr_term,
    match_status,
    conflict_flag,
    n_candidate_rcmr_terms,
    candidate_rcmr_terms,
    assigned_column,
    assigned_value,
    assigned_rule_label,
    assigned_pattern,
    n_obs_rows,
    ROIGroup,
    ROIGroupCoarse,
    ROIGroupFine,
    dissection,
    candidate_evidence
  ) %>%
  arrange(
    is.na(rcmr_term),
    rcmr_term,
    conflict_flag,
    dissection
  )

write_csv(
  dissection_lookup_review,
  file.path(out_dir, "rcmr_dissection_lookup_REVIEW.csv")
)

# -----------------------------
# 7. Unmatched outputs at two levels
# -----------------------------
# unmatched_dissection_terms:
# one row per unique obs$dissection value that has no final rCMR assignment
# after rule matching and ambiguous-anatomy exclusions.

unmatched_obs_anatomy_keys <- obs_match_audit %>%
  filter(is.na(rcmr_term)) %>%
  arrange(ROIGroup, ROIGroupCoarse, ROIGroupFine, dissection)

write_csv(
  unmatched_obs_anatomy_keys,
  file.path(out_dir, "rcmr_unmatched_obs_anatomy_keys.csv")
)

unmatched_dissection_terms <- unmatched_obs_anatomy_keys %>%
  filter(!is.na(dissection), dissection != "") %>%
  group_by(dissection) %>%
  summarise(
    n_unmatched_roi_contexts = n(),
    n_obs_rows_total = sum(n_obs_rows),
    roi_contexts = paste(
      sort(unique(paste(ROIGroup, ROIGroupCoarse, ROIGroupFine, sep = " / "))),
      collapse = " || "
    ),
    match_statuses = paste(sort(unique(match_status)), collapse = " | "),
    .groups = "drop"
  ) %>%
  arrange(dissection)

write_csv(
  unmatched_dissection_terms,
  file.path(out_dir, "rcmr_unmatched_unique_dissection_terms.csv")
)

# -----------------------------
# Filter-friendly unmatched table
# -----------------------------

unmatched_filterable <- unmatched_obs_anatomy_keys %>%
  select(
    dissection,
    ROIGroup,
    ROIGroupCoarse,
    ROIGroupFine,
    n_obs_rows,
    match_status,
    conflict_flag,
    candidate_rcmr_terms,
    candidate_evidence
  ) %>%
  pivot_longer(
    cols = c(ROIGroup, ROIGroupCoarse, ROIGroupFine),
    names_to = "roi_column",
    values_to = "roi_value"
  ) %>%
  filter(!is.na(roi_value), roi_value != "") %>%
  arrange(roi_column, roi_value, dissection)

write_csv(
  unmatched_filterable,
  file.path(out_dir, "rcmr_unmatched_filterable_by_roi_term.csv")
)

# -----------------------------
# 8. Summary by rCMR term
# -----------------------------

summary_by_rcmr <- dissection_lookup_review %>%
  filter(!is.na(rcmr_term)) %>%
  group_by(rcmr_term) %>%
  summarise(
    n_dissection_terms = n_distinct(dissection),
    n_obs_rows_total = sum(n_obs_rows),
    dissection_terms = paste(sort(unique(dissection)), collapse = " || "),
    .groups = "drop"
  ) %>%
  arrange(rcmr_term)

write_csv(
  summary_by_rcmr,
  file.path(out_dir, "rcmr_dissection_terms_by_rcmr_term.csv")
)

message("Wrote files to: ", out_dir)
message("Candidate matches: ", nrow(candidate_matches))
message("Assigned unique obs rows: ", sum(!is.na(obs_match_audit$rcmr_term)))
message("Unmatched unique obs rows: ", sum(is.na(obs_match_audit$rcmr_term)))
message("Possible conflicts: ", sum(obs_match_audit$conflict_flag))

# How many rows were matched or not at the unique dissection level?
tibble(
  total = n_distinct(obs$dissection[!is.na(obs$dissection) & obs$dissection != ""]),
  matched = n_distinct(obs_match_audit$dissection[
    !is.na(obs_match_audit$rcmr_term) &
      !is.na(obs_match_audit$dissection) &
      obs_match_audit$dissection != ""
  ]),
  unmatched = n_distinct(obs_match_audit$dissection[
    is.na(obs_match_audit$rcmr_term) &
      !is.na(obs_match_audit$dissection) &
      obs_match_audit$dissection != ""
  ])
)