# ============================================================
# Improve anatomy rules: map obs rows to rcmr_term equivalents
# ============================================================

setwd("~/Library/CloudStorage/Dropbox/COLLABORATIVE/Do expensive brain regions increase less in humans/analyses_metabol_rate_structure")

library(tidyverse)
library(stringr)

dir.create("data", showWarnings = FALSE, recursive = TRUE)
dir.create("tables", showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Inputs
# Edit these if your file names differ
# ============================================================

############################
## Load saved obs metadata
############################
obs <- readRDS("data/linnarsson_adult_human_brain_obs_metadata_nonneuronal.rds")

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

#####################################
# Anatomical grouping of dissections
#####################################
# Required columns:
# obs:  dissection, ROIGroupFine, ROIGroupCoarse, ROIGroup
# rcmr: rcmr_term

stopifnot(all(c("dissection", "ROIGroupFine", "ROIGroupCoarse", "ROIGroup") %in% names(obs)))
stopifnot("rcmr_term" %in% names(rcmr))

obs <- obs %>%
  mutate(
    dissection = as.character(dissection),
    ROIGroupFine = as.character(ROIGroupFine),
    ROIGroupCoarse = as.character(ROIGroupCoarse),
    ROIGroup = as.character(ROIGroup)
  )

rcmr_terms <- sort(unique(as.character(rcmr$rcmr_term)))

# ============================================================
# Improved anatomy rules
# anatomy_group should be one of rcmr$rcmr_term
# First match wins.
# dissection_pattern is primary.
# roi_pattern is fallback / supporting evidence.
# ============================================================

anatomy_rules <- tribble(
  ~rule_id, ~anatomy_group, ~dissection_pattern, ~roi_pattern,
  
  # ===================== CORTEX: OCCIPITAL =====================
  "occipital_cortex",
  "Occipital lobe",
  "^Cerebral cortex \\(Cx\\).*(Cuneus|Lingual|occipital|Peristriate|Prostriata|Primary Visual|Visual Cortex|\\bV1C?\\b|\\bV2\\b|\\bA19\\b|Areas 19 and MT)",
  #"Cerebral cortex",
  
  # ===================== CORTEX: TEMPORAL =====================
  "temporal_cortex",
  "Temporal lobe",
  "^Cerebral cortex \\(Cx\\).*(Temporal|temporal|parahippocampal|entorhinal|perirhinal|fusiform|occipitotemporal|Transverse temporal|auditory|\\bA1C\\b|Temporal pole|\\bA38\\b|\\bA35\\b|\\bA36\\b)",
  #"Cerebral cortex",
  
  # ===================== CORTEX: INSULAR =====================
  "insular_cortex",
  "Insular lobe",
  "^Cerebral cortex \\(Cx\\).*(insular|Insular|insula|Short insular|Long insular|agranular insular|granular insular|dysgranular insular)",
  #"Cerebral cortex",
  
  # ===================== CORTEX: PARIETAL =====================
  "parietal_cortex",
  "Parietal lobe",
  "^Cerebral cortex \\(Cx\\).*(Parietal|parietal|Postcentral|postcentral|somatosensory|Supramarginal|supramarginal|Supraparietal|SPL|operculum|A40|A43|A5|A7|S1C?)",
  #"Cerebral cortex",
  
  # ===================== CORTEX: FRONTAL =====================
  "frontal_cortex",
  "Frontal lobe",
  "^Cerebral cortex \\(Cx\\).*(Frontal|frontal|Precentral|precentral|motor cortex|M1C|orbital|Orbitofrontal|gyrus rectus|Gyrus rectus|Subcallosal|subcallosal|Middle frontal|Inferior frontal|Rostral gyrus|MFC|OFC|A13|A14|A24|A25|A32|A44|A45|A46)",
  #"Cerebral cortex",
  
  # ===================== HIPPOCAMPUS =====================
  # Avoid broad 'parahippocampal' here because those are cortical temporal rows.
  "hippocampus",
  "Hippocampus",
  "^(Head|Body|Tail) of [Hh]ippocampus|Hippocampus \\(Hi|Uncal CA|Rostral CA|Caudal Hippocampus|DG-CA4|CA1|CA2|CA3|CA4-DGC|Subicular cortex",
  "Hippocampus",
  
  # ===================== AMYGDALA =====================
  "amygdala",
  "Corpus amygdaloideum",
  "Amygdaloid complex|\\bAMY\\b|Basolateral nuclear group|Central nuclear group|Corticomedial nuclear group|lateral nucleus|basolateral nucleus|basomedial nucleus|anterior cortical nucleus",
  "Amygdaloid complex",
  
  # ===================== BASAL FOREBRAIN =====================
  "basal_forebrain",
  "Basal forebrain",
  "Basal forebrain|septal nuclei|substantia innominata|\\bBF\\b|\\bSEP\\b|\\bSI\\b",
  "Basal forebrain",
  
  # ===================== BASAL GANGLIA =====================
  "caudatum",
  "Caudatum",
  "Caudate|Body of the Caudate|\\bCaB\\b",
  #"Basal nuclei",
  
  "putamen",
  "Putamen",
  "Putamen|\\bPu\\b",
  #"Basal nuclei",
  
  "nucleus_accumbens",
  "Nucleus accumbens",
  "Accumbens|Nucleus Accumbens|\\bNAC\\b",
  #"Basal nuclei",
  
  "pallidum",
  "Pallidum",
  "Pallidum|Globus pallidus|External segment of globus pallidus|Internal segment of globus pallidus|\\bGPe\\b|\\bGPi\\b",
  #"Basal nuclei",
  
  # ===================== THALAMUS =====================
  "nucleus_medial_thalami",
  "Nucleus medial thalami",
  "medial nuclear complex of thalamus|mediodorsal nucleus|\\bMNC\\b|\\bMD\\b|reuniens nucleus",
  "Thalamus",
  
  "corpus_geniculatum_laterale",
  "Corpus geniculatum laterale",
  "lateral geniculate nucleus|lateral geniculate|\\bLG\\b",
  "Thalamus",
  
  "corpus_geniculatum_mediale",
  "Corpus geniculatum mediale",
  "medial geniculate nuclei|medial geniculate|\\bMG\\b",
  "Thalamus",
  
  "nucleus_subthalamicus",
  "Nucleus subthalamicus",
  "Subthalamic nucleus|Subthalamic|\\bSTH\\b",
  "Thalamus",
  
  # ===================== MIDBRAIN =====================
  "substantia_nigra",
  "Substantia nigra",
  "Substantia Nigra|Substantia nigra|\\bSN\\b",
  "Midbrain",
  
  "nucleus_ruber",
  "Nucleus ruber",
  "Red Nucleus|Red nucleus|\\bRN\\b",
  "Midbrain",
  
  "colliculus_inferior",
  "Colliculus inferior",
  "Inferior colliculus|\\bIC\\b",
  "Midbrain",
  
  "colliculus_superior",
  "Colliculus superior",
  "Superior colliculus|\\bSC\\b",
  "Midbrain",
  
  # ===================== CEREBELLUM =====================
  "vermis",
  "Vermis",
  "Vermis|Cerebellar Vermis|\\bCBV\\b",
  #"Cerebellum",
  
  "nucleus_dentatus_cerebelli",
  "Nucleus dentatus cerebelli",
  "Dentate nucleus|Cerebellar deep nuclei|\\bCbDN\\b",
  #"Cerebellum",
  
  "cerebellar_cortex",
  "Cerebellar cortex",
  "Cerebellar cortex|Lateral hemisphere of cerebellum|\\bCBL\\b",
  #"Cerebellum",
  
  # ===================== WHITE MATTER =====================
  "centrum_semiovale",
  "Centrum semiovale",
  "Centrum semiovale",
  NA_character_,
  
  "capsula_interna",
  "Capsula interna",
  "Internal capsule|Capsula interna",
  NA_character_
)

# Keep only rules whose output exists in rcmr
bad_rules <- anatomy_rules %>%
  filter(!anatomy_group %in% rcmr_terms)

if (nrow(bad_rules) > 0) {
  stop(
    "These anatomy_group values are not present in rcmr$rcmr_term:\n",
    paste(unique(bad_rules$anatomy_group), collapse = "\n")
  )
}

write.csv(
  anatomy_rules,
  "tables/anatomy_rules_improved.csv",
  row.names = FALSE
)

# ============================================================
# Apply rules
# Priority:
# 1. dissection_pattern against obs$dissection
# 2. roi_pattern against ROIGroupFine, ROIGroupCoarse, ROIGroup
# First match wins.
# ============================================================

# ============================================================
# Apply rules
# IMPORTANT:
# - anatomy_group assignment is based ONLY on dissection
# - ROIGroupFine / ROIGroupCoarse / ROIGroup are diagnostic only
# - this prevents broad labels like "Basal nuclei" from assigning
#   all basal nuclei rows to the first basal-ganglia rule
# ============================================================

obs_mapped <- obs %>%
  mutate(
    anatomy_group = NA_character_,
    anatomy_rule_id = NA_character_,
    anatomy_match_source = NA_character_,
    anatomy_match_pattern = NA_character_
  )

for (i in seq_len(nrow(anatomy_rules))) {
  
  rule <- anatomy_rules[i, ]
  
  dissection_hit <- !is.na(obs_mapped$dissection) &
    str_detect(
      obs_mapped$dissection,
      regex(rule$dissection_pattern, ignore_case = TRUE)
    )
  
  idx <- is.na(obs_mapped$anatomy_group) & dissection_hit
  
  obs_mapped$anatomy_group[idx] <- rule$anatomy_group
  obs_mapped$anatomy_rule_id[idx] <- rule$rule_id
  obs_mapped$anatomy_match_source[idx] <- "dissection"
  obs_mapped$anatomy_match_pattern[idx] <- rule$dissection_pattern
}

obs_mapped <- obs_mapped %>%
  mutate(
    anatomy_group = if_else(is.na(anatomy_group), "Unmapped", anatomy_group)
  )

# ============================================================
# Build dissection-level equivalent table
# This is the main inspection output:
# for each rcmr_term, which obs dissection terms got mapped?
# ============================================================

rcmr_obs_dissection_equivalents <- obs_mapped %>%
  filter(anatomy_group != "Unmapped") %>%
  count(
    anatomy_group,
    anatomy_match_source,
    anatomy_rule_id,
    ROIGroup,
    ROIGroupCoarse,
    ROIGroupFine,
    dissection,
    name = "n_obs_rows"
  ) %>%
  arrange(anatomy_group, anatomy_match_source, desc(n_obs_rows), dissection)

write.csv(
  rcmr_obs_dissection_equivalents,
  "tables/rcmr_obs_dissection_equivalents.csv",
  row.names = FALSE
)

# ============================================================
# Candidate search table
# This helps detect overlooked terms.
# It searches each rcmr_term literally across the four obs columns.
# ============================================================

candidate_matches <- map_dfr(rcmr_terms, function(term) {
  
  term_regex <- regex(str_replace_all(term, "\\s+", ".*"), ignore_case = TRUE)
  
  obs %>%
    mutate(
      rcmr_term = term,
      hit_dissection = str_detect(coalesce(dissection, ""), term_regex),
      hit_ROIGroupFine = str_detect(coalesce(ROIGroupFine, ""), term_regex),
      hit_ROIGroupCoarse = str_detect(coalesce(ROIGroupCoarse, ""), term_regex),
      hit_ROIGroup = str_detect(coalesce(ROIGroup, ""), term_regex)
    ) %>%
    filter(hit_dissection | hit_ROIGroupFine | hit_ROIGroupCoarse | hit_ROIGroup) %>%
    select(
      rcmr_term,
      hit_dissection,
      hit_ROIGroupFine,
      hit_ROIGroupCoarse,
      hit_ROIGroup,
      ROIGroup,
      ROIGroupCoarse,
      ROIGroupFine,
      dissection
    )
}) %>%
  distinct() %>%
  arrange(rcmr_term, desc(hit_dissection), desc(hit_ROIGroupFine), dissection)

write.csv(
  candidate_matches,
  "tables/rcmr_literal_candidate_matches.csv",
  row.names = FALSE
)

# ============================================================
# Summary outputs
# ============================================================

mapping_summary <- obs_mapped %>%
  count(anatomy_group, anatomy_match_source, name = "n_obs_rows") %>%
  arrange(anatomy_group, anatomy_match_source)

unmapped_summary <- obs_mapped %>%
  filter(anatomy_group == "Unmapped") %>%
  count(ROIGroup, ROIGroupCoarse, ROIGroupFine, dissection, name = "n_obs_rows") %>%
  arrange(desc(n_obs_rows), ROIGroup, ROIGroupCoarse, ROIGroupFine, dissection)

rcmr_terms_without_obs_match <- tibble(rcmr_term = rcmr_terms) %>%
  anti_join(
    obs_mapped %>%
      filter(anatomy_group != "Unmapped") %>%
      distinct(rcmr_term = anatomy_group),
    by = "rcmr_term"
  )

write.csv(
  obs_mapped,
  "data/obs_with_improved_anatomy_group.csv",
  row.names = FALSE
)

write.csv(
  mapping_summary,
  "tables/anatomy_mapping_summary.csv",
  row.names = FALSE
)

write.csv(
  unmapped_summary,
  "tables/unmapped_obs_dissection_summary.csv",
  row.names = FALSE
)

write.csv(
  rcmr_terms_without_obs_match,
  "tables/rcmr_terms_without_obs_match.csv",
  row.names = FALSE
)

# ============================================================
# Console diagnostics
# ============================================================

cat("\nMapped obs rows by anatomy_group:\n")
print(mapping_summary)

cat("\nrcmr terms without obs matches:\n")
print(rcmr_terms_without_obs_match)

cat("\nOutputs written:\n")
cat("  data/obs_with_improved_anatomy_group.csv\n")
cat("  tables/anatomy_rules_improved.csv\n")
cat("  tables/rcmr_obs_dissection_equivalents.csv\n")
cat("  tables/rcmr_literal_candidate_matches.csv\n")
cat("  tables/anatomy_mapping_summary.csv\n")
cat("  tables/unmapped_obs_dissection_summary.csv\n")
cat("  tables/rcmr_terms_without_obs_match.csv\n")