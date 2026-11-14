# Source-row decisions and independent-region counts under the primary policy.
source("helpers/read_heiss_rates.R")
paths <- list.files("metadata/anatomy", pattern = "_.*map[.]csv$", full.names = TRUE)
paths <- paths[!grepl("heiss_2004_region_map", paths)]
rows <- lapply(paths, function(p) {
 m <- read_study_anatomy_map(p)
 m$study_map <- basename(p)
 m
})
# Different maps have extra fields; use the shared audit columns only.
cols <- c("study_map", "source_region", "analysis_anatomy_id", "heiss_anatomy_id",
          "relationship", "primary_included", "primary_match_note")
audit <- do.call(rbind, lapply(rows, function(d) d[cols]))
dir.create("tables/anatomy", recursive = TRUE, showWarnings = FALSE)
write.csv(audit, "tables/anatomy/primary_mapping_decisions.csv", row.names = FALSE, na = "")
counts <- do.call(rbind, lapply(rows, function(m) {
 data.frame(study_map = m$study_map[1],
  n_previous_mapped_regions = length(unique(m$analysis_anatomy_id[!is.na(m$heiss_anatomy_id)])),
  n_primary_regions = length(unique(m$analysis_anatomy_id[m$primary_included])),
  excluded_source_rows = sum(!m$primary_included),
  note = "Mapping coverage ceiling; actual n depends on source and outcome completeness. BEIP map spans two waves; s5 spans volume sets.")
}))
write.csv(counts, "tables/anatomy/primary_mapping_counts.csv", row.names = FALSE)
write.csv(data.frame(study = "s4a", n_primary_metabolic_regions = 0,
 status = "No compartment-matched endocast metabolic analysis; existing calculations are proxy sensitivity only",
 location = "figs/s4a_proxy_sensitivity"), "tables/anatomy/s4a_primary_status.csv", row.names = FALSE)
print(counts)
# Earlier primary counts are the documented thalamus-only policy from this session.
# These count regional observations, never participants. s1a detailed glial
# counts/subtypes have fewer complete regions than its neuron/glia densities.
region_counts <- data.frame(
 analysis = c("s1a neuron/glia density", "s1b mapped parent groups", "s1c Johansen",
              "s1c Holmes", "s1d proteins", "s2 ERA", "s2 BEIP childhood",
              "s2 BEIP adolescence (each contrast)", "s5a regional scaling", "s4a primary metabolic"),
 earlier_primary_n = c(9, 23, 13, 12, 10, 7, 8, 7, 18, 6),
 strict_primary_n = c(
  counts$n_primary_regions[match(c("s1a_stereology_region_map.csv", "s1b_linnarsson_roi_map.csv",
   "s1c_johansen_region_map.csv", "s1c_holmes_region_map.csv", "s1d_sjolin_region_map.csv",
   "s2_stress_volume_region_map.csv"), counts$study_map)],
  7, 6, counts$n_primary_regions[match("s5_heiss_region_volume_map.csv", counts$study_map)], 0),
 note = "Regional points; before = earlier 12-Nov primary policy. Complete-case n is outcome-specific.")
write.csv(region_counts, "tables/anatomy/analysis_region_counts.csv", row.names = FALSE)
