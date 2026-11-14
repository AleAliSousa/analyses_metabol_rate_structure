# Matched-region comparisons within research questions and thalamus sensitivities.
# Read outputs after s1c, s1d, s2 and both replication scripts have run.
setwd(local({ d <- normalizePath(getwd()); while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d); d }))
library(dplyr)
source("helpers/read_heiss_rates.R")
source("helpers/regional_robustness.R")
rd <- function(p) read.csv(p, check.names = FALSE, stringsAsFactors = FALSE)
j <- rd("tables/s1c/synaptic_density_region_audit.csv") %>%
  rename(rcmr_value = rcmrglc_mean)
h <- rd("tables/s1c/holmes_rcmrglc_matched_regions.csv")
p <- rd("data_intermediate/s1d_neuroglial_protein_regional_summary.csv") %>%
  filter(is.finite(rcmr_value))
era <- rd("data_intermediate/s2_stress_volume_comparison.csv") %>%
  left_join(derive_study_heiss_rates("metadata/anatomy/s2_stress_volume_region_map.csv"),
            by = "anatomy_id") %>% filter(is.finite(rcmr_value))
beip <- rd("tables/s2/beip_volume_contrasts.csv")
protein_vars <- c("gfap_ug_per_g", "mbp_ug_per_g", "nfl_ug_per_g", "tau_ug_per_g",
                  "uchl1_ug_per_g", "total_protein_mg_per_g", "gfap_to_nfl_ratio",
                  "gfap_ug_per_mg_protein", "mbp_ug_per_mg_protein", "nfl_ug_per_mg_protein",
                  "tau_ug_per_mg_protein", "uchl1_ug_per_mg_protein")
blocks <- list(Johansen = list(d = j, vars = "synaptic_density_mean"),
               Holmes = list(d = h, vars = "mean"),
               Proteins = list(d = p, vars = protein_vars),
               ERA = list(d = era, vars = c("norm_all", "norm_english", "norm_romanian")))
for (nm in unique(beip$contrast)) blocks[[nm]] <- list(
  d = beip[beip$contrast == nm, ], vars = "volume_deficit_pct")
map_paths <- c(Johansen = "s1c_johansen_region_map.csv", Holmes = "s1c_holmes_region_map.csv",
               Proteins = "s1d_sjolin_region_map.csv", ERA = "s2_stress_volume_region_map.csv")
coverage <- bind_rows(lapply(names(blocks), function(nm) {
  d <- blocks[[nm]]$d
  path <- if (nm %in% names(map_paths)) map_paths[[nm]] else "s2_beip_region_map.csv"
  m <- read_study_anatomy_map(file.path("metadata/anatomy", path))
  rel <- m %>% group_by(analysis_anatomy_id) %>%
    summarise(relationship = paste(sort(unique(relationship)), collapse = " | "),
              mapping_note = paste(unique(mapping_note), collapse = " | "), .groups = "drop")
  data.frame(analysis = nm, anatomy_id = d$anatomy_id,
             primary_included = is.finite(d$rcmr_value)) %>%
    left_join(rel, by = c("anatomy_id" = "analysis_anatomy_id"))
}))
write_analysis_csv(coverage, "tables/s1_synthesis/regional_comparison_coverage.csv")
fits <- list(); counts <- list()
fit_set <- function(name, d, scope, comparison) {
  for (v in blocks[[name]]$vars) {
    fits[[length(fits) + 1L]] <<- cbind(analysis = name, variable = v,
      scope = scope, comparison = comparison, regional_fit(d, "rcmr_value", v))
  }
}
for (nm in names(blocks)) {
  d <- blocks[[nm]]$d
  cv <- coverage[coverage$analysis == nm, ]
  strict_ids <- cv$anatomy_id[!grepl("proxy|composite", cv$relationship)]
  fit_set(nm, primary_region_rows(d), "primary", "within_analysis")
}
pairs <- list(c("Johansen", "Proteins"), c("Johansen", "Holmes"))
for (nm in setdiff(names(blocks), names(map_paths))) pairs[[length(pairs) + 1L]] <- c("ERA", nm)
# s2 tests volume response, not tissue composition. Do not enforce an
# intersection between its regions and the s1 tissue-measurement studies.
for (pair in pairs) {
  comparison <- paste(pair, collapse = " / ")
  # Complete cases for EVERY outcome in each comparison; thus every fit in
  # that comparison uses exactly the same region set.
  eligible <- lapply(pair, function(nm) {
    b <- blocks[[nm]]; d <- primary_region_rows(b$d)
    d$anatomy_id[apply(d[, c("rcmr_value", b$vars), drop = FALSE], 1,
                       function(r) all(is.finite(r)))]
  })
  common <- Reduce(intersect, eligible)
  proxy_ids <- coverage$anatomy_id[coverage$analysis %in% pair &
                                   grepl("proxy|composite", coverage$relationship)]
  shared_sets <- list(shared_regions = common)
  for (scope in names(shared_sets)) {
    ids <- shared_sets[[scope]]
    counts[[length(counts) + 1L]] <- data.frame(comparison = comparison, scope = scope,
      n_shared_regions = length(ids), region_ids = paste(sort(ids), collapse = " | "),
      status = if (length(ids) >= 4L) "fits reported" else "too few regions; no inference")
    for (nm in pair) {
      d <- blocks[[nm]]$d
      fit_set(nm, d[d$anatomy_id %in% ids, ], scope, comparison)
    }
  }
}
write_analysis_csv(bind_rows(counts), "tables/s1_synthesis/shared_region_counts.csv")
write_analysis_csv(bind_rows(fits), "tables/s1_synthesis/regional_sensitivity_fits.csv")
print(bind_rows(counts))
