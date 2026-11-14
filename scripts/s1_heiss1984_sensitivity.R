# Alternate published x variable on identical regions: Heiss 1984 versus 2004.
# No imputation and no reconstruction of composites from missing components.
setwd(local({ d <- normalizePath(getwd()); while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d); d }))
library(dplyr)
source("helpers/read_heiss_rates.R")
source("helpers/regional_robustness.R")
rd <- function(p) read.csv(p, check.names = FALSE, stringsAsFactors = FALSE)
blocks <- list(
  s1a = list(path = "data_intermediate/s1a_stereology_comparison.csv",
             vars = c("NeurDensity", "GliaDensity", "GliaNeurDensityRatio")),
  s1b_EI = list(path = "data_analysis/original_broad_EI_ratio_all_regions_with_rcmr.csv",
                vars = c("EI_ratio_no_MSN", "EI_ratio_with_MSN")),
  s1b_astrocytes = list(path = "data_analysis/type1_type2_astrocyte_composition_with_rcmr_by_region.csv",
                        vars = "comp_Astrocyte_Type1"),
  s1b_OPC = list(path = "data_analysis/opc_oligo_lineage_composition_by_region_with_rcmr.csv",
                 vars = "opc_frac_lineage_mean"),
  s1c = list(path = "tables/s1c/synaptic_density_rcmrglc_matched_regions.csv", vars = "synaptic_density_mean"),
  s1d = list(path = "tables/s1d/neuroglial_protein_rcmrglc_matched_regions.csv",
             vars = c("mbp_ug_per_g", "tau_ug_per_g", "mbp_ug_per_mg_protein", "tau_ug_per_mg_protein"))
)
rates <- read_heiss_regions() %>% select(anatomy_id,
  rate_2004 = rcmrglc_mean_umol_100g_min, rate_1984 = heiss_1984_mean_umol_100g_min)
fits <- list(); matched <- list()
for (nm in names(blocks)) {
  b <- blocks[[nm]]
  d <- rd(b$path) %>% left_join(rates, by = "anatomy_id")
  if (nm == "s1a") d <- primary_region_rows(d, "metadata/anatomy/s1a_stereology_region_map.csv")
  # The older broad thalamus and mixed midbrain entries were placed next to
  # narrower 2004 nuclei in the secondary table: they are not matched anatomy.
  d <- d[!d$anatomy_id %in% c("medial_thalamic_nuclei", "substantia_nigra", "red_nucleus"), ]
  for (v in b$vars) {
    z <- d[is.finite(d[[v]]) & is.finite(d$rate_2004) & is.finite(d$rate_1984), ]
    matched[[paste(nm, v)]] <- data.frame(study = nm, variable = v,
      anatomy_id = z$anatomy_id, outcome = z[[v]], rate_2004 = z$rate_2004, rate_1984 = z$rate_1984)
    for (scope in c("matched_regions", "matched_without_white_matter")) {
      a <- if (scope == "matched_regions") z else
        z[!z$anatomy_id %in% c("centrum_semiovale", "internal_capsule"), ]
      for (year in c("2004", "1984")) {
        fits[[length(fits) + 1L]] <- cbind(study = nm, variable = v, scope = scope,
          rate_source = paste("Heiss", year), source_file = b$path,
          regional_fit(a, paste0("rate_", year), v))
      }
    }
  }
}
write_analysis_csv(bind_rows(matched), "tables/s1_synthesis/heiss1984_matched_inputs.csv")
write_analysis_csv(bind_rows(fits), "tables/s1_synthesis/heiss1984_sensitivity_fits.csv")
print(bind_rows(fits)[, c("study", "variable", "scope", "rate_source", "n", "p_ols", "p_spearman")])
