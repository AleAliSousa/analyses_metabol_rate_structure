# Scientific policy checks: rejected anatomy cannot silently recover a rate.
source("helpers/read_heiss_rates.R")
source("helpers/regional_robustness.R")
paths <- list.files("metadata/anatomy", pattern = "_.*map[.]csv$", full.names = TRUE)
paths <- paths[!grepl("heiss_2004_region_map", paths)]
for (p in paths) {
 m <- read_study_anatomy_map(p)
 z <- derive_study_heiss_rates(m)
 stopifnot(all(is.na(z$rcmr_value[!z$primary_included])),
           all(is.finite(z$rcmr_value[z$primary_included])))
}
m <- read_study_anatomy_map("metadata/anatomy/s2_stress_volume_region_map.csv")
z <- derive_study_heiss_rates(m)
old <- derive_study_heiss_rates(m, policy = "legacy_proxy")
stopifnot(is.na(z$rcmr_value[z$anatomy_id == "thalamus"]),
          old$rcmr_value[old$anatomy_id == "thalamus"] == 26.6)
rois <- read_s1b_roi_map()
stopifnot(!any(c("Human BNST", "Human MD", "Human MD-Re", "Human CbDN",
                "Human SEP", "Human SI", "Human CBL") %in% rois$roi),
          all(c("Human LG", "Human MG") %in% rois$roi),
          length(unique(rois$anatomy_id)) == 19L)
rd <- function(p) read.csv(p, check.names = FALSE)
j <- rd("tables/s1c/synaptic_density_rcmrglc_matched_regions.csv")
p <- rd("tables/s1d/neuroglial_protein_rcmrglc_matched_regions.csv")
stopifnot(nrow(j) == 12L, nrow(p) == 7L,
          !any(c("thalamus", "centrum_semiovale") %in% j$anatomy_id),
          "cerebellar_cortex" %in% j$anatomy_id,
          !any(c("thalamus", "centrum_semiovale", "midbrain", "cerebellar_cortex") %in% p$anatomy_id))
b <- rd("tables/s2/beip_replication_fits.csv")
stopifnot(identical(sort(b$n), c(6L, 6L, 7L)))
v <- rd("tables/s5/s5a_region_budgets.csv")
stopifnot(sum(v$primary_included) == 8L,
          all(is.na(v$glucose_umol_min[!v$primary_included])),
          all(v$volume_basis[v$primary_included] %in% c("grey", "white")))
wb <- rd("tables/s5/s5b_relative_brain_metabolism.csv")
stopifnot(wb$wholebrain_CMRgl_umol_100g_min[1] == 31,
          grepl("Clarke", wb$rate_source[1]))
s5c <- rd("tables/s5/s5c_regression_volume_change_vs_rate.csv")
stopifnot(all(s5c$n_regions == 7L))
a <- rd("tables/s1c/andersen_paired_regions.csv")
stopifnot(nrow(a) == 26L, all(table(a$group_as_published) == 13L),
          !anyDuplicated(a[c("group_as_published", "anatomy_id")]),
          sum(a$anatomy_id == "Thal") == 2L)
f <- rd("tables/s1_synthesis/heiss1984_matched_inputs.csv")
stopifnot(!any(f$anatomy_id %in% c("medial_thalamic_nuclei", "substantia_nigra", "red_nucleus")))
cat("Strict anatomical-matching scientific invariants passed.\n")
