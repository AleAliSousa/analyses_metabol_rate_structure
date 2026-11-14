# Scientific invariants for the 12 Nov analysis expansion. Run from repo root.
source("helpers/regional_robustness.R")
rd <- function(p) read.csv(p, check.names = FALSE, stringsAsFactors = FALSE)
st <- rd("data_intermediate/s1a_stereology_comparison.csv")
fr <- st[st$source_region == "Frontal lobe", ]
stopifnot(nrow(st) == 9L, nrow(fr) == 1L,
          abs(fr$Volume - 86.2428) < 1e-8,
          abs(fr$Neuron_N - 6064700000) < 1,
          abs(fr$NeurDensity - 6064700000 / 86.2428) < 1e-5)
au <- rd("data_intermediate/s1a_volume_basis_audit.csv")
stopifnot(nrow(au) == 5L, all(au$density_multiplier_vs_previous > 1.7),
          all(au$density_multiplier_vs_previous < 1.91))
h <- rd("tables/s1c/holmes_rcmrglc_matched_regions.csv")
stopifnot(nrow(h) == 13L, sum(h$primary_included) == 11L,
          all(h$is_healthy_control == 1), all(h$measure == "BPND"),
          !anyDuplicated(h$anatomy_id),
          !any(h$region_as_published == "Ventral striatum"))
b <- rd("tables/s2/beip_volume_contrasts.csv")
stopifnot(nrow(b) == 25L, length(unique(b$contrast)) == 3L,
          !anyDuplicated(b[c("contrast", "anatomy_id")]))
ca <- b[b$timepoint == "BEIP age 8-10" & b$anatomy_id == "caudate_nucleus", ]
stopifnot(nrow(ca) == 1L,
          abs(ca$volume_deficit_pct - 100 * (8693.7 - 8338.2) / 8693.7) < 1e-10)
am <- b[b$timepoint == "BEIP age 16" & b$group_code == "FCG" & b$anatomy_id == "amygdala", ]
stopifnot(nrow(am) == 1L, am$volume == 3, am$control_volume == 3,
          am$volume_deficit_pct == 0)
f <- rd("tables/s1_synthesis/heiss1984_sensitivity_fits.csv")
keys <- paste(f$study, f$variable, f$scope)
for (k in unique(keys)) {
  z <- f[keys == k, ]
  stopifnot(nrow(z) == 2L, length(unique(z$n)) == 1L)
}
counts <- rd("tables/s1_synthesis/shared_region_counts.csv")
stopifnot(!any(grepl("Johansen / ERA|Proteins / ERA", counts$comparison)),
          counts$n_shared_regions[counts$comparison == "Johansen / Proteins" &
                                     counts$scope == "shared_regions"] == 6L)
tiny <- data.frame(anatomy_id = c("a", "b"), x = 1:2, y = 2:3)
stopifnot(is.na(regional_fit(tiny, "x", "y")$p_ols))
constant <- data.frame(anatomy_id = letters[1:5], x = 1:5, y = rep(2, 5))
stopifnot(is.na(regional_fit(constant, "x", "y")$p_ols))
duplicate <- data.frame(anatomy_id = rep("a", 5), x = 1:5, y = 2:6)
stopifnot(inherits(try(regional_fit(duplicate, "x", "y"), silent = TRUE), "try-error"))
cat("All 12 Nov scientific invariants passed.\n")
