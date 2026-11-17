setwd(local({ d <- normalizePath(getwd()); while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d); d }))  # repo root

## ============================================================
## Study s1e_1 (PROTOTYPE) — receptor densities x parcellated cortical metabolism
##
## Joins the two Evo-M1 receptor items to the Julich-Brain-parcellated
## metabolic maps produced by scripts/s1e_0_extract_metabolic_maps_julich.R and
## reports DESCRIPTIVE Spearman correlations only (no spatial-null correction,
## no model). Not part of run_all.R / output_manifest.csv.
##
## Inputs (synced by scripts/00_sync_evo_m1_inputs.R; see metadata/evo_m1_input_manifest.csv):
##   data_raw/Zilles_PalomeroGallagher_2017_SupplementaryTable2.csv
##     Zilles K, Palomero-Gallagher N (2017) Front Neuroanat 11:78,
##     doi 10.3389/fnana.2017.00078 — 15 receptors x 44 areas x 4 strata,
##     fmol/mg protein, 3 post-mortem brains.
##   data_raw/Zilles_PalomeroGallagher_2017_EBRAINSFingerprints.csv
##     EBRAINS receptor-density fingerprints (Palomero-Gallagher, Zilles,
##     Amunts et al.; per-dataset DOIs 10.25493/...) anchored to Julich-Brain areas.
##   tables/s1e/s1e_zilles_area_metabolism.csv   (s1e_0; Zilles area x hemisphere)
##   data_intermediate/s1e_julich_area_metabolism.csv (s1e_0; Julich region x hemisphere)
##   metadata/anatomy/s1e_zilles2017_julich_area_map.csv (crosswalk)
##
## Metabolic measures are RELATIVE (Raichle maps as distributed by neuromaps are
## not in umol/100 g/min; see metadata/anatomy/S1E_PARCELLATION_PROTOTYPE.md).
## Outputs: tables/s1e/s1e_receptor_metabolism_join.csv
##          tables/s1e/s1e_receptor_metabolism_spearman.csv
##          tables/s1e/s1e_ebrains_metabolism_join.csv
## ============================================================

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })
dir.create("tables/s1e", showWarnings = FALSE, recursive = TRUE)

## ---- 1. Load --------------------------------------------------------------
zil <- read_csv("data_raw/Zilles_PalomeroGallagher_2017_SupplementaryTable2.csv",
                show_col_types = FALSE, col_types = cols(area = col_character(), .default = col_guess()))
ebr <- read_csv("data_raw/Zilles_PalomeroGallagher_2017_EBRAINSFingerprints.csv", show_col_types = FALSE)
met_z <- read_csv("tables/s1e/s1e_zilles_area_metabolism.csv", show_col_types = FALSE,
                  col_types = cols(zilles_area = col_character(), .default = col_guess()))
met_j <- read_csv("data_intermediate/s1e_julich_area_metabolism.csv", show_col_types = FALSE)
xw    <- read_csv("metadata/anatomy/s1e_zilles2017_julich_area_map.csv", show_col_types = FALSE,
                  col_types = cols(zilles_area = col_character(), .default = col_guess()))

stopifnot(nrow(zil) == 2640, length(unique(zil$area)) == 44, nrow(xw) == 44)

## ---- 2. Receptor table: all-layers densities, wide; class aggregates ------
# Receptor class membership (pharmacology, not data-driven):
#   iGlu = AMPA + NMDA + kainate; GABA = GABAA + GABAA/BZ + GABAB;
#   ACh = M1 + M2 + M3 + alpha4beta2; monoamine GPCR = alpha1 + alpha2 + 5-HT1A + 5-HT2 + D1
rec_names <- sort(unique(zil$receptor))
classes <- list(
  iGlu          = c("AMPA", "NMDA", "kainate"),
  GABA          = grep("^GABA", rec_names, value = TRUE),
  ACh           = c("M1", "M2", "M3", grep("alpha4beta2|a4b2", rec_names, value = TRUE)),
  monoamineGPCR = c(grep("alpha1|^a1$", rec_names, value = TRUE), grep("alpha2|^a2$", rec_names, value = TRUE),
                    grep("5-HT1A|5-HT1a", rec_names, value = TRUE), grep("5-HT2", rec_names, value = TRUE), "D1")
)
classes <- lapply(classes, function(v) intersect(v, rec_names))
missing_members <- setdiff(c("AMPA","NMDA","kainate","M1","M2","M3","D1"), rec_names)
if (length(missing_members)) stop("Receptor names not found in TableS2 CSV: ", paste(missing_members, collapse = ", "))

zil_all <- zil %>% filter(stratum == "all_layers") %>%
  select(area, receptor, density_fmol_mg_protein, is_primary_sensory, is_M1)
zil_wide <- zil_all %>% select(area, receptor, density_fmol_mg_protein) %>%
  pivot_wider(names_from = receptor, values_from = density_fmol_mg_protein)
for (cl in names(classes)) zil_wide[[paste0("class_", cl)]] <- rowSums(zil_wide[, classes[[cl]], drop = FALSE])
zil_wide$EI_receptor_index <- zil_wide$class_iGlu / zil_wide$class_GABA
zil_wide$total_receptor_density <- rowSums(zil_wide[, rec_names, drop = FALSE])
flags <- zil_all %>% distinct(area, is_primary_sensory, is_M1)
zil_wide <- left_join(zil_wide, flags, by = "area")

## ---- 3. Join to bilateral Zilles-area metabolism --------------------------
met_b <- met_z %>% filter(hemisphere == "B") %>%
  select(zilles_area, lobe, match_relationship, julich_regions, n_julich_regions,
         cmrglc, cmr02, cbf, cbv, gi_proxy, castrillon_cmrglc, cmrglc_rel_cortex,
         castrillon_rel_cortex, n_vertices, few_vertices_flag)
join <- zil_wide %>% rename(zilles_area = area) %>%
  left_join(met_b, by = "zilles_area") %>%
  mutate(match_relationship = ifelse(is.na(match_relationship), "none", match_relationship),
         match_strict = match_relationship %in% c("exact", "part_of", "union"),
         match_any    = match_relationship != "none")
write_csv(join, "tables/s1e/s1e_receptor_metabolism_join.csv", na = "")

n_strict <- sum(join$match_strict & !is.na(join$cmrglc)); n_any <- sum(join$match_any & !is.na(join$cmrglc))
message("Zilles areas with metabolism: strict = ", n_strict, ", any match = ", n_any, " of 44")

## ---- 4. Descriptive Spearman correlations ---------------------------------
predictors <- c(rec_names, paste0("class_", names(classes)), "EI_receptor_index", "total_receptor_density")
outcomes   <- c("cmrglc", "cmr02", "cbf", "cbv", "gi_proxy", "castrillon_cmrglc")
spear <- function(d, basis) {
  bind_rows(lapply(predictors, function(p) bind_rows(lapply(outcomes, function(o) {
    ok <- is.finite(d[[p]]) & is.finite(d[[o]])
    if (sum(ok) < 5) return(NULL)
    ct <- suppressWarnings(cor.test(d[[p]][ok], d[[o]][ok], method = "spearman", exact = FALSE))
    tibble(basis = basis, predictor = p, outcome = o, n = sum(ok),
           rho = unname(ct$estimate), p_naive = ct$p.value)
  }))))
}
sp <- bind_rows(spear(filter(join, match_strict), "strict_matches_exact_partof_union"),
                spear(filter(join, match_any),    "any_match_incl_approximate")) %>%
  mutate(note = "descriptive; no spatial-autocorrelation null; metabolic values relative (neuromaps scale)")
write_csv(sp, "tables/s1e/s1e_receptor_metabolism_spearman.csv")

## ---- 5. EBRAINS fingerprints joined directly on Julich region -------------
met_jb <- met_j %>% filter(hemisphere == "B") %>%
  select(julich_region, cmrglc, cmr02, cbf, cbv, gi_proxy, castrillon_cmrglc, cmrglc_rel_cortex, n_vertices)
ebr_join <- ebr %>% select(julich_area, dataset_locus, dataset_name, dataset_version, receptor,
                           density_fmol_mg_protein, sd, dataset_doi) %>%
  left_join(met_jb, by = c("julich_area" = "julich_region"))
write_csv(ebr_join, "tables/s1e/s1e_ebrains_metabolism_join.csv", na = "")
ebr_cov <- ebr_join %>% distinct(julich_area, dataset_doi, matched = !is.na(cmrglc))
message("EBRAINS datasets: ", n_distinct(ebr_cov$dataset_doi), "; Julich anchors: ", n_distinct(ebr_cov$julich_area),
        "; anchors with cortical metabolism: ", n_distinct(ebr_cov$julich_area[ebr_cov$matched]))

## ---- 6. Console summary (strict basis, all-layers cmrglc) -----------------
top <- sp %>% filter(basis == "strict_matches_exact_partof_union", outcome == "cmrglc") %>%
  arrange(desc(abs(rho))) %>% select(predictor, n, rho, p_naive)
print(as.data.frame(top), digits = 3)
