setwd(local({ d <- normalizePath(getwd()); while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d); d }))  # repo root

## ============================================================
## Study s1e_2 (PROTOTYPE) — middle tier: Zilles & Palomero-Gallagher 2017
## receptor areas aggregated to the gyral-group regions of Heiss et al. 1985
## (absolute rCMRglc, FDG-PET, umol/100 g/min), with Chugani et al. 1987 adult
## lCMRGlc as a second outcome where a Heiss group has a Chugani counterpart.
##
## DESCRIPTIVE only: unweighted means over member areas, Spearman rho with
## naive p-values, no spatial-autocorrelation null, no model.
## Not part of run_all.R / output_manifest.csv.
##
## Inputs (synced by scripts/00_sync_evo_m1_inputs.R; metadata/evo_m1_input_manifest.csv):
##   data_raw/Zilles_PalomeroGallagher_2017_SupplementaryTable2.csv
##     15 receptors x 44 areas x 4 strata, fmol/mg protein (3 post-mortem brains).
##   data_raw/Heiss_etal_1985_TABLEI.csv
##     Heiss W-D, Pawlik G, Herholz K, Wagner R, Wienhard K (1985) Brain Res 327:362-366,
##     doi 10.1016/0006-8993(85)91537-9. TABLE I, long format: one row per
##     structure x arm. Arms: nondreamers_awake (n = 3, mean +/- SD),
##     nondreamers_pct_change_asleep, dreamer_awake (n = 1), dreamer_pct_change_asleep.
##     4 healthy male volunteers in total (3 non-dreamers + 1 dreamer); sleep-deprived
##     waking condition (see Heiss_etal_1985/ANATOMICAL_MATCHING.md in Evo-M1).
##   data_raw/Chugani_etal_1987_Table2.csv
##     Chugani HT, Phelps ME, Mazziotta JC (1987) Ann Neurol 22:487-497,
##     doi 10.1002/ana.410220408. Table 2, long format; adult reference group
##     (19-30 y, n = 7) flagged by is_adult_reference_group; mean +/- SEM, umol/min/100 g.
##   metadata/anatomy/s1e_zilles2017_heiss1985_group_map.csv  (44 rows; crosswalk)
##
## Heiss 1985 column used: the AWAKE arms only. Primary outcome =
##   cmrglc_heiss1985_awake_nd3  (non-dreamers, n = 3, the only awake arm with a
##   printed SD). Secondary = cmrglc_heiss1985_awake_pooled4, the n-weighted mean
##   of all 4 awake subjects, (3 x non-dreamer mean + dreamer value) / 4 (the pool
##   Karbowski 2007 used; no SD derivable). Sleep arms (% change) are NOT used.
##
## Outputs: tables/s1e/s1e_middle_tier_group_receptors.csv
##          tables/s1e/s1e_middle_tier_spearman.csv
##          tables/s1e/s1e_middle_tier_group_coverage.csv
## ============================================================

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })
dir.create("tables/s1e", showWarnings = FALSE, recursive = TRUE)

## ---- 1. Load --------------------------------------------------------------
zil <- read_csv("data_raw/Zilles_PalomeroGallagher_2017_SupplementaryTable2.csv",
                show_col_types = FALSE, col_types = cols(area = col_character(), .default = col_guess()))
gmap <- read_csv("metadata/anatomy/s1e_zilles2017_heiss1985_group_map.csv",
                 show_col_types = FALSE, col_types = cols(.default = col_character()))
heiss <- read_csv("data_raw/Heiss_etal_1985_TABLEI.csv", show_col_types = FALSE,
                  col_types = cols(.default = col_character()))
chug <- read_csv("data_raw/Chugani_etal_1987_Table2.csv", show_col_types = FALSE)

stopifnot(nrow(zil) == 2640, length(unique(zil$area)) == 44, nrow(gmap) == 44,
          setequal(gmap$zilles_area, unique(zil$area)),
          all(gmap$assignment_relationship %in% c("contained", "partial", "approximate", "none")),
          nrow(heiss) == 96)

## ---- 2. Heiss 1985 awake rCMRglc per structure ----------------------------
heiss_num <- heiss %>%
  mutate(value = as.numeric(CMRglc.umol_100g_min), sd = suppressWarnings(as.numeric(dispersion_numeric)),
         n = as.integer(n))
nd <- heiss_num %>% filter(arm_id == "nondreamers_awake") %>%
  select(heiss1985_structure_as_published = structure_as_published,
         cmrglc_heiss1985_awake_nd3 = value, cmrglc_heiss1985_awake_nd3_sd = sd, heiss1985_n_nd = n)
dr <- heiss_num %>% filter(arm_id == "dreamer_awake") %>%
  select(heiss1985_structure_as_published = structure_as_published,
         cmrglc_heiss1985_awake_dreamer1 = value, heiss1985_n_dr = n)
stopifnot(nrow(nd) == 24, nrow(dr) == 24, all(nd$heiss1985_n_nd == 3), all(dr$heiss1985_n_dr == 1))
heiss_awake <- full_join(nd, dr, by = "heiss1985_structure_as_published") %>%
  mutate(cmrglc_heiss1985_awake_pooled4 = (heiss1985_n_nd * cmrglc_heiss1985_awake_nd3 +
                                           heiss1985_n_dr * cmrglc_heiss1985_awake_dreamer1) /
                                          (heiss1985_n_nd + heiss1985_n_dr),
         heiss1985_n_subjects_total = heiss1985_n_nd + heiss1985_n_dr) %>%
  select(-heiss1985_n_dr)
stopifnot(all(heiss_awake$heiss1985_n_subjects_total == 4))

## ---- 3. Chugani 1987 adult lCMRGlc per Heiss group ------------------------
# Group-level counterpart (only where a Chugani structure corresponds to a Heiss gyral group;
# Chugani's lobe-level rows 'Frontal/Parietal/Temporal cortex' have no single Heiss group).
chug_adult <- chug %>% filter(is_adult_reference_group) %>%
  select(chugani1987_structure_as_printed = structure_as_printed,
         lcmrglc_chugani1987_adult = lCMRGlc_mean.umol_min_100g,
         lcmrglc_chugani1987_adult_sem = lCMRGlc_SEM.umol_min_100g,
         chugani1987_n_adult = n_subjects)
stopifnot(nrow(chug_adult) == 14, all(chug_adult$chugani1987_n_adult == 7))
chug_groups <- tribble(
  ~heiss1985_group_id, ~chugani1987_structure_as_printed, ~chugani_match,
  "sensorimotor", "Sensorimotor cortex",        "match",
  "prim_visual",  "Calcarine cortex",           "match",
  "cingulate",    "Anterior cingulate cortex",  "partial (anterior cingulate only; Heiss group includes 23/31)",
  "sup_temporal", "Transverse temporal cortex", "partial (Heschl's gyrus only; Heiss group includes 22)",
  "vis_assoc",    "Occipital cortex",           "approximate (Chugani extrastriate 'Occipital cortex' vs Heiss 'Vis. assoc. cortex')"
) %>% left_join(chug_adult, by = "chugani1987_structure_as_printed")
stopifnot(!anyNA(chug_groups$lcmrglc_chugani1987_adult))
# consistency with the per-area crosswalk column
xw_chug <- gmap %>% filter(!is.na(chugani1987_structure_as_printed), chugani1987_structure_as_printed != "NA") %>%
  distinct(heiss1985_group_id, chugani1987_structure_as_printed)
stopifnot(nrow(anti_join(xw_chug, chug_groups, by = c("heiss1985_group_id", "chugani1987_structure_as_printed"))) == 0)

## ---- 4. Receptors: unweighted mean over member areas, per group x stratum --
rec_names <- sort(unique(zil$receptor))
classes <- list(
  iGlu          = c("AMPA", "NMDA", "kainate"),
  GABA          = grep("^GABA", rec_names, value = TRUE),
  ACh           = c("M1", "M2", "M3", grep("alpha4beta2|a4b2", rec_names, value = TRUE)),
  monoamineGPCR = c(grep("alpha1|^a1$", rec_names, value = TRUE), grep("alpha2|^a2$", rec_names, value = TRUE),
                    grep("5-HT1A|5-HT1a", rec_names, value = TRUE), grep("5-HT2", rec_names, value = TRUE), "D1")
)
classes <- lapply(classes, function(v) intersect(v, rec_names))
stopifnot(length(classes$iGlu) == 3, length(classes$GABA) == 3, length(classes$ACh) == 4, length(classes$monoamineGPCR) == 5)

members <- gmap %>% filter(assignment_relationship != "none") %>%
  select(zilles_area, heiss1985_group_id, heiss1985_structure_as_published, assignment_relationship)

long <- zil %>% inner_join(members, by = c("area" = "zilles_area")) %>%
  group_by(heiss1985_group_id, heiss1985_structure_as_published, stratum, receptor) %>%
  summarise(mean_density = mean(density_fmol_mg_protein, na.rm = TRUE), n_areas = n(),
            n_areas_nonmissing = sum(!is.na(density_fmol_mg_protein)),
            member_areas = paste(sort(area), collapse = ";"), .groups = "drop") %>%
  mutate(mean_density = ifelse(is.nan(mean_density), NA_real_, mean_density))
# Area 24 is agranular and has no granular-stratum values in the source table; with na.rm the
# cingulate granular mean is over 23/31/32 only. Minimum n_areas_nonmissing per group x stratum
# is recorded in the output (n_areas_nonmissing_min).

wide <- long %>% select(-member_areas, -n_areas_nonmissing) %>%
  pivot_wider(names_from = receptor, values_from = mean_density)
# n_areas and member list are receptor-invariant: keep one copy
meta <- long %>% group_by(heiss1985_group_id, stratum, n_areas, member_areas) %>%
  summarise(n_areas_nonmissing_min = min(n_areas_nonmissing), .groups = "drop")
stopifnot(nrow(meta) == nrow(wide))
wide <- wide %>% select(-n_areas) %>% left_join(meta, by = c("heiss1985_group_id", "stratum"))
for (cl in names(classes)) wide[[paste0("class_", cl)]] <- rowSums(wide[, classes[[cl]], drop = FALSE])
wide$EI_receptor_index <- wide$class_iGlu / wide$class_GABA
wide$total_receptor_density <- rowSums(wide[, rec_names, drop = FALSE])

cov_tbl <- members %>% group_by(heiss1985_group_id) %>%
  summarise(n_contained = sum(assignment_relationship == "contained"),
            n_partial = sum(assignment_relationship == "partial"),
            n_approximate = sum(assignment_relationship == "approximate"), .groups = "drop")

## ---- 5. Join outcomes ------------------------------------------------------
grp <- wide %>%
  left_join(cov_tbl, by = "heiss1985_group_id") %>%
  left_join(heiss_awake, by = "heiss1985_structure_as_published") %>%
  left_join(chug_groups, by = "heiss1985_group_id") %>%
  arrange(factor(stratum, levels = c("all_layers", "supragranular", "granular", "infragranular")),
          heiss1985_group_id) %>%
  relocate(heiss1985_group_id, heiss1985_structure_as_published, stratum, n_areas, n_areas_nonmissing_min,
           n_contained, n_partial, n_approximate, member_areas)
stopifnot(!anyNA(grp$cmrglc_heiss1985_awake_nd3))
write_csv(grp, "tables/s1e/s1e_middle_tier_group_receptors.csv", na = "")

## ---- 6. Coverage table (all 24 Heiss 1985 structures) ---------------------
coverage <- heiss_awake %>% select(heiss1985_structure_as_published) %>%
  left_join(gmap %>% filter(assignment_relationship != "none") %>%
              group_by(heiss1985_structure_as_published, heiss1985_group_id) %>%
              summarise(n_zilles_areas = n(),
                        n_contained = sum(assignment_relationship == "contained"),
                        n_partial = sum(assignment_relationship == "partial"),
                        n_approximate = sum(assignment_relationship == "approximate"),
                        zilles_areas = paste(zilles_area, collapse = ";"), .groups = "drop"),
            by = "heiss1985_structure_as_published") %>%
  mutate(across(c(n_zilles_areas, n_contained, n_partial, n_approximate), ~ replace_na(.x, 0L)),
         zilles_areas = replace_na(zilles_areas, ""),
         has_chugani1987_counterpart = heiss1985_group_id %in% chug_groups$heiss1985_group_id) %>%
  left_join(heiss_awake %>% select(heiss1985_structure_as_published, cmrglc_heiss1985_awake_nd3,
                                   cmrglc_heiss1985_awake_pooled4), by = "heiss1985_structure_as_published")
write_csv(coverage, "tables/s1e/s1e_middle_tier_group_coverage.csv", na = "")

## ---- 7. Descriptive Spearman correlations ---------------------------------
predictors <- c(rec_names, paste0("class_", names(classes)), "EI_receptor_index", "total_receptor_density")
outcomes   <- c("cmrglc_heiss1985_awake_nd3", "cmrglc_heiss1985_awake_pooled4", "lcmrglc_chugani1987_adult")
sp <- bind_rows(lapply(unique(grp$stratum), function(st) {
  d <- grp %>% filter(stratum == st)
  bind_rows(lapply(predictors, function(p) bind_rows(lapply(outcomes, function(o) {
    ok <- is.finite(d[[p]]) & is.finite(d[[o]])
    if (sum(ok) < 5) return(NULL)
    ct <- suppressWarnings(cor.test(d[[p]][ok], d[[o]][ok], method = "spearman", exact = FALSE))
    tibble(stratum = st, predictor = p, outcome = o, n = sum(ok),
           rho = unname(ct$estimate),
           p_naive = if (sum(ok) >= 8) ct$p.value else NA_real_)  # n < 8: rho only
  }))))
})) %>%
  mutate(unit_of_analysis = "Heiss 1985 gyral group (unweighted mean of member Zilles areas)",
         note = "descriptive; no spatial-autocorrelation null; p_naive omitted when n < 8 groups; outcomes absolute umol/100 g/min; Heiss n = 3 (nd3) / 4 (pooled4) subjects, Chugani adult n = 7")
write_csv(sp, "tables/s1e/s1e_middle_tier_spearman.csv")

## ---- 8. Console summary ----------------------------------------------------
message("Groups with receptor data: ", n_distinct(grp$heiss1985_group_id), " of ", nrow(heiss_awake),
        " Heiss 1985 structures; Zilles areas used: ", sum(coverage$n_zilles_areas), " of 44")
print(as.data.frame(coverage %>% filter(n_zilles_areas > 0) %>%
                      select(heiss1985_group_id, n_zilles_areas, n_contained, has_chugani1987_counterpart,
                             cmrglc_heiss1985_awake_nd3)), digits = 3)
top <- sp %>% filter(stratum == "all_layers", outcome == "cmrglc_heiss1985_awake_nd3") %>%
  arrange(desc(abs(rho))) %>% select(predictor, n, rho, p_naive)
print(as.data.frame(top), digits = 3)
