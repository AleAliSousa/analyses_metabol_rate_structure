setwd(local({ d <- normalizePath(getwd()); while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d); d }))  # repo root

## ============================================================
## Study s1e_3 (PROTOTYPE) — human <-> macaque (<-> rat) region-group comparison:
## receptor classes vs regional CMRglc, aggregated to the shared region groups of
## metadata/anatomy/human_macaque_region_group_homology.csv (66 groups; names only;
## 63 original + 3 rat-only groups added 2026-11-17: olfactory_cortex, thalamus_lateral_nuclei,
## thalamus_ventral_nuclear_group).
##
## DESCRIPTIVE only: unweighted means over member areas / structures, within-species
## z-scores, Spearman rho with naive p-values, no spatial-autocorrelation null, no model.
## Not part of run_all.R / output_manifest.csv.
##
## Inputs (Evo-M1 items synced by scripts/00_sync_evo_m1_inputs.R; metadata/evo_m1_input_manifest.csv):
##   HUMAN receptors : data_raw/Zilles_PalomeroGallagher_2017_SupplementaryTable2.csv
##                     15 receptors x 44 areas x 4 strata, fmol/mg PROTEIN; stratum all_layers used.
##   HUMAN CMRglc    : data_raw/Heiss_etal_1985_TABLEI.csv  (FDG PET, awake non-dreamer arm, n = 3;
##                     column choice as in s1e_2: cmrglc_heiss1985_awake_nd3)
##                     data_intermediate/heiss_2004_regions.csv (Heiss 2004 Table 1; lobes + subcortex)
##                     data_raw/Chugani_etal_1987_Table2.csv (adult reference group; optional fill)
##   MACAQUE receptors: data_raw/FroudistWalsh_etal_2023_SupplementaryTable1.csv
##                     14 receptors x 109 areas, fmol per NEURON, Macaca fascicularis.
##   RAT receptors   : data_raw/PalomeroGallagher_etal_2018_rat*.csv — 10 EBRAINS items (fmol/mg PROTEIN, Wistar,
##                     n = 5 brains): Oc1B, Oc1M (primary visual), Oc2L, Oc2M (secondary visual), CA1, CA2, CA3, DG,
##                     and CA stratum cellulare / stratum moleculare. The two stratum items are laminar parts of CA and
##                     are EXCLUDED from group means (they would double-count CA1-3 at a different grain). The four
##                     hippocampal items carry 10 receptors only (no D1, GABAA/BZ, GABAB, M3, alpha4beta2), so for rat
##                     hippocampus_dentate only class_iGlu is complete; the other classes are NA there.
##   RAT CMRglc      : data_raw/Sokoloff_etal_1977_Table3.csv (2-DG, conscious Sprague-Dawley rats, 30-min arm, n = 10,
##                     33 structures; Table 4 (thiopental) is synced but not used here)
##   MACAQUE CMRglc  : data_raw/Kennedy_etal_1978_Table.csv (2-DG, conscious M. mulatta, n = 7)
##                     data_raw/Palombo_etal_1990_Table1..5.csv (2-DG, CONTROL animals n = 4, right+left
##                     hemisphere arms averaged; MPTP arms NOT used)
##                     data_raw/Noda_etal_2002_Table3.csv (FDG PET, conscious young rhesus; average arm)
##                     data_raw/Meguro_etal_1999_Table2.csv (FDG PET, ANAESTHETISED BABOON, baseline arms
##                     of both surgery groups averaged; sensitivity only)
##   metadata/anatomy/human_macaque_region_group_homology.csv (crosswalk, verbatim names per dataset)
##
## Receptor classes (pharmacology, as in s1e_1/s1e_2):
##   iGlu = AMPA + NMDA + kainate; GABA = GABAA + GABAA/BZ + GABAB; monoamineGPCR = alpha1 + alpha2 +
##   5-HT1A + 5-HT2 + D1; EI = iGlu / GABA. ACh: FW2023 has NO alpha4beta2, so the cross-species
##   class_ACh is harmonised to M1 + M2 + M3 in BOTH species; the human 4-receptor version
##   (M1+M2+M3+alpha4beta2) is kept as class_ACh_incl_a4b2 (human only). total_shared14 = sum of the
##   14 receptors both datasets share; human total_15 adds alpha4beta2.
##
## UNITS: human receptor values are fmol/mg protein, macaque fmol/neuron; human CMRglc is FDG-PET,
## macaque 2-DG autoradiography (Noda/Meguro PET). Raw values must NEVER be compared across species;
## only within-species z-scores / ranks are compared.
##
## Outputs: tables/s1e/species_group_panel.csv      (+ copy under the old name human_macaque_group_panel.csv)
##          tables/s1e/species_group_spearman.csv   (+ copy human_macaque_spearman.csv)
##          tables/s1e/species_group_coverage.csv   (+ copy human_macaque_group_coverage.csv)
##          The old filenames are kept as byte-identical copies so existing links keep working.
## Figure: scripts/s1e_3_fig_human_macaque.py (reads the panel CSV) -> figs/s1e/human_macaque_*.png/.pdf and
##         figs/s1e/cmrglc_concordance_three_species.png/.pdf
## ============================================================

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })
dir.create("tables/s1e", showWarnings = FALSE, recursive = TRUE)
chr <- cols(.default = col_character())
split_names <- function(x) if (is.na(x) || x == "") character(0) else strsplit(x, ";", fixed = TRUE)[[1]]

## ---- 1. Load --------------------------------------------------------------
hom <- read_csv("metadata/anatomy/human_macaque_region_group_homology.csv", show_col_types = FALSE, col_types = chr)
zil <- read_csv("data_raw/Zilles_PalomeroGallagher_2017_SupplementaryTable2.csv", show_col_types = FALSE,
                col_types = cols(area = col_character(), .default = col_guess()))
fw  <- read_csv("data_raw/FroudistWalsh_etal_2023_SupplementaryTable1.csv", show_col_types = FALSE,
                col_types = cols(area = col_character(), .default = col_guess()))
heiss <- read_csv("data_raw/Heiss_etal_1985_TABLEI.csv", show_col_types = FALSE, col_types = chr)
h04 <- read_csv("data_intermediate/heiss_2004_regions.csv", show_col_types = FALSE, col_types = chr)
chug <- read_csv("data_raw/Chugani_etal_1987_Table2.csv", show_col_types = FALSE)
ken <- read_csv("data_raw/Kennedy_etal_1978_Table.csv", show_col_types = FALSE, col_types = chr)
pal <- bind_rows(lapply(1:5, function(i) read_csv(sprintf("data_raw/Palombo_etal_1990_Table%d.csv", i),
                                                   show_col_types = FALSE, col_types = chr) %>% mutate(palombo_table = i)))
noda <- read_csv("data_raw/Noda_etal_2002_Table3.csv", show_col_types = FALSE, col_types = chr)
meg <- read_csv("data_raw/Meguro_etal_1999_Table2.csv", show_col_types = FALSE, col_types = chr)
sok <- read_csv("data_raw/Sokoloff_etal_1977_Table3.csv", show_col_types = FALSE, col_types = chr)
rat_files <- sort(list.files("data_raw", pattern = "^PalomeroGallagher_etal_2018_rat.*\\.csv$", full.names = TRUE))
rat <- bind_rows(lapply(rat_files, read_csv, show_col_types = FALSE,
                        col_types = cols(density_fmol_mg_protein = col_double(), sd = col_double(), n_brains = col_integer(),
                                         .default = col_character())))

stopifnot(nrow(hom) == 66, nrow(zil) == 2640, n_distinct(zil$area) == 44,
          n_distinct(fw$area) == 109, nrow(fw) == 109 * 14, all(fw$Species == "Macaca fascicularis"),
          nrow(heiss) == 96, nrow(h04) == 26, nrow(ken) == 39, nrow(noda) == 57, nrow(meg) == 152,
          nrow(sok) == 66, length(rat_files) == 10, nrow(rat) == 4 * 15 + 6 * 10,
          all(rat$Species == "Rattus norvegicus"), all(rat$unit_as_published == "fmol/mg"), all(rat$n_brains == 5L))

## ---- 2. Receptor classes on the shared receptor vocabulary ----------------
rec_h <- sort(unique(zil$receptor)); rec_m <- sort(unique(fw$receptor))
shared_rec <- intersect(rec_h, rec_m)
stopifnot(length(shared_rec) == 14, setdiff(rec_h, rec_m) == "alpha4beta2", length(setdiff(rec_m, rec_h)) == 0)
classes <- list(
  iGlu          = c("AMPA", "NMDA", "kainate"),
  GABA          = grep("^GABA", shared_rec, value = TRUE),
  ACh           = c("M1", "M2", "M3"),                                   # harmonised (no alpha4beta2 in FW2023)
  monoamineGPCR = c("alpha1", "alpha2", "5-HT1A", "5-HT2", "D1")
)
stopifnot(all(unlist(classes) %in% shared_rec), length(classes$GABA) == 3)
class_names <- c(paste0("class_", names(classes)), "EI_receptor_index", "total_shared14")

area_classes <- function(wide, recs) {
  for (cl in names(classes)) wide[[paste0("class_", cl)]] <- rowSums(wide[, classes[[cl]], drop = FALSE])
  wide$EI_receptor_index <- wide$class_iGlu / wide$class_GABA
  wide$total_shared14 <- rowSums(wide[, shared_rec, drop = FALSE])
  wide
}
zil_wide <- zil %>% filter(stratum == "all_layers") %>% select(area, receptor, density_fmol_mg_protein) %>%
  pivot_wider(names_from = receptor, values_from = density_fmol_mg_protein) %>% area_classes()
zil_wide$class_ACh_incl_a4b2 <- zil_wide$class_ACh + zil_wide$alpha4beta2
zil_wide$total_15 <- zil_wide$total_shared14 + zil_wide$alpha4beta2
fw_wide <- fw %>% select(area, receptor, value_fmol_per_neuron) %>%
  pivot_wider(names_from = receptor, values_from = value_fmol_per_neuron) %>% area_classes()
stopifnot(nrow(zil_wide) == 44, nrow(fw_wide) == 109, !anyNA(zil_wide[, class_names]), !anyNA(fw_wide[, class_names]))
## Rat: 10 items -> 10 structures (strata included in the wide table; the homology column lists only the 8 whole structures).
## Receptors absent from an item give NA in the wide table; a class is NA wherever any member receptor is NA (no partial sums).
rat_wide <- rat %>% select(area = structure_as_published, receptor, density_fmol_mg_protein) %>%
  pivot_wider(names_from = receptor, values_from = density_fmol_mg_protein)
stopifnot(nrow(rat_wide) == 10, setequal(setdiff(names(rat_wide), "area"), rec_h))   # rat vocabulary = the human 15 receptors
rat_wide <- area_classes(rat_wide)
rat_missing <- rat %>% group_by(area = structure_as_published) %>%
  summarise(missing = paste(setdiff(rec_h, receptor), collapse = ";"), n_receptors = n(), .groups = "drop")
stopifnot(all(rat_missing$n_receptors[grepl("^Area Oc", rat_missing$area)] == 15),
          all(rat_missing$n_receptors[!grepl("^Area Oc", rat_missing$area)] == 10))

## Unweighted mean over the member areas listed verbatim in the homology table
group_receptors <- function(name_col, wide, vars) {
  bind_rows(lapply(seq_len(nrow(hom)), function(i) {
    areas <- split_names(hom[[name_col]][i]); if (!length(areas)) return(NULL)
    miss <- setdiff(areas, wide$area); if (length(miss)) stop(name_col, ": areas not in receptor table: ", paste(miss, collapse = ", "))
    d <- wide[match(areas, wide$area), ]
    out <- tibble(group_id = hom$group_id[i], n_member_areas = length(areas), member_areas = paste(areas, collapse = ";"))
    for (v in vars) out[[v]] <- mean(d[[v]])
    out
  }))
}
rec_human <- group_receptors("Zilles_PalomeroGallagher_2017_areas", zil_wide, c(class_names, "class_ACh_incl_a4b2", "total_15"))
rec_macaque <- group_receptors("FroudistWalsh_etal_2023_SuppTable1_areas", fw_wide, class_names)
stopifnot(sum(rec_human$n_member_areas) == 44, sum(rec_macaque$n_member_areas) == 109)   # every area used exactly once
rec_rat <- group_receptors("EBRAINS_rat_receptors", rat_wide, class_names)
stopifnot(nrow(rec_rat) == 3, sum(rec_rat$n_member_areas) == 8)   # 10 items minus the 2 laminar (stratum) items
rec_rat$receptors_missing <- sapply(rec_rat$member_areas, function(m) {
  a <- split_names(m); paste(unique(unlist(lapply(rat_missing$missing[match(a, rat_missing$area)], split_names))), collapse = ";") })

## ---- 3. CMRglc per dataset, then per group (mean over member structures) --
num <- function(x) suppressWarnings(as.numeric(x))
# Human — Heiss 1985 awake non-dreamers (n = 3), column choice as in s1e_2
h85 <- heiss %>% filter(arm_id == "nondreamers_awake") %>%
  transmute(structure = structure_as_published, value = num(CMRglc.umol_100g_min), n = as.integer(n))
stopifnot(nrow(h85) == 24, all(h85$n == 3), !anyNA(h85$value))
# Human — Heiss 2004 Table 1 (canonical labels of heiss_2004_regions.csv)
h04v <- h04 %>% transmute(structure = canonical_label, value = num(rcmrglc_mean_umol_100g_min))
stopifnot(!anyNA(h04v$value))
# Human — Chugani 1987 adult reference group (n = 7)
chugv <- chug %>% filter(is_adult_reference_group) %>%
  transmute(structure = structure_as_printed, value = lCMRGlc_mean.umol_min_100g, n = n_subjects)
stopifnot(nrow(chugv) == 14, all(chugv$n == 7))
# Macaque — Kennedy 1978 conscious rhesus (n = 7)
kenv <- ken %>% filter(arm_id == "conscious_normal") %>%
  transmute(structure = structure_as_published, value = num(CMRglc.umol_100g_min), n = as.integer(n))
stopifnot(nrow(kenv) == 39, !anyNA(kenv$value), all(kenv$n %in% c(6L, 7L)))  # 3 rows (spinal gray/white, optic chiasm) n = 6
# Macaque — Palombo 1990 CONTROL animals only: mean of control_right and control_left per structure
palv <- pal %>% filter(arm_id %in% c("control_right", "control_left")) %>%
  mutate(value = num(CMRglc.umol_100g_min)) %>%
  group_by(palombo_table, structure = structure_as_published) %>%
  summarise(value = mean(value), n_arms = n(), .groups = "drop")
stopifnot(all(palv$n_arms == 2), !anyNA(palv$value), nrow(palv) == 18 + 15 + 14 + 15 + 20,
          all(pal$group_as_published[pal$arm_id %in% c("control_right", "control_left")] == "Controls (n = 4)"))
# Macaque — Noda 2002 young adults, hemisphere-average arm
nodav <- noda %>% filter(arm_id == "young_average", side_or_part == "average") %>%
  transmute(structure = structure_as_published, value = num(CMRglc.umol_100g_min))
stopifnot(nrow(nodav) == 7, !anyNA(nodav$value))
# Rat — Sokoloff 1977 Table 3, conscious 30-min arm (n = 10); the 45-min arm is a replication, not used
sokv <- sok %>% filter(arm_id == "conscious_30min") %>%
  transmute(structure = structure_as_published, value = num(CMRglc.umol_100g_min), n = as.integer(n))
stopifnot(nrow(sokv) == 33, !anyNA(sokv$value), all(sokv$n == 10L))
# Baboon — Meguro 1999 pre-surgery baseline of both groups (RH, SH; n = 3 each), anaesthetised
megv <- meg %>% filter(pre_or_post_surgery == "pre-surgery") %>%
  mutate(value = num(CMRglc.umol_100g_min)) %>% group_by(structure = structure_as_published) %>%
  summarise(value = mean(value), n_arms = n(), .groups = "drop")
stopifnot(nrow(megv) == 19, all(megv$n_arms == 2), all(meg$species_as_published == "Papio anubis"))

group_cmrglc <- function(name_cols, src, label) {
  bind_rows(lapply(seq_len(nrow(hom)), function(i) {
    structs <- unlist(lapply(name_cols, function(cc) split_names(hom[[cc]][i]))); if (!length(structs)) return(NULL)
    miss <- setdiff(structs, src$structure); if (length(miss)) stop(label, ": structures not in table: ", paste(miss, collapse = ", "))
    tibble(group_id = hom$group_id[i], value = mean(src$value[match(structs, src$structure)]),
           n_structures = length(structs), structures = paste(structs, collapse = ";"))
  })) %>% rename_with(~ paste0(.x, "_", label), -group_id)
}
cm_h85  <- group_cmrglc("Heiss_etal_1985_TABLEI", h85, "heiss1985_nd3")
cm_h04  <- group_cmrglc("Heiss_etal_2004", h04v, "heiss2004")
cm_chug <- group_cmrglc("Chugani_etal_1987_Table2", chugv, "chugani1987_adult")
cm_ken  <- group_cmrglc("Kennedy_etal_1978_Table", kenv, "kennedy1978")
cm_pal  <- group_cmrglc(paste0("Palombo_etal_1990_Table", 1:5), palv, "palombo1990_ctrl")
cm_noda <- group_cmrglc("Noda_etal_2002_Table3", nodav, "noda2002_young")
cm_meg  <- group_cmrglc("Meguro_etal_1999_Table2_baboon", megv, "meguro1999_baboon")
cm_sok  <- group_cmrglc("Sokoloff_etal_1977_Table3", sokv, "sokoloff1977_30min")
stopifnot(sum(cm_h85$n_structures_heiss1985_nd3) == 24, sum(cm_h04$n_structures_heiss2004) == 26,
          sum(cm_ken$n_structures_kennedy1978) == 39, sum(cm_pal$n_structures_palombo1990_ctrl) == 82,
          sum(cm_noda$n_structures_noda2002_young) == 7, sum(cm_meg$n_structures_meguro1999_baboon) == 19,
          sum(cm_sok$n_structures_sokoloff1977_30min) == 33, nrow(cm_sok) == 26)

base <- hom %>% transmute(group_id, group_label, tissue_class,
                          homology_confidence = sub(":.*$", "", homology_confidence),
                          homology_confidence_reason = homology_confidence)
cm_human <- base %>% left_join(cm_h85, "group_id") %>% left_join(cm_h04, "group_id") %>% left_join(cm_chug, "group_id") %>%
  mutate(
    # primary: cortex -> Heiss 1985 gyral group (fallback Heiss 2004 lobe); non-cortex -> Heiss 2004 (fallback Heiss 1985)
    cmrglc = case_when(tissue_class == "cortex" & !is.na(value_heiss1985_nd3) ~ value_heiss1985_nd3,
                       tissue_class == "cortex" ~ value_heiss2004,
                       !is.na(value_heiss2004) ~ value_heiss2004,
                       TRUE ~ value_heiss1985_nd3),
    cmrglc_source = case_when(is.na(cmrglc) ~ NA_character_,
                              tissue_class == "cortex" & !is.na(value_heiss1985_nd3) ~ "Heiss1985_nondreamers_awake_n3",
                              tissue_class == "cortex" ~ "Heiss2004_Table1",
                              !is.na(value_heiss2004) ~ "Heiss2004_Table1",
                              TRUE ~ "Heiss1985_nondreamers_awake_n3"),
    cmrglc_n_datasets = as.integer(!is.na(cmrglc)),
    species = "human")
## Extended (SENSITIVITY) fill: a secondary dataset fills groups the primary lacks, after a multiplicative
## scale match k = mean(primary) / mean(secondary) over the groups BOTH have (any tissue class). Chugani adult
## values run ~25-30 % below Heiss for the same structures, so a raw fill would be mis-ranked.
scale_k <- function(primary, secondary) {
  ok <- !is.na(primary) & !is.na(secondary); stopifnot(sum(ok) >= 3)
  c(k = mean(primary[ok]) / mean(secondary[ok]), n = sum(ok))
}
k_chug <- scale_k(cm_human$cmrglc, cm_human$value_chugani1987_adult)
cm_human <- cm_human %>% mutate(
  cmrglc_extended = coalesce(cmrglc, k_chug[["k"]] * value_chugani1987_adult),
  cmrglc_extended_source = case_when(!is.na(cmrglc) ~ cmrglc_source,
                                     !is.na(value_chugani1987_adult) ~ sprintf("Chugani1987_adult_n7 x %.3f (scale-matched over %d groups)", k_chug[["k"]], k_chug[["n"]]),
                                     TRUE ~ NA_character_))
cm_macaque <- base %>% left_join(cm_ken, "group_id") %>% left_join(cm_pal, "group_id") %>%
  left_join(cm_noda, "group_id") %>% left_join(cm_meg, "group_id") %>%
  rowwise() %>%
  mutate(cmrglc_n_datasets = sum(!is.na(c(value_kennedy1978, value_palombo1990_ctrl))),
         cmrglc = if (cmrglc_n_datasets > 0) mean(c(value_kennedy1978, value_palombo1990_ctrl), na.rm = TRUE) else NA_real_,
         cmrglc_source = paste(c("Kennedy1978_conscious_n7", "Palombo1990_controls_n4")[!is.na(c(value_kennedy1978, value_palombo1990_ctrl))], collapse = "+"),
         cmrglc_source = ifelse(cmrglc_source == "", NA_character_, cmrglc_source)) %>%
  ungroup() %>% mutate(species = "macaque")
# extended (SENSITIVITY): primary (2-DG conscious), else Noda (PET conscious), else Meguro (PET, anaesthetised BABOON),
# each scale-matched to the primary over overlapping groups (see scale_k above)
k_noda <- scale_k(cm_macaque$cmrglc, cm_macaque$value_noda2002_young)
k_meg  <- scale_k(cm_macaque$cmrglc, cm_macaque$value_meguro1999_baboon)
cm_macaque <- cm_macaque %>% mutate(
  cmrglc_extended = coalesce(cmrglc, k_noda[["k"]] * value_noda2002_young, k_meg[["k"]] * value_meguro1999_baboon),
  cmrglc_extended_source = case_when(!is.na(cmrglc) ~ cmrglc_source,
                                     !is.na(value_noda2002_young) ~ sprintf("Noda2002_young_conscious_PET x %.3f (scale-matched over %d groups)", k_noda[["k"]], k_noda[["n"]]),
                                     !is.na(value_meguro1999_baboon) ~ sprintf("Meguro1999_BABOON_anaesthetised_PET x %.3f (scale-matched over %d groups)", k_meg[["k"]], k_meg[["n"]]),
                                     TRUE ~ NA_character_),
  cmrglc_extended_is_meguro = is.na(cmrglc) & is.na(value_noda2002_young) & !is.na(value_meguro1999_baboon))

cm_rat <- base %>% left_join(cm_sok, "group_id") %>%
  mutate(cmrglc = value_sokoloff1977_30min,
         cmrglc_source = ifelse(is.na(cmrglc), NA_character_, "Sokoloff1977_Table3_conscious_30min_n10"),
         cmrglc_n_datasets = as.integer(!is.na(cmrglc)),
         cmrglc_extended = cmrglc, cmrglc_extended_source = cmrglc_source,   # no secondary rat dataset: extended = primary
         species = "rat")

## ---- 4. Panel: one row per group x species, with within-species z-scores ---
zs <- function(x) if (sum(!is.na(x)) >= 3) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE) else NA_real_
zvars <- c(class_names, "cmrglc", "cmrglc_extended")
panel_species <- function(cm, rec, units) {
  p <- cm %>% left_join(rec, "group_id") %>%
    mutate(receptor_units = ifelse(is.na(n_member_areas), NA_character_, units),
           n_member_areas = replace_na(n_member_areas, 0L),
           has_receptors = n_member_areas > 0, has_cmrglc = !is.na(cmrglc), has_cmrglc_extended = !is.na(cmrglc_extended))
  for (v in zvars) p[[paste0(v, "_z")]] <- zs(p[[v]])                     # z over all groups of the species that have the variable
  cz <- p$tissue_class == "cortex"
  p$cmrglc_z_cortex <- NA_real_; p$cmrglc_z_cortex[cz] <- zs(p$cmrglc[cz])  # z over cortical groups only
  p
}
ph <- panel_species(cm_human, rec_human, "fmol/mg protein (Zilles & Palomero-Gallagher 2017, all layers)")
ph$class_ACh_incl_a4b2_z <- zs(ph$class_ACh_incl_a4b2); ph$total_15_z <- zs(ph$total_15)
pm <- panel_species(cm_macaque, rec_macaque, "fmol/neuron (Froudist-Walsh et al. 2023)")
pr <- panel_species(cm_rat, rec_rat, "fmol/mg protein (EBRAINS rat items, Palomero-Gallagher & Zilles; Wistar, n = 5)")
## z over the groups having each variable: receptor classes exist in 3 rat groups (2 for every class but iGlu), so the
## receptor z-scores are computed over n = 3 (iGlu) / n = 2 -> NA (zs() needs >= 3); only cmrglc_z (26 groups) is informative.

shared <- intersect(ph$group_id[ph$has_receptors & ph$has_cmrglc], pm$group_id[pm$has_receptors & pm$has_cmrglc])
shared_ext <- intersect(ph$group_id[ph$has_receptors & ph$has_cmrglc_extended], pm$group_id[pm$has_receptors & pm$has_cmrglc_extended])
shared_rec_only <- intersect(ph$group_id[ph$has_receptors], pm$group_id[pm$has_receptors])
## Rat sharing (CMRglc only; the rat receptor set covers 3 groups)
shared_cm_hr <- intersect(ph$group_id[ph$has_cmrglc], pr$group_id[pr$has_cmrglc])
shared_cm_mr <- intersect(pm$group_id[pm$has_cmrglc], pr$group_id[pr$has_cmrglc])
shared_cm_hm <- intersect(ph$group_id[ph$has_cmrglc], pm$group_id[pm$has_cmrglc])
shared_cm_hmr <- Reduce(intersect, list(shared_cm_hr, shared_cm_mr))
shared_rec_hr <- intersect(ph$group_id[ph$has_receptors], pr$group_id[pr$has_receptors])
panel <- bind_rows(ph, pm, pr) %>%
  mutate(shared_group_primary = group_id %in% shared,
         shared_group_extended = group_id %in% shared_ext,
         shared_receptors_both_species = group_id %in% shared_rec_only,
         shared_cmrglc_human_macaque = group_id %in% shared_cm_hm,
         shared_cmrglc_human_rat = group_id %in% shared_cm_hr,
         shared_cmrglc_macaque_rat = group_id %in% shared_cm_mr,
         shared_cmrglc_all_three = group_id %in% shared_cm_hmr) %>%
  arrange(factor(species, levels = c("human", "macaque", "rat")), factor(tissue_class, levels = c("cortex", "subcortex", "cerebellum", "white_matter")), group_id) %>%
  relocate(group_id, group_label, tissue_class, homology_confidence, species, n_member_areas, member_areas, receptor_units,
           all_of(class_names), class_ACh_incl_a4b2, total_15, all_of(paste0(class_names, "_z")), class_ACh_incl_a4b2_z, total_15_z,
           cmrglc, cmrglc_z, cmrglc_z_cortex, cmrglc_source, cmrglc_n_datasets,
           cmrglc_extended, cmrglc_extended_z, cmrglc_extended_source,
           has_receptors, has_cmrglc, has_cmrglc_extended, shared_group_primary, shared_group_extended, shared_receptors_both_species,
           shared_cmrglc_human_macaque, shared_cmrglc_human_rat, shared_cmrglc_macaque_rat, shared_cmrglc_all_three)
write_csv(panel, "tables/s1e/species_group_panel.csv", na = "")
invisible(file.copy("tables/s1e/species_group_panel.csv", "tables/s1e/human_macaque_group_panel.csv", overwrite = TRUE))  # legacy name

## Coverage table (one row per group)
coverage <- base %>% select(group_id, tissue_class, homology_confidence) %>%
  left_join(ph %>% transmute(group_id, human_n_zilles_areas = n_member_areas, human_cmrglc_source = cmrglc_source,
                             human_cmrglc_extended_source = cmrglc_extended_source), "group_id") %>%
  left_join(pm %>% transmute(group_id, macaque_n_fw2023_areas = n_member_areas, macaque_cmrglc_source = cmrglc_source,
                             macaque_cmrglc_n_datasets = cmrglc_n_datasets, macaque_cmrglc_extended_source = cmrglc_extended_source), "group_id") %>%
  left_join(pr %>% transmute(group_id, rat_n_ebrains_items = n_member_areas, rat_receptors_missing = receptors_missing,
                             rat_cmrglc_source = cmrglc_source), "group_id") %>%
  mutate(shared_group_primary = group_id %in% shared, shared_group_extended = group_id %in% shared_ext,
         shared_cmrglc_human_rat = group_id %in% shared_cm_hr, shared_cmrglc_macaque_rat = group_id %in% shared_cm_mr,
         shared_cmrglc_all_three = group_id %in% shared_cm_hmr,
         why_not_shared = case_when(
           shared_group_primary ~ "",
           human_n_zilles_areas == 0 & macaque_n_fw2023_areas == 0 ~ "no receptor data in either species",
           human_n_zilles_areas == 0 ~ "no human (Zilles 2017) receptor areas",
           macaque_n_fw2023_areas == 0 ~ "no macaque (FW2023) receptor areas",
           is.na(human_cmrglc_source) & is.na(macaque_cmrglc_source) ~ "no CMRglc in either species (Heiss / Kennedy+Palombo)",
           is.na(human_cmrglc_source) ~ "no human Heiss CMRglc group",
           is.na(macaque_cmrglc_source) ~ "no macaque Kennedy/Palombo CMRglc structure",
           TRUE ~ "other"))
write_csv(coverage, "tables/s1e/species_group_coverage.csv", na = "")
invisible(file.copy("tables/s1e/species_group_coverage.csv", "tables/s1e/human_macaque_group_coverage.csv", overwrite = TRUE))  # legacy name

## ---- 5. Spearman: within species, shared groups, cross-species concordance --
sp_row <- function(x, y, min_n = 5, min_n_rho = 3) {
  ## rho is reported from n >= min_n_rho (3) so that tiny sets (the rat receptor groups) are shown as what they are;
  ## p_naive only from n >= min_n (5). Below min_n_rho both are NA.
  ok <- is.finite(x) & is.finite(y); n <- sum(ok)
  if (n < min_n_rho) return(tibble(n = n, rho = NA_real_, p_naive = NA_real_))
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  tibble(n = n, rho = unname(ct$estimate), p_naive = if (n < min_n) NA_real_ else ct$p.value)
}
pred_h <- c(class_names, "class_ACh_incl_a4b2", "total_15"); pred_m <- class_names
within <- function(p, species, basis, preds, outcome = "cmrglc") {
  bind_rows(lapply(preds, function(v) bind_cols(tibble(analysis = "within_species", basis = basis, species = species,
                                                        variable = v, outcome = outcome), sp_row(p[[v]], p[[outcome]]))))
}
meguro_groups <- pm$group_id[pm$cmrglc_extended_is_meguro]
shared_ext_nomeg <- setdiff(shared_ext, meguro_groups)
subsets <- list(
  all_groups_with_data        = function(p) p,
  shared_groups_primary       = function(p) filter(p, group_id %in% shared),             # Kennedy/Palombo + Heiss only (= excluding Meguro)
  shared_primary_excl_low_confidence = function(p) filter(p, group_id %in% shared, homology_confidence != "low"),
  shared_extended_incl_noda_meguro_chugani = function(p) filter(p, group_id %in% shared_ext),
  shared_extended_excl_meguro = function(p) filter(p, group_id %in% shared_ext_nomeg)
)
sp_within <- bind_rows(lapply(names(subsets), function(b) {
  oc <- if (startsWith(b, "shared_extended")) "cmrglc_extended" else "cmrglc"
  bind_rows(within(subsets[[b]](ph), "human", b, pred_h, oc), within(subsets[[b]](pm), "macaque", b, pred_m, oc))
}))
## Rat within-species: all rat groups with both variables (receptors exist in 3 groups; classes other than iGlu in 2 -> NA)
sp_within_rat <- within(pr, "rat", "all_groups_with_data", pred_m, "cmrglc")
sp_within <- bind_rows(sp_within, sp_within_rat)
cross_vars <- c(class_names, "cmrglc")
cross <- function(gids, basis, outcome_var = "cmrglc") {
  bind_rows(lapply(cross_vars, function(v) {
    vv <- if (v == "cmrglc") outcome_var else v
    x <- ph[[paste0(vv, "_z")]][match(gids, ph$group_id)]; y <- pm[[paste0(vv, "_z")]][match(gids, pm$group_id)]
    bind_cols(tibble(analysis = "cross_species_concordance", basis = basis, species = "human_vs_macaque",
                     variable = paste0(vv, "_z"), outcome = "same variable, other species"), sp_row(x, y))
  }))
}
## Three-species CMRglc concordance (and rat receptor classes vs human over the 3 rat receptor groups)
cross2 <- function(pa, pb, gids, basis, label, vars) {
  bind_rows(lapply(vars, function(v) {
    x <- pa[[paste0(v, "_z")]][match(gids, pa$group_id)]; y <- pb[[paste0(v, "_z")]][match(gids, pb$group_id)]
    bind_cols(tibble(analysis = "cross_species_concordance", basis = basis, species = label,
                     variable = paste0(v, "_z"), outcome = "same variable, other species"), sp_row(x, y))
  }))
}
low <- base$group_id[base$homology_confidence == "low"]
sp_cross_rat <- bind_rows(
  cross2(ph, pm, shared_cm_hm, "shared_cmrglc_groups", "human_vs_macaque", "cmrglc"),
  cross2(ph, pr, shared_cm_hr, "shared_cmrglc_groups", "human_vs_rat", "cmrglc"),
  cross2(pm, pr, shared_cm_mr, "shared_cmrglc_groups", "macaque_vs_rat", "cmrglc"),
  cross2(ph, pm, setdiff(shared_cm_hm, low), "shared_cmrglc_groups_excl_low_confidence", "human_vs_macaque", "cmrglc"),
  cross2(ph, pr, setdiff(shared_cm_hr, low), "shared_cmrglc_groups_excl_low_confidence", "human_vs_rat", "cmrglc"),
  cross2(pm, pr, setdiff(shared_cm_mr, low), "shared_cmrglc_groups_excl_low_confidence", "macaque_vs_rat", "cmrglc"),
  cross2(ph, pm, shared_cm_hmr, "shared_cmrglc_groups_all_three_species", "human_vs_macaque", "cmrglc"),
  cross2(ph, pr, shared_cm_hmr, "shared_cmrglc_groups_all_three_species", "human_vs_rat", "cmrglc"),
  cross2(pm, pr, shared_cm_hmr, "shared_cmrglc_groups_all_three_species", "macaque_vs_rat", "cmrglc"),
  cross2(ph, pr, shared_rec_hr, "all_groups_with_receptors_in_both_species", "human_vs_rat", class_names)
)
sp_cross <- bind_rows(
  cross(shared, "shared_groups_primary"),
  cross(setdiff(shared, low), "shared_primary_excl_low_confidence"),
  cross(shared_ext, "shared_extended_incl_noda_meguro_chugani", "cmrglc_extended"),
  cross(shared_ext_nomeg, "shared_extended_excl_meguro", "cmrglc_extended"),
  cross(shared_rec_only, "all_groups_with_receptors_in_both_species"),   # receptor classes only; cmrglc row will have fewer n
  sp_cross_rat
)
spearman <- bind_rows(sp_within, sp_cross) %>%
  mutate(groups = case_when(basis == "shared_groups_primary" ~ paste(shared, collapse = ";"),
                            basis == "shared_primary_excl_low_confidence" ~ paste(setdiff(shared, low), collapse = ";"),
                            basis == "shared_extended_incl_noda_meguro_chugani" ~ paste(shared_ext, collapse = ";"),
                            basis == "shared_extended_excl_meguro" ~ paste(shared_ext_nomeg, collapse = ";"),
                            basis == "all_groups_with_receptors_in_both_species" & species == "human_vs_rat" ~ paste(shared_rec_hr, collapse = ";"),
                            basis == "all_groups_with_receptors_in_both_species" ~ paste(shared_rec_only, collapse = ";"),
                            basis == "shared_cmrglc_groups" & species == "human_vs_macaque" ~ paste(shared_cm_hm, collapse = ";"),
                            basis == "shared_cmrglc_groups" & species == "human_vs_rat" ~ paste(shared_cm_hr, collapse = ";"),
                            basis == "shared_cmrglc_groups" & species == "macaque_vs_rat" ~ paste(shared_cm_mr, collapse = ";"),
                            basis == "shared_cmrglc_groups_excl_low_confidence" & species == "human_vs_macaque" ~ paste(setdiff(shared_cm_hm, low), collapse = ";"),
                            basis == "shared_cmrglc_groups_excl_low_confidence" & species == "human_vs_rat" ~ paste(setdiff(shared_cm_hr, low), collapse = ";"),
                            basis == "shared_cmrglc_groups_excl_low_confidence" & species == "macaque_vs_rat" ~ paste(setdiff(shared_cm_mr, low), collapse = ";"),
                            basis == "shared_cmrglc_groups_all_three_species" ~ paste(shared_cm_hmr, collapse = ";"),
                            TRUE ~ "all groups of the species with both variables"),
         note = paste("descriptive; Spearman on group-level values; no spatial-autocorrelation null;",
                      "within-species: raw class value vs raw CMRglc of that species (rho is scale-free);",
                      "cross-species: human z vs macaque z of the same variable over the listed groups;",
                      "extended bases fill missing CMRglc from Chugani 1987 (human) / Noda 2002, Meguro 1999 baboon (macaque) after a",
                      sprintf("multiplicative scale match (k_chugani = %.3f over %d groups; k_noda = %.3f over %d; k_meguro = %.3f over %d);",
                              k_chug[["k"]], k_chug[["n"]], k_noda[["k"]], k_noda[["n"]], k_meg[["k"]], k_meg[["n"]]),
                      "rho reported from n >= 3 (rat receptor rows have n = 2-3 and are descriptive only), p_naive NA when n < 5;",
                      "receptor units differ between species (human and rat fmol/mg protein, macaque fmol/neuron); rat CMRglc = Sokoloff 1977 Table 3 30-min conscious arm;",
                      "rat hippocampal receptor items lack D1, GABAA/BZ, GABAB, M3, alpha4beta2 so every rat class but iGlu is NA for hippocampus_dentate"))
write_csv(spearman, "tables/s1e/species_group_spearman.csv", na = "")
invisible(file.copy("tables/s1e/species_group_spearman.csv", "tables/s1e/human_macaque_spearman.csv", overwrite = TRUE))  # legacy name

## ---- 6. Console summary ----------------------------------------------------
message("Shared groups (receptors + primary CMRglc in both species): ", length(shared), " -> ", paste(shared, collapse = ", "))
message("Shared groups, extended CMRglc (Noda/Meguro/Chugani fills): ", length(shared_ext), " -> ", paste(setdiff(shared_ext, shared), collapse = ", "), " added")
message("Groups with receptors in both species: ", length(shared_rec_only))
message("Human groups with receptors + CMRglc: ", sum(ph$has_receptors & ph$has_cmrglc), "; macaque: ", sum(pm$has_receptors & pm$has_cmrglc),
        "; rat: ", sum(pr$has_receptors & pr$has_cmrglc))
message("Rat: CMRglc groups ", sum(pr$has_cmrglc), "; receptor groups ", sum(pr$has_receptors), " (", paste(pr$group_id[pr$has_receptors], collapse = ", "), ")")
message("Rat receptors missing per receptor group: ", paste(rec_rat$group_id, "=", ifelse(rec_rat$receptors_missing == "", "none", rec_rat$receptors_missing), collapse = "; "))
message("CMRglc shared groups: human-macaque ", length(shared_cm_hm), "; human-rat ", length(shared_cm_hr), "; macaque-rat ", length(shared_cm_mr), "; all three ", length(shared_cm_hmr))
print(as.data.frame(sp_cross_rat %>% filter(variable == "cmrglc_z") %>% select(basis, species, n, rho, p_naive)), digits = 3)
print(as.data.frame(sp_cross %>% select(basis, variable, n, rho, p_naive)), digits = 3)
print(as.data.frame(sp_within %>% filter(basis %in% c("all_groups_with_data", "shared_groups_primary"), variable %in% paste0("class_", names(classes))) %>%
                      select(basis, species, variable, n, rho, p_naive)), digits = 3)
