setwd(local({ d <- normalizePath(getwd()); while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d); d }))  # repo root

## ============================================================
## Study s1e_4 (PROTOTYPE) — laminar prototype: receptor strata x BigBrain
## cell-body density x GLI, human (cross-species GLI where data exist).
##
## DESCRIPTIVE only (within-area contrasts, Spearman on a handful of cells; no
## model, no spatial null, no multiple-comparison control). Not part of
## run_all.R / output_manifest.csv.
##
## Inputs (synced by scripts/00_sync_evo_m1_inputs.R; metadata/evo_m1_input_manifest.csv):
##   data_raw/Zilles_PalomeroGallagher_2017_SupplementaryTable2.csv
##     15 receptors x 44 areas x 4 strata (supragranular I-III, granular IV,
##     infragranular V-VI, all layers), fmol/mg protein, means over 3 brains.
##     Strata were delimited cytoarchitectonically on neighbouring cell-body
##     sections (no depth percentages published); layer IV of areas 4 and 6 was
##     defined as a 3 % depth stripe below layer III (Zilles & Palomero-Gallagher 2017, Methods).
##   data_raw/Dickscheid_etal_2021_Area<X>.csv (12 areas)
##     BigBrain layer-wise cell-body density (# detected cells / 0.1 mm^3, all cell
##     types) and layer thickness (um, shrunken tissue) per layer I-VI, one donor.
##   data_raw/PalomeroGallagher_Zilles_2018_TableS3.csv  GLI % by stratum, areas 44/45, 7 primate species
##   data_raw/Zilles_etal_1986_Table1.csv / Table2.csv    posterior cingulate laminar % thickness and GLI, 17 primates
##   data_raw/Preuss_etal_1999_Fig2.csv, _MainText.csv, Preuss_Coleman_2002_Fig5.csv, _Fig7.csv,
##   data_raw/WongRiley_etal_1993_Text.csv, Hevner_WongRiley_1990_Table1.csv   (Preuss alignment, step 5)
##   metadata/anatomy/s1e_laminar_area_overlap.csv       curated area crosswalk (step 1)
##   Laminar receptor PROFILES (NOT items; read in place from the Evo-M1 root, sha256-verified):
##     <Evo-M1>/PalomeroGallagher_etal_2018/ebrains_downloads/profiles/*.tsv (+ profiles_manifest.csv)
##     8 human areas x 16 receptors, "example" profiles, cortical depth 0-100 % (0 = pial surface) vs fmol/mg.
##
## Outputs: tables/s1e/s1e_laminar_strata_panel.csv        area x stratum panel (primary + sensitivity sets)
##          tables/s1e/s1e_laminar_stratum_contrasts.csv   within-area granular/supra and infra/supra ratios
##          tables/s1e/s1e_laminar_spearman.csv            receptor class vs cell density across area x stratum cells
##          tables/s1e/s1e_laminar_profile_rebinned.csv    profile means per stratum / BigBrain layer
##          tables/s1e/s1e_laminar_profile_vs_tableS2.csv  re-binned profile vs Table S2
##          tables/s1e/s1e_laminar_gli_species.csv         cross-species GLI by stratum (44/45; posterior cingulate)
##          tables/s1e/s1e_laminar_gli_rank_concordance.csv Spearman of species rank orders between strata
##          tables/s1e/s1e_laminar_preuss_alignment.csv    human V1 layer 4A consistency table
## ============================================================

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })
dir.create("tables/s1e", showWarnings = FALSE, recursive = TRUE)

## ---- 0. Evo-M1 root (for the laminar profiles only) ----------------------
configured_root <- trimws(Sys.getenv("EVO_M1_TRAIT_DATA", unset = ""))
standard_root <- path.expand("~/Library/CloudStorage/OneDrive-AllenInstitute/Species/Evo-M1-Trait-Data")
evo_m1_root <- if (nzchar(configured_root)) configured_root else standard_root
profile_dir <- file.path(evo_m1_root, "PalomeroGallagher_etal_2018", "ebrains_downloads", "profiles")
if (!dir.exists(profile_dir)) stop("Laminar profile directory not found: ", profile_dir, call. = FALSE)

## ---- 1. Load --------------------------------------------------------------
zil <- read_csv("data_raw/Zilles_PalomeroGallagher_2017_SupplementaryTable2.csv", show_col_types = FALSE,
                col_types = cols(area = col_character(), .default = col_guess()))
stopifnot(nrow(zil) == 2640, length(unique(zil$area)) == 44)

xw <- read_csv("metadata/anatomy/s1e_laminar_area_overlap.csv", show_col_types = FALSE,
               col_types = cols(.default = col_character()))

bb_files <- list.files("data_raw", pattern = "^Dickscheid_etal_2021_Area.*\\.csv$", full.names = TRUE)
stopifnot(length(bb_files) == 12)
bb <- bind_rows(lapply(bb_files, function(f) {
  read_csv(f, show_col_types = FALSE) %>% mutate(bigbrain_area = sub("^Dickscheid_etal_2021_(Area.*)\\.csv$", "\\1", basename(f)))
}))
stopifnot(all(c("I","II","III","IV","V","VI","I-VI") %in% bb$layer), nrow(bb) == 12 * 7)

## ---- 2. Receptor classes per area x stratum ------------------------------
# Pharmacological classes as in s1e_1/s1e_3: iGlu = AMPA+NMDA+kainate; GABA = GABAA+GABAA/BZ+GABAB;
# ACh3 = M1+M2+M3 (the harmonised 3-receptor cholinergic class of s1e_3; alpha4beta2 kept separately);
# monoamineGPCR = alpha1+alpha2+5-HT1A+5-HT2+D1; total_15 = all 15; EI = iGlu/GABA.
classes <- list(iGlu = c("AMPA","NMDA","kainate"), GABA = c("GABAA","GABAA/BZ","GABAB"),
                ACh3 = c("M1","M2","M3"), monoamineGPCR = c("alpha1","alpha2","5-HT1A","5-HT2","D1"))
stopifnot(all(unlist(classes) %in% zil$receptor))
class_cols <- paste0("class_", names(classes))
strata3 <- c("supragranular", "granular", "infragranular")

zil_wide <- zil %>% select(area, stratum, receptor, density_fmol_mg_protein) %>%
  pivot_wider(names_from = receptor, values_from = density_fmol_mg_protein)
for (cl in names(classes)) zil_wide[[paste0("class_", cl)]] <- rowSums(zil_wide[, classes[[cl]], drop = FALSE])
zil_wide <- zil_wide %>% mutate(total_15 = rowSums(across(all_of(unique(zil$receptor)))),
                                EI_receptor_index = class_iGlu / class_GABA)
# area 24 granular is printed "--" (agranular; see S1E_RECEPTOR_PROTOTYPE.md) -> NA classes there; not in this panel.

## ---- 3. BigBrain layers -> Zilles strata (thickness-weighted) --------------
stratum_of_layer <- c(I = "supragranular", II = "supragranular", III = "supragranular", IV = "granular", V = "infragranular", VI = "infragranular")
bb_layers <- bb %>% filter(layer != "I-VI") %>%
  mutate(stratum = unname(stratum_of_layer[layer])) %>%
  group_by(bigbrain_area) %>%
  mutate(depth_frac_top = (cumsum(layer_thickness) - layer_thickness) / sum(layer_thickness),
         depth_frac_bottom = cumsum(layer_thickness) / sum(layer_thickness)) %>% ungroup()
bb_strata <- bb_layers %>% group_by(bigbrain_area) %>% mutate(cortex_um = sum(layer_thickness)) %>%
  group_by(bigbrain_area, stratum) %>%
  summarise(cell_density_per_0.1mm3 = sum(cell_density * layer_thickness) / sum(layer_thickness),
            stratum_thickness_um = sum(layer_thickness),
            stratum_depth_frac = sum(layer_thickness) / first(cortex_um),
            n_cells_detected = sum(n_cells), .groups = "drop")
bb_all <- bb %>% filter(layer == "I-VI") %>% transmute(bigbrain_area, cell_density_all_layers = cell_density, cortex_thickness_um = layer_thickness)
# check: thickness-weighted mean over the six layers reproduces the published I-VI value (same patches)
chk <- bb_layers %>% group_by(bigbrain_area) %>% summarise(w = sum(cell_density * layer_thickness) / sum(layer_thickness)) %>%
  inner_join(bb_all, by = "bigbrain_area") %>% mutate(rel = w / cell_density_all_layers - 1)
message("Thickness-weighted layer mean vs published I-VI density: max |rel diff| = ", signif(max(abs(chk$rel)), 3))

## ---- 4. Strata panel ------------------------------------------------------
# Panel sets: primary = exact/part_of BigBrain matches (V2d+V2v pooled to hOc2; area 4 agranular flagged);
#             unpooled = same but V2d and V2v as separate cells; sensitivity = + approximate matches (V3A, 37B, 10L).
rec_cols <- c(class_cols, "total_15", "EI_receptor_index", "alpha4beta2")
pairs <- xw %>% filter(bigbrain_match != "none", zilles2017_area != "") %>%
  select(zilles2017_area, bigbrain_dickscheid_area, bigbrain_match, in_strata_panel) %>%
  separate_rows(bigbrain_dickscheid_area, sep = ";")

make_panel <- function(pairs_sub, pool_v2 = TRUE, set_name) {
  ps <- pairs_sub %>% mutate(area_key = if (pool_v2) ifelse(zilles2017_area %in% c("V2d", "V2v"), "V2d+V2v", zilles2017_area) else zilles2017_area)
  # receptor side: unweighted mean of the Zilles areas pooled into one key (only V2d+V2v)
  r <- zil_wide %>% filter(stratum %in% strata3) %>% select(area, stratum, all_of(rec_cols)) %>%
    inner_join(ps %>% distinct(zilles2017_area, area_key, bigbrain_match, in_strata_panel), by = c(area = "zilles2017_area")) %>%
    group_by(area_key, stratum) %>%
    summarise(zilles_areas_pooled = paste(sort(unique(area)), collapse = "+"), bigbrain_match = first(bigbrain_match),
              in_strata_panel = first(in_strata_panel), across(all_of(rec_cols), mean), .groups = "drop")
  # cell side: thickness-weighted mean when a key maps to two BigBrain areas (37B -> FG3;FG4)
  b <- bb_strata %>% inner_join(ps %>% distinct(area_key, bigbrain_dickscheid_area), by = c(bigbrain_area = "bigbrain_dickscheid_area"),
                                relationship = "many-to-many") %>%   # hOc2 serves V2d and V2v in the unpooled set
    group_by(area_key, stratum) %>%
    summarise(bigbrain_dickscheid_area = paste(sort(unique(bigbrain_area)), collapse = ";"),
              cell_density_per_0.1mm3 = sum(cell_density_per_0.1mm3 * stratum_thickness_um) / sum(stratum_thickness_um),
              stratum_thickness_um = mean(stratum_thickness_um), stratum_depth_frac = mean(stratum_depth_frac), .groups = "drop")
  r %>% inner_join(b, by = c("area_key", "stratum")) %>% rename(area = area_key) %>% mutate(panel_set = set_name)
}
pairs_primary <- pairs %>% filter(bigbrain_match %in% c("exact", "part_of"))
pairs_sens    <- pairs
panel <- bind_rows(make_panel(pairs_primary, TRUE,  "primary"),
                   make_panel(pairs_primary, FALSE, "unpooled_V2"),
                   make_panel(pairs_sens,    TRUE,  "sensitivity_incl_approximate"))

# GLI (human, Table S3 mean over brains) for area 44 / 45; posterior cingulate GLI has no BigBrain partner (kept in gli tables only)
gli_s3 <- read_csv("data_raw/PalomeroGallagher_Zilles_2018_TableS3.csv", show_col_types = FALSE)
gli_h <- gli_s3 %>% filter(Species == "Homo sapiens") %>%
  group_by(area_as_published) %>%
  summarise(supragranular = mean(GLI_pct_mean_supragranular), granular = mean(GLI_pct_mean_granular),
            infragranular = mean(GLI_pct_mean_infragranular), n_brains = n(), .groups = "drop") %>%
  pivot_longer(all_of(strata3), names_to = "stratum", values_to = "GLI_pct_human_mean") %>%
  mutate(area = sub("^area", "", area_as_published), GLI_n_brains = n_brains) %>% select(area, stratum, GLI_pct_human_mean, GLI_n_brains)

panel <- panel %>% left_join(gli_h, by = c("area", "stratum")) %>%
  mutate(across(all_of(c(class_cols, "total_15")), ~ .x / cell_density_per_0.1mm3, .names = "{.col}_per_celldensity_ratio"),
         agranular_flag = area == "4",
         stratum = factor(stratum, levels = strata3)) %>%
  group_by(panel_set, area) %>%
  mutate(across(all_of(c(class_cols, "total_15", "EI_receptor_index", "cell_density_per_0.1mm3")),
                ~ if (sd(.x) > 0) (.x - mean(.x)) / sd(.x) else NA_real_, .names = "{.col}_z_within_area")) %>%
  ungroup() %>% arrange(panel_set, area, stratum) %>%
  relocate(panel_set, area, bigbrain_dickscheid_area, bigbrain_match, in_strata_panel, agranular_flag, stratum,
           all_of(class_cols), total_15, EI_receptor_index, alpha4beta2,
           cell_density_per_0.1mm3, stratum_thickness_um, stratum_depth_frac, GLI_pct_human_mean, GLI_n_brains)
# Units note (carried in the definitions block of the MD note): *_per_celldensity_ratio = (fmol receptor / mg protein) / (cells / 0.1 mm^3);
# a ratio of two differently normalised quantities (protein mass vs tissue volume, different brains, different shrinkage), NOT receptors per cell.
write_csv(panel, "tables/s1e/s1e_laminar_strata_panel.csv", na = "")
message("Panel rows: ", nrow(panel), " (primary cells n = ", sum(panel$panel_set == "primary"), ")")

## ---- 5. Within-area stratum contrasts -------------------------------------
contrast_vars <- c(class_cols, "total_15", "EI_receptor_index", "cell_density_per_0.1mm3")
contrasts <- panel %>% select(panel_set, area, agranular_flag, stratum, all_of(contrast_vars)) %>%
  pivot_longer(all_of(contrast_vars), names_to = "variable") %>%
  pivot_wider(names_from = stratum, values_from = value) %>%
  mutate(granular_over_supra = granular / supragranular, infra_over_supra = infragranular / supragranular,
         granular_over_infra = granular / infragranular) %>%
  rowwise() %>%
  mutate(rank_order = paste(c("supra", "gran", "infra")[order(-c(supragranular, granular, infragranular))], collapse = " > ")) %>% ungroup()
write_csv(contrasts, "tables/s1e/s1e_laminar_stratum_contrasts.csv", na = "")

## ---- 6. Spearman: receptor class vs cell density over area x stratum cells --
sp_row <- function(d, var, label, set) {
  d <- d[is.finite(d[[var]]) & is.finite(d$cell_density_per_0.1mm3), ]
  ct <- suppressWarnings(cor.test(d[[var]], d$cell_density_per_0.1mm3, method = "spearman", exact = FALSE))
  tibble(panel_set = set, cells = label, variable = var, n = nrow(d), rho = unname(ct$estimate), p_uncorrected = ct$p.value)
}
sp <- list()
for (set in unique(panel$panel_set)) {
  p <- panel %>% filter(panel_set == set)
  for (v in c(class_cols, "total_15", "EI_receptor_index")) {
    sp[[length(sp) + 1]] <- sp_row(p, v, "all area x stratum cells", set)
    sp[[length(sp) + 1]] <- sp_row(p %>% filter(!(agranular_flag & stratum == "granular")), v, "excluding agranular area 4 granular cell", set)
    for (s in strata3) sp[[length(sp) + 1]] <- sp_row(p %>% filter(stratum == s), v, paste("within stratum:", s, "(across areas)"), set)
    # within-area (across the 3 strata) rank agreement: mean of per-area Spearman
    wa <- p %>% group_by(area) %>% summarise(r = suppressWarnings(cor(.data[[v]], cell_density_per_0.1mm3, method = "spearman")), .groups = "drop")
    sp[[length(sp) + 1]] <- tibble(panel_set = set, cells = "within area (3 strata), mean of per-area rho", variable = v,
                                   n = nrow(wa), rho = mean(wa$r, na.rm = TRUE), p_uncorrected = NA_real_)
  }
}
spearman <- bind_rows(sp)
write_csv(spearman, "tables/s1e/s1e_laminar_spearman.csv", na = "")

## ---- 7. Laminar receptor profiles: parse, verify, re-bin -------------------
pm <- read_csv(file.path(evo_m1_root, "PalomeroGallagher_etal_2018", "ebrains_downloads", "profiles_manifest.csv"), show_col_types = FALSE)
stopifnot(nrow(pm) == 128)
sha_local <- vapply(file.path(profile_dir, pm$local_file), function(f) strsplit(system2("shasum", c("-a", "256", shQuote(f)), stdout = TRUE), " ")[[1]][1], character(1))
if (!all(sha_local == pm$sha256)) stop("sha256 mismatch for laminar profile file(s): ", paste(pm$local_file[sha_local != pm$sha256], collapse = ", "), call. = FALSE)
message("128 laminar profile TSVs verified against profiles_manifest.csv (sha256).")

profile_area <- c("PalomeroGallagher_etal_2018_Area3bv1.0" = "3b", "PalomeroGallagher_etal_2018_Area4pv1.0" = "4p",
                  "PalomeroGallagher_etal_2018_Area7Av1.0" = "7A", "PalomeroGallagher_etal_2018_dorsalpartofArea44v1.0" = "44d",
                  "PalomeroGallagher_etal_2018_Area45v1.0" = "45", "PalomeroGallagher_etal_2018_Area47v1.0" = "47",
                  "PalomeroGallagher_etal_2018_AreaPFmv1.0" = "PFm", "PalomeroGallagher_etal_2018_AreahOc1v1.0" = "hOc1")
# receptor names in the TSV folder -> Table S2 names
rec_map <- c("5-HT1A" = "5-HT1A", "5-HT2" = "5-HT2", AMPA = "AMPA", BZ = "GABAA/BZ", D1 = "D1", GABAA = "GABAA", GABAB = "GABAB",
             M1 = "M1", M2 = "M2", M3 = "M3", NMDA = "NMDA", alpha1 = "alpha1", alpha2 = "alpha2", alpha4beta2 = "alpha4beta2",
             kainate = "kainate", mGlu2_3 = NA_character_)
# Zilles 2017 area partner of each profile area (overlap table)
prof_zilles <- c("3b" = "3b", "4p" = "4", "7A" = NA, "44d" = "44", "45" = "45", "47" = "47", "PFm" = "PFm", "hOc1" = "V1")
# BigBrain boundary source for depth re-binning: hOc1 -> AreahOc1, 44d -> Area44, 4p -> Area4a (different subdivision; flagged)
prof_bb <- c("hOc1" = "AreahOc1", "44d" = "Area44", "4p" = "Area4a")

prof <- bind_rows(lapply(seq_len(nrow(pm)), function(i) {
  t <- read_tsv(file.path(profile_dir, pm$local_file[i]), show_col_types = FALSE, col_types = cols(.default = col_character()))
  stopifnot(names(t)[1] == "cortical depth (value)", names(t)[3] == "receptor density (value)")
  t <- t[!is.na(t[[1]]) & t[[1]] != "", ]
  tibble(profile_area = unname(profile_area[pm$item_name[i]]), receptor_folder = pm$receptor_from_url[i],
         depth_pct = as.numeric(t[[1]]), density = as.numeric(t[[3]]), depth_unit = t[[2]][1], density_unit = t[[4]][1],
         local_file = pm$local_file[i])
}))
stopifnot(all(prof$depth_unit == "%"), all(prof$density_unit == "fmol/mg"), !anyNA(prof$density))
n_samples <- prof %>% count(profile_area, receptor_folder)
stopifnot(all(n_samples$n == 101))
# the 3b "M3" file is a byte copy of the 3b M2 file (EBRAINS upload duplicate): drop it from the analysis
dup_3b_m3 <- pm$local_file[grepl("^3b_pr_examples_M3__", pm$local_file)]
prof <- prof %>% mutate(receptor = unname(rec_map[receptor_folder]),
                        receptor = ifelse(local_file %in% dup_3b_m3, NA_character_, receptor))

# (a) all-layers mean of each profile (101 equidistant samples 0..100 %)
prof_all <- prof %>% group_by(profile_area, receptor_folder, receptor, local_file) %>%
  summarise(profile_mean_all = mean(density), .groups = "drop")

# (b) stratum means using BigBrain cumulative layer-depth fractions (hOc1, 44d; 4p via 4a flagged)
bb_bounds <- bb_layers %>% select(bigbrain_area, layer, stratum, depth_frac_top, depth_frac_bottom)
bin_profile <- function(area_key) {
  bbk <- prof_bb[[area_key]]
  bounds <- bb_bounds %>% filter(bigbrain_area == bbk)
  p <- prof %>% filter(profile_area == area_key, !is.na(receptor))
  # assign each depth sample to a BigBrain layer (sample at depth d % lies in layer with top <= d/100 < bottom; d = 100 -> VI)
  p$layer <- sapply(p$depth_pct / 100, function(d) { l <- bounds$layer[d >= bounds$depth_frac_top & d < bounds$depth_frac_bottom]; if (length(l) == 0) "VI" else l })
  p$stratum <- unname(stratum_of_layer[p$layer])
  by_layer <- p %>% group_by(profile_area, receptor, layer) %>% summarise(mean_density = mean(density), n_samples = n(), .groups = "drop") %>%
    mutate(bin_type = "bigbrain_layer", bin = layer) %>% select(-layer)
  by_str <- p %>% group_by(profile_area, receptor, stratum) %>% summarise(mean_density = mean(density), n_samples = n(), .groups = "drop") %>%
    mutate(bin_type = "stratum_bigbrain_bounds", bin = stratum) %>% select(-stratum)
  out <- bind_rows(by_layer, by_str) %>% mutate(boundary_source = bbk)
  if (area_key == "4p") {
    # Zilles 2017 rule for agranular areas 4/6: layer IV = 3 % of depth below the layer III border
    top4 <- bounds$depth_frac_bottom[bounds$layer == "III"]
    p2 <- p %>% mutate(d = depth_pct / 100,
                       stratum = case_when(d < top4 ~ "supragranular", d < top4 + 0.03 ~ "granular", TRUE ~ "infragranular"))
    out <- bind_rows(out, p2 %>% group_by(profile_area, receptor, stratum) %>%
                       summarise(mean_density = mean(density), n_samples = n(), .groups = "drop") %>%
                       mutate(bin_type = "stratum_zilles2017_3pct_rule", bin = stratum, boundary_source = paste0(bbk, " layer III border + 3 % stripe")) %>% select(-stratum))
  }
  out
}
prof_bins <- bind_rows(lapply(names(prof_bb), bin_profile))
prof_rebinned <- bind_rows(prof_all %>% transmute(profile_area, receptor, bin_type = "all_layers_mean", bin = "all_layers", mean_density = profile_mean_all,
                                                  n_samples = 101L, boundary_source = NA_character_),
                           prof_bins) %>% arrange(profile_area, receptor, bin_type, bin)
write_csv(prof_rebinned, "tables/s1e/s1e_laminar_profile_rebinned.csv", na = "")

# (c) compare with Table S2 (same lab, 3-brain means; profiles are single "example" profiles)
s2_long <- zil %>% transmute(zilles_area = area, receptor, stratum, tableS2_density = density_fmol_mg_protein, tableS2_sd = sd)
# The EBRAINS fingerprint datasets (3-brain means, e.g. 44d AMPA 492 / kainate 498 / NMDA 964) equal Table S2, so the
# profiles ("pr_examples": one brain, one section) are compared with the 3-brain mean and its between-brain SD.
cmp <- prof_rebinned %>% filter(bin_type %in% c("all_layers_mean", "stratum_bigbrain_bounds", "stratum_zilles2017_3pct_rule"), !is.na(receptor)) %>%
  mutate(zilles_area = unname(prof_zilles[profile_area]), stratum = ifelse(bin == "all_layers", "all_layers", bin)) %>%
  filter(!is.na(zilles_area)) %>%
  inner_join(s2_long, by = c("zilles_area", "receptor", "stratum")) %>%
  mutate(ratio_profile_over_S2 = mean_density / tableS2_density, pct_diff = 100 * (mean_density / tableS2_density - 1),
         z_vs_S2_sd = (mean_density - tableS2_density) / tableS2_sd) %>%
  arrange(profile_area, bin_type, stratum, receptor)
write_csv(cmp, "tables/s1e/s1e_laminar_profile_vs_tableS2.csv", na = "")
cmp_summary <- cmp %>% group_by(profile_area, zilles_area, bin_type, stratum) %>%
  summarise(n_receptors = n(), median_pct_diff = median(pct_diff), median_abs_pct_diff = median(abs(pct_diff)),
            n_within_10pct = sum(abs(pct_diff) <= 10), n_within_25pct = sum(abs(pct_diff) <= 25),
            n_within_1sd = sum(abs(z_vs_S2_sd) <= 1, na.rm = TRUE), n_within_2sd = sum(abs(z_vs_S2_sd) <= 2, na.rm = TRUE),
            spearman_rho_across_receptors = suppressWarnings(cor(mean_density, tableS2_density, method = "spearman")), .groups = "drop")
write_csv(cmp_summary, "tables/s1e/s1e_laminar_profile_vs_tableS2_summary.csv", na = "")
print(as.data.frame(cmp_summary), digits = 3)

## ---- 8. Cross-species GLI ---------------------------------------------------
# 8a. Areas 44/45, 7 species (Table S3): species means over brains (human n=4; macaque n=2/3), then rank concordance between strata
gli_sp <- gli_s3 %>% group_by(area = sub("^area", "", area_as_published), Species) %>%
  summarise(n_brains = n(), supragranular = mean(GLI_pct_mean_supragranular), granular = mean(GLI_pct_mean_granular),
            infragranular = mean(GLI_pct_mean_infragranular), all_layers = mean(GLI_pct_mean_all_layers), .groups = "drop") %>%
  mutate(source = "PalomeroGallagher_Zilles_2018_TableS3", granular_over_supra = granular / supragranular,
         infra_over_supra = infragranular / supragranular, region = "ventrolateral prefrontal (Broca)")
# 8b. Posterior cingulate, 17 species (Zilles 1986 Table 2: O = outer main lamina (II-III), G = granular (areas 23, 31 only), I = inner main (V-VI))
z86 <- read_csv("data_raw/Zilles_etal_1986_Table2.csv", show_col_types = FALSE)
stopifnot(nrow(z86) == 17)
z86_long <- z86 %>% pivot_longer(-Species, names_to = "key", values_to = "GLI") %>%
  extract(key, into = c("area", "sublayer"), regex = "^area(\\d+)_(outer_main|granular|inner_main)_GLI$") %>%
  mutate(stratum = recode(sublayer, outer_main = "supragranular", granular = "granular", inner_main = "infragranular")) %>%
  select(-sublayer) %>% pivot_wider(names_from = stratum, values_from = GLI) %>%
  mutate(n_brains = NA_integer_, all_layers = NA_real_, source = "Zilles_etal_1986_Table2",
         granular_over_supra = granular / supragranular, infra_over_supra = infragranular / supragranular, region = "posterior cingulate")
gli_species <- bind_rows(gli_sp, z86_long) %>% relocate(source, region, area, Species, n_brains, supragranular, granular, infragranular, all_layers)
write_csv(gli_species, "tables/s1e/s1e_laminar_gli_species.csv", na = "")

rank_conc <- function(d, a, src) {
  d <- d %>% filter(area == a)
  prs <- list(c("supragranular", "granular"), c("supragranular", "infragranular"), c("granular", "infragranular"))
  bind_rows(lapply(prs, function(p) {
    ok <- is.finite(d[[p[1]]]) & is.finite(d[[p[2]]])
    if (sum(ok) < 4) return(NULL)
    ct <- suppressWarnings(cor.test(d[[p[1]]][ok], d[[p[2]]][ok], method = "spearman", exact = FALSE))
    tibble(source = src, area = a, stratum_a = p[1], stratum_b = p[2], n_species = sum(ok), rho = unname(ct$estimate), p_uncorrected = ct$p.value)
  }))
}
rank_tab <- bind_rows(rank_conc(gli_sp, "44", "PalomeroGallagher_Zilles_2018_TableS3"), rank_conc(gli_sp, "45", "PalomeroGallagher_Zilles_2018_TableS3"),
                      lapply(c("23", "31", "29", "30"), function(a) rank_conc(z86_long, a, "Zilles_etal_1986_Table2")))
# laminar pattern conservation: per species rank of strata, and how many species share the modal order
pattern <- gli_species %>% filter(is.finite(granular)) %>% rowwise() %>%
  mutate(rank_order = paste(c("supra", "gran", "infra")[order(-c(supragranular, granular, infragranular))], collapse = " > ")) %>% ungroup() %>%
  count(source, area, rank_order) %>% group_by(source, area) %>% mutate(n_species_with_granular = sum(n), share = n / sum(n)) %>% ungroup() %>% arrange(source, area, -n)
write_csv(rank_tab, "tables/s1e/s1e_laminar_gli_rank_concordance.csv", na = "")
write_csv(pattern, "tables/s1e/s1e_laminar_gli_stratum_order_by_species.csv", na = "")
print(as.data.frame(rank_tab), digits = 3); print(as.data.frame(pattern), digits = 3)

## ---- 9. Preuss alignment: human V1 layer 4A ---------------------------------
p99f2 <- read_csv("data_raw/Preuss_etal_1999_Fig2.csv", show_col_types = FALSE, col_types = cols(.default = col_character()))
p99mt <- read_csv("data_raw/Preuss_etal_1999_MainText.csv", show_col_types = FALSE, col_types = cols(.default = col_character()))
pc5   <- read_csv("data_raw/Preuss_Coleman_2002_Fig5.csv", show_col_types = FALSE, col_types = cols(.default = col_character()))
pc7   <- read_csv("data_raw/Preuss_Coleman_2002_Fig7.csv", show_col_types = FALSE, col_types = cols(.default = col_character()))
wr93  <- read_csv("data_raw/WongRiley_etal_1993_Text.csv", show_col_types = FALSE, col_types = cols(.default = col_character()))
h90   <- read_csv("data_raw/Hevner_WongRiley_1990_Table1.csv", show_col_types = FALSE, col_types = cols(.default = col_character()))

pick <- function(d, src, ...) d %>% filter(...) %>% transmute(source = src, Species, layer, marker, trait, state, printed_wording, source_locator)
al <- bind_rows(
  pick(p99f2, "Preuss_etal_1999_Fig2", trait == "L4A_CO_dense_band"),
  pick(p99mt, "Preuss_etal_1999_MainText", trait %in% c("L4A_CO_dense_band", "L4B_CO_band")),
  pick(pc5,  "Preuss_Coleman_2002_Fig5", trait == "L4A_Cat301_mesh"),
  pick(pc7,  "Preuss_Coleman_2002_Fig7", trait == "L4A_CB_pattern"),
  wr93 %>% filter(text_row %in% c("40", "58")) %>%
    transmute(source = "WongRiley_etal_1993_Text", Species, layer = layer_as_printed, marker = "CO histochemistry",
              trait = "IVa_demarcable_by_CO", state = "not discernable / not distinct", printed_wording = passage_as_printed,
              source_locator = paste0("text row ", text_row, ", p. ", page, " (age group: ", age_group_as_printed, ")")),
  wr93 %>% filter(text_row == "61") %>%
    transmute(source = "WongRiley_etal_1993_Text", Species, layer = layer_as_printed, marker = "CO histochemistry",
              trait = "IVcalpha_vs_IVcbeta_CO_OD", state = paste(value_printed, units_printed), printed_wording = passage_as_printed,
              source_locator = paste0("text row ", text_row, ", p. ", page)),
  bb_layers %>% filter(bigbrain_area == "AreahOc1", layer == "IV") %>%
    transmute(source = "Dickscheid_etal_2021_AreahOc1", Species, layer = "IV (single layer; IVa/IVb/IVc not subdivided)", marker = "cell bodies (BigBrain, Merker stain)",
              trait = "layer_IV_cell_density_per_0.1mm3", state = as.character(round(cell_density, 1)),
              printed_wording = paste0("layer IV thickness ", round(layer_thickness), " um = ", round(100 * (depth_frac_bottom - depth_frac_top), 1), " % of cortical depth; n_cells = ", n_cells),
              source_locator = "data_raw/Dickscheid_etal_2021_AreahOc1.csv"),
  h90 %>% transmute(source = "Hevner_WongRiley_1990_Table1", Species, layer = cortical_layer_as_published, marker = stain_method,
                    trait = paste0("CO_relOD_", compartment), state = paste(mean_printed, "+/-", dispersion_printed),
                    printed_wording = paste0(region_as_published, " (", condition, "); n_readings=", n_readings, ", n_animals=", n_animals),
                    source_locator = paste0("printed row ", printed_row, "; Table 1 quantifies puff/interpuff (layers 2-3) and layer 4C columns only - no layer 4A value"))
) %>% mutate(alignment_theme = case_when(
  grepl("CO_dense|IVa_demarcable", trait) ~ "human V1 4A: CO-dense band absent (Preuss 1999) vs IVa not demarcable by CO (Wong-Riley 1993)",
  grepl("Cat301|CB_pattern", trait) ~ "human-specific 4A compartments (Cat-301 mesh / calbindin territories; Preuss & Coleman 2002)",
  grepl("layer_IV_cell_density", trait) ~ "BigBrain hOc1 layer IV is one layer: 4A cannot be isolated in the cell-density item",
  grepl("L4B_CO", trait) ~ "hominoid 4B CO band (complement of the 4A loss)",
  TRUE ~ "macaque CO compartments (Hevner 1990: puffs/4C columns; no 4A row)"))
write_csv(al, "tables/s1e/s1e_laminar_preuss_alignment.csv", na = "")
message("Preuss alignment rows: ", nrow(al))
message("s1e_4 done.")
