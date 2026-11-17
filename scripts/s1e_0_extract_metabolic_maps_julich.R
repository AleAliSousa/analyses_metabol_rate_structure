# s1e_0 (PROTOTYPE, not in run_all.R / output_manifest.csv):
# parcellate cortical metabolic maps by Julich-Brain v3.1 areas, join the
# Zilles & Palomero-Gallagher 2017 receptor-area crosswalk, and sanity-check
# lobe aggregates against Heiss et al. 2004 lobe rCMRglc.
#
# Python step: helpers/parcellate_metabolic_maps_julich.py (numpy, scipy,
# pandas, nibabel, neuromaps). Interpreter selection, in order:
#   1. S1E_PYTHON environment variable
#   2. .venv-atlas/bin/python (the s1c virtual environment), if it imports neuromaps
#   3. ~/.claude-science/conda/envs/neuro-atlas/bin/python, if present
#   4. python3 on PATH
# Set S1E_SKIP_PYTHON=1 to reuse an existing data_intermediate/s1e_julich_area_metabolism.csv.
# See metadata/anatomy/S1E_PARCELLATION_PROTOTYPE.md for methods and limitations.
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })
setwd(local({ d <- normalizePath(getwd())
  while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d); d }))
stopifnot(file.exists("helpers/parcellate_metabolic_maps_julich.py"))
source("helpers/regional_robustness.R")  # write_analysis_csv()

python_imports_ok <- function(py) {
  if (!nzchar(py) || !file.exists(py)) return(FALSE)
  st <- suppressWarnings(system2(py, c("-c", shQuote("import neuromaps, nibabel, pandas, scipy, numpy")),
                                 stdout = FALSE, stderr = FALSE))
  identical(st, 0L)
}
find_python <- function() {
  cands <- c(Sys.getenv("S1E_PYTHON"),
             file.path(getwd(), ".venv-atlas", "bin", "python"),
             path.expand("~/.claude-science/conda/envs/neuro-atlas/bin/python"),
             Sys.which("python3"))
  cands <- cands[nzchar(cands)]
  for (p in cands) if (python_imports_ok(p)) return(p)
  stop("No Python interpreter with neuromaps + nibabel + pandas + scipy found. ",
       "Set S1E_PYTHON to one (e.g. a conda env with `pip install neuromaps nibabel pandas scipy`).")
}

area_csv <- "data_intermediate/s1e_julich_area_metabolism.csv"
if (Sys.getenv("S1E_SKIP_PYTHON", "0") != "1") {
  python <- find_python()
  message("s1e_0: running helper with ", python)
  status <- system2(python, shQuote("helpers/parcellate_metabolic_maps_julich.py"))
  if (status != 0L) stop("parcellate_metabolic_maps_julich.py failed (exit ", status, ")")
}
stopifnot(file.exists(area_csv))

# ---------------------------------------------------------------------------
area <- read.csv(area_csv, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(c("julich_region", "hemisphere", "n_vertices", "n_vertices_valid", "cmrglc", "cmr02",
                "cbf", "cbv", "gi_proxy", "castrillon_cmrglc", "n_voxels", "n_voxels_valid") %in% names(area)),
          !anyDuplicated(area[c("julich_region", "hemisphere")]))
cw <- read.csv("metadata/anatomy/s1e_zilles2017_julich_area_map.csv", stringsAsFactors = FALSE)
stopifnot(nrow(cw) == 44L, !anyDuplicated(cw$zilles_area),
          all(cw$match_relationship %in% c("exact", "part_of", "union", "approximate", "none")))
heiss <- read.csv("data_intermediate/heiss_2004_regions.csv", stringsAsFactors = FALSE) %>%
  select(anatomy_id, heiss_rcmrglc = rcmrglc_mean_umol_100g_min)
heiss_cortex <- heiss$heiss_rcmrglc[heiss$anatomy_id == "cerebral_cortex"]
stopifnot(length(heiss_cortex) == 1L)

measures <- c("cmrglc", "cmr02", "cbf", "cbv", "gi_proxy", "castrillon_cmrglc")
weight_col <- function(m) if (m == "castrillon_cmrglc") "n_voxels_valid" else "n_vertices_valid"

# Long crosswalk: one row per zilles_area x julich_region; every name must exist in the parcellation.
cw_long <- cw %>% filter(match_relationship != "none") %>%
  separate_rows(julich_regions, sep = ";") %>% rename(julich_region = julich_regions) %>%
  mutate(julich_region = trimws(julich_region))
missing <- setdiff(cw_long$julich_region, area$julich_region)
if (length(missing)) stop("Crosswalk names absent from parcellation: ", paste(missing, collapse = "; "))

# Weighted aggregation of Julich regions into Zilles areas (unions pooled by valid vertices / voxels).
wmean <- function(x, w) { ok <- is.finite(x) & is.finite(w) & w > 0; if (!any(ok)) NA_real_ else sum(x[ok] * w[ok]) / sum(w[ok]) }
stopifnot_no_na_unions <- function(tbl) {
  bad <- tbl %>% filter(n_vertices_valid > 0, !is.finite(cmrglc))
  if (nrow(bad)) stop("Aggregation produced NA cmrglc despite valid vertices: ", paste(bad$zilles_area, collapse = ", "))
  tbl
}
# .wv/.wx copies avoid the summarised totals shadowing the per-region weights inside summarise().
zilles_area_tbl <- cw_long %>% inner_join(area, by = "julich_region", relationship = "many-to-many") %>%  # hOc2 -> V2d and V2v; 3 hemisphere rows per region
  mutate(.wv = n_vertices_valid, .wx = n_voxels_valid) %>%
  group_by(zilles_area, lobe, heiss_lobe_anatomy_id, match_relationship, hemisphere) %>%
  summarise(julich_regions = paste(julich_region, collapse = ";"), n_julich_regions = n(),
            across(all_of(setdiff(measures, "castrillon_cmrglc")), ~ wmean(.x, .wv)),
            castrillon_cmrglc = wmean(castrillon_cmrglc, .wx),
            n_vertices = sum(n_vertices), n_vertices_valid = sum(n_vertices_valid),
            n_voxels = sum(n_voxels), n_voxels_valid = sum(n_voxels_valid), .groups = "drop") %>%
  stopifnot_no_na_unions() %>%
  mutate(frac_vertices_excluded = 1 - n_vertices_valid / n_vertices,
         few_vertices_flag = n_vertices < 200, few_voxels_flag = n_voxels < 20,
         cmrglc_rel_cortex = cmrglc / area$cmrglc_cortex_mean[match(hemisphere, area$hemisphere)],
         castrillon_rel_cortex = castrillon_cmrglc / area$castrillon_cmrglc_cortex_mean[match(hemisphere, area$hemisphere)]) %>%
  arrange(match(zilles_area, cw$zilles_area), match(hemisphere, c("L", "R", "B")))
unmatched <- cw %>% filter(match_relationship == "none") %>%
  select(zilles_area, lobe, heiss_lobe_anatomy_id, match_relationship, note)
write_analysis_csv(zilles_area_tbl, "tables/s1e/s1e_zilles_area_metabolism.csv")
write_analysis_csv(unmatched, "tables/s1e/s1e_zilles_area_unmatched.csv")

# ---------------------------------------------------------------------------
# Lobe aggregation. Basis A: the Zilles-2017 areas with exact/part_of/union matches
# (approximate matches added in basis A2). Basis B: all Julich v3.1 cortical areas with
# an unambiguous lobe from the gyrus token in the Julich name (allocortex, cingulate,
# parieto-occipital sulcus and cross-lobe GapMaps excluded).
gyrus_lobe <- c(PostCG = "parietal", PostCS = "parietal", IPL = "parietal", SPL = "parietal", IPS = "parietal",
                POperc = "parietal",
                IFG = "frontal", OFC = "frontal", FPole = "frontal", PreCG = "frontal", SFS = "frontal",
                MFG = "frontal", SFG = "frontal", IFS = "frontal", "IFS, PreCS" = "frontal",
                "SMA, mesial SFG" = "frontal", "preSMA, mesial SFG" = "frontal", "Frontal Operculum" = "frontal",
                "MFG, IFG" = "frontal",
                LOC = "occipital", "V1, 17, CalcS" = "occipital", "V2, 18" = "occipital", LingG = "occipital",
                Cuneus = "occipital",
                FusG = "temporal", HESCHL = "temporal", STG = "temporal", "STG, SMG" = "temporal", STS = "temporal",
                PhG = "temporal", CoS = "temporal", OTS = "temporal",
                Insula = "insular")
gapmap_lobe <- c("GapMap-Frontal-I.1" = "frontal", "GapMap-Frontal-I.2" = "frontal", "GapMap-Frontal-II" = "frontal")
julich_lobe <- area %>% distinct(julich_label, julich_region) %>%
  mutate(token = sub("^.*\\(([^()]*)\\)\\s*$", "\\1", julich_region),
         lobe = ifelse(grepl("^GapMap", julich_region), unname(gapmap_lobe[julich_region]),
                       unname(gyrus_lobe[token])),
         lobe = ifelse(grepl("Hippocamp", julich_region) | grepl("Cingulate|ACC", julich_region) | token %in% c("POS", "ACC", "sACC", "pACC"), NA, lobe))
write_analysis_csv(julich_lobe %>% select(-token), "tables/s1e/s1e_julich_region_lobe_assignment.csv")

lobe_ids <- c(frontal = "frontal_lobe", parietal = "parietal_lobe", temporal = "temporal_lobe",
              occipital = "occipital_lobe", insular = "insular_lobe")
aggregate_lobes <- function(tbl, basis) {
  tbl <- tbl %>% filter(hemisphere == "B", lobe %in% names(lobe_ids))
  bind_rows(lapply(measures, function(m) {
    w <- tbl[[weight_col(m)]]
    tbl %>% mutate(.w = w) %>% group_by(lobe) %>%
      summarise(n_areas = n(), weight_sum = sum(.w), value_native = wmean(.data[[m]], .w), .groups = "drop") %>%
      mutate(basis = basis, measure = m)
  })) %>%
    mutate(heiss_anatomy_id = unname(lobe_ids[lobe])) %>%
    left_join(heiss, by = c("heiss_anatomy_id" = "anatomy_id"))
}
basisA  <- aggregate_lobes(zilles_area_tbl %>% filter(match_relationship %in% c("exact", "part_of", "union")), "zilles44_strict_matches")
basisA2 <- aggregate_lobes(zilles_area_tbl, "zilles44_incl_approximate")
basisB  <- aggregate_lobes(area %>% select(-any_of("lobe")) %>% inner_join(julich_lobe %>% select(julich_label, lobe), by = "julich_label") %>%
                             filter(!is.na(lobe)), "all_julich_by_gyrus")
lobes <- bind_rows(basisA, basisA2, basisB) %>%
  group_by(basis, measure) %>%
  mutate(cortex_mean_native = sum(value_native * weight_sum) / sum(weight_sum),
         value_rel_lobes_mean = value_native / cortex_mean_native,
         # Raichle maps are in neuromaps' normalised units: express relative to the Heiss cortex mean
         value_scaled_to_heiss_cortex = ifelse(measure %in% c("castrillon_cmrglc", "gi_proxy"), NA_real_,
                                               value_rel_lobes_mean * heiss_cortex),
         heiss_rel_cortex = heiss_rcmrglc / heiss_cortex) %>% ungroup() %>%
  select(basis, measure, lobe, heiss_anatomy_id, n_areas, weight_sum, value_native, value_rel_lobes_mean,
         value_scaled_to_heiss_cortex, heiss_rcmrglc, heiss_rel_cortex)
summary_tbl <- lobes %>% group_by(basis, measure) %>%
  summarise(n_lobes = n(),
            pearson_r = if (n() >= 3) cor(value_native, heiss_rcmrglc) else NA_real_,
            spearman_rho = if (n() >= 3) cor(value_native, heiss_rcmrglc, method = "spearman") else NA_real_,
            rank_order_maps = paste(lobe[order(-value_native)], collapse = " > "),
            rank_order_heiss = paste(lobe[order(-heiss_rcmrglc)], collapse = " > "), .groups = "drop") %>%
  mutate(interpretation = "Descriptive sanity check on 4-5 lobe means; not a statistical test. Raichle values are relative (neuromaps units); Castrillon values are umol/100 g/min from a different cohort.")
write_analysis_csv(lobes, "tables/s1e/s1e_lobe_check_vs_heiss.csv")
write_analysis_csv(summary_tbl, "tables/s1e/s1e_lobe_check_vs_heiss_summary.csv")

coverage <- cw %>% count(match_relationship) %>% arrange(match(match_relationship, c("exact", "part_of", "union", "approximate", "none")))
write_analysis_csv(coverage, "tables/s1e/s1e_crosswalk_coverage.csv")
print(coverage)
print(as.data.frame(summary_tbl %>% select(basis, measure, n_lobes, pearson_r, spearman_rho)))
print(as.data.frame(lobes %>% filter(measure %in% c("cmrglc", "castrillon_cmrglc")) %>%
        select(basis, measure, lobe, n_areas, value_native, value_scaled_to_heiss_cortex, heiss_rcmrglc)))
