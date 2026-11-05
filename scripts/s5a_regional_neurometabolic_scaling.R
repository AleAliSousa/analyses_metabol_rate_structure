# ============================================================
# s5a_regional_neurometabolic_scaling.R   (Study s5: Coto-style analyses)
#
# QUESTION  Within the human brain, how does regional glucose use scale with
#           regional volume?  Coto et al. 2026 report a cross-species slope of
#           brain metabolic rate on brain mass of 0.87 (ants + 7 mammals);
#           Karbowski 2007 reports regional CMRgl exponents on brain mass of
#           about -0.15 across species.  Here the analysis unit is the Heiss
#           2004 region: budget_r = rCMRGlc_r x volume_r.
#
# INPUTS
#   data_intermediate/heiss_2004_regions.csv          (rCMRGlc, umol/100 g/min)
#   metadata/anatomy/s5_heiss_region_volume_map.csv   (volume term per Heiss region;
#        reuses the Table 1 pairings and the s4a Kochiyama lobe volumes)
#   data_raw/volumes_wide.csv (Homo sapiens row; Stephan/Zilles-collection volumes, mm3,
#        both hemispheres) and data_intermediate/s4a_endocranial_region_budgets.csv
#   Evo-M1 Karbowski__2007/Karbowski__2007_Table1.csv (published regional exponents)
#
# OUTPUTS  tables/s5/s5a_region_budgets.csv, s5a_scaling_fits.csv,
#          s5a_excluded_regions.csv, s5a_karbowski_comparison.csv
#          figs/s5/s5a_budget_vs_volume.{png,pdf}, s5a_rate_vs_volume.{png,pdf}
# ============================================================
setwd(local({ d <- normalizePath(getwd()); while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d); d }))  # repo root (portable; replaces hardcoded path -- see helpers/project_root.R)
source("helpers/model_diagnostics.R"); diag_stage("s5")  # model diagnostics -> figs/s5/diagnostics/

suppressPackageStartupMessages({ library(tidyverse); library(smatr); library(ggrepel) })
source("helpers/plot_settings.R")
source("helpers/read_heiss_rates.R")
source("helpers/s5_helpers.R")
dir.create("tables/s5", showWarnings = FALSE, recursive = TRUE)
dir.create("figs/s5", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 1. Rates (Heiss) and the s5 volume map
# ------------------------------------------------------------
heiss <- read_heiss_regions()
map   <- read_study_anatomy_map("metadata/anatomy/s5_heiss_region_volume_map.csv")
map$allocation_weight <- as.numeric(map$allocation_weight)
stopifnot(all(map$volume_set %in% c("primary", "kochiyama_lobes")))

# Resolve one rate per (volume term, analysis id); composites follow the map weights.
rates <- derive_study_heiss_rates(map) %>%
  select(source_region, anatomy_id, anatomy_group, rcmr_value, rcmr_sd, n_heiss_components) %>%
  # display labels follow the shared palette (helpers/plot_settings.R)
  mutate(anatomy_group = recode(anatomy_group,
                                "Midbrain" = "Mesencephalon",
                                "Neocortical grey matter" = "Cerebral cortex (global average)"))

# ------------------------------------------------------------
# 2. Volumes (mm3, both hemispheres)
# ------------------------------------------------------------
vw   <- read.csv("data_raw/volumes_wide.csv", check.names = FALSE, stringsAsFactors = FALSE)
homo <- vw[vw$Species == "Homo sapiens", , drop = FALSE]
stopifnot(nrow(homo) == 1)
s4a  <- read.csv("data_intermediate/s4a_endocranial_region_budgets.csv", stringsAsFactors = FALSE)
TOTAL_BRAIN_MM3 <- as.numeric(homo[["Total_brain_net_volume_Vol.mm3"]])

lookup_volume <- function(file, column, subtract, anatomy_id) {
  if (file == "data_raw/volumes_wide.csv") {
    v <- as.numeric(homo[[column]])
    if (!is.na(subtract) && nzchar(subtract)) {
      for (s in strsplit(subtract, ";")[[1]]) v <- v - as.numeric(homo[[trimws(s)]])
    }
    v
  } else if (file == "data_intermediate/s4a_endocranial_region_budgets.csv") {
    as.numeric(s4a[[column]][s4a$anatomy_id == anatomy_id]) * 1000  # cc -> mm3
  } else stop("unknown volume source: ", file)
}

vol_rows <- map %>%
  distinct(source_region, analysis_anatomy_id, volume_set, volume_source_file, volume_source_column,
           volume_subtract_columns, allocation_weight, volume_basis, volume_source_citation, relationship)
vol_rows$parent_volume_mm3 <- mapply(lookup_volume, vol_rows$volume_source_file, vol_rows$volume_source_column,
                                     vol_rows$volume_subtract_columns, vol_rows$analysis_anatomy_id)
# normalise split weights within a shared parent volume (de Jong proportions sum to 1.01)
vol_rows <- vol_rows %>% group_by(volume_set, source_region) %>%
  mutate(allocation_weight_norm = allocation_weight / sum(allocation_weight)) %>% ungroup() %>%
  mutate(volume_mm3 = parent_volume_mm3 * allocation_weight_norm)

budgets <- vol_rows %>%
  inner_join(rates, by = c("source_region", "analysis_anatomy_id" = "anatomy_id")) %>%
  mutate(
    mass_g             = volume_mm3 * TISSUE_DENSITY_MG_PER_MM3 / 1000,
    glucose_umol_min   = rcmr_value * mass_g / 100,
    energy_W           = glucose_umol_min * J_PER_UMOL_GLUCOSE / 60,
    fraction_of_brain_volume = volume_mm3 / TOTAL_BRAIN_MM3,
    is_white_matter    = volume_basis == "white",
    grey_only_volume   = volume_basis == "grey"
  ) %>%
  arrange(volume_set, desc(volume_mm3))

check_region_palette(budgets, region_col = "anatomy_group")
write.csv(budgets %>% select(volume_set, source_region, analysis_anatomy_id, anatomy_group, relationship,
                             n_heiss_components, rcmr_value, rcmr_sd, parent_volume_mm3, allocation_weight_norm,
                             volume_mm3, volume_basis, mass_g, glucose_umol_min, energy_W,
                             fraction_of_brain_volume, volume_source_citation),
          "tables/s5/s5a_region_budgets.csv", row.names = FALSE)

# Heiss regions that have no volume in the primary set
primary_ids <- unique(map$heiss_anatomy_id[map$volume_set == "primary"])
excluded <- heiss %>% filter(!anatomy_id %in% primary_ids) %>%
  transmute(heiss_region, anatomy_id, rcmrglc_mean_umol_100g_min,
            reason = case_when(
              anatomy_id %in% c("frontal_lobe", "parietal_lobe", "temporal_lobe") ~
                "No non-overlapping Stephan-basis lobe volume; covered inside 'Neocortex grey (remaining)'. Kochiyama grey+white lobe volumes used only in the kochiyama_lobes sensitivity set.",
              anatomy_id == "basal_forebrain" ~ "No defensible human volume (Stephan septum is only a partial proxy).",
              TRUE ~ "No volume term in Table 1 / volumes_wide."))
# composite members are represented, not excluded
members <- map %>% filter(volume_set == "primary", n_distinct(heiss_anatomy_id) > 0) %>%
  group_by(source_region, analysis_anatomy_id) %>% filter(n() > 1) %>% ungroup()
excluded <- excluded %>% mutate(
  represented_as_composite = anatomy_id %in% members$heiss_anatomy_id)
# add rows for Heiss regions carried only inside a composite
composite_only <- heiss %>% filter(anatomy_id %in% members$heiss_anatomy_id) %>%
  transmute(heiss_region, anatomy_id, rcmrglc_mean_umol_100g_min,
            reason = paste0("Carried inside composite '", members$source_region[match(anatomy_id, members$heiss_anatomy_id)],
                            "' (no separate human volume)"), represented_as_composite = TRUE)
write.csv(bind_rows(excluded, composite_only), "tables/s5/s5a_excluded_regions.csv", row.names = FALSE)

# ------------------------------------------------------------
# 3. Scaling fits: log10(budget) ~ log10(volume) and log10(rate) ~ log10(volume)
# ------------------------------------------------------------
prim <- budgets %>% filter(volume_set == "primary")
koch_ids <- unique(budgets$analysis_anatomy_id[budgets$volume_set == "kochiyama_lobes"])
# Kochiyama lobe volumes are grey + white, so the Stephan neocortical white matter
# (centrum semiovale row) is dropped from this set to avoid double counting.
sens <- bind_rows(
  prim %>% filter(!analysis_anatomy_id %in% c(koch_ids, "neocortex_grey_matter", "centrum_semiovale")),
  budgets %>% filter(volume_set == "kochiyama_lobes"))

subsets <- list(
  all                = prim,
  excl_white_matter  = prim %>% filter(!is_white_matter),
  grey_only_volumes  = prim %>% filter(grey_only_volume),
  kochiyama_lobes    = sens)

fits <- bind_rows(lapply(names(subsets), function(nm) {
  d <- subsets[[nm]]
  bind_rows(
    fit_sma_ols(d, "volume_mm3", "glucose_umol_min", "Homo sapiens", "Heiss 2004 x Table 1 volumes",
                "log10 regional glucose use (umol/min)", "log10 regional volume (mm3)", nm,
                note = paste0("within-brain, ", nrow(d), " regions; both hemispheres; density ", TISSUE_DENSITY_MG_PER_MM3, " mg/mm3;"),
                test_slope = 1),
    fit_sma_ols(d, "volume_mm3", "rcmr_value", "Homo sapiens", "Heiss 2004 x Table 1 volumes",
                "log10 rCMRGlc (umol/100g/min)", "log10 regional volume (mm3)", nm,
                note = paste0("within-brain, ", nrow(d), " regions; Karbowski-style mass-specific exponent; use OLS - SMA is undefined in meaning when r2 ~ 0;"),
                test_slope = 0))
}))
fits$coverage_fraction_of_brain_volume <- sapply(fits$subset, function(nm) sum(subsets[[nm]]$volume_mm3) / TOTAL_BRAIN_MM3)
write.csv(fits, "tables/s5/s5a_scaling_fits.csv", row.names = FALSE)

# ------------------------------------------------------------
# 4. Comparison with Coto (0.87) and Karbowski 2007 Table 1
# ------------------------------------------------------------
karb_path <- file.path("/Users/crossmodal/Library/CloudStorage/OneDrive-AllenInstitute/Species/Evo-M1-Trait-Data",
                       "Karbowski__2007", "Karbowski__2007_Table1.csv")
karb <- read.csv(karb_path, stringsAsFactors = FALSE)
karb_summary <- karb %>% summarise(n_structures = n(), median_exponent = median(scaling_exponent),
                                   min_exponent = min(scaling_exponent), max_exponent = max(scaling_exponent))
main_sma   <- fits %>% filter(subset == "all", method == "SMA", grepl("glucose", response))
main_rate  <- fits %>% filter(subset == "all", method == "OLS", grepl("rCMRGlc", response))
comparison <- tibble(
  statistic = c("Within-brain SMA slope, log10 budget ~ log10 volume (all regions)",
                "Within-brain OLS slope, log10 rCMRGlc ~ log10 volume (all regions)",
                "Coto 2026 cross-species SMA slope, brain MR ~ brain mass",
                "Karbowski 2007 regional CMRgl exponents on brain mass (median; range)"),
  value = c(sprintf("%.3f [%.3f, %.3f], n = %d", main_sma$slope, main_sma$ci_low, main_sma$ci_high, main_sma$n),
            sprintf("%.3f [%.3f, %.3f], n = %d", main_rate$slope, main_rate$ci_low, main_rate$ci_high, main_rate$n),
            "0.87 +/- 0.07 (r2 0.99, n = 8 incl. ant)",
            sprintf("%.2f; %.2f to %.2f (n = %d structures)", karb_summary$median_exponent, karb_summary$min_exponent,
                    karb_summary$max_exponent, karb_summary$n_structures)),
  interpretation = c("Across regions within one brain: slope > 1 means larger regions use disproportionately more glucose per unit volume.",
                     "Equivalent to (budget slope - 1) under OLS; the within-brain analogue of Karbowski's exponent, but the predictor is regional volume, not species brain mass.",
                     "Across species: larger brains have lower mass-specific metabolic rate (hypometric).",
                     "Across species: mass-specific regional CMRgl falls with brain mass to about the -0.15 power for nearly every structure."))
write.csv(comparison, "tables/s5/s5a_karbowski_comparison.csv", row.names = FALSE)

# ------------------------------------------------------------
# 5. Figures
# ------------------------------------------------------------
plot_df <- bind_rows(prim %>% mutate(panel = "Primary set (Table 1 volume terms)"),
                     sens %>% mutate(panel = "Sensitivity: Kochiyama lobes (s4a)"))
labs_df <- bind_rows(
  fits %>% filter(subset == "all", method == "OLS", grepl("glucose", response)) %>% mutate(panel = "Primary set (Table 1 volume terms)"),
  fits %>% filter(subset == "kochiyama_lobes", method == "OLS", grepl("glucose", response)) %>% mutate(panel = "Sensitivity: Kochiyama lobes (s4a)")) %>%
  mutate(label = sprintf("italic(R)^2~'='~'%.2f'*','~italic(p)~'='~'%s'*','~SMA~slope~'='~'%.2f'", r2, fmt_p(p),
                         fits$slope[fits$subset == subset & fits$method == "SMA" & grepl("glucose", fits$response)]))
pal <- make_region_palette(plot_df$anatomy_group)

p1 <- ggplot(plot_df, aes(x = volume_mm3, y = glucose_umol_min)) +
  geom_smooth(method = "lm", se = TRUE, colour = "steelblue", fill = "steelblue", alpha = 0.15, formula = y ~ x) +
  geom_abline(slope = 1, intercept = log10(mean(prim$rcmr_value) / 100 / 1000), linetype = "dashed", colour = "grey50") +
  geom_point(aes(colour = anatomy_group, shape = volume_basis), size = 3) +
  geom_text_repel(aes(label = anatomy_group), size = 3, max.overlaps = 30, show.legend = FALSE) +
  geom_text(data = labs_df, aes(x = min(plot_df$volume_mm3), y = max(plot_df$glucose_umol_min), label = label),
            inherit.aes = FALSE, hjust = 0, vjust = 1, parse = TRUE, size = 3.6) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = pal, guide = "none") +
  scale_shape_manual(values = c(grey = 16, grey_white = 17, mixed = 15, white = 1), name = "Volume basis") +
  facet_wrap(~ panel) +
  labs(x = expression(Regional~volume~(mm^3*", both hemispheres")),
       y = expression(Regional~glucose~use~(mu*mol~min^-1)),
       caption = "Dashed line: isometry (slope 1) at the mean regional rCMRGlc. Blue: OLS fit with 95% CI.") +
  theme_classic(base_size = 14) + theme(strip.text = element_text(face = "bold"), legend.position = "bottom")
save_png_pdf(p1, "figs/s5/s5a_budget_vs_volume", width = 13, height = 7)

labs2 <- fits %>% filter(subset %in% c("all", "kochiyama_lobes"), method == "OLS", grepl("rCMRGlc", response)) %>%
  mutate(panel = ifelse(subset == "all", "Primary set (Table 1 volume terms)", "Sensitivity: Kochiyama lobes (s4a)"),
         label = sprintf("slope~'='~'%.3f'*','~italic(R)^2~'='~'%.2f'*','~italic(p)~'='~'%s'", slope, r2, fmt_p(p)))
p2 <- ggplot(plot_df, aes(x = volume_mm3, y = rcmr_value)) +
  geom_smooth(method = "lm", se = TRUE, colour = "steelblue", fill = "steelblue", alpha = 0.15, formula = y ~ x) +
  geom_point(aes(colour = anatomy_group, shape = volume_basis), size = 3) +
  geom_text_repel(aes(label = anatomy_group), size = 3, max.overlaps = 30, show.legend = FALSE) +
  geom_text(data = labs2, aes(x = Inf, y = Inf, label = label), inherit.aes = FALSE, hjust = 1.05, vjust = 1.4, parse = TRUE, size = 3.6) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = pal, guide = "none") +
  scale_shape_manual(values = c(grey = 16, grey_white = 17, mixed = 15, white = 1), name = "Volume basis") +
  facet_wrap(~ panel) +
  labs(x = expression(Regional~volume~(mm^3*", both hemispheres")),
       y = expression(rCMRGlc~(mu*mol~100*g^-1~min^-1))) +
  theme_classic(base_size = 14) + theme(strip.text = element_text(face = "bold"), legend.position = "bottom")
save_png_pdf(p2, "figs/s5/s5a_rate_vs_volume", width = 13, height = 7)

# ------------------------------------------------------------
# 6. Console summary
# ------------------------------------------------------------
cat("\n== s5a region budgets (primary set) ==\n")
print(prim %>% select(anatomy_group, rcmr_value, volume_mm3, volume_basis, glucose_umol_min) %>% as.data.frame(), digits = 4)
cat(sprintf("\nPrimary set: %d regions, covering %.1f%% of total brain volume (%.0f of %.0f mm3); summed glucose use %.1f umol/min = %.2f W\n",
            nrow(prim), 100 * sum(prim$volume_mm3) / TOTAL_BRAIN_MM3, sum(prim$volume_mm3), TOTAL_BRAIN_MM3,
            sum(prim$glucose_umol_min), sum(prim$energy_W)))
cat("\n== scaling fits ==\n")
print(fits %>% select(response, subset, method, slope, ci_low, ci_high, r2, p, n) %>% as.data.frame(), digits = 3)
cat("\n== comparison ==\n"); print(as.data.frame(comparison[, 1:2]))
