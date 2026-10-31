# ============================================================
# s5e_export_comparison_inputs.R   (Study s5)
# Writes the standardized human rows for the humans-mice-ants comparison
# (Phase 3) from the s5 tables:
#   ../analyses_metabol_rate_structure_nonhuman/comparison/inputs/
#       human_comparison_rows.csv, human_regional_change_vs_rate.csv, human_scaling_fits.csv
# Run after s5a-s5d.
# ============================================================
setwd(local({ d <- normalizePath(getwd()); while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d); d }))  # repo root (portable; replaces hardcoded path -- see helpers/project_root.R)
suppressPackageStartupMessages({ library(tidyverse) })
source("helpers/read_heiss_rates.R"); source("helpers/s5_helpers.R")
OUT <- normalizePath(file.path("..", "analyses_metabol_rate_structure_nonhuman", "comparison", "inputs"), mustWork = FALSE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
SP <- "Homo sapiens"

s5b <- read.csv("tables/s5/s5b_relative_brain_metabolism.csv", stringsAsFactors = FALSE)
s5c <- read.csv("tables/s5/s5c_budget_vs_volume_change.csv", stringsAsFactors = FALSE)
s5c_reg <- read.csv("tables/s5/s5c_regression_volume_change_vs_rate.csv", stringsAsFactors = FALSE)
s5c_bud <- read.csv("tables/s5/s5c_region_budgets_by_group.csv", stringsAsFactors = FALSE)
s5a_fits <- read.csv("tables/s5/s5a_scaling_fits.csv", stringsAsFactors = FALSE)
s5d_fits <- read.csv("tables/s5/s5d_scaling_fits.csv", stringsAsFactors = FALSE)
s2 <- read.csv("data_intermediate/s2_stress_volume_comparison.csv", stringsAsFactors = FALSE, check.names = FALSE)

# ---- (a) comparison rows
hs <- s5b[s5b$group == "Heiss 2004 (volume-weighted) + Karbowski O2", ]
wb <- tibble(
  species = SP, dataset = "Heiss 2004 + Evo-M1 (Karbowski 2007; body_ecology)", stressor = "none", stress_timing = NA,
  level = c("brain_vs_body", "brain_vs_body", "brain_vs_body", "whole_brain", "whole_brain", "whole_brain", "whole_brain"),
  variable = c("relative_brain_mass", "relative_brain_MR_glucose", "relative_brain_MR_O2",
               "brain_MR_W_glucose", "brain_MR_W_O2", "wholebrain_CMRgl_umol_100g_min", "brain_mass_specific_MR_mlO2_min_g"),
  control_value = c(hs$relative_brain_mass, hs$relative_brain_MR_glucose, hs$relative_brain_MR_O2,
                    hs$brain_MR_W_glucose, hs$brain_MR_W_O2, hs$wholebrain_CMRgl_umol_100g_min, hs$brain_mass_specific_MR_mlO2_min_g),
  stressed_value = NA, pct_change = NA,
  units = c("g/g", "W/W", "W/W", "W", "W", "umol/100 g/min", "ml O2/min/g"),
  statistic_type = "point estimate", statistic = NA, p_value = NA, n_control = NA, n_stressed = NA,
  note = c("brain 1300.4 g / body 61535 g (Evo-M1)",
           "Heiss volume-weighted rCMRGlc x brain mass x 2.87 J/umol over BMR 12925 ml O2/h x 20.1 J/ml",
           "Karbowski S1 48.6 ml O2/min x 20.1 J/ml over BMR",
           "glucose route (Heiss regions weighted by s5a volumes, 96.5% of brain volume)",
           "O2 route (Karbowski 2007 S1)",
           "volume-weighted mean of 18 Heiss regions (incl. white matter); cortex global average 33.5; Clarke & Sokoloff 31",
           "48.6 ml O2/min / 1300.4 g; Coto quotes ~0.035"))
reg_rows <- s5c %>% transmute(
  species = SP, dataset = "Mackes 2020 adoptees vs Walhovd 2011 normative (s2 inputs) x Heiss 2004",
  stressor = "early institutional deprivation (adoptees)", stress_timing = "early childhood; measured in adulthood",
  level = "regional", variable = paste0("summed_", quantity, "_", group, "_adoptees"),
  control_value, stressed_value, pct_change, units,
  statistic_type = "95% CI of % change (parametric bootstrap, 2000 draws)",
  statistic = sprintf("[%.2f, %.2f]", ci_low, ci_high), p_value = boot_p_two_sided,
  n_control, n_stressed,
  note = paste0(contrast, "; 8 regions with a Heiss rate (Total Brain excluded); per hemisphere"))
# per-region rows (volume per hemisphere), English and Romanian
per_region <- s2 %>% filter(anatomy_id != "brain") %>%
  pivot_longer(c(english_mean, romanian_mean), names_to = "group", values_to = "stressed_value") %>%
  mutate(group = sub("_mean", "", group),
         species = SP, dataset = "Mackes 2020 adoptees vs Walhovd 2011 normative (s2 inputs)",
         stressor = "early institutional deprivation (adoptees)", stress_timing = "early childhood; measured in adulthood",
         level = "regional", variable = paste0("volume_", anatomy_id, "_", group, "_adoptees"),
         control_value = norm_mean, pct_change = 100 * (stressed_value - norm_mean) / norm_mean,
         units = "cm3 per hemisphere", statistic_type = "group means (SD in source)",
         statistic = sprintf("SD control %.3f; SD stressed %.3f", norm_sd, ifelse(group == "english", english_sd, romanian_sd)),
         p_value = NA, n_control = 262L, n_stressed = ifelse(group == "english", 21L, 67L),
         note = paste0(source_region, "; sign = (adoptee - normative)/normative")) %>%
  select(species, dataset, stressor, stress_timing, level, variable, control_value, stressed_value, pct_change, units,
         statistic_type, statistic, p_value, n_control, n_stressed, note)
slope_rows <- s5c_reg %>% transmute(
  species = SP, dataset = "Mackes 2020 adoptees vs Walhovd 2011 normative x Heiss 2004 (s2 regression)",
  stressor = "early institutional deprivation (adoptees)", stress_timing = "early childhood; measured in adulthood",
  level = "regional", variable = paste0("lm_slope_pct_volume_change_on_rCMRGlc_", group, "_adoptees"),
  control_value = NA, stressed_value = NA, pct_change = NA, units = "% change per umol/100 g/min",
  statistic_type = "lm slope [95% CI]; R2; Spearman rho",
  statistic = sprintf("slope %.3f [%.3f, %.3f]; R2 %.3f; rho %.3f (p %.3f)", lm_slope, lm_ci_low, lm_ci_high, lm_r2, spearman_rho, spearman_p),
  p_value = lm_p, n_control = n_regions, n_stressed = n_regions,
  note = paste0(contrast, "; response = (adoptee - normative)/normative x 100; n = regions"))
rows <- bind_rows(wb, reg_rows, per_region, slope_rows)
write.csv(rows, file.path(OUT, "human_comparison_rows.csv"), row.names = FALSE)

# ---- (b) regional change vs rate
rates <- derive_study_heiss_rates("metadata/anatomy/s2_stress_volume_region_map.csv") %>% select(anatomy_id, rcmr_value) %>% distinct()
rcr <- s2 %>% left_join(rates, by = "anatomy_id") %>%
  pivot_longer(c(english_mean, romanian_mean), names_to = "group", values_to = "stressed_volume") %>%
  mutate(group = sub("_mean", "", group),
         species = SP, contrast = paste0("normative -> ", group, " adoptees"),
         region_label = source_region, analysis_anatomy_id = anatomy_id,
         metabolic_rate_umol_100g_min = rcmr_value, control_volume = norm_mean,
         pct_change_volume = 100 * (stressed_volume - control_volume) / control_volume,
         control_budget = rcmr_value * control_volume / 100, stressed_budget = rcmr_value * stressed_volume / 100,
         pct_change_budget = 100 * (stressed_budget - control_budget) / control_budget,
         included_in_regression = !is.na(rcmr_value) & anatomy_id != "brain",
         note = ifelse(anatomy_id == "brain", "Total Brain: no Heiss rate; volume only",
                       ifelse(anatomy_id == "thalamus", "rate = unweighted mean of 3 Heiss thalamic terms (s2 rule); budget = umol/min per hemisphere (volume cm3 at 1 g/cm3)",
                              "budget = umol glucose/min per hemisphere (rCMRGlc x volume cm3 at 1 g/cm3 / 100)"))) %>%
  select(species, contrast, region_label, analysis_anatomy_id, metabolic_rate_umol_100g_min, control_volume, stressed_volume,
         pct_change_volume, control_budget, stressed_budget, pct_change_budget, included_in_regression, note)
write.csv(rcr, file.path(OUT, "human_regional_change_vs_rate.csv"), row.names = FALSE)

# ---- (c) scaling fits
sf <- bind_rows(s5a_fits %>% select(-coverage_fraction_of_brain_volume), s5d_fits) %>%
  select(species, dataset, response, predictor, subset, method, slope, ci_low, ci_high, intercept, r2, p, n, note)
write.csv(sf, file.path(OUT, "human_scaling_fits.csv"), row.names = FALSE)

cat("wrote", nrow(rows), "comparison rows,", nrow(rcr), "regional rows,", nrow(sf), "scaling fits to", OUT, "\n")
