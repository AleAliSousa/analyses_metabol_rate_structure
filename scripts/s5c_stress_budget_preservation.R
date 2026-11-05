# ============================================================
# s5c_stress_budget_preservation.R   (Study s5: Coto-style analyses)
#
# Is regional energy preferentially preserved under deprivation?  Coto 2026:
# ant brain MR was maintained while body MR fell ~30% under food deprivation.
# Human analogue (s2 inputs): Mackes et al. 2020 adult adoptees (early
# institutional deprivation; English and Romanian adoptees) vs normative
# volumes (Walhovd et al. 2011).  Per region: budget = rCMRGlc x volume.
# Summed over the regions with a Heiss rate, compare % change of budget with
# % change of volume (normative -> adoptee); the difference is the signature
# of expensive regions being spared (budget falls less than volume).
#
# Uncertainty: parametric bootstrap of the group means from the reported
# mean / SD / n per region (2000 draws).  n: normative 262 (Walhovd 2011
# Overall_N262), English adoptees 21, Romanian adoptees 67 (Mackes 2020
# source files data_raw/Mackes_etal_2020_*_derived.csv).
#
# CAVEAT  The stressor (early institutional deprivation measured in adulthood,
#         volume outcome) differs from acute food deprivation in ants
#         (metabolic-rate outcome) and from perinatal protein restriction in
#         mice; the comparison is of the direction of budget preservation only.
#
# INPUTS  data_intermediate/s2_stress_volume_comparison.csv (volumes per hemisphere, cm3)
#         metadata/anatomy/s2_stress_volume_region_map.csv
# OUTPUTS tables/s5/s5c_region_budgets_by_group.csv, s5c_budget_vs_volume_change.csv,
#         s5c_regression_volume_change_vs_rate.csv
#         figs/s5/s5c_budget_vs_volume_change.{png,pdf}
# ============================================================
setwd(local({ d <- normalizePath(getwd()); while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d); d }))  # repo root (portable; replaces hardcoded path -- see helpers/project_root.R)
source("helpers/model_diagnostics.R"); diag_stage("s5")  # model diagnostics -> figs/s5/diagnostics/
suppressPackageStartupMessages({ library(tidyverse) })
source("helpers/plot_settings.R"); source("helpers/read_heiss_rates.R"); source("helpers/s5_helpers.R")
set.seed(20261031)
N_BOOT <- 2000

data <- read.csv("data_intermediate/s2_stress_volume_comparison.csv", stringsAsFactors = FALSE, check.names = FALSE)
rates <- derive_study_heiss_rates("metadata/anatomy/s2_stress_volume_region_map.csv") %>%
  select(anatomy_id, anatomy_group, rcmr_value, n_heiss_components) %>% distinct()
data <- data %>% left_join(rates, by = "anatomy_id")

# group sample sizes from the source files
mack <- read.csv("data_raw/Mackes_etal_2020_subcortical_derived.csv", stringsAsFactors = FALSE)
n_group <- c(norm = 262L,
             english = unique(mack$n[mack$Group == "UK"]),
             romanian = unique(mack$n[mack$Group == "Romanian"]))
stopifnot(length(n_group) == 3, all(is.finite(n_group)))
walh <- read.csv("data_raw/Walhovd_etal_2011_derived.csv", stringsAsFactors = FALSE)
stopifnot(grepl("N262", names(walh)[2]))

long <- data %>% filter(anatomy_id != "brain", !is.na(rcmr_value)) %>%
  select(source_region, anatomy_id, anatomy_group, rcmr_value, ends_with("_mean"), ends_with("_sd")) %>%
  pivot_longer(c(ends_with("_mean"), ends_with("_sd")), names_to = c("group", ".value"), names_sep = "_") %>%
  mutate(n = n_group[group],
         # volume per hemisphere in cm3 -> g at 1 g/cm3; rCMRGlc per 100 g
         budget_umol_min = rcmr_value * mean * (TISSUE_DENSITY_MG_PER_MM3 * 1000 / 1000) / 100,
         budget_W = budget_umol_min * J_PER_UMOL_GLUCOSE / 60)
write.csv(long %>% arrange(group, desc(mean)), "tables/s5/s5c_region_budgets_by_group.csv", row.names = FALSE)

# ---- point estimates and bootstrap of summed volume / budget / weighted rate
summ <- function(vols, rate) {
  v <- sum(vols); b <- sum(rate * vols / 100)
  c(volume_cm3 = v, budget_umol_min = b, weighted_rate = 100 * b / v)
}
groups <- c("norm", "english", "romanian")
wide <- long %>% select(anatomy_id, rcmr_value, group, mean, sd, n) %>%
  pivot_wider(names_from = group, values_from = c(mean, sd, n))
point <- sapply(groups, function(g) summ(wide[[paste0("mean_", g)]], wide$rcmr_value))
boot <- array(NA_real_, dim = c(N_BOOT, 3, length(groups)), dimnames = list(NULL, rownames(point), groups))
for (g in groups) {
  mu <- wide[[paste0("mean_", g)]]; se <- wide[[paste0("sd_", g)]] / sqrt(wide[[paste0("n_", g)]])
  for (i in seq_len(N_BOOT)) boot[i, , g] <- summ(rnorm(length(mu), mu, se), wide$rcmr_value)
}
pct <- function(a, b) 100 * (a - b) / b   # (stressed - control) / control
res <- bind_rows(lapply(c("english", "romanian"), function(g) {
  bind_rows(lapply(c("volume_cm3", "budget_umol_min", "weighted_rate"), function(q) {
    draws <- pct(boot[, q, g], boot[, q, "norm"])
    tibble(contrast = paste0("normative -> ", g, " adoptees"), group = g, quantity = q,
           control_value = point[q, "norm"], stressed_value = point[q, g],
           pct_change = pct(point[q, g], point[q, "norm"]),
           ci_low = unname(quantile(draws, 0.025)), ci_high = unname(quantile(draws, 0.975)),
           boot_p_two_sided = 2 * min(mean(draws <= 0), mean(draws >= 0)),
           n_control = n_group[["norm"]], n_stressed = n_group[[g]], n_regions = nrow(wide), n_boot = N_BOOT)
  })) %>% bind_rows(
    { d <- pct(boot[, "budget_umol_min", g], boot[, "budget_umol_min", "norm"]) - pct(boot[, "volume_cm3", g], boot[, "volume_cm3", "norm"])
      tibble(contrast = paste0("normative -> ", g, " adoptees"), group = g, quantity = "budget_minus_volume_pct_change",
             control_value = NA, stressed_value = NA,
             pct_change = pct(point["budget_umol_min", g], point["budget_umol_min", "norm"]) - pct(point["volume_cm3", g], point["volume_cm3", "norm"]),
             ci_low = unname(quantile(d, 0.025)), ci_high = unname(quantile(d, 0.975)),
             boot_p_two_sided = 2 * min(mean(d <= 0), mean(d >= 0)),
             n_control = n_group[["norm"]], n_stressed = n_group[[g]], n_regions = nrow(wide), n_boot = N_BOOT) })
})) %>% mutate(units = recode(quantity, volume_cm3 = "cm3 per hemisphere, summed over regions with a Heiss rate",
                              budget_umol_min = "umol glucose/min per hemisphere (rCMRGlc x volume at 1 g/cm3)",
                              weighted_rate = "umol/100 g/min (volume-weighted mean rCMRGlc)",
                              budget_minus_volume_pct_change = "percentage points"))
write.csv(res, "tables/s5/s5c_budget_vs_volume_change.csv", row.names = FALSE)

# ---- s2 regression in standardized form: % volume change (stressed - control) on rCMRGlc
reg_df <- data %>% filter(anatomy_id != "brain", !is.na(rcmr_value)) %>%
  transmute(anatomy_id, anatomy_group, rcmr_value,
            english = -norm_english, romanian = -norm_romanian, all = -norm_all) %>%  # s2 sign is (norm - adoptee); flip
  pivot_longer(c(english, romanian, all), names_to = "group", values_to = "pct_change_volume")
reg <- reg_df %>% group_by(group) %>% group_modify(function(d, k) {
  m <- lm(pct_change_volume ~ rcmr_value, data = d); s <- summary(m); ci <- confint(m)
  diag_lm(m, paste("pct volume change ~ rCMRGlc,", k$group), labels = d$anatomy_group)
  sp <- suppressWarnings(cor.test(d$rcmr_value, d$pct_change_volume, method = "spearman", exact = FALSE))
  tibble(contrast = paste0("normative -> ", k$group, " adoptees"), n_regions = nrow(d),
         lm_slope = coef(m)[[2]], lm_ci_low = ci[2, 1], lm_ci_high = ci[2, 2], lm_intercept = coef(m)[[1]],
         lm_r2 = s$r.squared, lm_p = s$coefficients[2, 4], spearman_rho = unname(sp$estimate), spearman_p = sp$p.value)
}) %>% ungroup()
write.csv(reg, "tables/s5/s5c_regression_volume_change_vs_rate.csv", row.names = FALSE)

# ---- figure: paired bars
fig_df <- res %>% filter(quantity %in% c("volume_cm3", "budget_umol_min")) %>%
  mutate(group = recode(group, english = "English adoptees", romanian = "Romanian adoptees"),
         quantity = recode(quantity, volume_cm3 = "% change in summed volume", budget_umol_min = "% change in summed energy budget"),
         quantity = factor(quantity, levels = c("% change in summed volume", "% change in summed energy budget")))
p <- ggplot(fig_df, aes(x = group, y = pct_change, fill = quantity)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, colour = "black") +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), position = position_dodge(width = 0.7), width = 0.2) +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values = c("% change in summed volume" = "#999999", "% change in summed energy budget" = "#D55E00"), name = NULL) +
  labs(x = NULL, y = "% change, normative -> adoptees",
       subtitle = "Summed over the 8 regions with a Heiss rate; budget = rCMRGlc x volume",
       caption = sprintf("Whiskers: 95%% CI from %d parametric bootstrap draws of the group means.", N_BOOT)) +
  theme_classic(base_size = 14) + theme(legend.position = "bottom")
save_png_pdf(p, "figs/s5/s5c_budget_vs_volume_change", width = 8, height = 6)

cat("\n== n per group ==\n"); print(n_group)
cat("\n== budget vs volume change ==\n"); print(res %>% select(group, quantity, control_value, stressed_value, pct_change, ci_low, ci_high, boot_p_two_sided) %>% as.data.frame(), digits = 4)
cat("\n== regression ==\n"); print(as.data.frame(reg), digits = 3)
