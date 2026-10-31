# ============================================================
# s5d_cross_species_brain_mr_scaling.R   (Study s5: Coto-style analyses)
#
# Coto et al. 2026 Fig. 7: log10 whole-brain O2 consumption on log10 brain
# mass across 7 mammals + T. immigrans; published SMA slope 0.87 +/- 0.07,
# r2 0.99.  Replicated with (a) Coto's own table and (b) the Evo-M1 species
# table (Karbowski 2007 S1 O2 where available; brain MR in uW from glucose
# where O2 is absent), each with and without the ant point.  Human, mouse
# and ant are tested against the 95% prediction band of the mammal-only OLS.
#
# INPUTS  ../analyses_metabol_rate_structure_nonhuman/data_intermediate/coto_fig7_points.csv
#         ../analyses_metabol_rate_structure_nonhuman/data_intermediate/species_brain_body_metabolism.csv
# OUTPUTS tables/s5/s5d_scaling_fits.csv, s5d_points.csv, s5d_prediction_band_checks.csv
#         figs/s5/s5d_cross_species_brain_mr_scaling.{png,pdf}
# ============================================================
setwd(local({ d <- normalizePath(getwd()); while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d); d }))  # repo root (portable; replaces hardcoded path -- see helpers/project_root.R)
suppressPackageStartupMessages({ library(tidyverse); library(smatr); library(ggrepel) })
source("helpers/plot_settings.R"); source("helpers/s5_helpers.R")

NONHUMAN <- normalizePath(file.path("..", "analyses_metabol_rate_structure_nonhuman"), mustWork = TRUE)
coto <- read.csv(file.path(NONHUMAN, "data_intermediate", "coto_fig7_points.csv"), stringsAsFactors = FALSE)
sp   <- read.csv(file.path(NONHUMAN, "data_intermediate", "species_brain_body_metabolism.csv"), stringsAsFactors = FALSE)

# ---- (a) Coto's own points
a <- coto %>% transmute(dataset = "Coto 2026 Fig. 7 table", species = Species, brain_mass_g = brain_mass_grams,
                        brain_MR_uW = brain_oxygen_consumption_rate_milliliters_per_min * J_PER_ML_O2 / 60 * 1e6,
                        brain_O2_ml_min = brain_oxygen_consumption_rate_milliliters_per_min,
                        is_ant = grepl("Tetramorium", Species), route = "O2 (Coto; ant Q10-corrected to 37 C)")
# ---- (b) Evo-M1 species table
b <- sp %>% transmute(dataset = "Evo-M1 species table", species,
                      brain_mass_g, brain_O2_ml_min,
                      brain_MR_uW = coalesce(brain_MR_uW_from_O2, brain_MR_uW_from_glucose),
                      route = ifelse(!is.na(brain_MR_uW_from_O2), "O2 (Karbowski S1)", "glucose (Karbowski S2 / Evo-M1 merge)"),
                      is_ant = grepl("Tetramorium", species)) %>%
  filter(is.finite(brain_mass_g), is.finite(brain_MR_uW))
# ant: keep the Queenright row (the group Coto plotted) for the "with ant" fits; 31 C, not Q10-corrected
b <- b %>% filter(!is_ant | species == "Tetramorium immigrans (Queenright)")
pts <- bind_rows(a, b) %>% mutate(highlight = case_when(is_ant ~ "Tetramorium immigrans",
                                                        species == "Homo sapiens" ~ "Homo sapiens",
                                                        species == "Mus musculus" ~ "Mus musculus",
                                                        TRUE ~ "Other mammal"))
write.csv(pts, "tables/s5/s5d_points.csv", row.names = FALSE)

# ---- fits
fits <- bind_rows(
  fit_sma_ols(a %>% filter(!is_ant), "brain_mass_g", "brain_O2_ml_min", "Homo sapiens", "Coto 2026 Fig. 7 table",
              "log10 brain O2 consumption (ml/min)", "log10 brain mass (g)", "mammals_only", note = "cross-species; Coto's 7 mammals;"),
  fit_sma_ols(a, "brain_mass_g", "brain_O2_ml_min", "Homo sapiens", "Coto 2026 Fig. 7 table",
              "log10 brain O2 consumption (ml/min)", "log10 brain mass (g)", "mammals_plus_ant", note = "cross-species; published 0.87 +/- 0.07, r2 0.99;"),
  fit_sma_ols(b %>% filter(!is_ant, !is.na(brain_O2_ml_min)), "brain_mass_g", "brain_O2_ml_min", "Homo sapiens", "Evo-M1 species table",
              "log10 brain O2 consumption (ml/min)", "log10 brain mass (g)", "mammals_only_O2", note = "cross-species; Karbowski S1 O2 species;"),
  fit_sma_ols(b %>% filter(!is.na(brain_O2_ml_min)), "brain_mass_g", "brain_O2_ml_min", "Homo sapiens", "Evo-M1 species table",
              "log10 brain O2 consumption (ml/min)", "log10 brain mass (g)", "mammals_plus_ant_O2", note = "cross-species; ant Queenright at 31 C;"),
  fit_sma_ols(b %>% filter(!is_ant), "brain_mass_g", "brain_MR_uW", "Homo sapiens", "Evo-M1 species table",
              "log10 brain MR (uW)", "log10 brain mass (g)", "mammals_only_all", note = "cross-species; O2 route where available else glucose;"),
  fit_sma_ols(b, "brain_mass_g", "brain_MR_uW", "Homo sapiens", "Evo-M1 species table",
              "log10 brain MR (uW)", "log10 brain mass (g)", "mammals_plus_ant_all", note = "cross-species; O2 route where available else glucose; ant Queenright at 31 C;"))
fits$dataset <- paste0(fits$dataset, " (species: all)")
fits$dataset <- sub(" \\(species: all\\)", "", fits$dataset)
write.csv(fits, "tables/s5/s5d_scaling_fits.csv", row.names = FALSE)

# ---- prediction-band checks (mammal-only OLS, log10 brain MR uW ~ log10 brain mass)
band_check <- function(train, test, label) {
  m <- lm(log10(brain_MR_uW) ~ log10(brain_mass_g), data = train)
  pr <- predict(m, newdata = test, interval = "prediction", level = 0.95)
  tibble(fit_dataset = label, species = test$species, brain_mass_g = test$brain_mass_g,
         observed_log10_uW = log10(test$brain_MR_uW), predicted_log10_uW = pr[, "fit"],
         pred_low = pr[, "lwr"], pred_high = pr[, "upr"],
         within_95_prediction_band = log10(test$brain_MR_uW) >= pr[, "lwr"] & log10(test$brain_MR_uW) <= pr[, "upr"],
         observed_over_predicted = 10^(log10(test$brain_MR_uW) - pr[, "fit"]),
         n_train = nrow(train), in_training_set = test$species %in% train$species)
}
checks <- bind_rows(
  band_check(a %>% filter(!is_ant), a %>% filter(species %in% c("Homo sapiens", "Tetramorium immigrans (Queenright)")), "Coto table, 7 mammals"),
  band_check(b %>% filter(!is_ant), b %>% filter(species %in% c("Homo sapiens", "Mus musculus", "Tetramorium immigrans (Queenright)")), sprintf("Evo-M1 mammals (n = %d)", sum(!b$is_ant))),
  # leave-one-out versions for human and mouse (they are in the training set)
  band_check(b %>% filter(!is_ant, species != "Homo sapiens"), b %>% filter(species == "Homo sapiens"), "Evo-M1 mammals excl. Homo (LOO)"),
  band_check(b %>% filter(!is_ant, species != "Mus musculus"), b %>% filter(species == "Mus musculus"), "Evo-M1 mammals excl. Mus (LOO)"))
write.csv(checks, "tables/s5/s5d_prediction_band_checks.csv", row.names = FALSE)

# ---- figure
band_df <- bind_rows(lapply(split(pts, pts$dataset), function(d) {
  tr <- d %>% filter(!is_ant)
  m <- lm(log10(brain_MR_uW) ~ log10(brain_mass_g), data = tr)
  grid <- tibble(brain_mass_g = 10^seq(log10(min(d$brain_mass_g)) - 0.3, log10(max(d$brain_mass_g)) + 0.3, length.out = 100))
  pr <- predict(m, newdata = grid, interval = "prediction", level = 0.95)
  grid %>% mutate(dataset = d$dataset[1], fit = 10^pr[, "fit"], lwr = 10^pr[, "lwr"], upr = 10^pr[, "upr"])
}))
pos <- pts %>% group_by(dataset) %>% summarise(x = min(brain_mass_g), y = max(brain_MR_uW), .groups = "drop")
lab_df <- fits %>% filter(method == "SMA", subset %in% c("mammals_only", "mammals_only_all")) %>%
  mutate(label = sprintf("Mammal-only SMA slope = %.2f [%.2f, %.2f], r2 = %.2f, n = %d", slope, ci_low, ci_high, r2, n)) %>%
  left_join(pos, by = "dataset")
lab_df2 <- fits %>% filter(method == "SMA", subset %in% c("mammals_plus_ant", "mammals_plus_ant_all")) %>%
  mutate(label = sprintf("With ant: SMA slope = %.2f [%.2f, %.2f], r2 = %.2f, n = %d", slope, ci_low, ci_high, r2, n)) %>%
  left_join(pos, by = "dataset")
p <- ggplot(pts, aes(x = brain_mass_g, y = brain_MR_uW)) +
  geom_ribbon(data = band_df, aes(x = brain_mass_g, ymin = lwr, ymax = upr), inherit.aes = FALSE, fill = "steelblue", alpha = 0.15) +
  geom_line(data = band_df, aes(x = brain_mass_g, y = fit), inherit.aes = FALSE, colour = "steelblue", linewidth = 0.9) +
  geom_point(aes(colour = highlight, shape = route), size = 3.2) +
  geom_text_repel(aes(label = species), size = 3, max.overlaps = 30) +
  geom_text(data = lab_df, aes(x = x, y = y, label = label), inherit.aes = FALSE, hjust = 0, vjust = 0.5, size = 3.4) +
  geom_text(data = lab_df2, aes(x = x, y = y, label = label), inherit.aes = FALSE, hjust = 0, vjust = 2.2, size = 3.4) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = c("Tetramorium immigrans" = "#D55E00", "Homo sapiens" = "#0072B2", "Mus musculus" = "#009E73", "Other mammal" = "grey40"), name = NULL) +
  scale_shape_manual(values = c(16, 17, 15), name = "Rate route") +
  facet_wrap(~ dataset, scales = "free") +
  labs(x = "Brain mass (g)", y = expression(Whole-brain~metabolic~rate~(mu*W)),
       caption = "Line and band: mammal-only OLS fit with 95% prediction interval (ant excluded from the fit).") +
  theme_classic(base_size = 14) + theme(strip.text = element_text(face = "bold"), legend.position = "bottom", legend.box = "vertical")
save_png_pdf(p, "figs/s5/s5d_cross_species_brain_mr_scaling", width = 14, height = 7.5)

cat("\n== s5d fits ==\n"); print(fits %>% select(dataset, subset, method, slope, ci_low, ci_high, r2, p, n) %>% as.data.frame(), digits = 3)
cat("\n== prediction band checks ==\n"); print(checks %>% select(fit_dataset, species, observed_log10_uW, predicted_log10_uW, pred_low, pred_high, within_95_prediction_band, observed_over_predicted) %>% as.data.frame(), digits = 3)
