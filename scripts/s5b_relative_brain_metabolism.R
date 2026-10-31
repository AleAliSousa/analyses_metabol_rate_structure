# ============================================================
# s5b_relative_brain_metabolism.R   (Study s5: Coto-style analyses)
#
# Relative brain mass (brain / body mass) and relative brain metabolic rate
# (brain MR / whole-body MR) for humans, with mouse and the three Coto 2026
# ant treatment groups side by side (Coto Fig. 2 / Table S1 analogue).
#
# Human brain MR, glucose route: volume-weighted whole-brain rCMRGlc from the
#   Heiss 2004 regions (s5a primary volume set), scaled to total brain mass.
#   Cross-checks: Heiss cortex global average (33.5); Clarke & Sokoloff 1994
#   whole-brain CMRgl (0.31 umol/g/min) from the Evo-M1 species table.
# Human brain MR, O2 route: Karbowski 2007 S1 whole-brain O2 (ml/min).
# Body MR: Evo-M1 body_ecology BMR (ml O2/h) via species_brain_body_metabolism.csv.
#
# INPUTS  tables/s5/s5a_region_budgets.csv, data_intermediate/heiss_2004_regions.csv,
#         ../analyses_metabol_rate_structure_nonhuman/data_intermediate/species_brain_body_metabolism.csv
# OUTPUTS tables/s5/s5b_relative_brain_metabolism.csv, tables/s5/s5b_human_wholebrain_rate_checks.csv
#         figs/s5/s5b_relative_mass_vs_relative_MR.{png,pdf}
# ============================================================
setwd(local({ d <- normalizePath(getwd()); while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d); d }))  # repo root (portable; replaces hardcoded path -- see helpers/project_root.R)
suppressPackageStartupMessages({ library(tidyverse); library(ggrepel) })
source("helpers/plot_settings.R"); source("helpers/read_heiss_rates.R"); source("helpers/s5_helpers.R")

NONHUMAN <- normalizePath(file.path("..", "analyses_metabol_rate_structure_nonhuman"), mustWork = TRUE)
sp <- read.csv(file.path(NONHUMAN, "data_intermediate", "species_brain_body_metabolism.csv"), stringsAsFactors = FALSE)
bud <- read.csv("tables/s5/s5a_region_budgets.csv", stringsAsFactors = FALSE) %>% filter(volume_set == "primary")
heiss <- read_heiss_regions()
vw <- read.csv("data_raw/volumes_wide.csv", check.names = FALSE, stringsAsFactors = FALSE)
TOTAL_BRAIN_MM3 <- as.numeric(vw[vw$Species == "Homo sapiens", "Total_brain_net_volume_Vol.mm3"])

# ---- 1. Human whole-brain rCMRGlc from Heiss regions
wb_rate_weighted <- sum(bud$rcmr_value * bud$volume_mm3) / sum(bud$volume_mm3)
coverage <- sum(bud$volume_mm3) / TOTAL_BRAIN_MM3
cortex_avg <- heiss$rcmrglc_mean_umol_100g_min[heiss$anatomy_id == "cerebral_cortex"]
grey_rate_weighted <- with(bud %>% filter(volume_basis != "white"), sum(rcmr_value * volume_mm3) / sum(volume_mm3))
hs <- sp[sp$species == "Homo sapiens", ]
checks <- tibble(
  quantity = c("Volume-weighted whole-brain rCMRGlc (s5a primary set)",
               "Volume-weighted rCMRGlc, non-white-matter regions only",
               "Heiss cortex global average",
               "Clarke & Sokoloff 1994 whole-brain CMRgl (Evo-M1 Karbowski S2)",
               "Evo-M1 merge whole-brain CMRgl, human (primary-only mean)",
               "Fraction of total brain volume covered by the s5a primary set"),
  value = c(wb_rate_weighted, grey_rate_weighted, cortex_avg, hs$brain_CMRgl_umol_g_min * 100,
            hs$merge_wholebrain_CMRgl_umol_100g_min, coverage),
  units = c(rep("umol/100 g/min", 5), "fraction"))
write.csv(checks, "tables/s5/s5b_human_wholebrain_rate_checks.csv", row.names = FALSE)

# ---- 2. Side-by-side table with every conversion step
mouse <- sp[sp$species == "Mus musculus", ]
ants  <- sp[grepl("Tetramorium", sp$species), ]
MOUSE_CMRGL <- mouse$merge_wholebrain_CMRgl_umol_100g_min  # Evo-M1 whole-brain CMRgl 74.8

row_from <- function(species, group, brain_g, body_g, bmr_ml_h, cmrgl_100g = NA, o2_ml_min = NA, body_W_direct = NA,
                     rate_source, note) {
  glucose_umol_min <- cmrgl_100g * brain_g / 100
  brain_W_glc <- glucose_umol_min * J_PER_UMOL_GLUCOSE / 60
  brain_W_o2  <- o2_ml_min * J_PER_ML_O2 / 60
  body_W <- if (!is.na(body_W_direct)) body_W_direct else bmr_ml_h * J_PER_ML_O2 / 3600
  tibble(species, group, brain_mass_g = brain_g, body_mass_g = body_g,
         relative_brain_mass = brain_g / body_g,
         wholebrain_CMRgl_umol_100g_min = cmrgl_100g, glucose_umol_min,
         brain_MR_W_glucose = brain_W_glc,
         brain_O2_ml_min = o2_ml_min, brain_MR_W_O2 = brain_W_o2,
         BMR_mlO2_h = bmr_ml_h, body_MR_W = body_W,
         relative_brain_MR_glucose = brain_W_glc / body_W,
         relative_brain_MR_O2 = brain_W_o2 / body_W,
         brain_mass_specific_MR_mlO2_min_g = o2_ml_min / brain_g,
         brain_mass_specific_MR_W_g = coalesce(brain_W_o2, brain_W_glc) / brain_g,
         rate_source, note)
}
tab <- bind_rows(
  row_from("Homo sapiens", "Heiss 2004 (volume-weighted) + Karbowski O2", hs$brain_mass_g, hs$body_mass_g, hs$BMR_mlO2_h,
           cmrgl_100g = wb_rate_weighted, o2_ml_min = hs$brain_O2_ml_min,
           rate_source = "glucose: Heiss 2004 regions weighted by s5a volumes; O2: Karbowski 2007 S1",
           note = sprintf("volume coverage %.3f of brain; brain mass %.1f g (Evo-M1 brain_mass merge)", coverage, hs$brain_mass_g)),
  row_from("Homo sapiens", "Clarke & Sokoloff whole-brain CMRgl (check)", hs$brain_mass_g, hs$body_mass_g, hs$BMR_mlO2_h,
           cmrgl_100g = hs$brain_CMRgl_umol_g_min * 100, o2_ml_min = hs$brain_O2_ml_min,
           rate_source = "glucose: Karbowski 2007 S2 (Clarke & Sokoloff 1994); O2: Karbowski S1",
           note = "cross-check row; same body MR"),
  row_from("Mus musculus", "Evo-M1 whole-brain CMRgl", mouse$brain_mass_g, mouse$body_mass_g, mouse$BMR_mlO2_h,
           cmrgl_100g = MOUSE_CMRGL, o2_ml_min = NA,
           rate_source = "glucose: Evo-M1 merge whole-brain CMRgl (primary-only mean)",
           note = sprintf("species table column brain_MR_uW_from_glucose = %.0f uW uses %.2f umol/min (Karbowski S2 total); here %.2f umol/min from 74.8 x brain mass",
                          mouse$brain_MR_uW_from_glucose, mouse$brain_glucose_umol_min, MOUSE_CMRGL * mouse$brain_mass_g / 100)))
for (i in seq_len(nrow(ants))) {
  a <- ants[i, ]
  tab <- bind_rows(tab, row_from("Tetramorium immigrans", sub(".*\\((.*)\\)", "\\1", a$species), a$brain_mass_g, a$body_mass_g, NA,
                                 cmrgl_100g = NA, o2_ml_min = a$brain_O2_ml_min, body_W_direct = a$body_MR_uW / 1e6,
                                 rate_source = "Coto 2026 Dataset S1 treatment means (31 C, not Q10-corrected); body MR = whole-body O2 at 31 C",
                                 note = sprintf("n = %d workers; relative_brain_MR_O2 in species table = %.3f", a$n, a$relative_brain_MR_O2)))
}
# Q10 = 2 correction of ant brain OCR to 37 C (Coto SI: queenright 0.102 ml/min/g at 37 C)
tab <- tab %>% mutate(
  brain_mass_specific_MR_mlO2_min_g_37C = ifelse(species == "Tetramorium immigrans",
                                                 brain_mass_specific_MR_mlO2_min_g * 2^((37 - 31) / 10), brain_mass_specific_MR_mlO2_min_g),
  ratio_to_human_mass_specific_O2 = brain_mass_specific_MR_mlO2_min_g_37C / brain_mass_specific_MR_mlO2_min_g[1])
# Reference values carried in the shared species table (ant rows there are Coto's per-worker-mean ratios,
# which match the published 24% / 45% / 30%; the rows above are ratios of treatment means).
ref_key <- c("Heiss 2004 (volume-weighted) + Karbowski O2" = "Homo sapiens", "Clarke & Sokoloff whole-brain CMRgl (check)" = "Homo sapiens",
             "Evo-M1 whole-brain CMRgl" = "Mus musculus", Control = "Tetramorium immigrans (Control)",
             "Food Deprived" = "Tetramorium immigrans (Food Deprived)", Queenright = "Tetramorium immigrans (Queenright)")
tab$relative_brain_MR_O2_species_table <- sp$relative_brain_MR_O2[match(ref_key[tab$group], sp$species)]
tab$relative_brain_MR_glucose_species_table <- sp$relative_brain_MR_glucose[match(ref_key[tab$group], sp$species)]
tab$coto_published_relative_brain_MR <- c(NA, NA, NA, 0.24, 0.45, 0.30)
write.csv(tab, "tables/s5/s5b_relative_brain_metabolism.csv", row.names = FALSE)

# ---- 3. Figure: relative brain mass vs relative brain MR (all species with both; ants as separate points)
all_sp <- sp %>% transmute(species, relative_brain_mass,
                           relative_brain_MR = coalesce(relative_brain_MR_O2, relative_brain_MR_glucose),
                           route = ifelse(!is.na(relative_brain_MR_O2), "O2", "glucose")) %>%
  filter(is.finite(relative_brain_mass), is.finite(relative_brain_MR)) %>%
  mutate(label = species,
         highlight = case_when(grepl("Tetramorium", species) ~ "Ant (Coto treatment)",
                               species == "Homo sapiens" ~ "Human", species == "Mus musculus" ~ "Mouse",
                               TRUE ~ "Other mammal"))
# replace the human glucose-route point by the Heiss-derived one as an extra point
heiss_pt <- tibble(species = "Homo sapiens (Heiss glucose route)", relative_brain_mass = tab$relative_brain_mass[1],
                   relative_brain_MR = tab$relative_brain_MR_glucose[1], route = "glucose",
                   label = "Homo sapiens (Heiss glucose)", highlight = "Human")
mouse_pt <- tibble(species = "Mus musculus (Evo-M1 74.8)", relative_brain_mass = tab$relative_brain_mass[3],
                   relative_brain_MR = tab$relative_brain_MR_glucose[3], route = "glucose",
                   label = "Mus musculus (CMRgl 74.8)", highlight = "Mouse")
plot_df <- bind_rows(all_sp, heiss_pt, mouse_pt)
write.csv(plot_df, "tables/s5/s5b_figure_points.csv", row.names = FALSE)
lab <- annotate_lm(plot_df %>% filter(highlight != "Ant (Coto treatment)"), "relative_brain_mass", "relative_brain_MR")

p <- ggplot(plot_df, aes(x = relative_brain_mass, y = relative_brain_MR)) +
  geom_smooth(data = plot_df %>% filter(highlight != "Ant (Coto treatment)"), method = "lm", se = TRUE,
              colour = "steelblue", fill = "steelblue", alpha = 0.15, formula = y ~ x) +
  geom_point(aes(colour = highlight, shape = route), size = 3.2) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 30) +
  geom_text(data = lab, aes(x = Inf, y = -Inf, label = label), inherit.aes = FALSE, hjust = 1.05, vjust = -0.6, parse = TRUE, size = 3.6) +
  scale_colour_manual(values = c("Ant (Coto treatment)" = "#D55E00", Human = "#0072B2", Mouse = "#009E73", "Other mammal" = "grey40"), name = NULL) +
  scale_shape_manual(values = c(O2 = 16, glucose = 17), name = "Brain MR route") +
  labs(x = "Relative brain mass (brain / body mass)", y = "Relative brain metabolic rate (brain MR / whole-body MR)",
       caption = "Blue: OLS fit across mammals only (ant points excluded). Ant values at 31 C, not Q10-corrected.") +
  theme_classic(base_size = 14) + theme(legend.position = "bottom")
save_png_pdf(p, "figs/s5/s5b_relative_mass_vs_relative_MR", width = 9, height = 7)

cat("\n== s5b checks ==\n"); print(as.data.frame(checks), digits = 4)
cat("\n== s5b table ==\n")
print(tab %>% select(species, group, relative_brain_mass, wholebrain_CMRgl_umol_100g_min, brain_MR_W_glucose, brain_MR_W_O2,
                     body_MR_W, relative_brain_MR_glucose, relative_brain_MR_O2, brain_mass_specific_MR_mlO2_min_g,
                     brain_mass_specific_MR_mlO2_min_g_37C, ratio_to_human_mass_specific_O2, relative_brain_MR_O2_species_table, coto_published_relative_brain_MR) %>% as.data.frame(), digits = 4)
