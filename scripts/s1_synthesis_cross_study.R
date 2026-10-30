setwd(local({ d <- normalizePath(getwd()); while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d); d }))  # repo root

## ============================================================
## Study 1 synthesis — one view across s1a, s1b, s1c and s1d
##
## Every Study 1 comparison asks the same question of a different
## tissue property: does it co-vary with regional glucose metabolic
## rate (Heiss et al. 2004 rCMRGlc) across human brain regions?
##
##   s1a  stereology           cell numbers / densities / glia:neuron
##   s1b  snRNA-seq (Siletti)  cell-class proportions, E:I, astrocyte
##                             Type1/Type2, OPC:oligodendrocyte lineage
##   s1c  SV2A PET (Johansen)  synaptic density
##   s1d  ELISA (Sjölin)       GFAP, MBP, NFL, tau, UCHL1, total protein
##
## This script does NOT rerun those studies. It reads the REGION-LEVEL
## tables they write (one row per anatomy_id with rcmr_value) and
## recomputes one uniform statistic for every predictor — Spearman rho
## vs rCMRGlc (cor.test, exact = FALSE) with BH adjustment within each
## study x measurement block, which is the convention each study script
## uses — plus the per-predictor linear fit (R^2, p) that the study
## figures annotate. Where a study writes its own Spearman table the
## recomputed values are checked against it and the script stops on a
## mismatch, so the synthesis cannot silently drift from the studies.
##
## Outputs
##   tables/s1_synthesis/s1_cross_study_rcmrglc_correlations.csv
##   tables/s1_synthesis/s1_cross_study_region_coverage.csv
##   tables/s1_synthesis/s1_cross_study_predictor_inventory.csv
##   figs/s1_synthesis/s1_cross_study_rho_forest.{png,pdf}
##   figs/s1_synthesis/s1_cross_study_region_coverage.{png,pdf}
## ============================================================

library(tidyverse)
source("helpers/plot_settings.R")
source("helpers/read_heiss_rates.R")

for (d in c("figs/s1_synthesis", "tables/s1_synthesis")) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

read_required <- function(path) {
  if (!file.exists(path)) {
    stop("Missing study output (run the producing study first): ", path,
         call. = FALSE)
  }
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
           na.strings = c("", "NA"))
}

heiss <- read_heiss_rcmr(include_global_average = TRUE)

## ----------------------------------------------------------------
## 1. Predictor blocks: one data frame per study x measurement block,
##    standardised to anatomy_id | rcmr_value | <predictors>
## ----------------------------------------------------------------

blocks <- list()
add_block <- function(study, block, source_file, data, predictors, labels,
                      note = "") {
  stopifnot(all(c("anatomy_id", "rcmr_value") %in% names(data)))
  missing_pred <- setdiff(predictors, names(data))
  if (length(missing_pred)) {
    stop(source_file, " lacks predictor column(s): ",
         paste(missing_pred, collapse = ", "), call. = FALSE)
  }
  if (anyDuplicated(data$anatomy_id)) {
    stop(source_file, " has more than one row per anatomy_id.", call. = FALSE)
  }
  blocks[[length(blocks) + 1]] <<- list(
    study = study, block = block, source_file = source_file,
    data = data[, c("anatomy_id", "rcmr_value", predictors)],
    predictors = predictors, labels = labels, note = note
  )
}

## --- s1a stereology (same join as s1a_1) -------------------------------
s1a_file <- "data_intermediate/s1a_stereology_comparison.csv"
s1a <- read_required(s1a_file) %>%
  inner_join(as_tibble(heiss), by = "anatomy_id")
add_block(
  "s1a", "Cell numbers", s1a_file, s1a,
  c("Neuron_N", "Glia_N", "Astro_N", "Oligo_N", "Microglia_N"),
  c("Neurons (N)", "Glia (N)", "Astrocytes (N)", "Oligodendrocytes (N)",
    "Microglia (N)"),
  "Karlsen & Pakkenberg 2011 lobes; Morgan 2014 amygdala; counts per region"
)
add_block(
  "s1a", "Cell densities", s1a_file, s1a,
  c("NeurDensity", "GliaDensity", "AstroDensity", "OligoDensity",
    "MicroDensity", "GliaNeurDensityRatio"),
  c("Neuron density", "Glia density", "Astrocyte density",
    "Oligodendrocyte density", "Microglia density", "Glia : neuron ratio"),
  "cells per cm3; Hanson 2018 accumbens adds density-only rows"
)

## --- s1b transcriptomic cell composition --------------------------------
bind_scopes <- function(stem) {
  files <- paste0("data_analysis/", stem, "_by_region_",
                  c("telencephalon", "nontelencephalon"), ".csv")
  bind_rows(lapply(files, read_required))
}
s1b_neuronal <- bind_scopes("neuronal_celltype_proportions")
add_block(
  "s1b", "Neuronal class proportions",
  "data_analysis/neuronal_celltype_proportions_by_region_{telencephalon,nontelencephalon}.csv",
  s1b_neuronal,
  c("p_cells_Excitatory_projection", "p_cells_Inhibitory_interneuron",
    "p_cells_Inhibitory_projection_MSN"),
  c("Excitatory projection (prop.)", "Inhibitory interneuron (prop.)",
    "Inhibitory projection / MSN (prop.)"),
  "proportion of neuronal nuclei, donor-averaged (s1b_3_n)"
)
s1b_nonneuronal <- bind_scopes("nonneuronal_celltype_proportions")
add_block(
  "s1b", "Non-neuronal class proportions",
  "data_analysis/nonneuronal_celltype_proportions_by_region_{telencephalon,nontelencephalon}.csv",
  s1b_nonneuronal,
  c("p_cells_Astrocyte", "p_cells_Oligodendrocyte",
    "p_cells_Oligodendrocyte precursor", "p_cells_Microglia",
    "p_cells_Vascular"),
  c("Astrocyte (prop.)", "Oligodendrocyte (prop.)", "OPC (prop.)",
    "Microglia (prop.)", "Vascular (prop.)"),
  "proportion of non-neuronal nuclei, donor-averaged (s1b_3_nn)"
)
s1b_ei_file <- "data_analysis/original_broad_EI_ratio_all_regions_with_rcmr.csv"
add_block(
  "s1b", "Excitation : inhibition", s1b_ei_file, read_required(s1b_ei_file),
  c("EI_ratio_no_MSN", "EI_ratio_with_MSN"),
  c("E:I ratio (MSN excluded)", "E:I ratio (MSN included)"),
  "broad Siletti classes, all mapped regions (s1b_5)"
)
s1b_astro_file <- "data_analysis/type1_type2_astrocyte_composition_with_rcmr_by_region.csv"
add_block(
  "s1b", "Astrocyte types", s1b_astro_file, read_required(s1b_astro_file),
  c("comp_Astrocyte_Type1"),
  c("Type 1 astrocyte share"),
  "Type1 / (Type1 + Type2) astrocytes (s1b_6)"
)
s1b_opc_file <- "data_analysis/opc_oligo_lineage_composition_by_region_with_rcmr.csv"
add_block(
  "s1b", "Oligodendrocyte lineage", s1b_opc_file, read_required(s1b_opc_file),
  c("opc_frac_lineage_mean"),
  c("OPC / (OPC + oligodendrocyte)"),
  "primary s1b_7 metric"
)

## --- s1c synaptic density -------------------------------------------------
s1c_file <- "tables/s1c/synaptic_density_rcmrglc_matched_regions.csv"
s1c <- read_required(s1c_file) %>% rename(rcmr_value = rcmrglc_mean)
add_block(
  "s1c", "Synaptic density", s1c_file, s1c,
  c("synaptic_density_mean"),
  c("SV2A synaptic density"),
  "Johansen 2024 [11C]UCB-J Bmax, lobe totals + subcortical"
)

## --- s1d neuroglial proteins ---------------------------------------------
s1d_file <- "tables/s1d/neuroglial_protein_rcmrglc_matched_regions.csv"
s1d <- read_required(s1d_file) %>% rename(rcmr_value = rcmrglc_mean)
add_block(
  "s1d", "Tissue concentration", s1d_file, s1d,
  c("gfap_ug_per_g", "mbp_ug_per_g", "nfl_ug_per_g", "tau_ug_per_g",
    "uchl1_ug_per_g", "total_protein_mg_per_g", "gfap_to_nfl_ratio"),
  c("GFAP (µg/g)", "MBP (µg/g)", "NFL (µg/g)", "Tau (µg/g)", "UCHL1 (µg/g)",
    "Total protein (mg/g)", "GFAP : NFL ratio"),
  "Sjölin 2022 ELISA, regional geometric means"
)
add_block(
  "s1d", "Protein-normalized", s1d_file, s1d,
  c("gfap_ug_per_mg_protein", "mbp_ug_per_mg_protein",
    "nfl_ug_per_mg_protein", "tau_ug_per_mg_protein",
    "uchl1_ug_per_mg_protein"),
  c("GFAP (µg/mg protein)", "MBP (µg/mg protein)", "NFL (µg/mg protein)",
    "Tau (µg/mg protein)", "UCHL1 (µg/mg protein)"),
  "Sjölin 2022 ELISA, µg per mg total protein"
)

## ----------------------------------------------------------------
## 2. Uniform statistics per predictor
## ----------------------------------------------------------------

stat_one <- function(df, v) {
  ok <- complete.cases(df$rcmr_value, df[[v]]) & is.finite(df[[v]])
  x <- df[[v]][ok]; y <- df$rcmr_value[ok]
  if (length(x) < 4) {
    return(data.frame(n = length(x), rho = NA_real_, p_spearman = NA_real_,
                      lm_slope = NA_real_, lm_r_squared = NA_real_,
                      lm_p = NA_real_))
  }
  ct <- suppressWarnings(cor.test(y, x, method = "spearman", exact = FALSE))
  fit <- lm(y ~ x)
  data.frame(
    n = length(x),
    rho = unname(ct$estimate),
    p_spearman = ct$p.value,
    lm_slope = unname(coef(fit)[["x"]]),
    lm_r_squared = summary(fit)$r.squared,
    lm_p = summary(fit)$coefficients["x", "Pr(>|t|)"]
  )
}

results <- bind_rows(lapply(blocks, function(b) {
  bind_rows(lapply(seq_along(b$predictors), function(i) {
    v <- b$predictors[i]
    cbind(
      data.frame(study = b$study, block = b$block, predictor = v,
                 label = b$labels[i], stringsAsFactors = FALSE),
      stat_one(b$data, v)
    )
  })) %>%
    mutate(p_adj_BH_within_block = p.adjust(p_spearman, method = "BH"))
})) %>%
  mutate(
    p_adj_BH_all_s1 = p.adjust(p_spearman, method = "BH"),
    direction = case_when(is.na(rho) ~ "n/a", rho > 0 ~ "positive",
                          rho < 0 ~ "negative", TRUE ~ "zero"),
    sig_within_block = !is.na(p_adj_BH_within_block) & p_adj_BH_within_block < 0.05,
    sig_nominal = !is.na(p_spearman) & p_spearman < 0.05
  )

## ----------------------------------------------------------------
## 3. Self-checks against the studies' own Spearman tables
## ----------------------------------------------------------------

check_against <- function(label, recomputed_rho, reference_rho, tol = 1e-6) {
  if (is.na(recomputed_rho) || is.na(reference_rho) ||
      abs(recomputed_rho - reference_rho) > tol) {
    stop("Synthesis disagrees with the study's own statistic for ", label,
         ": recomputed rho = ", recomputed_rho, ", study rho = ",
         reference_rho, call. = FALSE)
  }
  message("  check OK: ", label, " rho = ", signif(recomputed_rho, 4))
}

message("Self-checks against study-written statistics:")
ref_s1b_opc <- read_required("data_analysis/opc_oligo_lineage_rcmr_correlations.csv") %>%
  filter(subset == "All regions", metric == "opc_frac_lineage_mean")
check_against("s1b OPC fraction (s1b_7)",
              results$rho[results$predictor == "opc_frac_lineage_mean"],
              ref_s1b_opc$spearman_rho)

ref_s1b_astro <- read_required("data_analysis/type1_type2_astrocyte_composition_rcmr_correlations.csv")
check_against("s1b Type 1 astrocyte share (s1b_6)",
              results$rho[results$predictor == "comp_Astrocyte_Type1"],
              ref_s1b_astro$spearman_rho_type1)

ref_s1d <- read_required("tables/s1d/spearman_rcmrglc_tissue_concentration.csv")
for (v in c("gfap_ug_per_g", "mbp_ug_per_g", "tau_ug_per_g")) {
  check_against(paste0("s1d ", v, " (s1d_1)"),
                results$rho[results$predictor == v],
                ref_s1d$rho[ref_s1d$variable == v])
}

## ----------------------------------------------------------------
## 4. Region coverage: which Heiss regions each study pairs
## ----------------------------------------------------------------

coverage_long <- bind_rows(lapply(blocks, function(b) {
  b$data %>%
    filter(!is.na(rcmr_value)) %>%
    distinct(anatomy_id) %>%
    mutate(study = b$study)
})) %>% distinct()

study_levels <- c("s1a", "s1b", "s1c", "s1d")
study_names <- c(s1a = "s1a stereology", s1b = "s1b transcriptomics",
                 s1c = "s1c synaptic density", s1d = "s1d neuroglial proteins")

# Composite analysis regions (s1d thalamus / midbrain) are not native Heiss
# rows; label them from the registry.
registry <- read_anatomy_registry()
region_label <- function(id) {
  lab <- heiss$rcmr_term[match(id, heiss$anatomy_id)]
  lab[is.na(lab)] <- paste0(
    registry$canonical_label[match(id[is.na(lab)], registry$anatomy_id)],
    " (composite)"
  )
  lab
}
region_rate <- function(id) {
  r <- heiss$rcmr_value[match(id, heiss$anatomy_id)]
  missing <- is.na(r)
  if (any(missing)) {
    # take the composite rate from whichever block carries it
    all_rates <- bind_rows(lapply(blocks, function(b)
      b$data[, c("anatomy_id", "rcmr_value")])) %>%
      distinct() %>% filter(!is.na(rcmr_value))
    r[missing] <- all_rates$rcmr_value[match(id[missing], all_rates$anatomy_id)]
  }
  r
}

coverage <- coverage_long %>%
  mutate(covered = TRUE) %>%
  complete(anatomy_id, study = study_levels, fill = list(covered = FALSE)) %>%
  mutate(region = region_label(anatomy_id),
         rcmr_value = region_rate(anatomy_id)) %>%
  arrange(rcmr_value)

coverage_wide <- coverage %>%
  pivot_wider(names_from = study, values_from = covered) %>%
  mutate(n_studies = rowSums(across(all_of(study_levels)))) %>%
  arrange(desc(n_studies), rcmr_value)
write.csv(coverage_wide,
          "tables/s1_synthesis/s1_cross_study_region_coverage.csv",
          row.names = FALSE)

shared_all <- coverage_wide %>% filter(n_studies == length(study_levels))
message(nrow(coverage_wide), " regions appear in at least one Study 1 ",
        "comparison; ", nrow(shared_all), " appear in all four: ",
        paste(shared_all$region, collapse = ", "))

## ----------------------------------------------------------------
## 5. Write results + predictor inventory
## ----------------------------------------------------------------

results_out <- results %>%
  mutate(study_name = unname(study_names[study])) %>%
  select(study, study_name, block, predictor, label, n, rho, p_spearman,
         p_adj_BH_within_block, p_adj_BH_all_s1, direction,
         sig_nominal, sig_within_block, lm_slope, lm_r_squared, lm_p) %>%
  arrange(study, block, p_spearman)
write.csv(results_out,
          "tables/s1_synthesis/s1_cross_study_rcmrglc_correlations.csv",
          row.names = FALSE)

inventory <- bind_rows(lapply(blocks, function(b) {
  data.frame(study = b$study, block = b$block, source_file = b$source_file,
             n_regions = sum(!is.na(b$data$rcmr_value)),
             n_predictors = length(b$predictors), note = b$note,
             stringsAsFactors = FALSE)
}))
write.csv(inventory,
          "tables/s1_synthesis/s1_cross_study_predictor_inventory.csv",
          row.names = FALSE)

print(as_tibble(results_out) %>%
        select(study, block, label, n, rho, p_spearman,
               p_adj_BH_within_block) %>%
        mutate(across(where(is.numeric), ~ signif(.x, 3))),
      n = Inf, width = 200)

## ----------------------------------------------------------------
## 6. Figure: Spearman rho forest, grouped by study / block
## ----------------------------------------------------------------

study_palette <- c(s1a = "#D55E00", s1b = "#0072B2", s1c = "#009E73",
                   s1d = "#CC79A7")

forest_df <- results_out %>%
  filter(!is.na(rho)) %>%
  mutate(
    panel = factor(paste0(study_name, "\n", block),
                   levels = unique(paste0(study_name, "\n", block))),
    label = factor(label, levels = rev(unique(label))),
    sig_mark = case_when(sig_within_block ~ "BH q < 0.05 (within block)",
                         sig_nominal ~ "p < 0.05 (nominal)",
                         TRUE ~ "not significant"),
    # Fisher-z approximate 95% interval for display only (n small)
    z = atanh(pmin(pmax(rho, -0.999), 0.999)),
    se = 1 / sqrt(pmax(n - 3, 1)),
    lo = tanh(z - 1.96 * se), hi = tanh(z + 1.96 * se),
    n_label = paste0("n=", n)
  )

p_forest <- ggplot(forest_df, aes(x = rho, y = label, color = study)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0,
                alpha = 0.5, linewidth = 0.5) +
  geom_point(aes(shape = sig_mark), size = 3) +
  geom_text(aes(x = 1.05, label = n_label), hjust = 0, size = 2.6,
            color = "grey30", show.legend = FALSE) +
  facet_grid(panel ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_color_manual(values = study_palette, labels = study_names,
                     name = "Study") +
  scale_shape_manual(values = c("BH q < 0.05 (within block)" = 16,
                                "p < 0.05 (nominal)" = 1,
                                "not significant" = 4),
                     name = "Significance") +
  scale_x_continuous(limits = c(-1, 1.3), breaks = seq(-1, 1, 0.5)) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Study 1: which regional tissue properties track glucose metabolic rate?",
    subtitle = "Spearman rho vs. Heiss et al. 2004 rCMRGlc, recomputed uniformly from each study's region-level table",
    x = "Spearman rho (predictor vs. rCMRGlc)",
    y = NULL,
    caption = paste0(
      "Bars = approximate 95% interval (Fisher z); n = regions paired with a ",
      "Heiss rate.\nBH adjustment within each study x block, as in the study ",
      "scripts.\nRegion sets differ between studies (see region-coverage ",
      "figure): rho values are comparable in sign and size\nbut are not ",
      "estimated on the same regions."
    )
  ) +
  theme_project(base_size = 11) +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = 8.5),
    strip.background = element_rect(fill = "grey95", color = NA),
    panel.spacing.y = unit(0.4, "lines"),
    legend.position = "right",
    plot.caption = element_text(size = 8, hjust = 0, color = "grey40"),
    plot.subtitle = element_text(size = 9.5, color = "grey40")
  )

ggsave("figs/s1_synthesis/s1_cross_study_rho_forest.png", p_forest,
       width = 11, height = 11.5, units = "in", dpi = 300, bg = "white")
ggsave("figs/s1_synthesis/s1_cross_study_rho_forest.pdf", p_forest,
       width = 11, height = 11.5, units = "in", device = cairo_pdf)

## ----------------------------------------------------------------
## 7. Figure: region coverage heatmap
## ----------------------------------------------------------------

cov_df <- coverage %>%
  mutate(
    region = factor(region, levels = coverage_wide$region),
    study_name = factor(unname(study_names[study]), levels = study_names),
    rate_label = sprintf("%.1f", rcmr_value)
  )

p_cov <- ggplot(cov_df, aes(x = study_name, y = region)) +
  geom_tile(aes(fill = ifelse(covered, study, NA_character_)),
            color = "white", linewidth = 0.8) +
  geom_text(data = cov_df %>% distinct(region, rate_label),
            aes(x = 0.35, y = region, label = rate_label),
            hjust = 1, size = 2.8, color = "grey30", inherit.aes = FALSE) +
  scale_fill_manual(values = study_palette, na.value = "grey93",
                    guide = "none") +
  scale_x_discrete(position = "top", expand = expansion(add = c(0.9, 0.5))) +
  scale_y_discrete(limits = rev) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Study 1 region coverage",
    subtitle = "Heiss regions paired by each study (left: rCMRGlc, µmol/100 g/min); sorted by number of studies, then rate",
    x = NULL, y = NULL
  ) +
  theme_project(base_size = 11) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x.top = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, color = "grey40"))

ggsave("figs/s1_synthesis/s1_cross_study_region_coverage.png", p_cov,
       width = 8, height = 8.5, units = "in", dpi = 300, bg = "white")
ggsave("figs/s1_synthesis/s1_cross_study_region_coverage.pdf", p_cov,
       width = 8, height = 8.5, units = "in", device = cairo_pdf)

message("Study 1 synthesis done: tables/s1_synthesis/, figs/s1_synthesis/.")
