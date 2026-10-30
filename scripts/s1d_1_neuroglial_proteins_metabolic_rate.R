setwd(local({ d <- normalizePath(getwd()); while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d); d }))  # repo root

## ============================================================
## Study s1d_1 — Neuroglial protein concentrations vs. rCMRGlc
##
## Replicates the design of s1a_1 (stereological cell counts and
## densities vs. rCMRGlc) with regional protein concentrations in
## place of cell densities.
##
## Sjölin et al. 2022, Additional file 2: ELISA concentrations of five
##   clinically important neuroglial proteins — GFAP, MBP, NFL, tau
##   (total) and UCHL1 — plus total protein, in post-mortem tissue from
##   17 CNS regions x 10 donors (168 samples; hippocampus missing for
##   two donors). Donors 24-50 y, no neurological disease.
##   Source: Mol Brain 15:52. DOI: 10.1186/s13041-022-00935-6
##   Local source copy: data_raw/Sjolin_etal_2022_SupplementaryData.csv
##   (synced from Evo-M1-Trait-Data/Sjölin_etal_2022/). The *_ug_per_g
##   and total_protein_mg_per_g columns are the Evo-M1 item's conversion
##   of the published supernatant concentrations to wet-tissue units
##   (10 mL lysing buffer per g tissue), i.e. the units of the paper's
##   Table 2.
##
## Heiss et al. 2004 Table 1: regional rCMRGlc (µmol/100g/min),
##   18F-FDG HRRT PET, N = 9 healthy adults.
##   Analysis input: data_intermediate/heiss_2004_regions.csv
##
## Design (parallel to s1a_1):
##   s1a_1 counts    -> here "tissue concentration": µg analyte per g wet
##                      tissue (plus total protein, mg/g)
##   s1a_1 densities -> here "protein-normalized": µg analyte per mg total
##                      protein (composition of the soluble proteome)
##   glia:neuron ratio -> here GFAP:NFL ratio (astroglial : axonal marker)
##
## Observation level: the source is one row per donor x region. Regions
## are summarised as the geometric mean across donors (distributions are
## log-normal; the paper itself reports geometric means), one value per
## analysis region. Source regions sharing one Heiss rate (the four lobar
## white-matter samples -> centrum semiovale) are pooled into one analysis
## region, so every point is one independent Heiss rate, as in s1a_1.
##
## Region pairing is defined in metadata/anatomy/s1d_sjolin_region_map.csv
## (validated by 00_prepare_heiss_rates.R). Pons, medulla and cervical
## spinal cord have no Heiss rate and are excluded from the metabolic
## analyses; they remain in the regional summary table for s1d_2.
## ============================================================

library(tidyverse)
library(ggpmisc)

source("helpers/plot_settings.R")
source("helpers/read_heiss_rates.R")

for (d in c("figs/s1d", "tables/s1d", "data_intermediate")) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

## ----------------------------------------------------------------
## 1. Load data
## ----------------------------------------------------------------

sjolin_source <- read.csv(
  "data_raw/Sjolin_etal_2022_SupplementaryData.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA")
)

required_sjolin_columns <- c(
  "donor", "region_code", "region_as_published",
  "GFAP_ug_per_g", "MBP_ug_per_g", "NFL_ug_per_g", "TAU_ug_per_g",
  "UCHL1_ug_per_g", "total_protein_mg_per_g", "observation_level"
)
missing_sjolin_columns <- setdiff(required_sjolin_columns, names(sjolin_source))
if (length(missing_sjolin_columns)) {
  stop(
    "Sjölin source is missing required column(s): ",
    paste(missing_sjolin_columns, collapse = ", "),
    call. = FALSE
  )
}
if (!all(sjolin_source$observation_level == "sample (one donor x one CNS region)")) {
  stop("Sjölin source rows are not all donor x region samples.", call. = FALSE)
}

# Import boundary: descriptive analysis names, units in the name.
samples <- sjolin_source %>%
  transmute(
    donor = as.integer(donor),
    region_code,
    source_region = trimws(region_as_published),
    gfap_ug_per_g  = as.numeric(GFAP_ug_per_g),
    mbp_ug_per_g   = as.numeric(MBP_ug_per_g),
    nfl_ug_per_g   = as.numeric(NFL_ug_per_g),
    tau_ug_per_g   = as.numeric(TAU_ug_per_g),
    uchl1_ug_per_g = as.numeric(UCHL1_ug_per_g),
    total_protein_mg_per_g = as.numeric(total_protein_mg_per_g)
  ) %>%
  mutate(
    # Protein-normalized concentrations: µg analyte per mg total protein.
    gfap_ug_per_mg_protein  = gfap_ug_per_g  / total_protein_mg_per_g,
    mbp_ug_per_mg_protein   = mbp_ug_per_g   / total_protein_mg_per_g,
    nfl_ug_per_mg_protein   = nfl_ug_per_g   / total_protein_mg_per_g,
    tau_ug_per_mg_protein   = tau_ug_per_g   / total_protein_mg_per_g,
    uchl1_ug_per_mg_protein = uchl1_ug_per_g / total_protein_mg_per_g,
    # Astroglial : axonal marker ratio (analogue of the s1a glia:neuron ratio).
    gfap_to_nfl_ratio = gfap_ug_per_g / nfl_ug_per_g
  )

## ----------------------------------------------------------------
## 2. Region mapping: Sjölin region -> project anatomy_id -> Heiss rate
## ----------------------------------------------------------------

map_path <- "metadata/anatomy/s1d_sjolin_region_map.csv"
study_map <- read_study_anatomy_map(map_path)
region_map <- derive_study_heiss_rates(study_map)

# Every source region must be in the map and every mapped region in the data.
unmapped <- setdiff(unique(samples$source_region), region_map$source_region)
orphans  <- setdiff(region_map$source_region, unique(samples$source_region))
if (length(unmapped) || length(orphans)) {
  stop(
    "s1d map / source mismatch. Unmapped source region(s): ",
    paste(unmapped, collapse = ", "),
    ". Map rows without data: ", paste(orphans, collapse = ", "),
    call. = FALSE
  )
}

# One relationship label per source region (composites repeat it).
relationship_by_source <- study_map %>%
  distinct(source_region, relationship)
if (anyDuplicated(relationship_by_source$source_region)) {
  stop(map_path, " gives more than one relationship to a source region.",
       call. = FALSE)
}

samples_matched <- samples %>%
  inner_join(region_map, by = "source_region") %>%
  inner_join(relationship_by_source, by = "source_region")

## ----------------------------------------------------------------
## 3. Regional summary: geometric mean across donors (and across pooled
##    source regions) per analysis region
## ----------------------------------------------------------------

measure_columns <- c(
  "gfap_ug_per_g", "mbp_ug_per_g", "nfl_ug_per_g", "tau_ug_per_g",
  "uchl1_ug_per_g", "total_protein_mg_per_g",
  "gfap_ug_per_mg_protein", "mbp_ug_per_mg_protein", "nfl_ug_per_mg_protein",
  "tau_ug_per_mg_protein", "uchl1_ug_per_mg_protein",
  "gfap_to_nfl_ratio"
)

geo_mean <- function(x) exp(mean(log(x[!is.na(x) & x > 0])))
# Multiplicative (geometric) SD: region mean x/÷ this factor spans ±1 SD on
# the log scale.
geo_sd <- function(x) {
  x <- x[!is.na(x) & x > 0]
  if (length(x) < 2) return(NA_real_)
  exp(sd(log(x)))
}

regional_summary <- samples_matched %>%
  group_by(anatomy_id, anatomy_group, rcmr_value, rcmr_sd, n_heiss_components) %>%
  summarise(
    source_regions = paste(sort(unique(source_region)), collapse = " | "),
    relationship   = paste(sort(unique(relationship)), collapse = " | "),
    n_samples = n(),
    n_donors  = n_distinct(donor),
    across(all_of(measure_columns), list(gmean = geo_mean, gsd = geo_sd),
           .names = "{.col}__{.fn}"),
    .groups = "drop"
  ) %>%
  arrange(desc(!is.na(rcmr_value)), rcmr_value)

# Wide table with geometric means under the plain measure names (the s1a_1
# analysis_df analogue) plus the geometric SDs as *_gsd columns.
regional_means <- regional_summary %>%
  rename_with(~ sub("__gmean$", "", .x), ends_with("__gmean")) %>%
  rename_with(~ sub("__gsd$", "_gsd", .x), ends_with("__gsd"))

write.csv(
  regional_means,
  "data_intermediate/s1d_neuroglial_protein_regional_summary.csv",
  row.names = FALSE, na = ""
)

# Analysis set: regions with a Heiss rate (as in s1a_1, rows lacking the
# comparison variables are dropped; here the rows lacking a rate).
analysis_df <- regional_means %>% filter(!is.na(rcmr_value))

message("s1d_1: ", nrow(analysis_df), " analysis regions with a Heiss rate; ",
        sum(is.na(regional_means$rcmr_value)), " source regions without one.")

##########################################################################
# Correlations between rCMRGlc and (1) tissue concentrations and
# (2) protein-normalized concentrations
##########################################################################

tissue_concentration <- c(
  "gfap_ug_per_g",
  "mbp_ug_per_g",
  "nfl_ug_per_g",
  "tau_ug_per_g",
  "uchl1_ug_per_g",
  "total_protein_mg_per_g",
  "gfap_to_nfl_ratio"
)

protein_normalized <- c(
  "gfap_ug_per_mg_protein",
  "mbp_ug_per_mg_protein",
  "nfl_ug_per_mg_protein",
  "tau_ug_per_mg_protein",
  "uchl1_ug_per_mg_protein",
  "gfap_to_nfl_ratio"
)

# Display labels for facets / tables
measure_labels <- c(
  gfap_ug_per_g  = "GFAP (µg/g)",
  mbp_ug_per_g   = "MBP (µg/g)",
  nfl_ug_per_g   = "NFL (µg/g)",
  tau_ug_per_g   = "Tau (µg/g)",
  uchl1_ug_per_g = "UCHL1 (µg/g)",
  total_protein_mg_per_g = "Total protein (mg/g)",
  gfap_ug_per_mg_protein  = "GFAP (µg/mg protein)",
  mbp_ug_per_mg_protein   = "MBP (µg/mg protein)",
  nfl_ug_per_mg_protein   = "NFL (µg/mg protein)",
  tau_ug_per_mg_protein   = "Tau (µg/mg protein)",
  uchl1_ug_per_mg_protein = "UCHL1 (µg/mg protein)",
  gfap_to_nfl_ratio = "GFAP : NFL ratio"
)

###########################
## MULTIVARIATE REGRESSION
###########################
# Note (as in s1a_1): this is only to get an overall sense of whether the
# protein concentrations jointly predict rCMRGlc. With ~11 regions and five
# correlated predictors the fit is close to saturated and must be
# interpreted with great caution. Total protein and the ratio are left out
# so the two fits have the same five analyte predictors.

five_proteins_tissue <- tissue_concentration[1:5]
five_proteins_normalized <- protein_normalized[1:5]

fit_tissue <- lm(
  reformulate(five_proteins_tissue, response = "rcmr_value"),
  data = analysis_df
)
print(summary(fit_tissue))

fit_normalized <- lm(
  reformulate(five_proteins_normalized, response = "rcmr_value"),
  data = analysis_df
)
print(summary(fit_normalized))

lm_table <- function(fit, model_label) {
  co <- summary(fit)$coefficients
  data.frame(
    model = model_label,
    term = rownames(co),
    estimate = co[, "Estimate"],
    std_error = co[, "Std. Error"],
    t_value = co[, "t value"],
    p_value = co[, "Pr(>|t|)"],
    n = nobs(fit),
    residual_df = fit$df.residual,
    adj_r_squared = summary(fit)$adj.r.squared,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}
bind_rows(
  lm_table(fit_tissue, "rcmr ~ five proteins (µg/g tissue)"),
  lm_table(fit_normalized, "rcmr ~ five proteins (µg/mg total protein)")
) %>%
  write.csv("tables/s1d/multivariate_lm_rcmrglc_five_proteins.csv",
            row.names = FALSE)

###########################
## CORRELATIONS (LOOP)
###########################

spearman_table <- function(df, variables) {
  do.call(
    rbind,
    lapply(variables, function(v) {
      ct <- cor.test(
        df$rcmr_value,
        df[[v]],
        method = "spearman",
        exact = FALSE
      )
      data.frame(
        variable = v,
        label    = unname(measure_labels[v]),
        rho      = unname(ct$estimate),
        p_value  = ct$p.value,
        n        = sum(complete.cases(df$rcmr_value, df[[v]])),
        stringsAsFactors = FALSE
      )
    })
  ) %>%
    dplyr::arrange(p_value) %>%
    dplyr::mutate(p_adj_BH = p.adjust(p_value, method = "BH"))
}

cor_results_tissue <- spearman_table(analysis_df, tissue_concentration)
cor_results_normalized <- spearman_table(analysis_df, protein_normalized)

print(cor_results_tissue)
print(cor_results_normalized)

write.csv(cor_results_tissue,
          "tables/s1d/spearman_rcmrglc_tissue_concentration.csv",
          row.names = FALSE)
write.csv(cor_results_normalized,
          "tables/s1d/spearman_rcmrglc_protein_normalized.csv",
          row.names = FALSE)

## Sensitivity 1: direct anatomical matches only (drop proxy and composite
## pairings: lobar white matter -> centrum semiovale, cerebellum ->
## cerebellar cortex, thalamus and midbrain composites).
direct_df <- analysis_df %>%
  filter(relationship %in% c("exact", "source_synonym"))
message("s1d_1 sensitivity: ", nrow(direct_df), " direct-match regions.")

## Sensitivity 1b: grey matter only. The white-matter regions (centrum
## semiovale, internal capsule) have the lowest rates and by far the highest
## MBP, so they can dominate any fit; tissue_class comes from the registry.
registry_tissue <- read_anatomy_registry() %>%
  select(anatomy_id, tissue_class)
grey_df <- analysis_df %>%
  inner_join(registry_tissue, by = "anatomy_id") %>%
  filter(tissue_class != "white_matter")
message("s1d_1 sensitivity: ", nrow(grey_df), " grey-matter regions.")

bind_rows(
  spearman_table(direct_df, tissue_concentration) %>%
    mutate(measure_set = "tissue_concentration",
           region_scope = "direct matches only (exact / source_synonym)"),
  spearman_table(direct_df, protein_normalized) %>%
    mutate(measure_set = "protein_normalized",
           region_scope = "direct matches only (exact / source_synonym)"),
  spearman_table(grey_df, tissue_concentration) %>%
    mutate(measure_set = "tissue_concentration",
           region_scope = "grey matter only (white-matter regions excluded)"),
  spearman_table(grey_df, protein_normalized) %>%
    mutate(measure_set = "protein_normalized",
           region_scope = "grey matter only (white-matter regions excluded)")
) %>%
  write.csv("tables/s1d/spearman_rcmrglc_region_scope_sensitivity.csv",
            row.names = FALSE)

## Sensitivity 2: within-donor correlations. The regional geometric means
## discard the donor structure; here each donor's own regional values are
## correlated with rCMRGlc (lobar white-matter samples averaged within the
## donor first), and the distribution of per-donor rho is summarised.
donor_region <- samples_matched %>%
  filter(!is.na(rcmr_value)) %>%
  group_by(donor, anatomy_id, rcmr_value) %>%
  summarise(across(all_of(measure_columns), geo_mean), .groups = "drop")

within_donor <- bind_rows(lapply(
  unique(c(tissue_concentration, protein_normalized)), function(v) {
    per_donor <- donor_region %>%
      group_by(donor) %>%
      summarise(
        n_regions = sum(complete.cases(rcmr_value, .data[[v]])),
        rho = suppressWarnings(cor(rcmr_value, .data[[v]],
                                   method = "spearman",
                                   use = "complete.obs")),
        .groups = "drop"
      )
    wt <- suppressWarnings(wilcox.test(per_donor$rho, mu = 0, exact = FALSE))
    data.frame(
      variable = v,
      label = unname(measure_labels[v]),
      n_donors = nrow(per_donor),
      regions_per_donor = paste(range(per_donor$n_regions), collapse = "-"),
      median_rho = median(per_donor$rho),
      min_rho = min(per_donor$rho),
      max_rho = max(per_donor$rho),
      n_donors_rho_positive = sum(per_donor$rho > 0),
      wilcoxon_p_rho_vs_0 = wt$p.value,
      stringsAsFactors = FALSE
    )
  }
))
print(within_donor)
write.csv(within_donor, "tables/s1d/spearman_rcmrglc_within_donor.csv",
          row.names = FALSE)

###########################
## FACETED SCATTERPLOT
###########################

# Region palette from the anatomy registry (stable colours by anatomy_id).
registry <- read_anatomy_registry()
region_pal <- setNames(registry$color_hex, registry$anatomy_id)
region_pal <- region_pal[analysis_df$anatomy_id]
names(region_pal) <- analysis_df$anatomy_group
if (anyNA(region_pal)) {
  stop("Missing registry colour for region(s): ",
       paste(names(region_pal)[is.na(region_pal)], collapse = ", "),
       call. = FALSE)
}
region_levels <- analysis_df %>%
  arrange(rcmr_value) %>% pull(anatomy_group)

make_facet_plot <- function(variables, x_label, title, subtitle) {
  # Long format: one row per region x measure (the s1a_1 "prop" column).
  plot_means <- analysis_df %>%
    dplyr::select(anatomy_group, rcmr_value, all_of(variables)) %>%
    tidyr::pivot_longer(
      cols = all_of(variables),
      names_to = "measure",
      values_to = "prop"
    )
  plot_gsd <- analysis_df %>%
    dplyr::select(anatomy_group, all_of(paste0(variables, "_gsd"))) %>%
    tidyr::pivot_longer(
      cols = -anatomy_group,
      names_to = "measure",
      values_to = "gsd"
    ) %>%
    dplyr::mutate(measure = sub("_gsd$", "", measure))
  plot_df <- plot_means %>%
    dplyr::left_join(plot_gsd, by = c("anatomy_group", "measure")) %>%
    dplyr::mutate(
      measure = factor(unname(measure_labels[measure]),
                       levels = unname(measure_labels[variables])),
      anatomy_group = factor(anatomy_group, levels = region_levels)
    )

  # One lm per facet; BH-adjusted asterisks (as in s1a_1).
  model_results <- plot_df %>%
    dplyr::filter(!is.na(prop), !is.na(rcmr_value)) %>%
    dplyr::group_by(measure) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      fit = purrr::map(data, ~ lm(rcmr_value ~ prop, data = .x)),
      n = purrr::map_int(data, nrow),
      slope = purrr::map_dbl(fit, ~ coef(.x)[["prop"]]),
      r_squared = purrr::map_dbl(fit, ~ summary(.x)$r.squared),
      p_value = purrr::map_dbl(
        fit, ~ summary(.x)$coefficients["prop", "Pr(>|t|)"]
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      p_adj = p.adjust(p_value, method = "BH"),
      sig_label = dplyr::case_when(
        p_adj < 0.001 ~ "***",
        p_adj < 0.01  ~ "**",
        p_adj < 0.05  ~ "*",
        TRUE ~ ""
      )
    ) %>%
    dplyr::select(measure, n, slope, r_squared, p_value, p_adj, sig_label)

  sig_labels <- model_results %>% dplyr::filter(sig_label != "")

  dims <- facet_dims(length(variables))

  p <- ggplot(plot_df, aes(x = prop, y = rcmr_value, color = anatomy_group)) +
    geom_errorbar(
      aes(xmin = prop / gsd, xmax = prop * gsd),
      orientation = "y", width = 0, alpha = 0.35, linewidth = 0.5
    ) +
    geom_point(size = 2.8, alpha = 0.85) +
    geom_smooth(aes(group = 1), method = "lm", se = TRUE,
                color = "steelblue", linewidth = 0.8) +
    stat_poly_eq(
      aes(label = after_stat(paste(rr.label, p.value.label, sep = "*\", \"*")),
          group = 1),
      formula = y ~ x,
      parse = TRUE,
      label.x = "right",
      label.y = "top",
      size = 3.2,
      color = "black"
    ) +
    scale_x_continuous(n.breaks = 4) +
    facet_wrap(~ measure, scales = "free_x", ncol = dims$ncol) +
    scale_color_manual(values = region_pal) +
    labs(
      title = title,
      subtitle = subtitle,
      x = x_label,
      y = "rCMRGlc (µmol/100 g/min)",
      color = "Region",
      caption = paste0(
        "Points = regional geometric means across donors (Sjölin et al. 2022, ",
        "n = 8-40 samples per region); horizontal bars = x/÷ geometric SD.\n",
        "rCMRGlc from Heiss et al. 2004. Asterisks = BH-adjusted p of the per-panel ",
        "linear fit (* <0.05, ** <0.01, *** <0.001)."
      )
    ) +
    theme_facet_compact(base_size = 12) +
    theme(legend.position = "right")

  p <- p +
    geom_text(
      data = sig_labels,
      aes(x = -Inf, y = Inf, label = sig_label),
      inherit.aes = FALSE,
      hjust = -0.5,
      vjust = 1.2,
      size = 6
    )

  list(plot = p, model_results = model_results, dims = dims)
}

############# -------- -------- #############
# tissue concentrations (s1a_1 "counts" analogue)

res_tissue <- make_facet_plot(
  tissue_concentration,
  x_label = "Protein concentration in wet tissue",
  title = "Neuroglial protein concentration vs. glucose metabolic rate across brain regions",
  subtitle = "Sjölin et al. 2022 (ELISA, µg per g wet tissue) × Heiss et al. 2004 (FDG PET)"
)
p_tissue <- res_tissue$plot
p_tissue

ggsave("figs/s1d/p_tissue_concentration.pdf", plot = p_tissue,
       width = res_tissue$dims$width, height = res_tissue$dims$height,
       units = "in", device = cairo_pdf)
ggsave("figs/s1d/p_tissue_concentration.jpg", plot = p_tissue,
       width = res_tissue$dims$width, height = res_tissue$dims$height,
       units = "in", dpi = 300, device = "jpeg", bg = "white")

# protein-normalized concentrations (s1a_1 "densities" analogue)

res_normalized <- make_facet_plot(
  protein_normalized,
  x_label = "Protein concentration normalized to total protein",
  title = "Protein-normalized neuroglial protein content vs. glucose metabolic rate",
  subtitle = "Sjölin et al. 2022 (ELISA, µg per mg total protein) × Heiss et al. 2004 (FDG PET)"
)
p_normalized <- res_normalized$plot
p_normalized

ggsave("figs/s1d/p_protein_normalized.pdf", plot = p_normalized,
       width = res_normalized$dims$width, height = res_normalized$dims$height,
       units = "in", device = cairo_pdf)
ggsave("figs/s1d/p_protein_normalized.jpg", plot = p_normalized,
       width = res_normalized$dims$width, height = res_normalized$dims$height,
       units = "in", dpi = 300, device = "jpeg", bg = "white")

# Per-panel linear fits behind the asterisks
bind_rows(
  res_tissue$model_results %>% mutate(measure_set = "tissue_concentration"),
  res_normalized$model_results %>% mutate(measure_set = "protein_normalized")
) %>%
  mutate(measure = as.character(measure)) %>%
  write.csv("tables/s1d/lm_rcmrglc_per_measure.csv", row.names = FALSE)

## ----------------------------------------------------------------
## Matched-region table for inspection
## ----------------------------------------------------------------

analysis_df %>%
  dplyr::select(
    anatomy_id, anatomy_group, source_regions, relationship,
    n_heiss_components, n_samples, n_donors,
    rcmrglc_mean = rcmr_value, rcmrglc_sd = rcmr_sd,
    all_of(tissue_concentration), all_of(setdiff(protein_normalized, "gfap_to_nfl_ratio"))
  ) %>%
  arrange(rcmrglc_mean) %>%
  write.csv("tables/s1d/neuroglial_protein_rcmrglc_matched_regions.csv",
            row.names = FALSE)

message("Study s1d_1 done. Figures in figs/s1d/, tables in tables/s1d/, ",
        "regional summary in data_intermediate/.")
