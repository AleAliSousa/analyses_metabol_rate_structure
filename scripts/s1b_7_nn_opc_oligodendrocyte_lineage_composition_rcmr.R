setwd(local({ d <- normalizePath(getwd()); while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d); d }))  # repo root (portable; replaces hardcoded path -- see helpers/project_root.R)

###################################################################
# Study 1b (s1b_7): OPC share of the oligodendrocyte lineage vs rCMRGlc
###################################################################
# Question. Does the proportion of oligodendrocyte precursor cells (OPCs)
# among oligodendrocyte-lineage cells -- OPC / (OPC + Oligodendrocyte) --
# track regional cerebral glucose metabolic rate (rCMRGlc, Heiss et al.
# 2004)?
#
# Why OPCs and oligodendrocytes belong in one denominator.
#   OPCs and oligodendrocytes are the two ends of a single developmental
#   lineage (OPC -> committed precursor (COP) -> newly formed -> mature
#   myelinating oligodendrocyte; Marques et al. 2016 Science; Jakel et al.
#   2019 Nature; Siletti et al. 2023 Science, whose superclusters we use).
#   All members share the lineage transcription factors OLIG1/2 and SOX10;
#   OPCs are distinguished by PDGFRA / CSPG4 (NG2), mature cells by
#   PLP1 / MBP / MOG. Adult OPCs are not a vestigial population: they
#   remain proliferative throughout life and are the sole source of new
#   oligodendrocytes in the adult CNS (Young et al. 2013 Neuron; Hughes et
#   al. 2018 Nat Neurosci), so the OPC : oligodendrocyte ratio is routinely
#   used as an index of lineage maturity / standing precursor reserve in the
#   myelination and remyelination literature. Siletti et al. themselves
#   present OPC and oligodendrocyte types side by side when describing
#   regional variation (the panel that s1b_6 mirrors for astrocytes).
#
# Why a lineage denominator rather than all nonneuronal cells.
#   s1b_3_nn expresses OPCs as a fraction of ALL nonneuronal cells. That
#   fraction mixes two things: (i) how much of the sample is oligodendrocyte
#   lineage at all, which is dominated by white-matter admixture in the
#   dissection and by ventricle / vascular sampling, and (ii) the precursor
#   vs mature balance within the lineage. Restricting the denominator to
#   the lineage isolates (ii) and asks a narrower, biologically closed
#   question: of the cells that are oligodendrocyte lineage, what share is
#   still precursor? Because the two fractions sum to 1 the analysis is a
#   two-part composition, exactly like Type 1 / Type 2 astrocytes in s1b_6.
#
# Known limitation (tested below, see "White-matter admixture check").
#   The lineage denominator removes astrocytes / microglia / vascular cells
#   from the ratio, but it cannot remove white-matter admixture: a sample
#   that includes more adjacent white matter contributes many mature
#   oligodendrocytes and few OPCs (grey-matter OPC : OL ratios are far
#   higher than white-matter ones), so OPC share falls. The script therefore
#   also reports the lineage's share of all nonneuronal cells (a proxy for
#   white-matter content) and a partial correlation of rCMRGlc with OPC
#   share controlling for it. Interpret the primary correlation together
#   with that check.
#
# COPs. Committed oligodendrocyte precursors are a transitional state and
#   < 1 % of the lineage here; they are excluded from the primary ratio and
#   folded into the precursor pool in the sensitivity metric.
#
# Definitions (Linnarsson / Siletti et al. 2023 supercluster_term):
#   OPC   = "Oligodendrocyte precursor"
#   COP   = "Committed oligodendrocyte precursor"
#   Oligo = "Oligodendrocyte"
#
# PRIMARY metric   opc_frac_lineage      = OPC / (OPC + Oligo)      (COP excluded)
# SENSITIVITY      opc_cop_frac_lineage  = (OPC + COP) / (OPC + COP + Oligo)
# Also reported    opc_to_oligo_ratio    = OPC / Oligo  and its log2
# COMPARATOR       p_cells_OPC_nonneuronal = OPC / all nonneuronal cells
#                  (the s1b_3_nn definition, recomputed here so both
#                  metrics are compared on identical donor x region samples)
#
# Aggregation follows s1b_3_nn / s1b_6: compositions are computed per
# donor x anatomy_group, then averaged across donors (mean of donor-level
# proportions). Pooled-cell compositions are carried alongside for
# reference.
#
# Inputs
#   data_intermediate/linnarsson_adult_human_brain_obs_metadata_nonneuronal.rds
#   metadata/anatomy/s1b_linnarsson_roi_map.csv  (via attach_s1b_anatomy)
#   data_intermediate/heiss_2004_regions.csv      (via join_heiss_rates)
# Outputs
#   data_analysis/opc_oligo_lineage_composition_long_by_region.csv
#   data_analysis/opc_oligo_lineage_composition_by_region_with_rcmr.csv
#   data_analysis/opc_oligo_lineage_rcmr_correlations.csv
#   data_analysis/opc_oligo_lineage_rcmr_lm.csv
#   data_analysis/opc_oligo_lineage_admixture_check.csv
#   figs/s1b/p_opc_lineage_fraction_vs_rcmr.{pdf,jpg}
#   figs/s1b/p_opc_lineage_fraction_vs_rcmr_by_telencephalon.{pdf,jpg}
#   figs/s1b/p_opc_oligo_composition_signed_bar_rcmr_ordered.{pdf,jpg}
#   figs/s1b/p_opc_metric_comparison_vs_rcmr.{pdf,jpg}
###################################################################

source("helpers/plot_settings.R")

## Load packages
library(ggplot2)
library(tidyverse)
# NOTE: unlike the other s1b scripts this one does not use ggpmisc. The
# R^2 / p annotations are computed from the same lm() fit that geom_smooth
# draws and placed with geom_text (see lm_label_layer below), so the script
# runs with base ggplot2 + tidyverse only.

############################
## Load saved obs metadata
############################
obs <- readRDS("data_intermediate/linnarsson_adult_human_brain_obs_metadata_nonneuronal.rds")

######################################################
# Read rCMRGlc values from Heiss et al. 2004
######################################################

# Shared anatomy readers and the source-independent Heiss preparation.
source("helpers/read_heiss_rates.R")

#####################################
# Anatomical grouping of rois
#####################################

# The one authoritative s1b ROI crosswalk is
# metadata/anatomy/s1b_linnarsson_roi_map.csv.
obs <- attach_s1b_anatomy(obs)

###################################################################
# Telencephalon classification (same construction as s1b_3_nn / s1b_6;
# not re-written to disk here -- s1b_3_nn_telencephalon owns that file)
###################################################################
telencephalon_table <- obs %>%
  filter(!is.na(anatomy_id), anatomy_group != "", anatomy_group != "Unmapped") %>%
  count(anatomy_id, anatomy_group, ROIGroupCoarse, is_telencephalon, name = "n_cells") %>%
  group_by(anatomy_id, anatomy_group, is_telencephalon) %>%
  mutate(n_total = sum(n_cells), frac = n_cells / n_total) %>%
  slice_max(n_cells, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    anatomy_id, anatomy_group,
    dominant_ROIGroupCoarse = ROIGroupCoarse,
    n_cells_total           = n_total,
    is_telencephalon,
    division                = ifelse(is_telencephalon, "Telencephalon", "Non-telencephalon")
  ) %>%
  arrange(desc(is_telencephalon), anatomy_group)

###################################################################
# Lineage assignment
###################################################################
lineage_terms <- c(
  "Oligodendrocyte precursor"           = "OPC",
  "Committed oligodendrocyte precursor" = "COP",
  "Oligodendrocyte"                     = "Oligo"
)
missing_terms <- setdiff(names(lineage_terms), unique(obs$supercluster_term))
if (length(missing_terms) > 0) {
  stop("Expected supercluster_term(s) not found in nonneuronal obs: ",
       paste(missing_terms, collapse = ", "))
}

obs_mapped <- obs %>%
  filter(anatomy_group != "", anatomy_group != "Unmapped") %>%
  # as.character(): supercluster_term is a factor; indexing a named vector
  # with a factor would use the integer codes, not the labels.
  mutate(lineage_class = unname(lineage_terms[as.character(supercluster_term)]))

# Guard: every lineage cell must be classified, and the classes must be
# populated (a silent mis-assignment here would make every ratio ~1 or ~0).
lineage_check <- obs_mapped %>%
  filter(supercluster_term %in% names(lineage_terms)) %>%
  count(supercluster_term, lineage_class)
stopifnot(
  all(!is.na(lineage_check$lineage_class)),
  all(c("OPC", "COP", "Oligo") %in% lineage_check$lineage_class),
  identical(unname(lineage_terms[as.character(lineage_check$supercluster_term)]),
            lineage_check$lineage_class)
)

cat("\n================ Oligodendrocyte-lineage cells in mapped regions ================\n")
print(obs_mapped %>% count(supercluster_term, lineage_class) %>% arrange(desc(n)) %>% tibble::as_tibble(),
      n = Inf)

###################################################################
# Donor x region counts and compositions
###################################################################
donor_region <- obs_mapped %>% distinct(anatomy_group, donor_id)

donor_counts <- obs_mapped %>%
  group_by(anatomy_group, donor_id) %>%
  summarise(
    n_nonneuronal = n(),
    n_OPC   = sum(lineage_class == "OPC",   na.rm = TRUE),
    n_COP   = sum(lineage_class == "COP",   na.rm = TRUE),
    n_Oligo = sum(lineage_class == "Oligo", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  right_join(donor_region, by = c("anatomy_group", "donor_id")) %>%
  mutate(across(starts_with("n_"), ~ replace_na(.x, 0L))) %>%
  mutate(
    n_lineage        = n_OPC + n_Oligo,
    n_lineage_cop    = n_OPC + n_COP + n_Oligo,
    # PRIMARY: OPC share of OPC + mature oligodendrocytes
    opc_frac_lineage = if_else(n_lineage > 0, n_OPC / n_lineage, NA_real_),
    # SENSITIVITY: precursors (OPC + COP) share of the whole lineage
    opc_cop_frac_lineage = if_else(n_lineage_cop > 0, (n_OPC + n_COP) / n_lineage_cop, NA_real_),
    # Ratio forms
    opc_to_oligo_ratio   = if_else(n_Oligo > 0, n_OPC / n_Oligo, NA_real_),
    log2_opc_to_oligo    = log2(opc_to_oligo_ratio),
    # COMPARATOR: s1b_3_nn definition (OPC over all nonneuronal cells)
    p_cells_OPC_nonneuronal   = n_OPC / n_nonneuronal,
    p_cells_Oligo_nonneuronal = n_Oligo / n_nonneuronal,
    p_cells_lineage_nonneuronal = n_lineage / n_nonneuronal
  )

# Minimum lineage cells per donor x region sample for a composition to
# count. Samples below this are dropped from the donor mean (not set to 0).
min_lineage_cells <- 20
n_dropped <- sum(donor_counts$n_lineage < min_lineage_cells)
cat("\nDonor x region samples:", nrow(donor_counts),
    "| dropped for < ", min_lineage_cells, " lineage cells:", n_dropped, "\n")
if (n_dropped > 0) {
  print(donor_counts %>% filter(n_lineage < min_lineage_cells) %>%
          select(anatomy_group, donor_id, n_OPC, n_Oligo, n_lineage))
}

donor_counts <- donor_counts %>%
  mutate(across(c(opc_frac_lineage, opc_cop_frac_lineage, opc_to_oligo_ratio, log2_opc_to_oligo),
                ~ if_else(n_lineage >= min_lineage_cells, .x, NA_real_)))

metric_cols <- c("opc_frac_lineage", "opc_cop_frac_lineage",
                 "opc_to_oligo_ratio", "log2_opc_to_oligo",
                 "p_cells_OPC_nonneuronal", "p_cells_Oligo_nonneuronal",
                 "p_cells_lineage_nonneuronal")

# Region-level: mean (and sd) of donor-level compositions, plus pooled counts.
region_composition <- donor_counts %>%
  group_by(anatomy_group) %>%
  summarise(
    n_donors          = n_distinct(donor_id),
    n_donors_lineage  = sum(is.finite(opc_frac_lineage)),
    across(all_of(metric_cols), list(mean = ~ mean(.x, na.rm = TRUE),
                                     sd   = ~ sd(.x,   na.rm = TRUE)),
           .names = "{.col}_{.fn}"),
    pooled_n_OPC   = sum(n_OPC),
    pooled_n_COP   = sum(n_COP),
    pooled_n_Oligo = sum(n_Oligo),
    pooled_n_nonneuronal = sum(n_nonneuronal),
    .groups = "drop"
  ) %>%
  mutate(
    pooled_opc_frac_lineage     = pooled_n_OPC / (pooled_n_OPC + pooled_n_Oligo),
    pooled_opc_cop_frac_lineage = (pooled_n_OPC + pooled_n_COP) /
                                  (pooled_n_OPC + pooled_n_COP + pooled_n_Oligo),
    oligo_frac_lineage_mean     = 1 - opc_frac_lineage_mean
  )

# Long form (OPC vs Oligo shares; sum to 1 within region) for the bar figures.
composition_long <- region_composition %>%
  select(anatomy_group, n_donors_lineage,
         OPC = opc_frac_lineage_mean, Oligodendrocyte = oligo_frac_lineage_mean) %>%
  pivot_longer(c(OPC, Oligodendrocyte), names_to = "lineage_class", values_to = "comp_mean") %>%
  mutate(lineage_class = factor(lineage_class, levels = c("OPC", "Oligodendrocyte")))

###################################################################
# Bind to rCMRGlc and the telencephalon flag
###################################################################
analysis_df <- region_composition %>%
  join_heiss_rates() %>%
  left_join(
    telencephalon_table %>% select(anatomy_id, anatomy_group, is_telencephalon, division),
    by = c("anatomy_id", "anatomy_group")
  ) %>%
  arrange(desc(rcmr_value))

cat("\n================ Regions in analysis by division ================\n")
print(analysis_df %>% count(division, name = "n_regions"))

cat("\n================ OPC share of oligodendrocyte lineage by region ================\n")
print(
  analysis_df %>%
    transmute(anatomy_group, division, rcmr_value,
              n_donors_lineage,
              opc_frac_lineage = round(opc_frac_lineage_mean, 3),
              opc_cop_frac     = round(opc_cop_frac_lineage_mean, 3),
              opc_over_oligo   = round(opc_to_oligo_ratio_mean, 3),
              opc_over_allNN   = round(p_cells_OPC_nonneuronal_mean, 3),
              lineage_over_allNN = round(p_cells_lineage_nonneuronal_mean, 3)) %>%
    tibble::as_tibble(),
  n = Inf, width = Inf
)

write.csv(composition_long,
          "data_analysis/opc_oligo_lineage_composition_long_by_region.csv",
          row.names = FALSE)
write.csv(analysis_df,
          "data_analysis/opc_oligo_lineage_composition_by_region_with_rcmr.csv",
          row.names = FALSE)

###################################################################
# Correlations with rCMRGlc: all regions / telencephalon / non-telencephalon
###################################################################
metric_labels <- c(
  opc_frac_lineage_mean            = "OPC / (OPC + Oligo)  [primary]",
  opc_cop_frac_lineage_mean        = "(OPC + COP) / (OPC + COP + Oligo)",
  log2_opc_to_oligo_mean           = "log2(OPC / Oligo)",
  p_cells_OPC_nonneuronal_mean     = "OPC / all nonneuronal  [s1b_3_nn comparator]",
  p_cells_Oligo_nonneuronal_mean   = "Oligo / all nonneuronal",
  p_cells_lineage_nonneuronal_mean = "(OPC + Oligo) / all nonneuronal",
  pooled_opc_frac_lineage          = "OPC / (OPC + Oligo), pooled cells"
)

cor_one <- function(df, var, subset_label) {
  ok <- complete.cases(df$rcmr_value, df[[var]])
  n  <- sum(ok)
  if (n < 4) {
    return(data.frame(subset = subset_label, metric = var, metric_label = metric_labels[[var]],
                      n = n, spearman_rho = NA, spearman_p = NA,
                      pearson_r = NA, pearson_p = NA))
  }
  sp <- suppressWarnings(cor.test(df$rcmr_value[ok], df[[var]][ok], method = "spearman", exact = FALSE))
  pe <- suppressWarnings(cor.test(df$rcmr_value[ok], df[[var]][ok], method = "pearson"))
  data.frame(subset = subset_label, metric = var, metric_label = metric_labels[[var]],
             n = n,
             spearman_rho = unname(sp$estimate), spearman_p = sp$p.value,
             pearson_r    = unname(pe$estimate), pearson_p    = pe$p.value)
}

subsets <- list(
  "All regions"       = analysis_df,
  "Telencephalon"     = analysis_df %>% filter(is_telencephalon),
  "Non-telencephalon" = analysis_df %>% filter(!is_telencephalon)
)

cor_results <- bind_rows(lapply(names(subsets), function(s) {
  bind_rows(lapply(names(metric_labels), function(v) cor_one(subsets[[s]], v, s)))
})) %>%
  group_by(subset) %>%
  mutate(spearman_p_adj_BH = p.adjust(spearman_p, method = "BH")) %>%
  ungroup()

cat("\n================ Spearman / Pearson correlations with rCMRGlc ================\n")
print(cor_results %>% select(subset, metric_label, n, spearman_rho, spearman_p, pearson_r, pearson_p) %>%
        mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
        tibble::as_tibble(),   # tibble print: plain data.frame would partial-match n= to na.print
      n = Inf, width = Inf)

write.csv(cor_results, "data_analysis/opc_oligo_lineage_rcmr_correlations.csv", row.names = FALSE)

# Simple linear models, rCMRGlc ~ metric, per subset (primary + comparator).
lm_one <- function(df, var, subset_label) {
  df <- df[complete.cases(df$rcmr_value, df[[var]]), ]
  if (nrow(df) < 5) return(NULL)
  fit <- lm(reformulate(var, response = "rcmr_value"), data = df)
  s   <- summary(fit)
  data.frame(subset = subset_label, metric = var, metric_label = metric_labels[[var]],
             n = nrow(df),
             slope = coef(fit)[[2]], slope_se = s$coefficients[2, 2],
             slope_p = s$coefficients[2, 4],
             r_squared = s$r.squared, adj_r_squared = s$adj.r.squared)
}
lm_results <- bind_rows(lapply(names(subsets), function(s) {
  bind_rows(lapply(c("opc_frac_lineage_mean", "opc_cop_frac_lineage_mean",
                     "p_cells_OPC_nonneuronal_mean"),
                   function(v) lm_one(subsets[[s]], v, s)))
}))
cat("\n================ Linear models rCMRGlc ~ metric ================\n")
print(lm_results %>% mutate(across(where(is.numeric), ~ signif(.x, 3))) %>% tibble::as_tibble(),
      n = Inf, width = Inf)
write.csv(lm_results, "data_analysis/opc_oligo_lineage_rcmr_lm.csv", row.names = FALSE)

###################################################################
# White-matter admixture check
###################################################################
# The lineage's share of all nonneuronal cells, (OPC + Oligo) / all NN, is
# used as a proxy for how much white matter each dissection contains. If
# OPC share is driven by admixture it will be strongly negatively related
# to this proxy, and the rCMRGlc association should weaken once the proxy
# is partialled out. Partial Spearman = Pearson correlation of the
# rank-transformed variables after residualising both on the ranked proxy.
partial_spearman <- function(y, x, z) {
  ok <- complete.cases(y, x, z)
  ry <- rank(y[ok]); rx <- rank(x[ok]); rz <- rank(z[ok])
  ey <- resid(lm(ry ~ rz)); ex <- resid(lm(rx ~ rz))
  ct <- cor.test(ey, ex)
  c(rho = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}

admixture_check <- bind_rows(lapply(names(subsets), function(s) {
  d <- subsets[[s]]
  if (nrow(d) < 5) return(NULL)
  a <- suppressWarnings(cor.test(d$opc_frac_lineage_mean, d$p_cells_lineage_nonneuronal_mean,
                                 method = "spearman", exact = FALSE))
  b <- suppressWarnings(cor.test(d$rcmr_value, d$p_cells_lineage_nonneuronal_mean,
                                 method = "spearman", exact = FALSE))
  c0 <- suppressWarnings(cor.test(d$rcmr_value, d$opc_frac_lineage_mean,
                                  method = "spearman", exact = FALSE))
  pc <- partial_spearman(d$rcmr_value, d$opc_frac_lineage_mean, d$p_cells_lineage_nonneuronal_mean)
  data.frame(
    subset = s, n = nrow(d),
    rho_opcshare_vs_lineageshare = unname(a$estimate), p_opcshare_vs_lineageshare = a$p.value,
    rho_rcmr_vs_lineageshare     = unname(b$estimate), p_rcmr_vs_lineageshare     = b$p.value,
    rho_rcmr_vs_opcshare         = unname(c0$estimate), p_rcmr_vs_opcshare        = c0$p.value,
    partial_rho_rcmr_vs_opcshare_given_lineageshare = pc[["rho"]],
    partial_p_rcmr_vs_opcshare_given_lineageshare   = pc[["p"]]
  )
}))

cat("\n================ White-matter admixture check ================\n")
cat("Proxy for white-matter content: (OPC + Oligo) / all nonneuronal cells.\n")
print(admixture_check %>% mutate(across(where(is.numeric), ~ signif(.x, 3))) %>% tibble::as_tibble(),
      n = Inf, width = Inf)
write.csv(admixture_check, "data_analysis/opc_oligo_lineage_admixture_check.csv", row.names = FALSE)

###################################################################
# Region palette (project palette, unused levels dropped)
###################################################################
region_scale <- tryCatch({
  check_region_palette(analysis_df, region_col = "anatomy_group")
  present <- intersect(region_order, unique(as.character(analysis_df$anatomy_group)))
  analysis_df$anatomy_group <- factor(analysis_df$anatomy_group, levels = present)
  ggplot2::scale_color_manual(values = region_palette, drop = TRUE)
}, error = function(e) {
  message("Region palette unavailable (", e$message, "); using default ggplot palette.")
  NULL
})

division_colors <- c("Telencephalon" = "#D7263D", "Non-telencephalon" = "#1B998B")
rcmr_lab <- "rCMRGlc (\u00b5mol/100 g/min.)"

# ---- R^2 / p annotation without ggpmisc -----------------------------------
# lm_label_df(): one row per group with an R^2 / p / n label from lm(y ~ x)
# (the same fit geom_smooth(method = "lm") draws). lm_label_layer(): a
# geom_text layer that pins those labels to the top-right of each panel.
# Labels are plotmath expressions so R^2 renders with a superscript and
# P in italics, matching the ggpmisc style used in the other s1b figures.
fmt_p <- function(p) ifelse(p < 0.001, "< 0.001", sprintf("= %.3f", p))
lm_label_df <- function(df, x, y, by = NULL) {
  df %>%
    group_by(across(all_of(by))) %>%
    summarise(
      label = {
        ok  <- complete.cases(.data[[x]], .data[[y]])
        fit <- lm(.data[[y]][ok] ~ .data[[x]][ok])
        s   <- summary(fit)
        sprintf("italic(R)^2~'='~%.2f*','~italic(P)~'%s'*','~italic(n)~'='~%d",
                s$r.squared, fmt_p(s$coefficients[2, 4]), sum(ok))
      },
      .groups = "drop"
    )
}
lm_label_layer <- function(label_df, size = 4, color = "black") {
  geom_text(
    data = label_df,
    aes(x = Inf, y = Inf, label = label),
    hjust = 1.05, vjust = 1.6, size = size, color = color,
    parse = TRUE, inherit.aes = FALSE, show.legend = FALSE
  )
}

###################################################################
# Figure A: rCMRGlc vs OPC share of the lineage, all regions
###################################################################
p_frac_vs_rcmr <- ggplot(analysis_df,
                         aes(x = opc_frac_lineage_mean, y = rcmr_value, color = anatomy_group)) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black") +
  geom_point(size = 3, alpha = 0.9) +
  lm_label_layer(lm_label_df(analysis_df, "opc_frac_lineage_mean", "rcmr_value"), size = 4) +
  labs(
    title = "rCMRGlc versus OPC share of the oligodendrocyte lineage",
    subtitle = paste0("OPC / (OPC + oligodendrocytes), mean of donor-level compositions; n = ",
                      nrow(analysis_df), " regions"),
    x = "OPC fraction among OPC + oligodendrocytes",
    y = rcmr_lab, color = "Region"
  ) +
  theme_classic(base_size = 13)
if (!is.null(region_scale)) p_frac_vs_rcmr <- p_frac_vs_rcmr + region_scale

ggsave("figs/s1b/p_opc_lineage_fraction_vs_rcmr.pdf", p_frac_vs_rcmr,
       width = 9, height = 6, units = "in")
ggsave("figs/s1b/p_opc_lineage_fraction_vs_rcmr.jpg", p_frac_vs_rcmr,
       width = 9, height = 6, units = "in", dpi = 300)

###################################################################
# Figure B: telencephalon split
###################################################################
p_frac_vs_rcmr_division <- ggplot(analysis_df,
                                  aes(x = opc_frac_lineage_mean, y = rcmr_value,
                                      color = division, fill = division)) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15) +
  geom_point(size = 3, alpha = 0.9) +
  lm_label_layer(lm_label_df(analysis_df, "opc_frac_lineage_mean", "rcmr_value", by = "division"),
                 size = 3.5) +
  facet_wrap(~ division) +
  scale_color_manual(values = division_colors) +
  scale_fill_manual(values = division_colors) +
  labs(
    title = "rCMRGlc versus OPC share of the oligodendrocyte lineage by telencephalon status",
    subtitle = "OPC / (OPC + oligodendrocytes); separate fits per division",
    x = "OPC fraction among OPC + oligodendrocytes",
    y = rcmr_lab, color = "Division", fill = "Division"
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top")

ggsave("figs/s1b/p_opc_lineage_fraction_vs_rcmr_by_telencephalon.pdf", p_frac_vs_rcmr_division,
       width = 10, height = 5.8, units = "in")
ggsave("figs/s1b/p_opc_lineage_fraction_vs_rcmr_by_telencephalon.jpg", p_frac_vs_rcmr_division,
       width = 10, height = 5.8, units = "in", dpi = 300)

###################################################################
# Figure C: signed paired bars (OPC up, oligodendrocyte down), regions
# ordered by rCMRGlc, bars filled by rCMRGlc (matches s1b_6 Figure A)
###################################################################
composition_bar_df <- composition_long %>%
  inner_join(analysis_df %>% select(anatomy_group, rcmr_value, division), by = "anatomy_group") %>%
  mutate(
    signed_comp   = if_else(lineage_class == "OPC", comp_mean, -comp_mean),
    anatomy_group = forcats::fct_reorder(as.character(anatomy_group), rcmr_value)
  )

p_signed_bar <- ggplot(composition_bar_df,
                       aes(x = anatomy_group, y = signed_comp, fill = rcmr_value)) +
  geom_col(width = 0.85, color = NA) +
  geom_hline(yintercept = 0, linewidth = 0.35) +
  facet_grid(lineage_class ~ ., scales = "free_y") +
  scale_y_continuous(labels = function(x) abs(x), breaks = c(-1, -0.5, 0, 0.5, 1)) +
  scale_fill_viridis_c(option = "magma", direction = -1) +
  labs(
    title = "Relative composition of OPCs and oligodendrocytes by region",
    subtitle = "Within each region OPC + oligodendrocyte fractions sum to 1; regions ordered by rCMRGlc",
    x = NULL, y = "Fraction of OPC + oligodendrocytes",
    fill = "rCMRGlc\n(\u00b5mol/100 g/min.)"
  ) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.background = element_blank(),
        strip.text.y = element_text(angle = 0))

ggsave("figs/s1b/p_opc_oligo_composition_signed_bar_rcmr_ordered.pdf", p_signed_bar,
       width = 11, height = 5.5, units = "in")
ggsave("figs/s1b/p_opc_oligo_composition_signed_bar_rcmr_ordered.jpg", p_signed_bar,
       width = 11, height = 5.5, units = "in", dpi = 300)

###################################################################
# Figure D: metric comparison -- the s1b_3_nn definition (OPC / all
# nonneuronal) next to the lineage definition, same regions, same rCMRGlc
###################################################################
comparison_levels <- c(
  "OPC / all nonneuronal cells (s1b_3_nn)",
  "OPC / (OPC + oligodendrocytes)",
  "(OPC + COP) / (OPC + COP + oligodendrocytes)"
)
comparison_df <- analysis_df %>%
  select(anatomy_group, division, rcmr_value,
         `OPC / all nonneuronal cells (s1b_3_nn)`     = p_cells_OPC_nonneuronal_mean,
         `OPC / (OPC + oligodendrocytes)`              = opc_frac_lineage_mean,
         `(OPC + COP) / (OPC + COP + oligodendrocytes)` = opc_cop_frac_lineage_mean) %>%
  pivot_longer(all_of(comparison_levels), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = comparison_levels))

p_metric_comparison <- ggplot(comparison_df, aes(x = value, y = rcmr_value, color = anatomy_group)) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "steelblue") +
  geom_point(size = 2.6, alpha = 0.9) +
  lm_label_layer(lm_label_df(comparison_df, "value", "rcmr_value", by = "metric"), size = 3.2) +
  facet_wrap(~ metric, scales = "free_x", ncol = 3) +
  labs(
    title = "How the OPC denominator changes the association with rCMRGlc",
    subtitle = "Same regions and donor x region samples; only the denominator differs",
    x = "OPC proportion (definition per panel)", y = rcmr_lab, color = "Region"
  ) +
  theme_facet_compact(12)
if (!is.null(region_scale)) p_metric_comparison <- p_metric_comparison + region_scale

ggsave("figs/s1b/p_opc_metric_comparison_vs_rcmr.pdf", p_metric_comparison,
       width = 13, height = 5.2, units = "in")
ggsave("figs/s1b/p_opc_metric_comparison_vs_rcmr.jpg", p_metric_comparison,
       width = 13, height = 5.2, units = "in", dpi = 300)

cat("\nDone. Outputs:\n")
cat("  data_analysis/opc_oligo_lineage_composition_long_by_region.csv\n")
cat("  data_analysis/opc_oligo_lineage_composition_by_region_with_rcmr.csv\n")
cat("  data_analysis/opc_oligo_lineage_rcmr_correlations.csv\n")
cat("  data_analysis/opc_oligo_lineage_rcmr_lm.csv\n")
cat("  data_analysis/opc_oligo_lineage_admixture_check.csv\n")
cat("  figs/s1b/p_opc_lineage_fraction_vs_rcmr.{pdf,jpg}\n")
cat("  figs/s1b/p_opc_lineage_fraction_vs_rcmr_by_telencephalon.{pdf,jpg}\n")
cat("  figs/s1b/p_opc_oligo_composition_signed_bar_rcmr_ordered.{pdf,jpg}\n")
cat("  figs/s1b/p_opc_metric_comparison_vs_rcmr.{pdf,jpg}\n")
