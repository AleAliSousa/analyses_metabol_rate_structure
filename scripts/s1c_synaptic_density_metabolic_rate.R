setwd(local({ d <- normalizePath(getwd()); while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d); d }))  # repo root

## ============================================================
## Study 6 — Synaptic Density (SV2A) vs. rCMRGlc
##
## Johansen et al. 2024 Table 2: in vivo SV2A Bmax (pmol/mL)
##   [11C]UCB-J PET atlas, N = 33 healthy adults.
##   Source: J Neurosci 44(33):e1750232024.
##   DOI: 10.1523/JNEUROSCI.1750-23.2024
##
## Heiss et al. 2004 Table 1: regional rCMRGlc (µmol/100g/min)
##   18F-FDG HRRT PET, N = 9 healthy adults.
##   Source: J Nucl Med 45(11):1811-1815.
##   Local data: data_raw/Heiss_etal_2004_TABLE1.csv
##
## Overlap regions are joined by name; the mapping is documented
## explicitly in region_map below. Only lobe-level / homologous
## regions are paired to avoid comparing mismatched anatomy.
## ============================================================

library(tidyverse)
library(ggpmisc)

source("R/plot_settings.R")

if (!dir.exists("figs/s1c")) dir.create("figs/s1c", recursive = TRUE)

## ----------------------------------------------------------------
## 1. Load data
## ----------------------------------------------------------------

johansen <- read.csv("data_raw/Johansen_2024_Table2.csv", header = TRUE)

heiss_raw <- read.csv("data_raw/Heiss_etal_2004_TABLE1.csv", header = TRUE)
names(heiss_raw) <- c(
  "category", "heiss_region",
  "rcmrglc_mean", "rcmrglc_sd",
  "lr_diff_mean", "lr_diff_sd", "lr_pval",
  "heiss1991_mean"
)

## ----------------------------------------------------------------
## 2. Region mapping: Johansen (Desikan-Killiany + FreeSurfer
##    subcortical) -> canonical Heiss region name
##
## Rules:
##  - Use lobe totals for cortical lobes (Johansen provides them).
##  - Use directly matched subcortical labels where the anatomy
##    aligns well. Exclude regions with no reasonable Heiss match
##    (e.g. Pons, Ventral dc, cingulate sub-parcels, individual
##    gyri that Heiss aggregates into lobes).
## ----------------------------------------------------------------

region_map <- tribble(
  ~johansen_region,          ~heiss_region,
  "Frontal (total)",         "Frontal lobe",
  "Parietal (total)",        "Parietal lobe",
  "Temporal (total)",        "Temporal lobe",
  "Occipital (total)",       "Occipital lobe",
  "Insula",                  "Insular lobe",
  "Hippocampus",             "Hippocampus",
  "Amygdala",                "Corpus amygdaloideum",
  "Caudate",                 "Caudatum",
  "Putamen",                 "Putamen",
  "Accumbens",               "Nucleus accumbens",
  "Pallidum",                "Pallidum",
  # Thalamus excluded: Johansen measures the whole thalamus (FreeSurfer ROI)
  # but Heiss only sampled Nucleus medial thalami + geniculate bodies — no
  # composite spanning the full thalamus exists, consistent with s3 and s4a.
  "Cerebellum",              "Cerebellar cortex",
  "White matter",            "Centrum semiovale"
)

## Attach Johansen SV2A data
johansen_matched <- johansen %>%
  inner_join(region_map, by = c("region" = "johansen_region")) %>%
  select(heiss_region, lobe,
         sv2a_mean = sv2a_bmax_total_mean,
         sv2a_sd   = sv2a_bmax_total_sd)

## Attach Heiss rCMRGlc data
heiss_matched <- heiss_raw %>%
  select(heiss_region, rcmrglc_mean, rcmrglc_sd)

## Join on heiss_region
plot_df <- johansen_matched %>%
  inner_join(heiss_matched, by = "heiss_region") %>%
  mutate(
    anatomy_group = canonical_region(heiss_region)
  )

## Sanity check: stop if any region is missing a palette colour.
check_region_palette(plot_df, region_col = "anatomy_group")

## ----------------------------------------------------------------
## 3. Scatter: SV2A Bmax vs. rCMRGlc
## ----------------------------------------------------------------

p1 <- ggplot(plot_df,
             aes(x = rcmrglc_mean, y = sv2a_mean, color = anatomy_group)) +
  geom_errorbar(
    aes(ymin = sv2a_mean - sv2a_sd, ymax = sv2a_mean + sv2a_sd),
    width = 0, alpha = 0.4, linewidth = 0.5
  ) +
  geom_errorbarh(
    aes(xmin = rcmrglc_mean - rcmrglc_sd, xmax = rcmrglc_mean + rcmrglc_sd),
    height = 0, alpha = 0.4, linewidth = 0.5
  ) +
  geom_point(size = 3, alpha = 0.9) +
  geom_smooth(
    aes(group = 1),
    method  = "lm",
    se      = TRUE,
    color   = "steelblue",
    linewidth = 0.8
  ) +
  stat_poly_eq(
    aes(label = after_stat(paste(rr.label, p.value.label, sep = "*\", \"*")),
        group = 1),
    formula    = y ~ x,
    parse      = TRUE,
    label.x    = "right",
    label.y    = "top",
    size       = 4,
    color      = "black"
  ) +
  scale_color_manual(values = region_palette, drop = TRUE) +
  labs(
    title    = "Synaptic density vs. glucose metabolic rate across brain regions",
    subtitle = "Johansen et al. 2024 (SV2A PET) × Heiss et al. 2004 (FDG PET)",
    x        = "rCMRGlc (µmol / 100 g / min)",
    y        = "SV2A Bₘₐˣ (pmol / mL)",
    color    = "Region",
    caption  = "Points = region means; error bars = SD."
  ) +
  theme_project(base_size = 14) +
  theme(
    plot.subtitle = element_text(size = 11, color = "grey40"),
    plot.caption  = element_text(size = 9, hjust = 0, color = "grey40")
  )

p1

ggsave(
  filename = "figs/s1c/sv2a_vs_rcmrglc.pdf",
  plot     = p1,
  width    = 9,
  height   = 6.5,
  units    = "in",
  device   = cairo_pdf
)

ggsave(
  filename = "figs/s1c/sv2a_vs_rcmrglc.jpg",
  plot     = p1,
  width    = 9,
  height   = 6.5,
  units    = "in",
  dpi      = 300,
  device   = "jpeg",
  bg       = "white"
)

## ----------------------------------------------------------------
## 4. Scatter: cortex-only vs. subcortical highlighted
## ----------------------------------------------------------------

plot_df2 <- plot_df %>%
  mutate(
    compartment = case_when(
      lobe == "Subcortical"   ~ "Subcortical",
      lobe == "White matter"  ~ "White matter",
      TRUE                    ~ "Cortical"
    )
  )

p2 <- ggplot(plot_df2,
             aes(x = rcmrglc_mean, y = sv2a_mean, color = anatomy_group,
                 shape = compartment)) +
  geom_errorbar(
    aes(ymin = sv2a_mean - sv2a_sd, ymax = sv2a_mean + sv2a_sd),
    width = 0, alpha = 0.4, linewidth = 0.5
  ) +
  geom_errorbarh(
    aes(xmin = rcmrglc_mean - rcmrglc_sd, xmax = rcmrglc_mean + rcmrglc_sd),
    height = 0, alpha = 0.4, linewidth = 0.5
  ) +
  geom_point(size = 3.2, alpha = 0.9) +
  ggrepel::geom_text_repel(
    aes(label = anatomy_group),
    size = 3, max.overlaps = 20, color = "grey30",
    box.padding = 0.35, show.legend = FALSE
  ) +
  geom_smooth(
    aes(group = 1),
    method    = "lm",
    se        = TRUE,
    color     = "steelblue",
    linewidth = 0.8
  ) +
  stat_poly_eq(
    aes(label = after_stat(paste(rr.label, p.value.label, sep = "*\", \"*")),
        group = 1),
    formula    = y ~ x,
    parse      = TRUE,
    label.x    = "right",
    label.y    = "top",
    size       = 4,
    color      = "black"
  ) +
  scale_color_manual(values = region_palette, drop = TRUE) +
  scale_shape_manual(
    values = c("Cortical" = 16, "Subcortical" = 17, "White matter" = 15)
  ) +
  labs(
    title    = "Synaptic density vs. glucose metabolic rate — labelled",
    subtitle = "Johansen et al. 2024 (SV2A PET) × Heiss et al. 2004 (FDG PET)",
    x        = "rCMRGlc (µmol / 100 g / min)",
    y        = "SV2A Bₘₐˣ (pmol / mL)",
    color    = "Region",
    shape    = "Compartment"
  ) +
  theme_project(base_size = 13)

p2

ggsave(
  filename = "figs/s1c/sv2a_vs_rcmrglc_labelled.pdf",
  plot     = p2,
  width    = 10,
  height   = 7,
  units    = "in",
  device   = cairo_pdf
)

ggsave(
  filename = "figs/s1c/sv2a_vs_rcmrglc_labelled.jpg",
  plot     = p2,
  width    = 10,
  height   = 7,
  units    = "in",
  dpi      = 300,
  device   = "jpeg",
  bg       = "white"
)

## ----------------------------------------------------------------
## 5. Save merged table for inspection
## ----------------------------------------------------------------

plot_df %>%
  select(
    anatomy_group, lobe,
    rcmrglc_mean, rcmrglc_sd,
    sv2a_mean, sv2a_sd
  ) %>%
  arrange(rcmrglc_mean) %>%
  write.csv("tables/s1c/sv2a_rcmrglc_matched_regions.csv", row.names = FALSE)

message("Study 1c done. Figures in figs/s1c/, merged table in tables/.")
