setwd(local({ d <- normalizePath(getwd()); while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d); d }))  # repo root

## ============================================================
## Study s1c — Synaptic density vs. rCMRGlc
##
## Johansen et al. 2024 Table 2: in vivo 3D atlas of synaptic density,
##   measured using the synaptic marker Synaptic Vesicle glycoprotein 2A
##   (SV2A; [11C]UCB-J PET Bmax in pmol/mL), N = 33 healthy adults.
##   Source: J Neurosci 44(33):e1750232024.
##   DOI: 10.1523/JNEUROSCI.1750-23.2024
##   Local source copy: data_raw/Johansen_etal_2024_Table2.csv
##
## Heiss et al. 2004 Table 1: regional rCMRGlc (µmol/100g/min)
##   18F-FDG HRRT PET, N = 9 healthy adults.
##   Source: J Nucl Med 45(11):1811-1815.
##   Authoritative source: data_raw/Heiss_etal_2004_TABLE1.csv
##   Analysis input: data_intermediate/heiss_2004_regions.csv
##
## Overlap regions are joined by anatomy_id; the mapping is documented
## in metadata/anatomy/s1c_johansen_region_map.csv. Only lobe-level / homologous
## regions are paired to avoid comparing mismatched anatomy.
## ============================================================

library(tidyverse)
library(ggpmisc)

source("helpers/plot_settings.R")
source("helpers/read_heiss_rates.R")

if (!dir.exists("figs/s1c")) dir.create("figs/s1c", recursive = TRUE)

## ----------------------------------------------------------------
## 1. Load data
## ----------------------------------------------------------------

johansen_source <- read.csv(
  "data_raw/Johansen_etal_2024_Table2.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Preserve the paper's uppercase SV2A column names in the raw source. At this
# import boundary, give the measurements descriptive analysis names so "s" is
# reserved for project study identifiers such as s1c.
required_johansen_columns <- c(
  "Lobe", "Region", "SV2A_total.pmol_mL", "SV2A_total_SD"
)
missing_johansen_columns <- setdiff(
  required_johansen_columns, names(johansen_source)
)
if (length(missing_johansen_columns)) {
  stop(
    "Johansen source is missing required column(s): ",
    paste(missing_johansen_columns, collapse = ", "),
    call. = FALSE
  )
}

synaptic_density_atlas <- johansen_source %>%
  transmute(
    source_lobe = Lobe,
    region = Region,
    synaptic_density_mean = as.numeric(`SV2A_total.pmol_mL`),
    synaptic_density_sd = as.numeric(SV2A_total_SD)
  )

## ----------------------------------------------------------------
## 2. Region mapping: Johansen (Desikan-Killiany + FreeSurfer
##    subcortical) -> stable project anatomy_id -> Heiss rate
##
## Rules:
##  - Use lobe totals for cortical lobes (Johansen provides them).
##  - Use directly matched subcortical labels where the anatomy
##    aligns well. Exclude regions with no reasonable Heiss match
##    (e.g. Pons, Ventral dc, cingulate sub-parcels, individual
##    gyri that Heiss aggregates into lobes).
## ----------------------------------------------------------------

region_map <- derive_study_heiss_rates(
  "metadata/anatomy/s1c_johansen_region_map.csv"
)

## Attach Johansen synaptic-density data
synaptic_density_matched <- synaptic_density_atlas %>%
  inner_join(region_map, by = c("region" = "source_region")) %>%
  select(anatomy_id, anatomy_group, source_lobe,
         synaptic_density_mean, synaptic_density_sd,
         rcmrglc_mean = rcmr_value,
         rcmrglc_sd = rcmr_sd)

plot_df <- synaptic_density_matched %>%
  mutate(
    anatomical_scope = case_when(
      anatomy_id == "centrum_semiovale" ~ "White matter",
      anatomy_id == "cerebellar_cortex" ~ "Cerebellar cortex",
      source_lobe == "Subcortical" ~ "Subcortical grey matter",
      TRUE ~ "Cerebral cortex"
    )
  )

## Sanity check: stop if any region is missing a palette colour.
check_region_palette(plot_df, region_col = "anatomy_group")

## ----------------------------------------------------------------
## 3. Scatter: synaptic density vs. rCMRGlc
## ----------------------------------------------------------------

p1 <- ggplot(plot_df,
             aes(x = rcmrglc_mean, y = synaptic_density_mean,
                 color = anatomy_group)) +
  geom_errorbar(
    aes(
      ymin = synaptic_density_mean - synaptic_density_sd,
      ymax = synaptic_density_mean + synaptic_density_sd
    ),
    width = 0, alpha = 0.4, linewidth = 0.5
  ) +
  geom_errorbar(
    aes(xmin = rcmrglc_mean - rcmrglc_sd, xmax = rcmrglc_mean + rcmrglc_sd),
    orientation = "y", width = 0, alpha = 0.4, linewidth = 0.5
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
    subtitle = "Johansen et al. 2024 (synaptic-density atlas) × Heiss et al. 2004 (FDG PET)",
    x        = "rCMRGlc (µmol / 100 g / min)",
    y        = "Synaptic density (SV2A Bmax, pmol / mL)",
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
  filename = "figs/s1c/synaptic_density_vs_rcmrglc.pdf",
  plot     = p1,
  width    = 9,
  height   = 6.5,
  units    = "in",
  device   = cairo_pdf
)

ggsave(
  filename = "figs/s1c/synaptic_density_vs_rcmrglc.jpg",
  plot     = p1,
  width    = 9,
  height   = 6.5,
  units    = "in",
  dpi      = 300,
  device   = "jpeg",
  bg       = "white"
)

## ----------------------------------------------------------------
## 4. Scatter with explicit cerebral, cerebellar, subcortical, and white-matter scopes
## ----------------------------------------------------------------

p2 <- ggplot(plot_df,
             aes(x = rcmrglc_mean, y = synaptic_density_mean,
                 color = anatomy_group,
                 shape = anatomical_scope)) +
  geom_errorbar(
    aes(
      ymin = synaptic_density_mean - synaptic_density_sd,
      ymax = synaptic_density_mean + synaptic_density_sd
    ),
    width = 0, alpha = 0.4, linewidth = 0.5
  ) +
  geom_errorbar(
    aes(xmin = rcmrglc_mean - rcmrglc_sd, xmax = rcmrglc_mean + rcmrglc_sd),
    orientation = "y", width = 0, alpha = 0.4, linewidth = 0.5
  ) +
  geom_point(size = 3.2, alpha = 0.9) +
  ggrepel::geom_text_repel(
    aes(label = anatomy_group),
    size = 3, max.overlaps = 20, color = "grey30",
    box.padding = 0.35, show.legend = FALSE
  ) +
  geom_smooth(
    data = plot_df,
    aes(
      x = rcmrglc_mean,
      y = synaptic_density_mean,
      group = 1
    ),
    inherit.aes = FALSE,
    method    = "lm",
    se        = TRUE,
    color     = "steelblue",
    linewidth = 0.8
  ) +
  stat_poly_eq(
    data = plot_df,
    aes(
      x = rcmrglc_mean,
      y = synaptic_density_mean,
      label = after_stat(paste(rr.label, p.value.label, sep = "*\", \"*")),
      group = 1
    ),
    inherit.aes = FALSE,
    formula    = y ~ x,
    parse      = TRUE,
    label.x    = "right",
    label.y    = "top",
    size       = 4,
    color      = "black"
  ) +
  scale_color_manual(values = region_palette, drop = TRUE) +
  scale_shape_manual(
    values = c(
      "Cerebral cortex" = 16,
      "Cerebellar cortex" = 18,
      "Subcortical grey matter" = 17,
      "White matter" = 15
    )
  ) +
  labs(
    title    = "Synaptic density vs. glucose metabolic rate: labelled",
    subtitle = "Johansen et al. 2024 (synaptic-density atlas) × Heiss et al. 2004 (FDG PET)",
    x        = "rCMRGlc (µmol / 100 g / min)",
    y        = "Synaptic density (SV2A Bmax, pmol / mL)",
    color    = "Region",
    shape    = "Compartment"
  ) +
  theme_project(base_size = 13)

p2

ggsave(
  filename = "figs/s1c/synaptic_density_vs_rcmrglc_labelled.pdf",
  plot     = p2,
  width    = 10,
  height   = 7,
  units    = "in",
  device   = cairo_pdf
)

ggsave(
  filename = "figs/s1c/synaptic_density_vs_rcmrglc_labelled.jpg",
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
    anatomy_id, anatomy_group, source_lobe, anatomical_scope,
    rcmrglc_mean, rcmrglc_sd,
    synaptic_density_mean, synaptic_density_sd
  ) %>%
  arrange(rcmrglc_mean) %>%
  write.csv(
    "tables/s1c/synaptic_density_rcmrglc_matched_regions.csv",
    row.names = FALSE
  )

message("Study s1c done. Figures in figs/s1c/, merged table in tables/s1c/.")
