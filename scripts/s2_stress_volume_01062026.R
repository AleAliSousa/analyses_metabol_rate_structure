setwd("~/Library/CloudStorage/Dropbox/COLLABORATIVE/Do expensive brain regions increase less in humans/analyses_metabol_rate_structure")

## Install and Load up packages
library(tidyverse)  # loads ggplot2, dplyr, tidyr, etc.
library(ggpmisc)    # for stat_poly_eq

## ================================
## Volume Change vs rCMRGlc
## ================================

## ----------------
## Load data
## ----------------

data <- read.csv("data_raw/stress_volume.csv", header = TRUE)

# Expected columns (in order): Region, rCMRGlc,
# norm-vs-(all|english|romanian) %, then mean/sd/lci/uci for english, romanian, normative.
names(data) <- c(
  "anatomy_group",
  "rcmr_value",
  "norm_all",
  "norm_english",
  "norm_romanian",
  "english_mean",
  "english_sd",
  "english_lci",
  "english_uci",
  "romanian_mean",
  "romanian_sd",
  "romanian_lci",
  "romanian_uci",
  "norm_mean",
  "norm_sd",
  "norm_lci",
  "norm_uci"
)

# Tidy region labels for display
data$anatomy_group <- gsub("Accumbens", "Nucleus accumbens", data$anatomy_group)
data$anatomy_group <- gsub("_", " ", data$anatomy_group)

## ----------------
## Reshape for plotting
## ----------------

plot_df <- data %>%
  pivot_longer(
    cols = c(norm_all, norm_english, norm_romanian),
    names_to = "group",
    values_to = "volume_change"
  ) %>%
  mutate(
    group = recode(
      group,
      norm_all      = "Normative vs All adoptees",
      norm_english  = "Normative vs English adoptees",
      norm_romanian = "Normative vs Romanian adoptees"
    )
  )

## ----------------
## Color palette
## ----------------

regions <- sort(unique(plot_df$anatomy_group))
base_cols <- c(
  "#000000", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
)

pal <- grDevices::colorRampPalette(base_cols)(length(regions))
names(pal) <- regions

## ----------------
## Plot
## ----------------

# Keep only rows with data for this plot and remove Total Brain
plot_df_clean <- plot_df %>%
  filter(
    anatomy_group != "Total Brain",
    !is.na(rcmr_value),
    !is.na(volume_change)
  ) %>%
  droplevels()

p1 <- ggplot(plot_df_clean, aes(x = rcmr_value, y = volume_change, color = anatomy_group)) +
  geom_point(size = 2.8, alpha = 0.85) +
  geom_smooth(
    aes(group = 1),
    method = "lm",
    se = TRUE,
    color = "steelblue"
  ) +
  stat_poly_eq(
    aes(label = after_stat(paste(rr.label, p.value.label, sep = "*\", \"*"))),
    formula = y ~ x,
    parse = TRUE,
    label.x = "right",
    label.y = "top",
    size = 4,
    color = "black"
  ) +
  facet_wrap(~ group, scales = "free_y") +
  scale_color_manual(values = pal) +
  labs(
    x = "rCMRGlc (µmol / 100 g / min)",
    y = "% Volume change (Normative - Adoptee Group)",
    color = "Region"
  ) +
  theme_classic(base_size = 14)

## ================================
## Box plots: Volumes per region
## Box = mean ± SD, whiskers = 95% CI, midline = mean
## (geom_boxplot with stat = "identity" because only summary
##  statistics are available — no raw subject-level data)
## ================================

box_df <- data %>%
  select(
    anatomy_group,
    english_mean, english_sd, english_lci, english_uci,
    romanian_mean, romanian_sd, romanian_lci, romanian_uci,
    norm_mean, norm_sd, norm_lci, norm_uci
  ) %>%
  pivot_longer(
    cols = -anatomy_group,
    names_to = c("group", ".value"),
    names_sep = "_"
  ) %>%
  mutate(
    group = factor(
      recode(group,
        english  = "English adoptees",
        romanian = "Romanian adoptees",
        norm     = "Normative"
      ),
      levels = c("English adoptees", "Romanian adoptees", "Normative")
    )
  )

# Box = mean ± SD; whiskers = 95% CI; midline = mean.
p2 <- ggplot(box_df, aes(x = group, fill = group)) +
  geom_boxplot(
    aes(
      ymin = lci,
      lower = mean - sd,
      middle = mean,
      upper = mean + sd,
      ymax = uci
    ),
    stat = "identity",
    width = 0.55,
    color = "black",
    linewidth = 0.4,
    alpha = 0.85
  ) +
  facet_wrap(~ anatomy_group, scales = "free_y") +
  labs(
    x = NULL,
    y = "Volume per hemisphere (cm³)",
    fill = "Group",
    caption = "Box = mean \u00b1 SD\u00a0\u00a0\u2502\u00a0\u00a0Whiskers = 95% CI\u00a0\u00a0\u2502\u00a0\u00a0Midline = mean"
  ) +
  scale_fill_manual(
    values = c(
      "English adoptees"  = "#56B4E9",
      "Romanian adoptees" = "#D55E00",
      "Normative"         = "#999999"
    )
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(face = "bold"),
    strip.background = element_blank(),
    panel.spacing = unit(0.8, "lines"),
    plot.caption = element_text(size = 10, hjust = 0.5, margin = margin(t = 8))
  )

# Print plot
p1

# Save as PDF and JPG
ggsave(
  filename = "figs/s2/rCMRGlc_volume_change.pdf",
  plot = p1,
  width = 12,
  height = 6.5,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = "figs/s2/rCMRGlc_volume_change.jpg",
  plot = p1,
  width = 12,
  height = 6.5,
  units = "in",
  dpi = 300,
  device = "jpeg",
  bg = "white"
)

# Print plot
p2

# Save as PDF and JPG
ggsave(
  filename = "figs/s2/mean_brain_region_volumes_boxplot.pdf",
  plot = p2,
  width = 10,
  height = 7,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = "figs/s2/mean_brain_region_volumes_boxplot.jpg",
  plot = p2,
  width = 10,
  height = 7,
  units = "in",
  dpi = 300,
  device = "jpeg",
  bg = "white"
)
