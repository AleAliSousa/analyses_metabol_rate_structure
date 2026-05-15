setwd("~/Library/CloudStorage/Dropbox/COLLABORATIVE/Do expensive brain regions increase less in humans/analyses_metabol_rate_structure")

## Install and Load up packages
library(here)
library(ggplot2)
library(ggpmisc)  # for stat_poly_eq
library(tidyverse)

# Paths relative to project root:
r_dir      <- here("R")
data_dir   <- here("data")
output_dir <- here("output")

# Brain Structure Metabolic Rates (rCMRGlc) compared to Anatomical Variables

## Preparation

### Load data

######## To do - make a cleaner sheet which can link to references. Calculate ratios and densities in the sheet, so that we can just import the relevant columns here.

RawData <- read.csv(here("data", "stereology.csv"), header = TRUE)
#RawData <- read.csv("~/Library/CloudStorage/Dropbox/COLLABORATIVE/Do expensive brain regions increase less in humans/analyses_metabol_rate_structure/data/stereology.csv", header = TRUE)
# Select relevant columns and rename them for clarity
data <- RawData[, c(
  "Region",
  "rCMRGlc",
  "Volume",
  "VolumeSD",
  "Neuron_N",
  "NeuronSD",
  "NeurDensity",
  "Glia_N",
  "GliaSD",
  "GliaDensity",
  "Astro_N",
  "AstroSD",
  "AstroDensity",
  "Oligo_N",
  "OligoSD",
  "OligoDensity",
  "Microglia_N",
  "MicroSD",
  "MicroDensity",
  "GliaNeurDensityRatio"
)]

# Rename columns for clarity
names(data) <- c(
  "anatomy_group",
  "rcmr_value",
  "volume",
  "volume_sd",
  "neurons_n",
  "neurons_sd",
  "neuron_densities",
  "glia_n",
  "glia_sd",
  "glia_densities",
  "astrocytes_n",
  "astrocytes_sd",
  "astrocyte_densities",
  "oligodendrocytes_n",
  "oligodendrocytes_sd",
  "oligodendrocyte_densities",
  "microglia_n",
  "microglia_sd",
  "microglia_densities",
  "glia_to_neuron_ratio"
)

# Create a new dataframe that omits rows where there is only data for anatomy_group and rcmr_value, but no data for the other variables (i.e. where all the other variables are NA)
analysis_df <- data[!(is.na(data$neurons_n) & is.na(data$glia_n) & is.na(data$astrocytes_n) & is.na(data$oligodendrocytes_n) & is.na(data$microglia_n) & is.na(data$neuron_densities) & is.na(data$glia_densities) & is.na(data$astrocyte_densities) & is.na(data$oligodendrocyte_densities) & is.na(data$microglia_densities) & is.na(data$glia_to_neuron_ratio)), ]
analysis_df$anatomy_group <- gsub("_", " ", analysis_df$anatomy_group)
##########################################################################
# Correlations between rCMRGlc and (1) cell counts and (2) cell densities
##########################################################################

cell_counts <- c(
  "neurons_n",
  "glia_n",
  "astrocytes_n",
  "oligodendrocytes_n",
  "microglia_n",
  "glia_to_neuron_ratio"
)

cell_densities <- c(
  "neuron_densities",
  "glia_densities",
  "astrocyte_densities",
  "oligodendrocyte_densities",
  "microglia_densities",
  "glia_to_neuron_ratio"
)

###########################
## MULTIVARIATE REGRESSION
###########################
# Note: this is just to get an overall sense of whether the cell counts and densities are significant predictors of rCMRGlc, but the results should be interpreted with caution given the small sample size and potential multicollinearity between predictors.

fit_counts <- lm(
  reformulate(cell_counts, response = "rcmr_value"),
  data = analysis_df
)
summary(fit_counts)

fit_densities <- lm(
  reformulate(cell_densities, response = "rcmr_value"),
  data = analysis_df
)
summary(fit_densities)

###########################
## CORRELATIONS (LOOP)
###########################

# Correlation results for cell counts
cor_results_counts <- do.call(
  rbind,
  lapply(cell_counts, function(v) {
    ct <- cor.test(
      analysis_df$rcmr_value,
      analysis_df[[v]],
      method = "spearman",
      exact = FALSE
    )
    data.frame(
      variable = v,
      r        = unname(ct$estimate),
      p_value  = ct$p.value,
      n        = sum(complete.cases(analysis_df$rcmr_value, analysis_df[[v]]))
    )
  })
) %>%
  dplyr::arrange(p_value) %>%
  dplyr::mutate(p_adj_BH = p.adjust(p_value, method = "BH"))

cor_results_counts

# Correlation results for cell densities
cor_results_densities <- do.call(
  rbind,
  lapply(cell_densities, function(v) {
    ct <- cor.test(
      analysis_df$rcmr_value,
      analysis_df[[v]],
      method = "spearman",
      exact = FALSE
    )
    data.frame(
      variable = v,
      r        = unname(ct$estimate),
      p_value  = ct$p.value,
      n        = sum(complete.cases(analysis_df$rcmr_value, analysis_df[[v]]))
    )
  })
) %>%
  dplyr::arrange(p_value) %>%
  dplyr::mutate(p_adj_BH = p.adjust(p_value, method = "BH"))

cor_results_densities

###########################
## FACETED SCATTERPLOT
###########################

# counts

# Create a long-format dataframe for plotting
plot_df_counts <- analysis_df %>%
  dplyr::select(anatomy_group, rcmr_value, all_of(cell_counts)) %>%
  tidyr::pivot_longer(
    cols = all_of(cell_counts),
    names_to = "cell_counts",
    values_to = "prop"
  ) %>%
  dplyr::mutate(
    cell_counts = gsub("_", " ", cell_counts)  # remove underscores in facet labels
  )

# Build high-contrast palette
regions <- sort(unique(plot_df_counts$anatomy_group))
base_cols <- c(
  "#000000", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
)
pal <- grDevices::colorRampPalette(base_cols)(length(regions))
names(pal) <- regions

# Plot with regression line and equation
ggplot(plot_df_counts, aes(x = prop, y = rcmr_value, color = anatomy_group)) +
  geom_point(size = 2.8, alpha = 0.85) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "steelblue") +
  stat_poly_eq(
    aes(label = paste(..rr.label.., ..p.value.label.., sep = "*\", \"*")),
    formula = y ~ x,
    parse = TRUE,
    label.x = "right",
    label.y = "top",
    size = 4,
    color = "black"
  ) +
  facet_wrap(~ cell_counts, scales = "free_x") +
  scale_color_manual(values = pal) +
  labs(
    x = "Cell-type counts",
    y = "rCMRGlc (µmol/100 g/min.)",
    color = "Region"
  ) +
  theme_classic(base_size = 14)

# densities

# Create a long-format dataframe for plotting
plot_df_densities <- analysis_df %>%
  dplyr::select(anatomy_group, rcmr_value, all_of(cell_densities)) %>%
  tidyr::pivot_longer(
    cols = all_of(cell_densities),
    names_to = "cell_densities",
    values_to = "prop"
  ) %>%
  dplyr::mutate(
    cell_densities = gsub("_", " ", cell_densities)  # remove underscores in facet labels
  )

# Build high-contrast palette
regions <- sort(unique(plot_df_densities$anatomy_group))
base_cols <- c(
  "#000000", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
)
pal <- grDevices::colorRampPalette(base_cols)(length(regions))
names(pal) <- regions

# Plot with regression line and equation
ggplot(plot_df_densities, aes(x = prop, y = rcmr_value, color = anatomy_group)) +
  geom_point(size = 2.8, alpha = 0.85) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "steelblue") +
  stat_poly_eq(
    aes(label = paste(..rr.label.., ..p.value.label.., sep = "*\", \"*")),
    formula = y ~ x,
    parse = TRUE,
    label.x = "right",
    label.y = "top",
    size = 4,
    color = "black"
  ) +
  facet_wrap(~ cell_densities, scales = "free_x") +
  scale_color_manual(values = pal) +
  labs(
    x = "Cell-type densities",
    y = "rCMRGlc (µmol/100 g/min.)",
    color = "Region"
  ) +
  theme_classic(base_size = 14)

