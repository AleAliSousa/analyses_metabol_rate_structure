setwd("~/Library/CloudStorage/Dropbox/COLLABORATIVE/Do expensive brain regions increase less in humans/analyses_metabol_rate_structure")

library(tidyverse)
library(scales)

# ============================================================
# Inputs
# ============================================================

core_with_rCMRGlc  <- read.csv("data/core_with_rCMRGlc_predicted_volumes.csv", check.names = FALSE)
celltype_proportion_table <- read.csv("data/nonneuronal_celltype_proportions_by_region.csv", check.names = FALSE)
num <- function(x) readr::parse_number(as.character(x))
# ============================================================
# Extract p_astrocyte
# ============================================================

astro_cols <- names(celltype_proportion_table)[
  str_detect(
    str_to_lower(names(celltype_proportion_table)),
    "^p_.*astro"
  )
]

if (length(astro_cols) == 0) {
  stop(
    "No astrocyte-like p_ column found in celltype_proportion_table. ",
    "Look at names(celltype_proportion_table) and adjust astro_cols manually."
  )
}

message("Astrocyte columns used:")
print(astro_cols)

p_astrocyte_by_anatomy <- celltype_proportion_table %>%
  mutate(
    p_astrocyte = rowSums(across(all_of(astro_cols)), na.rm = TRUE)
  ) %>%
  select(
    anatomy_group,
    p_astrocyte,
    n_all_nonneurons,
    all_of(astro_cols)
  )

# ============================================================
# Map transcriptomic anatomy groups onto PGLS volume structures
# ============================================================

structure_anatomy_map <- tribble(
  ~Structure, ~anatomy_group,
  
  "Insular cortex (grey)", "Insular lobe",
  
  "Hippocampus", "Hippocampus",
  
  "Cerebellum", "Cerebellar cortex",
  "Cerebellum", "Vermis",
  "Cerebellum", "Nucleus dentatus cerebelli",
  
  "Mesencephalon", "Colliculus inferior",
  "Mesencephalon", "Colliculus superior",
  "Mesencephalon", "Substantia nigra",
  "Mesencephalon", "Nucleus ruber",
  
  "Neocortex grey", "Frontal lobe",
  "Neocortex grey", "Parietal lobe",
  "Neocortex grey", "Temporal lobe",
  "Neocortex grey", "Occipital lobe",
  "Neocortex grey", "Cingulate / limbic",
  
  "Area striata grey", "Occipital lobe",
  
  "Striatum", "Caudatum",
  "Striatum", "Putamen",
  "Striatum", "Nucleus accumbens",
  
  "Capsula interna", "Capsula interna",
  
  "Nucleus subthalamicus Luysi", "Nucleus subthalamicus",
  
  "Neocortex white", "Centrum semiovale",
  
  "Pallidum", "Pallidum",
  
  "Amygdala", "Corpus amygdaloideum",
  
  "Corpus geniculatum laterale", "Corpus geniculatum laterale"
)

missing_anatomy_groups <- structure_anatomy_map %>%
  anti_join(p_astrocyte_by_anatomy, by = "anatomy_group")

if (nrow(missing_anatomy_groups) > 0) {
  warning("Some anatomy groups in the map were not found:")
  print(missing_anatomy_groups)
}

p_astrocyte_by_structure <- structure_anatomy_map %>%
  left_join(p_astrocyte_by_anatomy, by = "anatomy_group") %>%
  group_by(Structure) %>%
  summarise(
    p_astrocyte = weighted.mean(
      p_astrocyte,
      w = n_all_nonneurons,
      na.rm = TRUE
    ),
    astrocyte_anatomy_groups = paste(anatomy_group, collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(
    p_astrocyte = if_else(is.nan(p_astrocyte), NA_real_, p_astrocyte)
  )

# ============================================================
# Append p_astrocyte to Brownian predicted volume table
# ============================================================

df_bm_astro <- core_with_rCMRGlc %>%
  filter(Model == "Brownian (BM)") %>%
  transmute(
    Structure,
    rCMRGlc = num(rCMRGlc),
    Predicted = num(Predicted),
    Observed = num(Observed),
    CI_min = num(`95% CI min`),
    CI_max = num(`95% CI max`),
    Diff.pre = num(Diff.pre),
    Diff.min = num(Diff.min),
    Diff.max = num(Diff.max),
    excluded = Structure == "Neocortex white"
  ) %>%
  left_join(p_astrocyte_by_structure, by = "Structure")

write_csv(
  df_bm_astro,
  "tables/BM_predicted_volume_rCMRGlc_p_astrocyte.csv"
)

# ============================================================
# Plot: Quadratic fit with 95% confidence interval
# Same as old plot, but point color/size reflect p_astrocyte
# ============================================================

df_plot <- df_bm_astro %>%
  drop_na(rCMRGlc, Diff.pre)

df_fit <- df_plot %>%
  filter(!excluded)

df_neo <- df_plot %>%
  filter(excluded)

xg <- seq(10, 40, length.out = 300)

quad_fit <- lm(
  Diff.pre ~ poly(rCMRGlc, 2, raw = TRUE),
  data = df_fit
)

quad_sum <- summary(quad_fit)
quad_adj_r2 <- quad_sum$adj.r.squared
quad_fstat <- quad_sum$fstatistic
quad_p <- pf(
  quad_fstat[1],
  quad_fstat[2],
  quad_fstat[3],
  lower.tail = FALSE
)

b <- coef(quad_fit)

quad_eqn <- paste0(
  "y = ", sprintf("%.4f", b[1]),
  ifelse(b[2] >= 0, " + ", " - "), sprintf("%.4f", abs(b[2])), "x",
  ifelse(b[3] >= 0, " + ", " - "), sprintf("%.4f", abs(b[3])), "x²"
)

quad_stats <- paste0(
  "Adj. R² = ", sprintf("%.3f", quad_adj_r2),
  "\nModel p = ", format.pval(quad_p, digits = 2, eps = 1e-4),
  "\n", quad_eqn
)

quad_pred <- data.frame(rCMRGlc = xg)

quad_pred <- cbind(
  quad_pred,
  as.data.frame(
    predict(
      quad_fit,
      newdata = quad_pred,
      interval = "confidence"
    )
  )
)

theme_paper <- theme_bw(base_size = 16) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 15),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 15),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.margin = margin(8, 8, 8, 8),
    legend.position = "right"
  )

x_lab <- "rCMRGlc (µmol/100 g/min)"
y_lab <- "Prediction error"

add_points_and_labels <- list(
  
  geom_point(
    data = df_fit,
    aes(
      x = rCMRGlc,
      y = Diff.pre,
      color = p_astrocyte,
      size = p_astrocyte
    ),
    inherit.aes = FALSE,
    alpha = 0.9
  ),
  
  geom_text(
    data = df_fit,
    aes(
      x = rCMRGlc,
      y = Diff.pre,
      label = Structure
    ),
    inherit.aes = FALSE,
    hjust = -0.08,
    size = 2.8,
    check_overlap = TRUE
  ),
  
  geom_point(
    data = df_neo,
    aes(
      x = rCMRGlc,
      y = Diff.pre,
      color = p_astrocyte,
      size = p_astrocyte
    ),
    inherit.aes = FALSE,
    shape = 16,
    alpha = 0.95
  ),
  
  geom_text(
    data = df_neo,
    aes(
      x = rCMRGlc,
      y = Diff.pre,
      label = "Neocortex white (excluded)"
    ),
    inherit.aes = FALSE,
    hjust = -0.08,
    size = 3.1
  ),
  
  scale_color_viridis_c(
    labels = percent_format(accuracy = 1),
    na.value = "grey70"
  ),
  
  scale_size_continuous(
    labels = percent_format(accuracy = 1),
    range = c(2, 6)
  ),
  
  labs(
    color = "p_astrocyte",
    size = "p_astrocyte"
  )
)

p_quad_ci_pastrocyte <- ggplot() +
  geom_ribbon(
    data = quad_pred,
    aes(x = rCMRGlc, ymin = lwr, ymax = upr),
    inherit.aes = FALSE,
    fill = "steelblue",
    alpha = 0.2
  ) +
  geom_line(
    data = quad_pred,
    aes(x = rCMRGlc, y = fit),
    inherit.aes = FALSE,
    color = "steelblue4",
    linewidth = 1.2
  ) +
  add_points_and_labels +
  annotate(
    "text",
    x = 10.6,
    y = 1.35,
    label = quad_stats,
    hjust = 0,
    vjust = 1,
    size = 3.2,
    color = "navy"
  ) +
  coord_cartesian(
    xlim = c(10, 40),
    ylim = c(-3.5, 1.5),
    clip = "on"
  ) +
  labs(
    x = x_lab,
    y = y_lab,
    color = "p_astrocyte"
  ) +
  guides(size = "none") +
  theme_paper

p_quad_ci_pastrocyte