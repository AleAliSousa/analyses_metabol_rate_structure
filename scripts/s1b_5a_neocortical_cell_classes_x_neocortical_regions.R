setwd(local({
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d)
  d
}))

# Jorstad-aligned analysis: the neocortical E:I cell-class crosswalk is applied
# only to s1b ROIs in the explicitly defined neocortex scope.

suppressPackageStartupMessages(library(tidyverse))
source("helpers/read_heiss_rates.R")
source("helpers/plot_settings.R")
source("helpers/s1b_EI.R")

run_s1b_neocortical_class_ei_scope(
  scope_flag = "is_neocortex",
  scope_label = "Neocortex",
  interpretation = "Jorstad-aligned anatomical scope",
  title = "Neocortical cell classes × neocortical regions",
  subtitle = "Jorstad-aligned anatomical scope",
  output_stem =
    "neocortical_cell_classes_x_neocortical_regions_EI_vs_rcmrglc"
)
