setwd(local({
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d)
  d
}))

# Scope-extrapolation analysis: the same neocortical E:I cell-class crosswalk
# is applied to the broader cerebral-cortex scope, which additionally includes
# Corpus amygdaloideum and Hippocampus. This is not a Jorstad anatomical scope.

suppressPackageStartupMessages(library(tidyverse))
source("helpers/read_heiss_rates.R")
source("helpers/plot_settings.R")
source("helpers/s1b_EI.R")

run_s1b_neocortical_class_ei_scope(
  scope_flag = "is_cerebral_cortex",
  scope_label = "Cerebral cortex (including amygdaloid complex and hippocampus)",
  interpretation = "Anatomical-scope extrapolation; not a Jorstad anatomical definition",
  title = "Neocortical cell classes × cerebral-cortex regions",
  subtitle = "Scope extrapolation: adds amygdaloid complex and hippocampus",
  output_stem =
    "neocortical_cell_classes_x_cerebral_cortex_regions_EI_vs_rcmrglc"
)
