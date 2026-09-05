setwd(local({
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d)
  d
}))

# Siletti cerebral-cortex cell classes in the project's cerebral-cortex scope.
# Superclusters associated with thalamus, hypothalamus, midbrain, hindbrain,
# and cerebellum are excluded even if rare cells occur in a cortical ROI.
# The primary E:I definition excludes medium spiny neurons (MSNs) from the I
# denominator. A separately labelled sensitivity figure includes MSNs.

suppressPackageStartupMessages(library(tidyverse))
source("helpers/read_heiss_rates.R")
source("helpers/plot_settings.R")
source("helpers/s1b_EI.R")

dir.create("data_analysis", showWarnings = FALSE, recursive = TRUE)
dir.create("figs/s1b", showWarnings = FALSE, recursive = TRUE)

output_stem <-
  "cerebral_cortex_cell_classes_x_cerebral_cortex_regions_EI_vs_rcmrglc"

comparison_data <- calculate_s1b_cerebral_cortex_class_ei()
comparison_stats <- comparison_data %>%
  group_by(MSN_included_in_denominator) %>%
  group_modify(~ summarize_s1b_ei(.x)) %>%
  ungroup()

write.csv(
  comparison_data,
  paste0("data_analysis/", output_stem, ".csv"),
  row.names = FALSE
)
write.csv(
  comparison_stats,
  paste0("data_analysis/", output_stem, "_statistics.csv"),
  row.names = FALSE
)

save_comparison_plot <- function(include_msn, filename_suffix, subtitle) {
  plot_data <- comparison_data %>%
    filter(MSN_included_in_denominator == include_msn)
  plot_stats <- comparison_stats %>%
    filter(MSN_included_in_denominator == include_msn)

  p <- plot_s1b_ei(
    data = plot_data,
    stats = plot_stats,
    title = "Cerebral-cortex cell classes × cerebral-cortex regions",
    subtitle = subtitle,
    x_label = "E:I ratio from Siletti cerebral-cortex cell classes"
  )

  figure_stem <- paste0("figs/s1b/", output_stem, filename_suffix)
  ggsave(paste0(figure_stem, ".pdf"), p, width = 8.5, height = 6.8, units = "in")
  ggsave(
    paste0(figure_stem, ".jpg"),
    p,
    width = 8.5,
    height = 6.8,
    units = "in",
    dpi = 300
  )
}

save_comparison_plot(
  include_msn = FALSE,
  filename_suffix = "",
  subtitle = "Cerebral-cortex classes; off-scope classes and MSNs excluded"
)
save_comparison_plot(
  include_msn = TRUE,
  filename_suffix = "_MSN_included_sensitivity",
  subtitle = "Off-scope classes excluded; MSNs included as a sensitivity"
)

print(comparison_stats)
