# Shared setup for the nonhuman (mouse) s1b scripts.
# Sets the working directory to the repo root and defines NH, the folder that
# holds the mouse pipeline, so scripts can be run from anywhere.
setwd(local({
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d)
  d
}))
NH <- "analyses_metabol_rate_structure_nonhuman"
nh <- function(...) file.path(NH, ...)
for (d in c("data_raw", "data_intermediate", "tables/mouse_s1b", "figs/mouse_s1b"))
  dir.create(nh(d), recursive = TRUE, showWarnings = FALSE)
