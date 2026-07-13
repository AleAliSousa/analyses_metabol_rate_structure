# =====================================================================
# R/project_root.R  --  portable project-root detection
# ---------------------------------------------------------------------
# Replaces the old machine-specific line
#     setwd("~/Library/CloudStorage/Dropbox/.../analyses_metabol_rate_structure")
# which only worked on one computer and was the source of output files
# landing in the wrong folder (e.g. stray PNGs in the repo root).
#
# Two ways to use it:
#
#   (a) Sourced helper (recommended for new scripts). From any script whose
#       working directory is somewhere inside the repo:
#           source("R/project_root.R"); setwd(repo_root())
#       or simply
#           source("R/project_root.R")        # also sets wd to the root
#
#   (b) Inline one-liner (no dependency, no source path needed) -- this is
#       what the existing analysis scripts now use in place of setwd():
#           setwd(local({ d <- normalizePath(getwd())
#                         while (!file.exists(file.path(d, ".git")) &&
#                                dirname(d) != d) d <- dirname(d); d }))
#
# RStudio users: an .Rproj file at the repo root also lets here::here()
# resolve the root from anywhere.
# =====================================================================

repo_root <- function(start = getwd(),
                       sentinels = c(".git", ".Rproj", "data_raw")) {
  d <- normalizePath(start, mustWork = FALSE)
  repeat {
    hit <- any(file.exists(file.path(d, sentinels))) ||
           length(list.files(d, pattern = "\\.Rproj$")) > 0
    if (hit) return(d)
    parent <- dirname(d)
    if (parent == d)
      stop("repo_root(): no sentinel (", paste(sentinels, collapse = ", "),
           ") found above ", start)
    d <- parent
  }
}

# Sourcing this file sets the working directory to the detected root.
setwd(repo_root())
