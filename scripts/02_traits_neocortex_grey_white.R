# ------------------------------------------------------------
# plot_neocortex_traits.R
# Comparative neuroanatomy background plot (clean version)
# ------------------------------------------------------------

# Set working directory
setwd(local({ d <- normalizePath(getwd()); while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d); d }))  # repo root (portable; replaces hardcoded path -- see helpers/project_root.R)

# The pipeline startup synchronizer manages the Evo-M1 root and refreshes this
# versioned local snapshot. Keeping that dependency logic out of this analysis
# ensures interactive and pipeline runs read the same file.
volume_path <- "data_raw/volumes_wide.csv"
if (!file.exists(volume_path)) {
  stop(
    "Missing synchronized Evo-M1 input: ", volume_path,
    ". Run scripts/00_sync_evo_m1_inputs.R first.",
    call. = FALSE
  )
}

df <- read.csv(volume_path, stringsAsFactors = FALSE, check.names = FALSE)

required_columns <- c(
  "Species",
  "Neocortex_Vol.mm3",
  "Neocortex_grey_matter_Vol.mm3",
  "Neocortex_white_matter_Vol.mm3"
)
missing_columns <- setdiff(required_columns, names(df))
if (length(missing_columns)) {
  stop(
    "volumes_wide.csv is missing required column(s): ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

# Keep the primate scope defined by the project's phylogeny. The alias helper
# only standardizes species labels; it does not supply anatomical measurements.
if (!requireNamespace("ape", quietly = TRUE)) {
  stop("Package 'ape' is required to read data_raw/species.nwk.", call. = FALSE)
}
source("helpers/species_aliases.R")
tree <- ape::read.tree("data_raw/species.nwk")
df$tree_label <- canon_species(df$Species)
df <- df[df$tree_label %in% tree$tip.label, , drop = FALSE]

# ---- Prepare variables ----
df$logNeo   <- log10(df$Neocortex_Vol.mm3)
df$logGray  <- log10(df$Neocortex_grey_matter_Vol.mm3)
df$logWhite <- log10(df$Neocortex_white_matter_Vol.mm3)

# Use a common finite, positive species set for a directly comparable pair of
# regressions. Duplicate taxon aliases with no measurements disappear here;
# multiple complete rows for one tree tip remain a hard error.
df <- subset(
  df,
  is.finite(logNeo) & is.finite(logGray) & is.finite(logWhite)
)

if (anyDuplicated(df$tree_label)) {
  duplicated_species <- unique(df$tree_label[duplicated(df$tree_label)])
  stop(
    "volumes_wide.csv has multiple complete rows resolving to tree tip(s): ",
    paste(duplicated_species, collapse = ", "),
    call. = FALSE
  )
}

dG <- df
dW <- df

if (nrow(dG) < 3L || nrow(dW) < 3L) {
  stop(
    "Too few complete positive primate observations for the neocortex plot ",
    "(grey n = ", nrow(dG), ", white n = ", nrow(dW), ").",
    call. = FALSE
  )
}

message(
  "Neocortex traits snapshot: ", normalizePath(volume_path),
  "; grey n = ", nrow(dG), ", white n = ", nrow(dW)
)

# ---- Fit regressions ----
fitG <- lm(logGray ~ logNeo, data = dG)
fitW <- lm(logWhite ~ logNeo, data = dW)

coefG <- coef(fitG)
coefW <- coef(fitW)

r2G <- summary(fitG)$r.squared
r2W <- summary(fitW)$r.squared

# ---- Labels (2 significant figures, slope-focused) ----
labelG <- paste0(
  "Gray: slope = ", signif(coefG[2], 3),
  "\nR² = ", signif(r2G, 4)
)

labelW <- paste0(
  "White: slope = ", signif(coefW[2], 3),
  "\nR² = ", signif(r2W, 4)
)

# ---- Output folder ----
dir.create("figs/traits", recursive = TRUE, showWarnings = FALSE)

# ---- Draw routine (base graphics; replayed once per device) ----
# Wrapped in a function so the identical plot can be written to every device
# without relying on dev.copy() of an interactive window (which is a fragile
# "snapshot" and fails when the script is run head-less via Rscript).
draw_neocortex_plot <- function() {
  plot(NULL,
       xlim = range(df$logNeo, na.rm = TRUE),
       ylim = range(c(df$logGray, df$logWhite), na.rm = TRUE),
       xlab = expression(log[10]*" neocortex volume (mm"^3*")"),
       ylab = expression(log[10]*" gray or white matter volume (mm"^3*")"))

  # Points
  points(dG$logNeo, dG$logGray,
         pch = 16, col = "grey40", cex = 1.3)
  points(dW$logNeo, dW$logWhite,
         pch = 1, col = "black", cex = 1.3)

  # Regression lines (match aesthetics)
  abline(fitG, col = "grey40", lwd = 2)
  abline(fitW, col = "black", lwd = 2)

  # Text labels (matched to data)
  text(3.5, 4.6, labelG, col = "grey40", adj = 0)
  text(3.5, 2.2, labelW, col = "black", adj = 0)

  # Legend
  legend(x = 3, y = 6,
         legend = c("Gray matter", "White matter"),
         pch = c(16, 1),
         col = c("grey40", "black"),
         bty = "n")
}

# ---- Optional interactive preview ----
if (interactive()) {
  dev.new()
  draw_neocortex_plot()
}

# ---- SAVE: raster for slides (PNG) + vector for print/PDF ----
png("figs/traits/neocortex_gray_white.png",
    width = 800, height = 600)
draw_neocortex_plot()
dev.off()

# cairo_pdf embeds the Unicode superscripts / R^2 glyph reliably.
cairo_pdf("figs/traits/neocortex_gray_white.pdf",
          width = 8, height = 6)
draw_neocortex_plot()
dev.off()
