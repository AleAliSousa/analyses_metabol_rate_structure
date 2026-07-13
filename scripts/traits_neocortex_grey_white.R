# ------------------------------------------------------------
# plot_neocortex_traits.R
# Comparative neuroanatomy background plot (clean version)
# ------------------------------------------------------------

# Set working directory
setwd(local({ d <- normalizePath(getwd()); while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d); d }))  # repo root (portable; replaces hardcoded path -- see R/project_root.R)

# Load data
df <- read.csv("data_raw/Stephan_primates.csv", stringsAsFactors = FALSE)

# ---- Prepare variables ----
df$logNeo   <- log10(df$NeoWG)
df$logGray  <- log10(df$NeoG_Frahm)
df$logWhite <- log10(df$NeoW_Frahm)

# Subset complete cases
dG <- subset(df, !is.na(logNeo) & !is.na(logGray))
dW <- subset(df, !is.na(logNeo) & !is.na(logWhite))

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