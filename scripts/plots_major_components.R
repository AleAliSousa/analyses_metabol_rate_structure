# ------------------------------------------------------------
# plot_brain_components_vs_total.R
# Major brain components vs total brain volume
# ------------------------------------------------------------

# Set working directory
setwd("~/Library/CloudStorage/Dropbox/COLLABORATIVE/Do expensive brain regions increase less in humans/analyses_metabol_rate_structure")

# Load data
df <- read.csv("data_raw/Stephan_primates.csv", stringsAsFactors = FALSE)

# ---- Define columns ----
# Total brain
df$logTotal <- log10(df$Total_brain_net_volume)

# Components
components <- c("Medulla_oblongata",
                "Cerebellum",
                "Mesencephalon",
                "Diencephalon",
                "Telencephalon")

# Colors (distinct but clean)
cols <- c("steelblue", "darkorange", "forestgreen", "purple", "firebrick")

# ---- Output folder ----
dir.create("figs/traits", recursive = TRUE, showWarnings = FALSE)

# ---- OPEN plotting window ----
dev.new()

# ---- Empty plot ----
plot(NULL,
     xlim = range(log10(df[,components]), na.rm = TRUE),
     ylim = range(df$logTotal, na.rm = TRUE),
     xlab = expression(log[10]*" component volume (mm"^3*")"),
     ylab = expression(log[10]*" total brain volume (mm"^3*")"))

# ---- Loop over components ----
fits <- list()

for (i in seq_along(components)) {
  
  comp <- components[i]
  
  # log transform
  x <- log10(df[[comp]])
  y <- df$logTotal
  
  # subset complete cases
  d <- data.frame(x, y)
  d <- d[complete.cases(d), ]
  
  # fit
  fit <- lm(y ~ x, data = d)
  fits[[comp]] <- fit
  
  # points
  points(d$x, d$y,
         pch = 16, col = cols[i], cex = 1.1)
  
  # line
  abline(fit, col = cols[i], lwd = 2)
}

# ---- Legend with slopes ----
labels <- sapply(seq_along(components), function(i) {
  fit <- fits[[components[i]]]
  slope <- signif(coef(fit)[2], 2)
  paste0(components[i], " (", slope, ")")
})

legend("topleft",
       legend = labels,
       col = cols,
       lwd = 2,
       pch = 16,
       bty = "n",
       title = "Component (slope)")

# ---- SAVE ----
dev.copy(png, "figs/traits/brain_components_vs_total.png",
         width = 900, height = 700)
dev.off()
