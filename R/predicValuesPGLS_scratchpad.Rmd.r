##################################################################################
# Make best plot for ppt - quadratic fit with on-plot stats,
# exclude Neocortex white from fit, then add Neocortex white on top
##################################################################################

quad_fit <- lm(Diff.pre ~ poly(rCMRGlc, 2, raw = TRUE), data = df)

neo <- core_with_rCMRGlc %>%
  filter(Model == "Brownian (BM)", Structure == "Neocortex white") %>%
  transmute(
    Structure,
    Diff.pre = as.numeric(Diff.pre),
    rCMRGlc  = as.numeric(rCMRGlc)
  ) %>%
  drop_na(Diff.pre, rCMRGlc)

# stats text
sm <- summary(quad_fit)
adj_r2 <- sm$adj.r.squared
f <- sm$fstatistic
model_p <- pf(f[1], f[2], f[3], lower.tail = FALSE)
b <- coef(quad_fit)
eqn <- paste0(
  "y = ", sprintf("%.4f", b[1]),
  ifelse(b[2] >= 0, " + ", " - "), sprintf("%.4f", abs(b[2])), "x",
  ifelse(b[3] >= 0, " + ", " - "), sprintf("%.4f", abs(b[3])), "x^2"
)

xg <- seq(10, 40, length.out = 300)
pr <- predict(quad_fit, newdata = data.frame(rCMRGlc = xg), interval = "confidence")

plot(df$rCMRGlc, df$Diff.pre,
     xlim = c(10, 40),
     ylim = c(-3.5, 1.5),
     pch = 16, cex = 1.3,
     xlab = "rCMRGlc (µmol/100 g/min)",
     ylab = "Difference from prediction (BM)")

polygon(c(xg, rev(xg)), c(pr[, "lwr"], rev(pr[, "upr"])),
        border = NA, col = adjustcolor("steelblue", 0.2))
lines(xg, pr[, "fit"], col = "steelblue4", lwd = 3)

# add excluded point after drawing line
points(neo$rCMRGlc, neo$Diff.pre, pch = 24, bg = "firebrick2", col = "firebrick4", cex = 1.8)
#text(neo$rCMRGlc, neo$Diff.pre, "Neocortex white (excluded)", pos = 4, cex = 0.95, col = "firebrick4")

# points
points(df$rCMRGlc, df$Diff.pre, pch = 16)

# labels for all points
text(df$rCMRGlc, df$Diff.pre,
     labels = df$Structure,
     pos = 4,        # 1=below, 2=left, 3=above, 4=right
     offset = 0.3,
     cex = 0.75)

points(neo$rCMRGlc, neo$Diff.pre, pch = 24, bg = "firebrick2")
text(neo$rCMRGlc, neo$Diff.pre,
     labels = "Neocortex white (excluded)",
     pos = 4, cex = 0.9, col = "firebrick4", xpd = NA)

# on-plot stats
usr <- par("usr")
text(usr[1] + 0.02*diff(usr[1:2]), usr[4] - 0.03*diff(usr[3:4]),
     labels = paste0("Adj. R² = ", sprintf("%.3f", adj_r2),
                     "\nModel p = ", format.pval(model_p, digits = 2, eps = 1e-4),
                     "\n", eqn),
     adj = c(0,1), cex = 0.95, col = "navy")






###################
# Polynomial comparison + on-plot stats (no console output)
###################

# ---------- Regions of interest ----------
regions_of_interest <- cor_df %>%
  filter(Structure != "Neocortex white") %>%
  pull(Structure) %>%
  unique()

# ---------- Prepare data ----------
df <- core_with_rCMRGlc %>%
  filter(
    Model == "Brownian (BM)",
    Structure %in% regions_of_interest
  ) %>%
  transmute(
    Structure,
    Diff.pre = as.numeric(Diff.pre),
    rCMRGlc  = as.numeric(rCMRGlc)
  ) %>%
  drop_na(Diff.pre, rCMRGlc)

# ---------- Fit polynomial models degree 1..5 ----------
fits <- lapply(1:5, function(d) {
  if (d == 1) {
    lm(Diff.pre ~ rCMRGlc, data = df)
  } else {
    lm(Diff.pre ~ poly(rCMRGlc, d, raw = TRUE), data = df)
  }
})

adj_r2 <- sapply(fits, function(f) summary(f)$adj.r.squared)
best_degree <- which.max(adj_r2)
best_fit <- fits[[best_degree]]
best_sum <- summary(best_fit)

# Model p-value (overall F-test)
fstat <- best_sum$fstatistic
model_p <- pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)

# ---------- Build equation string ----------
fmt_num <- function(x) formatC(x, digits = 4, format = "f")
coefs <- coef(best_fit)

eqn <- paste0("y = ", fmt_num(coefs[1]))
if (length(coefs) > 1) {
  for (i in 2:length(coefs)) {
    sign_txt <- if (coefs[i] >= 0) " + " else " - "
    pow <- i - 1
    term <- paste0(fmt_num(abs(coefs[i])), "·x", if (pow > 1) paste0("^", pow) else "")
    eqn <- paste0(eqn, sign_txt, term)
  }
}

# Wrap long equation so it stays in plot area
eqn_wrapped <- paste(strwrap(eqn, width = 60), collapse = "\n")

# ---------- Prediction grid ----------
x_axis <- seq(min(df$rCMRGlc), max(df$rCMRGlc), length.out = 300)
curve_df <- bind_rows(lapply(1:5, function(d) {
  tibble(
    degree = d,
    rCMRGlc = x_axis,
    y = predict(fits[[d]], newdata = data.frame(rCMRGlc = x_axis))
  )
}))

# ---------- Plot ----------
cols <- c("green4", "red3", "purple3", "dodgerblue3", "orange3")

plot(df$rCMRGlc, df$Diff.pre, pch = 16,
     xlab = "rCMRGlc (µmol/100 g/min)",
     ylab = "Diff.pre (BM)",
     main = "Polynomial fits (Degree 1–5)")

# all model curves
for (d in 1:5) {
  lines(
    curve_df$rCMRGlc[curve_df$degree == d],
    curve_df$y[curve_df$degree == d],
    col = cols[d], lwd = 2
  )
}

# best-fit curve highlighted
lines(
  curve_df$rCMRGlc[curve_df$degree == best_degree],
  curve_df$y[curve_df$degree == best_degree],
  col = "black", lwd = 3, lty = 2
)

# point labels
text(df$rCMRGlc, df$Diff.pre, labels = df$Structure, pos = 4, cex = 0.75)

# legend with Adj R² per degree
legend(
  "topright",
  legend = c(
    paste0("Degree ", 1:5, " (Adj R²=", sprintf("%.3f", adj_r2), ")"),
    paste0("Best fit: Degree ", best_degree)
  ),
  col = c(cols, "black"),
  lwd = c(rep(2, 5), 3),
  lty = c(rep(1, 5), 2),
  cex = 0.8,
  bty = "n"
)

# On-plot annotation (instead of console output)
usr <- par("usr")
x_annot <- usr[1] + 0.02 * (usr[2] - usr[1])
y_annot <- usr[4] - 0.02 * (usr[4] - usr[3])

stats_txt <- paste(
  paste0("Best degree: ", best_degree),
  paste0("Adj R²: ", sprintf("%.4f", adj_r2[best_degree])),
  paste0("Model p: ", format.pval(model_p, digits = 3, eps = 1e-4)),
  eqn_wrapped,
  sep = "\n"
)

text(x_annot, y_annot, labels = stats_txt, adj = c(0, 1), cex = 0.8, col = "blue4")





# Draw line without NeoW_Frahm (outlier) for visual clarity (but keep the point on the scatter)

# First filter for Brownian model
temp <- core_with_rCMRGlc %>% filter(Model == "Brownian (BM)")

# Data used for the loess fit: exclude "Neocortex white"
temp_no_nc_white <- temp %>% filter(Structure != "Neocortex white")  # exclude only for fit

# Option A — Exclude from the loess fit only (keep the point on the scatter)
plot_loess <- ggplot(temp, aes(x = rCMRGlc, y = `Diff.pre`)) +
  geom_point() +
  # loess line computed without "Neocortex white"
  geom_smooth(
    data = temp_no_nc_white,
    aes(x = rCMRGlc, y = `Diff.pre`),
    method = "loess",
    se = FALSE,
    color = "#377eb8",
    span = 0.75
  ) +
  coord_cartesian(ylim = c(-3.5, 0.6)) +
  theme_bw() +
  labs(
    title = "Correlation between rCMRGlc and Prediction Error (BM)",
    x = "rCMRGlc (µmol/100 g/min.)",
    y = "Prediction Error (Diff.pre) under Brownian Model"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.title = element_text(size = 10)
  ) +
  geom_text(aes(label = Structure), nudge_x = 0.5, size = 3, check_overlap = TRUE)

plot_loess

# Option B — Use weights to exclude from the fit (keep the point)
plot_loess_2 <- ggplot(temp, aes(x = rCMRGlc, y = `Diff.pre`)) +
  geom_point() +
  geom_smooth(
    aes(weight = ifelse(Structure == "Neocortex white", 0, 1)),
    method = "loess",
    se = FALSE,
    color = "#377eb8",
    span = 0.75
  ) +
  coord_cartesian(ylim = c(-3.5, 0.6)) +
  theme_bw() +
  labs(
    title = "Correlation between rCMRGlc and Prediction Error (BM)",
    x = "rCMRGlc (µmol/100 g/min.)",
    y = "Prediction Error (Diff.pre) under Brownian Model"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.title = element_text(size = 10)
  ) +
  geom_text(aes(label = Structure), nudge_x = 0.5, size = 3, check_overlap = TRUE)

plot_loess_2

####################
# Plot with loess line including all points (no exclusion) and showing error bars for Diff.min and Diff.max 
temp <- core_with_rCMRGlc %>%
  filter(Model == "Brownian (BM)") %>%
  mutate(
    rCMRGlc = readr::parse_number(as.character(rCMRGlc)),
  )
plot_loess_bars<-ggplot(temp, aes(x = rCMRGlc, y = Diff.pre)) +
  geom_errorbar(aes(ymin = Diff.min, ymax = Diff.max), width = 0) +
  geom_point(size = 2.5) +
  geom_smooth(method = "loess", se = FALSE, span = 0.75) +
  geom_text(aes(label = Structure), nudge_x = 0.5, size = 3, check_overlap = TRUE) +
  theme_bw() +
  labs(x = "rCMRGlc (µmol/100 g/min.)", y = "Diff.pre (BM)")
plot_loess_bars  

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ggrepel)

# ---------- Prepare data ----------
df <- core_with_rCMRGlc %>%
  filter(
    Model == "Brownian (BM)",
    Structure %in% setdiff(unique(cor_df$Structure), "Neocortex white")
  ) %>%
  transmute(
    Structure,
    Diff.pre = as.numeric(Diff.pre),
    rCMRGlc  = as.numeric(rCMRGlc)
  ) %>%
  drop_na(Diff.pre, rCMRGlc)

# ---------- Fit degree 1:5 polynomial models ----------
fit_poly <- function(d) {
  form <- if (d == 1) {
    Diff.pre ~ rCMRGlc
  } else {
    as.formula(sprintf("Diff.pre ~ poly(rCMRGlc, %d, raw = TRUE)", d))
  }
  lm(form, data = df)
}

fits <- map(1:5, fit_poly)

adj_r2 <- map_dbl(fits, ~ summary(.x)$adj.r.squared)
best_degree <- which.max(adj_r2)
best_fit <- fits[[best_degree]]
best_sum <- summary(best_fit)

model_p <- with(
  as.list(best_sum$fstatistic),
  pf(value, numdf, dendf, lower.tail = FALSE)
)

# ---------- Equation text ----------
fmt <- function(x) formatC(x, digits = 4, format = "f")

make_eqn <- function(coefs) {
  terms <- map_chr(seq_along(coefs), function(i) {
    if (i == 1) return(fmt(coefs[i]))
    
    pow <- i - 1
    paste0(
      ifelse(coefs[i] >= 0, " + ", " - "),
      fmt(abs(coefs[i])),
      "·x",
      ifelse(pow > 1, paste0("^", pow), "")
    )
  })
  
  paste0("y = ", paste0(terms, collapse = ""))
}

eqn <- make_eqn(coef(best_fit)) |>
  strwrap(width = 60) |>
  paste(collapse = "\n")

stats_txt <- paste(
  paste0("Best degree: ", best_degree),
  paste0("Adj R²: ", sprintf("%.4f", adj_r2[best_degree])),
  paste0("Model p: ", format.pval(model_p, digits = 3, eps = 1e-4)),
  eqn,
  sep = "\n"
)

# ---------- Prediction grid ----------
x_grid <- seq(min(df$rCMRGlc), max(df$rCMRGlc), length.out = 300)

curve_df <- map_dfr(1:5, function(d) {
  tibble(
    degree = d,
    degree_label = paste0("Degree ", d, " (Adj R²=", sprintf("%.3f", adj_r2[d]), ")"),
    rCMRGlc = x_grid,
    Diff.pre = predict(fits[[d]], newdata = data.frame(rCMRGlc = x_grid))
  )
})

degree_labels <- unique(curve_df$degree_label)

cols <- setNames(
  c("green4", "red3", "purple3", "dodgerblue3", "orange3"),
  degree_labels
)

# ---------- Plot ----------
ggplot(df, aes(rCMRGlc, Diff.pre)) +
  geom_point(size = 2) +
  geom_line(
    data = curve_df,
    aes(y = Diff.pre, color = degree_label),
    linewidth = 0.9
  ) +
  geom_line(
    data = filter(curve_df, degree == best_degree),
    aes(y = Diff.pre),
    color = "black",
    linewidth = 1.2,
    linetype = "dashed"
  ) +
  geom_text_repel(
    aes(label = Structure),
    size = 3,
    max.overlaps = Inf
  ) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = stats_txt,
    hjust = -0.05,
    vjust = 1.05,
    size = 3,
    color = "blue4"
  ) +
  scale_color_manual(values = cols, name = NULL) +
  labs(
    title = "Polynomial fits, degree 1–5",
    subtitle = paste0("Best fit: degree ", best_degree),
    x = "rCMRGlc (µmol/100 g/min)",
    y = "Diff.pre (BM)"
  ) +
  theme_classic()
##################################################################################
# ggplot version: p5
# quadratic fit with on-plot stats,
# exclude Neocortex white from fit, then add Neocortex white on top
##################################################################################

## Optional safeguard: make sure df excludes Neocortex white
## If df already excludes it, you can skip this.
df_fit <- df %>%
  filter(Structure != "Neocortex white") %>%
  drop_na(Diff.pre, rCMRGlc)

quad_fit <- lm(Diff.pre ~ poly(rCMRGlc, 2, raw = TRUE), data = df_fit)

neo <- core_with_rCMRGlc %>%
  filter(Model == "Brownian (BM)", Structure == "Neocortex white") %>%
  transmute(
    Structure,
    Diff.pre = as.numeric(Diff.pre),
    rCMRGlc  = as.numeric(rCMRGlc)
  ) %>%
  drop_na(Diff.pre, rCMRGlc)

## stats text
sm <- summary(quad_fit)
adj_r2 <- sm$adj.r.squared
f <- sm$fstatistic
model_p <- pf(f[1], f[2], f[3], lower.tail = FALSE)
b <- coef(quad_fit)

eqn <- paste0(
  "y = ", sprintf("%.4f", b[1]),
  ifelse(b[2] >= 0, " + ", " - "), sprintf("%.4f", abs(b[2])), "x",
  ifelse(b[3] >= 0, " + ", " - "), sprintf("%.4f", abs(b[3])), "x²"
)

## prediction grid
xg <- seq(10, 40, length.out = 300)

pr <- predict(
  quad_fit,
  newdata = data.frame(rCMRGlc = xg),
  interval = "confidence"
)

pred_df <- data.frame(
  rCMRGlc = xg,
  fit = pr[, "fit"],
  lwr = pr[, "lwr"],
  upr = pr[, "upr"]
)

stats_label <- paste0(
  "Adj. R² = ", sprintf("%.3f", adj_r2),
  "\nModel p = ", format.pval(model_p, digits = 2, eps = 1e-4),
  "\n", eqn
)

p5 <- ggplot(df_fit, aes(x = rCMRGlc, y = Diff.pre)) +
  geom_ribbon(
    data = pred_df,
    aes(x = rCMRGlc, ymin = lwr, ymax = upr),
    inherit.aes = FALSE,
    fill = "steelblue",
    alpha = 0.2
  ) +
  geom_line(
    data = pred_df,
    aes(x = rCMRGlc, y = fit),
    inherit.aes = FALSE,
    color = "steelblue4",
    linewidth = 1.2
  ) +
  geom_point(
    pch = 16,
    size = 3.2
  ) +
  geom_text(
    aes(label = Structure),
    hjust = -0.08,
    size = 4.2,
    check_overlap = TRUE
  ) +
  geom_point(
    data = neo,
    aes(x = rCMRGlc, y = Diff.pre),
    inherit.aes = FALSE,
    shape = 24,
    fill = "firebrick2",
    color = "firebrick4",
    size = 5
  ) +
  geom_text(
    data = neo,
    aes(x = rCMRGlc, y = Diff.pre, label = "Neocortex white (excluded)"),
    inherit.aes = FALSE,
    hjust = -0.08,
    size = 4.4,
    color = "firebrick4"
  ) +
  annotate(
    "text",
    x = 10.6,
    y = 1.35,
    label = stats_label,
    hjust = 0,
    vjust = 1,
    size = 4.6,
    color = "navy"
  ) +
  coord_cartesian(
    xlim = c(10, 40),
    ylim = c(-3.5, 1.5),
    clip = "on"
  ) +
  labs(
    x = "rCMRGlc (µmol/100 g/min)",
    y = "Difference from prediction (BM)"
  ) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 15),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 15),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.margin = margin(8, 8, 8, 8)
  )

p5

ggsave(
  filename = "plot5_quadratic_fit_BM_rCMRGlc.pdf",
  plot = p5,
  width = 7,
  height = 5
)

ggsave(
  filename = "plot5_quadratic_fit_BM_rCMRGlc.png",
  plot = p5,
  width = 10,
  height = 7,
  dpi = 300
)