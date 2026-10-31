# Study s5 helpers (Coto-style analyses on Heiss 2004). New file; mirrors the
# nonhuman repo helpers/plot_settings.R so the two repos annotate and export
# figures identically. ggpmisc is not available in the conda channel, so lm
# annotations are geom_text labels built by annotate_lm().
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })

# Energy constants (documented once in the nonhuman 00_prepare_species_metabolism.R)
J_PER_ML_O2       <- 20.1   # oxycaloric equivalent
J_PER_UMOL_GLUCOSE <- 2.87  # dG(ox) ~ -2870 kJ/mol
TISSUE_DENSITY_MG_PER_MM3 <- 1.0  # task convention (s4a used 1.036 g/cc; see METHODS_NOTE_s5.md)

fmt_p <- function(p) ifelse(is.na(p), "NA", ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3)))

# Fit y ~ x (optionally per group) and return a label data frame with R^2 and p.
annotate_lm <- function(df, x, y, group = NULL) {
  f <- function(d) {
    d <- d[is.finite(d[[x]]) & is.finite(d[[y]]), , drop = FALSE]
    if (nrow(d) < 3) return(data.frame(n = nrow(d), r2 = NA_real_, p = NA_real_, slope = NA_real_))
    m <- lm(reformulate(x, y), data = d)
    s <- summary(m)
    data.frame(n = nrow(d), r2 = s$r.squared, p = s$coefficients[2, 4], slope = coef(m)[[2]])
  }
  out <- if (is.null(group)) f(df) else df %>% group_by(across(all_of(group))) %>% group_modify(~ f(.x)) %>% ungroup()
  out$label <- ifelse(is.na(out$r2), "italic(n) < 3",
                      sprintf("italic(R)^2~'='~'%.2f'*','~italic(p)~'='~'%s'", out$r2, fmt_p(out$p)))
  out
}

# Fallback palette for labels that are not in the shared region_palette.
region_base_cols <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
                      "#D55E00", "#CC79A7", "#882255", "#44AA99", "#117733", "#88CCEE",
                      "#999933", "#DDCC77", "#6699CC", "#661100", "#332288", "#AA4499",
                      "#7570B3", "#66A61E")
make_region_palette <- function(regions, shared = if (exists("region_palette")) region_palette else NULL) {
  regions <- sort(unique(as.character(regions)))
  pal <- grDevices::colorRampPalette(region_base_cols)(length(regions))
  names(pal) <- regions
  if (!is.null(shared)) { hit <- intersect(regions, names(shared)); pal[hit] <- shared[hit] }
  pal
}

# PNG (300 dpi) + PDF (cairo) side by side.
save_png_pdf <- function(plot, stem, width = 10, height = 7, dpi = 300) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  ggsave(paste0(stem, ".png"), plot, width = width, height = height, units = "in", dpi = dpi, bg = "white")
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, units = "in", device = cairo_pdf)
  invisible(stem)
}

# SMA and OLS fits of log10(y) on log10(x) in the standardized scaling_fits schema.
# test_slope: null slope for the SMA slope test (smatr::sma slope.test).
fit_sma_ols <- function(df, x, y, species, dataset, response, predictor, subset, note = "",
                        test_slope = 1) {
  d <- df[is.finite(df[[x]]) & is.finite(df[[y]]) & df[[x]] > 0 & df[[y]] > 0, , drop = FALSE]
  n <- nrow(d)
  lx <- log10(d[[x]]); ly <- log10(d[[y]])
  base <- data.frame(species = species, dataset = dataset, response = response, predictor = predictor,
                     subset = subset, stringsAsFactors = FALSE)
  if (n < 3) {
    return(rbind(
      cbind(base, method = "SMA", slope = NA, ci_low = NA, ci_high = NA, intercept = NA, r2 = NA, p = NA, n = n, note = paste(note, "n < 3")),
      cbind(base, method = "OLS", slope = NA, ci_low = NA, ci_high = NA, intercept = NA, r2 = NA, p = NA, n = n, note = paste(note, "n < 3"))))
  }
  sm <- smatr::sma(ly ~ lx, slope.test = test_slope)
  sc <- sm$coef[[1]]
  sma_row <- cbind(base, method = "SMA", slope = sc["slope", "coef(SMA)"],
                   ci_low = sc["slope", "lower limit"], ci_high = sc["slope", "upper limit"],
                   intercept = sc["elevation", "coef(SMA)"], r2 = sm$r2[[1]], p = sm$pval[[1]],
                   n = n, note = paste0(note, " SMA p = correlation test; slope test vs ", test_slope,
                                        ": p = ", signif(sm$slopetest[[1]]$p, 3)))
  m <- lm(ly ~ lx); ci <- confint(m); s <- summary(m)
  ols_row <- cbind(base, method = "OLS", slope = coef(m)[[2]], ci_low = ci[2, 1], ci_high = ci[2, 2],
                   intercept = coef(m)[[1]], r2 = s$r.squared, p = s$coefficients[2, 4], n = n,
                   note = paste0(note, " OLS p = slope != 0; slope vs ", test_slope, ": t-test p = ",
                                 signif(2 * pt(-abs((coef(m)[[2]] - test_slope) / s$coefficients[2, 2]), df = n - 2), 3)))
  out <- rbind(sma_row, ols_row); rownames(out) <- NULL
  out
}
