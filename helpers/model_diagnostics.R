# =====================================================================
# helpers/model_diagnostics.R  --  regression diagnostics for lm() fits
# ---------------------------------------------------------------------
# Usage (after the repo-root setwd at the top of a script):
#     source("helpers/model_diagnostics.R"); diag_stage("s1a")
#     fit <- lm(y ~ x, data = d)
#     diag_lm(fit, "short label", labels = d$anatomy_group)   # one figure
#     diag_facets(df, "facet_col", y = "rcmr_value", x = "prop", tag = "counts")
#
# Each call writes a 6-panel figure to figs/<stage>/diagnostics/<label>.jpg
# AND .pdf, and upserts one row into tables/<stage>/diagnostics_summary.csv.
# Base graphics + stats only (no extra packages). Never stops the host
# script: failures become a message.
#
# Panels: (1) residuals vs fitted  (2) normal Q-Q of standardised residuals
# (3) scale-location  (4) Cook's distance (4/n guide)  (5) residuals vs
# leverage (2p/n guide, Cook's D contours)  (6) observed vs fitted.
# Tests in the header/summary: Shapiro-Wilk (normality of residuals),
# Breusch-Pagan, studentised/Koenker form (heteroscedasticity).
# NOTE: with the small n here (often 5-26 regions) the tests have little
# power -- read the plots first.
# =====================================================================
.diag_env <- new.env()
.diag_env$stage <- "misc"

diag_stage <- function(stage) { .diag_env$stage <- stage; invisible(stage) }

.diag_slug <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x); x <- gsub("^_+|_+$", "", x)
  substr(x, 1, 90)
}

.diag_bp <- function(fit) {            # Koenker (studentised) Breusch-Pagan
  e2 <- residuals(fit)^2
  aux <- stats::lm(e2 ~ fitted(fit))
  n <- length(e2); stat <- n * summary(aux)$r.squared
  stats::pchisq(stat, df = 1, lower.tail = FALSE)
}

diag_lm <- function(fit, label, stage = .diag_env$stage, labels = NULL,
                    dir = file.path("figs", stage, "diagnostics")) {
  tryCatch({
    n <- length(residuals(fit)); p <- length(coef(fit)); p <- sum(!is.na(coef(fit)))
    if (n - p < 1) { message("[diag] skip '", label, "': saturated (n=", n, ", p=", p, ")"); return(invisible(NULL)) }
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    base <- file.path(dir, .diag_slug(label))
    e <- residuals(fit); f <- fitted(fit)
    h <- hatvalues(fit)
    sres <- tryCatch(rstandard(fit), error = function(.) rep(NA_real_, n))
    stud <- tryCatch(rstudent(fit),  error = function(.) rep(NA_real_, n))
    cd <- tryCatch(cooks.distance(fit), error = function(.) rep(NA_real_, n))
    if (is.null(labels)) labels <- names(e) else labels <- as.character(labels)
    if (length(labels) != n) labels <- names(e)
    if (is.null(labels)) labels <- seq_len(n)
    sh <- if (n >= 3 && n <= 5000 && diff(range(e)) > 0) tryCatch(shapiro.test(e)$p.value, error = function(.) NA_real_) else NA_real_
    bp <- if (n - p >= 2) tryCatch(.diag_bp(fit), error = function(.) NA_real_) else NA_real_
    sm <- summary(fit)
    fmtp <- function(x) ifelse(is.na(x), "NA", ifelse(x < 0.001, "<0.001", sprintf("%.3f", x)))
    hdr <- sprintf("%s\n%s   n=%d, params=%d, R2=%.3f, adj.R2=%.3f | Shapiro p=%s, Breusch-Pagan p=%s",
                   label, paste(deparse(formula(fit)), collapse = ""), n, p,
                   sm$r.squared, sm$adj.r.squared, fmtp(sh), fmtp(bp))
    draw <- function() {
      op <- par(mfrow = c(2, 3), mar = c(4.2, 4.2, 2.4, 1), oma = c(0, 0, 4, 0), mgp = c(2.4, 0.7, 0), cex.main = 0.95)
      on.exit(par(op))
      lab_pts <- function(x, y, idx) if (length(idx)) text(x[idx], y[idx], labels[idx], pos = 4, cex = 0.6, xpd = NA)
      top <- function(v, k = 3) { v <- abs(v); v[!is.finite(v)] <- -Inf; head(order(v, decreasing = TRUE), min(k, n)) }
      # 1 resid vs fitted
      plot(f, e, pch = 19, col = "#00000099", xlab = "Fitted", ylab = "Residuals", main = "Residuals vs fitted")
      abline(h = 0, lty = 2, col = "grey50")
      if (n >= 5) lines(stats::lowess(f, e), col = "firebrick", lwd = 1.5)
      lab_pts(f, e, top(e))
      # 2 QQ
      ok <- is.finite(sres)
      if (sum(ok) >= 3) {
        qq <- qqnorm(sres[ok], plot.it = FALSE)
        plot(qq$x, qq$y, pch = 19, col = "#00000099", xlab = "Theoretical quantiles", ylab = "Standardised residuals", main = "Normal Q-Q")
        qqline(sres[ok], lty = 2, col = "grey50")
        lab_pts(qq$x, qq$y, top(qq$y))
      } else { plot.new(); title("Normal Q-Q (n/a)") }
      # 3 scale-location
      sl <- sqrt(abs(sres))
      plot(f, sl, pch = 19, col = "#00000099", xlab = "Fitted", ylab = expression(sqrt("|standardised residuals|")), main = "Scale-location")
      if (n >= 5 && all(is.finite(sl))) lines(stats::lowess(f, sl), col = "firebrick", lwd = 1.5)
      lab_pts(f, sl, top(sl))
      # 4 Cook's distance
      if (all(is.finite(cd))) {
        plot(seq_len(n), cd, type = "h", lwd = 3, col = "#00000099", xlab = "Observation index", ylab = "Cook's distance", main = "Cook's distance")
        abline(h = 4 / n, lty = 2, col = "firebrick")
        lab_pts(seq_len(n), cd, top(cd)[cd[top(cd)] > 4 / n])
      } else { plot.new(); title("Cook's distance (n/a)") }
      # 5 resid vs leverage
      plot(h, sres, pch = 19, col = "#00000099", xlab = "Leverage", ylab = "Standardised residuals", main = "Residuals vs leverage",
           xlim = c(0, max(h, 2 * p / n, na.rm = TRUE) * 1.05))
      abline(h = 0, lty = 2, col = "grey50"); abline(v = 2 * p / n, lty = 2, col = "firebrick")
      hh <- seq(0.001, max(h, na.rm = TRUE) * 1.05, length.out = 100)
      for (cdv in c(0.5, 1)) {
        r <- sqrt(cdv * p * (1 - hh) / hh); lines(hh, r, lty = 3, col = "orange"); lines(hh, -r, lty = 3, col = "orange")
      }
      lab_pts(h, sres, top(cd))
      # 6 observed vs fitted
      yobs <- f + e
      plot(f, yobs, pch = 19, col = "#00000099", xlab = "Fitted", ylab = "Observed", main = "Observed vs fitted")
      abline(0, 1, lty = 2, col = "grey50")
      mtext(hdr, outer = TRUE, line = 0.5, cex = 0.8, font = 2)
    }
    jpeg(paste0(base, ".jpg"), width = 11, height = 7.5, units = "in", res = 200, quality = 92, bg = "white"); draw(); dev.off()
    pdf(paste0(base, ".pdf"), width = 11, height = 7.5); draw(); dev.off()

    row <- data.frame(label = label, formula = paste(deparse(formula(fit)), collapse = ""),
      n = n, n_params = p, df_resid = n - p, r_squared = sm$r.squared, adj_r_squared = sm$adj.r.squared,
      resid_sigma = sm$sigma, shapiro_p = sh, breusch_pagan_p = bp,
      max_cooks_d = suppressWarnings(max(cd, na.rm = TRUE)), n_cooks_gt_4_over_n = sum(cd > 4 / n, na.rm = TRUE),
      max_leverage = max(h), n_leverage_gt_2p_over_n = sum(h > 2 * p / n),
      max_abs_stud_resid = suppressWarnings(max(abs(stud), na.rm = TRUE)),
      most_influential = labels[which.max(replace(cd, !is.finite(cd), -Inf))],
      figure = paste0(base, ".jpg"), stringsAsFactors = FALSE)
    tdir <- file.path("tables", stage); dir.create(tdir, recursive = TRUE, showWarnings = FALSE)
    sp <- file.path(tdir, "diagnostics_summary.csv")
    if (file.exists(sp)) {
      old <- read.csv(sp, stringsAsFactors = FALSE)
      old <- old[old$label != label, , drop = FALSE]
      row <- rbind(old[, names(row), drop = FALSE], row)
    }
    write.csv(row, sp, row.names = FALSE)
    invisible(row)
  }, error = function(err) { message("[diag] failed for '", label, "': ", conditionMessage(err)); invisible(NULL) })
}

# One simple lm (y ~ x) per level of a facet column (mirrors the nest()/map()
# fits behind the faceted scatterplots).
diag_facets <- function(df, facet_col, y = "rcmr_value", x = "prop", tag = "", label_col = "anatomy_group",
                        stage = .diag_env$stage) {
  for (lv in unique(df[[facet_col]])) {
    d <- df[df[[facet_col]] == lv & complete.cases(df[[x]], df[[y]]), , drop = FALSE]
    if (nrow(d) < 4) next
    fit <- stats::lm(reformulate(x, y), data = d)
    diag_lm(fit, paste(c(tag, lv), collapse = " - "), stage = stage,
            labels = if (label_col %in% names(d)) d[[label_col]] else NULL)
  }
  invisible(NULL)
}
