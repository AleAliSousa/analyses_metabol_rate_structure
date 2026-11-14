# Region-level OLS, rank correlation and leave-one-region-out checks.
# Inference is across region means, not subjects. P values are exploratory.
primary_region_rows <- function(data, map_path = NULL) {
  keep <- rep(TRUE, nrow(data))
  if (!is.null(map_path)) {
    m <- read_study_anatomy_map(map_path)
    ids <- unique(m$analysis_anatomy_id[m$primary_included])
    keep <- keep & data$anatomy_id %in% ids
  }
  if ("primary_included" %in% names(data))
    keep <- keep & !is.na(data$primary_included) & data$primary_included
  for (v in intersect(c("rcmr_value", "rcmrglc_mean"), names(data)))
    keep <- keep & is.finite(data[[v]])
  data[keep, , drop = FALSE]
}

regional_fit <- function(data, x, y, id = "anatomy_id") {
  stopifnot(all(c(x, y, id) %in% names(data)))
  d <- data.frame(x = data[[x]], y = data[[y]], id = data[[id]])
  d <- d[is.finite(d$x) & is.finite(d$y), , drop = FALSE]
  if (anyNA(d$id) || anyDuplicated(d$id)) stop("Expected one row per region")
  out <- data.frame(n = nrow(d), slope = NA_real_, slope_ci_low = NA_real_,
                    slope_ci_high = NA_real_, r_squared = NA_real_,
                    p_ols = NA_real_, rho = NA_real_, p_spearman = NA_real_,
                    loo_max_p = NA_real_, loo_max_p_dropped_region = NA_character_,
                    loo_slope_min = NA_real_, loo_slope_max = NA_real_,
                    status = "insufficient variation or fewer than four regions")
  if (nrow(d) < 4L || length(unique(d$x)) < 2L ||
      length(unique(d$y)) < 2L) return(out)
  fit <- lm(y ~ x, data = d)
  cf <- summary(fit)$coefficients
  ci <- confint(fit, "x")
  rank <- suppressWarnings(cor.test(d$x, d$y, method = "spearman", exact = FALSE))
  loo <- lapply(seq_len(nrow(d)), function(i) {
    z <- d[-i, , drop = FALSE]
    if (length(unique(z$x)) < 2L || length(unique(z$y)) < 2L)
      return(c(slope = NA_real_, p = NA_real_))
    m <- lm(y ~ x, data = z)
    c(slope = unname(coef(m)["x"]), p = summary(m)$coefficients["x", "Pr(>|t|)"])
  })
  loo <- do.call(rbind, loo)
  out$slope <- cf["x", "Estimate"]
  out$slope_ci_low <- ci[1]; out$slope_ci_high <- ci[2]
  out$r_squared <- summary(fit)$r.squared
  out$p_ols <- cf["x", "Pr(>|t|)"]
  out$rho <- unname(rank$estimate); out$p_spearman <- rank$p.value
  if (any(is.finite(loo[, "p"]))) {
    worst <- which.max(replace(loo[, "p"], !is.finite(loo[, "p"]), -Inf))
    out$loo_max_p <- loo[worst, "p"]
    out$loo_max_p_dropped_region <- d$id[worst]
  }
  if (any(is.finite(loo[, "slope"]))) {
    out$loo_slope_min <- min(loo[, "slope"], na.rm = TRUE)
    out$loo_slope_max <- max(loo[, "slope"], na.rm = TRUE)
  }
  out$status <- "ok"
  out
}

write_analysis_csv <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(data, path, row.names = FALSE, na = "")
}
