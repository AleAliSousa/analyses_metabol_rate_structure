# ============================================================
# Network residual autocorrelation analysis for regional human
# volume prediction errors and rCMRGlc
# ============================================================
# Purpose
# -------
# This script asks whether regions that are anatomically or functionally linked
# have more similar human prediction residuals than expected after accounting for
# the metabolic model: residual volume deviation ~ rCMRGlc + rCMRGlc^2.
#
# The analysis is deliberately low parameter because the present dataset contains
# few regions. The primary test is therefore not a sensory-system fixed effect.
# It is a permutation test for residual network autocorrelation using an a priori
# region-region adjacency matrix W.
#
# Expected input
# --------------
# Option A, recommended: run this script after your existing s3 prediction script,
# so that the object `core_with_rCMRGlc` already exists in the R session.
#
# Option B, standalone: place the exported Excel table
#   "Table 4 Predicted volumes for Brownian and Pagel's lambda models.xlsx"
# in the working directory, or keep it at tables/s3/all/ within your project.
#
# Main outputs
# ------------
# network_residual_autocorrelation/
#   analysis_df_Brownian_BM_log_ratio.csv
#   primary_edges_used.csv
#   moran_permutation_results.csv
#   edge_residual_products.csv
#   model_comparison_metabolic_vs_network_eigenvector.csv
#   leave_one_region_moran_sensitivity.csv
#   baseline_fit_network_residuals.pdf/png
#   edge_residual_products.pdf/png
#   adjacency_matrix_heatmap.pdf/png
#
# Statistical notes
# -----------------
# 1. The recommended response is log(Observed / Predicted), which is symmetric
#    around zero and avoids the lower tail problem of (Observed - Predicted) / Observed.
# 2. Moran's I is computed on residuals from y ~ x + x^2, where x is centered rCMRGlc.
# 3. The adjacency matrix W is a hypothesis matrix. It must be specified before
#    looking at residual signs if the p-values are to be interpreted inferentially.
# 4. The network eigenvector model is a sensitivity analysis, not a definitive
#    causal model.

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
})

# ============================================================
# User settings
# ============================================================

model_use <- "Brownian (BM)"     # Alternatives: "Pagel's lambda (ML)"
response_metric <- "log_ratio"   # Alternatives: "diff_pre"
exclude_structures <- c("Neocortex white")

B_perm <- 9999                    # Main permutation count
B_perm_loo <- 1999                # Smaller count for leave-one-region sensitivity
set_seed <- 20260530

out_dir <- "network_residual_autocorrelation"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Candidate Excel-table locations for standalone use.
predicted_table_candidates <- c(
  file.path("tables", "s3", "all", "Table 4 Predicted volumes for Brownian and Pagel's lambda models.xlsx"),
  file.path("Table 4 Predicted volumes for Brownian and Pagel's lambda models.xlsx"),
  file.path("/mnt/data", "Table 4 Predicted volumes for Brownian and Pagel's lambda models.xlsx")
)

# ============================================================
# Small utilities
# ============================================================

as_num <- function(x) {
  readr::parse_number(as.character(x))
}

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  gsub("_+$", "", x)
}

AICc_lm <- function(fit) {
  n <- stats::nobs(fit)
  k <- attr(stats::logLik(fit), "df")  # includes residual variance
  if (!is.finite(n) || !is.finite(k) || n <= k + 1) return(NA_real_)
  stats::AIC(fit) + (2 * k * (k + 1)) / (n - k - 1)
}

loocv_rmse <- function(formula, data) {
  n <- nrow(data)
  pred <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    fit_i <- try(stats::lm(formula, data = data[-i, , drop = FALSE]), silent = TRUE)
    if (!inherits(fit_i, "try-error")) {
      pred[i] <- tryCatch(
        as.numeric(stats::predict(fit_i, newdata = data[i, , drop = FALSE])),
        error = function(e) NA_real_
      )
    }
  }
  sqrt(mean((data$y - pred)^2, na.rm = TRUE))
}

# ============================================================
# Read or construct the analysis table
# ============================================================

read_predicted_table_from_xlsx <- function(path, model = c("Brownian (BM)", "Pagel's lambda (ML)")) {
  model <- match.arg(model)
  raw <- readxl::read_xlsx(path, col_names = FALSE)

  if (nrow(raw) < 3 || ncol(raw) < 19) {
    stop("The Excel table does not have the expected Table 4 layout.")
  }

  rows <- 3:nrow(raw)

  if (model == "Brownian (BM)") {
    cols <- list(lambda = 2, Lower = 3, Predicted = 4, Upper = 5,
                 Observed = 6, Diff.min = 7, Diff.pre = 8, Diff.max = 9)
  } else {
    cols <- list(lambda = 10, Lower = 11, Predicted = 12, Upper = 13,
                 Observed = 14, Diff.min = 15, Diff.pre = 16, Diff.max = 17)
  }

  tibble(
    Structure = as.character(raw[[1]][rows]),
    Model = model,
    lambda = as_num(raw[[cols$lambda]][rows]),
    Lower = as_num(raw[[cols$Lower]][rows]),
    Predicted = as_num(raw[[cols$Predicted]][rows]),
    Upper = as_num(raw[[cols$Upper]][rows]),
    Observed = as_num(raw[[cols$Observed]][rows]),
    Diff.min = as_num(raw[[cols$Diff.min]][rows]),
    Diff.pre = as_num(raw[[cols$Diff.pre]][rows]),
    Diff.max = as_num(raw[[cols$Diff.max]][rows]),
    N = as.integer(as_num(raw[[18]][rows])),
    rCMRGlc = as_num(raw[[19]][rows])
  ) %>%
    filter(!is.na(Structure), Structure != "")
}

construct_analysis_df <- function(model_use, response_metric, exclude_structures) {
  if (exists("core_with_rCMRGlc", inherits = TRUE)) {
    message("Using object `core_with_rCMRGlc` from the current R session.")
    dat0 <- get("core_with_rCMRGlc", inherits = TRUE) %>%
      filter(.data$Model == model_use) %>%
      transmute(
        Structure = as.character(.data$Structure),
        Model = as.character(.data$Model),
        lambda = as_num(.data$lambda),
        Lower = as_num(.data$`95% CI min`),
        Predicted = as_num(.data$Predicted),
        Upper = as_num(.data$`95% CI max`),
        Observed = as_num(.data$Observed),
        Diff.min = as_num(.data$Diff.min),
        Diff.pre = as_num(.data$Diff.pre),
        Diff.max = as_num(.data$Diff.max),
        N = as.integer(as_num(.data$N)),
        rCMRGlc = as_num(.data$rCMRGlc)
      )
  } else {
    predicted_table_xlsx <- predicted_table_candidates[file.exists(predicted_table_candidates)][1]
    if (is.na(predicted_table_xlsx)) {
      stop(
        "Could not find `core_with_rCMRGlc` or the exported Excel table.\n",
        "Either run this after your s3 prediction script, or place Table 4 in the working directory."
      )
    }
    message("Reading standalone Excel table: ", predicted_table_xlsx)
    dat0 <- read_predicted_table_from_xlsx(predicted_table_xlsx, model = model_use)
  }

  dat <- dat0 %>%
    mutate(
      excluded = .data$Structure %in% exclude_structures,
      log_ratio = log(.data$Observed / .data$Predicted),
      diff_pre = .data$Diff.pre
    ) %>%
    filter(!.data$excluded) %>%
    drop_na(rCMRGlc, Observed, Predicted)

  if (response_metric == "log_ratio") {
    dat <- dat %>% mutate(y = .data$log_ratio)
  } else if (response_metric == "diff_pre") {
    dat <- dat %>% mutate(y = .data$diff_pre)
  } else {
    stop("response_metric must be either 'log_ratio' or 'diff_pre'.")
  }

  dat %>%
    mutate(
      x_raw = .data$rCMRGlc,
      x = .data$rCMRGlc - mean(.data$rCMRGlc, na.rm = TRUE),
      Structure = as.character(.data$Structure)
    ) %>%
    arrange(.data$Structure)
}

analysis_df <- construct_analysis_df(model_use, response_metric, exclude_structures)

if (nrow(analysis_df) < 8) {
  warning("Very few regions remain after filtering. Treat all network tests as descriptive.")
}

readr::write_csv(
  analysis_df,
  file.path(out_dir, paste0("analysis_df_", safe_name(model_use), "_", response_metric, ".csv"))
)

# ============================================================
# Define an a priori anatomical or functional network
# ============================================================
# Edges are undirected. The edge list is intentionally sparse and conservative.
# You can edit this table to reflect a stricter anatomical hypothesis.
#
# Important: do not tune this network to maximize Moran's I. Use sensitivity
# analyses below to check whether any one module or edge drives the result.

edges_primary <- tribble(
  ~from, ~to, ~module, ~weight, ~rationale,
  "Corpus geniculatum laterale", "Area striata grey", "visual", 1,
    "LGN to primary visual cortex thalamocortical visual pathway",

  "Striatum", "Pallidum", "basal_ganglia", 1,
    "Canonical striatopallidal basal ganglia linkage",
  "Striatum", "Nucleus subthalamicus Luysi", "basal_ganglia", 1,
    "Basal ganglia indirect and hyperdirect circuit membership",
  "Pallidum", "Nucleus subthalamicus Luysi", "basal_ganglia", 1,
    "Reciprocal basal ganglia circuit coupling",
  "Capsula interna", "Striatum", "basal_ganglia_white_matter", 1,
    "Internal capsule courses adjacent to and through basal ganglia territory",
  "Capsula interna", "Pallidum", "basal_ganglia_white_matter", 1,
    "Projection fiber system coupled to basal ganglia territory",

  "Amygdala", "Hippocampus", "limbic", 1,
    "Medial temporal limbic system coupling",
  "Amygdala", "Insular cortex (grey)", "limbic_interoceptive", 1,
    "Salience, interoceptive, and affective coupling",
  "Hippocampus", "Insular cortex (grey)", "limbic_interoceptive", 1,
    "Broader limbic and interoceptive association",

  "Mesencephalon", "Cerebellum", "motor_hindbrain", 1,
    "Brainstem and cerebellar motor network coupling",
  "Mesencephalon", "Nucleus subthalamicus Luysi", "basal_ganglia_midbrain", 1,
    "Midbrain and basal ganglia motor circuit coupling",
  "Mesencephalon", "Pallidum", "basal_ganglia_midbrain", 1,
    "Midbrain and pallidal motor circuit coupling"
)

# Optional edge sets for sensitivity tests.
edge_sets <- list(
  primary = edges_primary,
  no_visual_LGN_V1 = edges_primary %>% filter(.data$module != "visual"),
  visual_only = edges_primary %>% filter(.data$module == "visual"),
  basal_ganglia_only = edges_primary %>% filter(grepl("basal_ganglia", .data$module)),
  limbic_only = edges_primary %>% filter(grepl("limbic", .data$module)),
  motor_hindbrain_only = edges_primary %>% filter(.data$module == "motor_hindbrain")
)

make_W <- function(structures, edges) {
  structures <- as.character(structures)
  n <- length(structures)
  W <- matrix(0, nrow = n, ncol = n, dimnames = list(structures, structures))

  edges_used <- edges %>%
    filter(.data$from %in% structures, .data$to %in% structures) %>%
    mutate(weight = ifelse(is.na(.data$weight), 1, .data$weight))

  if (nrow(edges_used) > 0) {
    for (i in seq_len(nrow(edges_used))) {
      a <- edges_used$from[i]
      b <- edges_used$to[i]
      w <- edges_used$weight[i]
      W[a, b] <- W[a, b] + w
      W[b, a] <- W[b, a] + w
    }
  }

  diag(W) <- 0

  dropped_edges <- edges %>%
    filter(!(.data$from %in% structures & .data$to %in% structures))

  list(W = W, edges_used = edges_used, dropped_edges = dropped_edges)
}

W_info <- make_W(analysis_df$Structure, edges_primary)
W <- W_info$W
edges_used <- W_info$edges_used

readr::write_csv(edges_used, file.path(out_dir, "primary_edges_used.csv"))
readr::write_csv(W_info$dropped_edges, file.path(out_dir, "primary_edges_dropped_because_region_missing.csv"))

message("Regions in analysis: ", nrow(analysis_df))
message("Primary network edges used: ", nrow(edges_used))
message("Primary network total symmetric weight S0: ", sum(W))

# ============================================================
# Moran permutation test on baseline metabolic residuals
# ============================================================

moran_I <- function(z, W) {
  z <- as.numeric(z)
  if (!is.matrix(W)) stop("W must be a matrix.")
  if (nrow(W) != length(z) || ncol(W) != length(z)) {
    stop("Dimensions of W do not match the residual vector.")
  }
  S0 <- sum(W)
  if (!is.finite(S0) || S0 <= 0) stop("W has no edges or positive weights.")

  z <- z - mean(z, na.rm = TRUE)
  denom <- sum(z^2, na.rm = TRUE)
  if (!is.finite(denom) || denom <= 0) return(NA_real_)

  as.numeric((length(z) / S0) * (t(z) %*% W %*% z) / denom)
}

moran_permutation <- function(z, W, B = 9999, seed = 1) {
  set.seed(seed)
  I_obs <- moran_I(z, W)
  I_perm <- replicate(B, moran_I(sample(z, length(z), replace = FALSE), W))

  # Null expectation for Moran's I is often close to -1/(n-1), not exactly zero.
  # The centered two-sided p-value compares deviations from the permutation mean.
  p_greater <- (sum(I_perm >= I_obs, na.rm = TRUE) + 1) / (B + 1)
  p_less <- (sum(I_perm <= I_obs, na.rm = TRUE) + 1) / (B + 1)
  p_two_centered <- (
    sum(abs(I_perm - mean(I_perm, na.rm = TRUE)) >=
          abs(I_obs - mean(I_perm, na.rm = TRUE)), na.rm = TRUE) + 1
  ) / (B + 1)

  list(
    I_obs = I_obs,
    I_perm = I_perm,
    expected_perm = mean(I_perm, na.rm = TRUE),
    sd_perm = stats::sd(I_perm, na.rm = TRUE),
    p_greater = p_greater,
    p_less = p_less,
    p_two_centered = p_two_centered
  )
}

# Baseline metabolic model: residual deviation as quadratic function of centered rCMRGlc.
fit_metabolic <- stats::lm(y ~ x + I(x^2), data = analysis_df)
analysis_df$baseline_fitted <- stats::fitted(fit_metabolic)
analysis_df$baseline_resid <- stats::resid(fit_metabolic)
analysis_df$baseline_std_resid <- as.numeric(scale(analysis_df$baseline_resid))

# Test the primary network.
moran_primary <- moran_permutation(analysis_df$baseline_resid, W, B = B_perm, seed = set_seed)

# Repeat across edge-set sensitivity definitions.
run_moran_for_edge_set <- function(edge_name, edge_tbl, dat, B, seed) {
  wi <- make_W(dat$Structure, edge_tbl)
  if (sum(wi$W) <= 0) {
    return(tibble(
      edge_set = edge_name,
      n_regions = nrow(dat),
      n_edges = 0,
      I_obs = NA_real_,
      expected_perm = NA_real_,
      sd_perm = NA_real_,
      p_greater = NA_real_,
      p_less = NA_real_,
      p_two_centered = NA_real_
    ))
  }

  mr <- moran_permutation(dat$baseline_resid, wi$W, B = B, seed = seed)
  tibble(
    edge_set = edge_name,
    n_regions = nrow(dat),
    n_edges = nrow(wi$edges_used),
    I_obs = mr$I_obs,
    expected_perm = mr$expected_perm,
    sd_perm = mr$sd_perm,
    p_greater = mr$p_greater,
    p_less = mr$p_less,
    p_two_centered = mr$p_two_centered
  )
}

moran_results <- purrr::imap_dfr(
  edge_sets,
  ~ run_moran_for_edge_set(.y, .x, analysis_df, B = B_perm, seed = set_seed + match(.y, names(edge_sets)))
)

readr::write_csv(moran_results, file.path(out_dir, "moran_permutation_results.csv"))

# ============================================================
# Edge-level residual products and concordance
# ============================================================

edge_residual_products <- edges_used %>%
  mutate(
    resid_from = analysis_df$baseline_resid[match(.data$from, analysis_df$Structure)],
    resid_to = analysis_df$baseline_resid[match(.data$to, analysis_df$Structure)],
    z_from = analysis_df$baseline_std_resid[match(.data$from, analysis_df$Structure)],
    z_to = analysis_df$baseline_std_resid[match(.data$to, analysis_df$Structure)],
    residual_product = .data$z_from * .data$z_to,
    same_sign = sign(.data$z_from) == sign(.data$z_to),
    edge = paste(.data$from, .data$to, sep = " -- ")
  ) %>%
  arrange(desc(.data$residual_product))

readr::write_csv(edge_residual_products, file.path(out_dir, "edge_residual_products.csv"))

edge_product_stat <- function(z, edges, structures) {
  z <- as.numeric(scale(z))
  names(z) <- structures
  mean(z[edges$from] * z[edges$to], na.rm = TRUE)
}

set.seed(set_seed)
edge_product_obs <- edge_product_stat(analysis_df$baseline_resid, edges_used, analysis_df$Structure)
edge_product_perm <- replicate(
  B_perm,
  edge_product_stat(sample(analysis_df$baseline_resid), edges_used, analysis_df$Structure)
)
edge_product_summary <- tibble(
  statistic = "mean standardized residual product across primary edges",
  observed = edge_product_obs,
  expected_perm = mean(edge_product_perm, na.rm = TRUE),
  sd_perm = sd(edge_product_perm, na.rm = TRUE),
  p_greater = (sum(edge_product_perm >= edge_product_obs, na.rm = TRUE) + 1) / (B_perm + 1),
  p_less = (sum(edge_product_perm <= edge_product_obs, na.rm = TRUE) + 1) / (B_perm + 1),
  p_two_centered = (
    sum(abs(edge_product_perm - mean(edge_product_perm, na.rm = TRUE)) >=
          abs(edge_product_obs - mean(edge_product_perm, na.rm = TRUE)), na.rm = TRUE) + 1
  ) / (B_perm + 1)
)
readr::write_csv(edge_product_summary, file.path(out_dir, "edge_product_permutation_summary.csv"))

# ============================================================
# Network eigenvector sensitivity model
# ============================================================
# This asks whether one low-dimensional network-structured axis improves the
# metabolic model. It is analogous to Moran eigenvector maps, but kept deliberately
# simple for the present small-n dataset.

network_eigenvectors <- function(W) {
  n <- nrow(W)
  H <- diag(n) - matrix(1 / n, nrow = n, ncol = n)
  M <- H %*% W %*% H
  ev <- eigen(M, symmetric = TRUE)
  keep <- which(ev$values > sqrt(.Machine$double.eps))
  if (length(keep) == 0) {
    out <- matrix(nrow = n, ncol = 0)
    rownames(out) <- rownames(W)
    return(out)
  }
  out <- ev$vectors[, keep, drop = FALSE]
  rownames(out) <- rownames(W)
  colnames(out) <- paste0("net", seq_len(ncol(out)))
  attr(out, "eigenvalues") <- ev$values[keep]
  out
}

mem <- network_eigenvectors(W)

model_comparison <- tibble()
fit_net1 <- NULL

model_row <- function(model_name, formula, dat) {
  fit <- stats::lm(formula, data = dat)
  s <- summary(fit)
  tibble(
    model = model_name,
    formula = paste(deparse(formula), collapse = " "),
    n = stats::nobs(fit),
    df_logLik = attr(stats::logLik(fit), "df"),
    adj_r2 = s$adj.r.squared,
    sigma = s$sigma,
    AIC = stats::AIC(fit),
    AICc = AICc_lm(fit),
    LOOCV_RMSE = loocv_rmse(formula, dat)
  )
}

model_comparison <- bind_rows(
  model_comparison,
  model_row("metabolic_quadratic", y ~ x + I(x^2), analysis_df)
)

if (ncol(mem) >= 1) {
  analysis_df$net1 <- mem[analysis_df$Structure, "net1"]
  fit_net1 <- stats::lm(y ~ x + I(x^2) + net1, data = analysis_df)
  model_comparison <- bind_rows(
    model_comparison,
    model_row("metabolic_quadratic_plus_network_eigenvector_1", y ~ x + I(x^2) + net1, analysis_df)
  )
}

# Diagnostic residual-lag association. This is not a proper predictive model,
# because neighbor_resid is computed from the baseline residuals. It is a compact
# descriptive version of the same question as Moran's I.
row_standardize <- function(W) {
  rs <- rowSums(W)
  Ws <- sweep(W, 1, rs, FUN = "/")
  Ws[!is.finite(Ws)] <- 0
  Ws
}

if (sum(W) > 0) {
  W_row <- row_standardize(W)
  analysis_df$neighbor_resid <- as.numeric(W_row %*% analysis_df$baseline_resid)
  fit_resid_lag <- stats::lm(baseline_resid ~ neighbor_resid, data = analysis_df)
  resid_lag_summary <- tibble(
    term = names(coef(fit_resid_lag)),
    estimate = unname(coef(fit_resid_lag)),
    p_value = coef(summary(fit_resid_lag))[, "Pr(>|t|)"]
  )
  readr::write_csv(resid_lag_summary, file.path(out_dir, "neighbor_residual_lag_diagnostic.csv"))
}

readr::write_csv(
  model_comparison,
  file.path(out_dir, "model_comparison_metabolic_vs_network_eigenvector.csv")
)

readr::write_csv(
  analysis_df,
  file.path(out_dir, paste0("analysis_df_with_residuals_", safe_name(model_use), "_", response_metric, ".csv"))
)

# ============================================================
# Leave-one-region sensitivity for the primary network test
# ============================================================

loo_moran <- purrr::map_dfr(seq_len(nrow(analysis_df)), function(i) {
  dat_i <- analysis_df[-i, , drop = FALSE]
  wi <- make_W(dat_i$Structure, edges_primary)

  if (sum(wi$W) <= 0 || nrow(dat_i) < 6) {
    return(tibble(
      removed_structure = analysis_df$Structure[i],
      n_regions = nrow(dat_i),
      n_edges = nrow(wi$edges_used),
      I_obs = NA_real_,
      expected_perm = NA_real_,
      p_greater = NA_real_,
      p_two_centered = NA_real_
    ))
  }

  fit_i <- stats::lm(y ~ x + I(x^2), data = dat_i)
  mr_i <- moran_permutation(stats::resid(fit_i), wi$W, B = B_perm_loo, seed = set_seed + i)

  tibble(
    removed_structure = analysis_df$Structure[i],
    n_regions = nrow(dat_i),
    n_edges = nrow(wi$edges_used),
    I_obs = mr_i$I_obs,
    expected_perm = mr_i$expected_perm,
    p_greater = mr_i$p_greater,
    p_two_centered = mr_i$p_two_centered
  )
})

readr::write_csv(loo_moran, file.path(out_dir, "leave_one_region_moran_sensitivity.csv"))

# ============================================================
# Plots
# ============================================================

plot_y_lab <- if (response_metric == "log_ratio") {
  "log(Observed / Predicted)"
} else {
  "Diff.pre = (Observed - Predicted) / Observed"
}

p_baseline <- ggplot(analysis_df, aes(x = rCMRGlc, y = y)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = TRUE,
              color = "steelblue4", fill = "steelblue", alpha = 0.20, linewidth = 1.1) +
  geom_point(size = 2.5) +
  geom_text(aes(label = Structure), hjust = -0.05, vjust = 0.4, size = 3, check_overlap = TRUE) +
  labs(
    title = "Baseline metabolic model for human volume residuals",
    subtitle = paste0("Model = ", model_use, "; response = ", response_metric,
                      "; neocortex white excluded"),
    x = "rCMRGlc (µmol/100 g/min)",
    y = plot_y_lab
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

p_edge <- ggplot(edge_residual_products, aes(x = reorder(edge, residual_product), y = residual_product, fill = same_sign)) +
  geom_hline(yintercept = 0, color = "grey50") +
  geom_col(width = 0.75) +
  coord_flip() +
  labs(
    title = "Residual concordance across a priori network edges",
    subtitle = "Positive values mean both connected regions deviate in the same direction after the metabolic model",
    x = NULL,
    y = "Product of standardized metabolic-model residuals",
    fill = "Same sign"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

W_long <- as.data.frame(as.table(W), stringsAsFactors = FALSE) %>%
  as_tibble() %>%
  rename(from = Var1, to = Var2, weight = Freq) %>%
  filter(.data$weight > 0)

p_W <- ggplot(W_long, aes(x = from, y = to, fill = weight)) +
  geom_tile() +
  coord_equal() +
  labs(
    title = "Primary a priori network adjacency matrix",
    x = NULL,
    y = NULL,
    fill = "Weight"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

ggsave(file.path(out_dir, "baseline_fit_network_residuals.pdf"), p_baseline, width = 7.5, height = 5)
ggsave(file.path(out_dir, "baseline_fit_network_residuals.png"), p_baseline, width = 7.5, height = 5, dpi = 300)

ggsave(file.path(out_dir, "edge_residual_products.pdf"), p_edge, width = 8, height = 5.5)
ggsave(file.path(out_dir, "edge_residual_products.png"), p_edge, width = 8, height = 5.5, dpi = 300)

ggsave(file.path(out_dir, "adjacency_matrix_heatmap.pdf"), p_W, width = 6.5, height = 5.5)
ggsave(file.path(out_dir, "adjacency_matrix_heatmap.png"), p_W, width = 6.5, height = 5.5, dpi = 300)

# ============================================================
# Console summary
# ============================================================

cat("\n================ Network residual autocorrelation summary ================\n")
cat("Model used: ", model_use, "\n", sep = "")
cat("Response metric: ", response_metric, "\n", sep = "")
cat("Regions analyzed: ", nrow(analysis_df), "\n", sep = "")
cat("Primary edges used: ", nrow(edges_used), "\n", sep = "")
cat("\nBaseline metabolic model:\n")
print(summary(fit_metabolic))

cat("\nPrimary Moran permutation test on baseline residuals:\n")
cat("Observed Moran's I: ", signif(moran_primary$I_obs, 4), "\n", sep = "")
cat("Permutation mean: ", signif(moran_primary$expected_perm, 4), "\n", sep = "")
cat("p greater, positive autocorrelation: ", signif(moran_primary$p_greater, 4), "\n", sep = "")
cat("p two-sided centered: ", signif(moran_primary$p_two_centered, 4), "\n", sep = "")

cat("\nEdge-set sensitivity results:\n")
print(moran_results)

cat("\nModel comparison:\n")
print(model_comparison)

cat("\nOutputs written to: ", normalizePath(out_dir), "\n", sep = "")
cat("========================================================================\n")
