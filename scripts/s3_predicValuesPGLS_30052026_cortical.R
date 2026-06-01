setwd("~/Library/CloudStorage/Dropbox/COLLABORATIVE/Do expensive brain regions increase less in humans/analyses_metabol_rate_structure")

library(here)
library(ape)
library(nlme)  
# Mac note: ML optimizer is different in Mac than PC so extra code was added
library(tidyverse)
library(dispRity)
library(scales)
library(readxl) 
library(writexl)

# --- SETUP ALL ---
# --- Read Table (for correlations)
heiss_stephan_tbl <- read.csv("data_intermediate/Heiss_Stephan_data.csv")

# --- Load tree and data
tr <- read.tree("data_raw/species.nwk")
Stephan_primates <- read.csv("data_raw/Stephan_primates.csv")

# --- Clean columns
data_clean <- subset(Stephan_primates, select = -c(X, order))

#  --- Create a unified “preferred brain volume” column
data <- data_clean %>%
  mutate(
    Preferred_brain_volume = coalesce(
      Brain_volume,           # 1st choice
      Brainvol,               # 2nd choice
      Total_brain_net_volume  # 3rd choice
    )
  )

#  --- Sanity check: did everyone get a value?
sum(is.na(data$Preferred_brain_volume))

# --- Parameters
options(scipen=999)

# --- Structures and display labels, in plotting/table order
region_labels <- c(
  "LGN_Sousa"             = "Corpus geniculatum laterale",
  "Amygdala"              = "Amygdala",
  "Pallidum"              = "Pallidum",
  "NeoW_Frahm"            = "Neocortex white",
  "Total_insula_volume_L" = "Insular cortex (grey)",
  "Nucleus_subthalamicus" = "Nucleus subthalamicus Luysi",
  "Capsula_interna"       = "Capsula interna",
  "Striatum"              = "Striatum",
  "ASG_Sousa"             = "Area striata grey",
  "NeoG_Frahm"            = "Neocortex grey",
  "Mesencephalon"         = "Mesencephalon",
  "Cerebellum"            = "Cerebellum",
  "Hippocampus"           = "Hippocampus"
)

target_cols <- names(region_labels)

# --- Cortical grouping schemes, in same raw-name space as region_labels ---
# Scheme 1: two-way split. Neocortex white is treated as cortical-associated.
region_groups_cortical2 <- c(
  "LGN_Sousa"             = "Non-cortical",
  "Amygdala"              = "Non-cortical",
  "Pallidum"              = "Non-cortical",
  "NeoW_Frahm"            = "Cortical",
  "Total_insula_volume_L" = "Cortical",
  "Nucleus_subthalamicus" = "Non-cortical",
  "Capsula_interna"       = "Non-cortical",
  "Striatum"              = "Non-cortical",
  "ASG_Sousa"             = "Cortical",
  "NeoG_Frahm"            = "Cortical",
  "Mesencephalon"         = "Non-cortical",
  "Cerebellum"            = "Non-cortical",
  "Hippocampus"           = "Cortical"
)

# Scheme 2: three-way split. White-matter structures are separated from cortical gray.
# This treats both Neocortex white and Capsula interna as white matter.
region_groups_cortical3 <- c(
  "LGN_Sousa"             = "Non-cortical",
  "Amygdala"              = "Non-cortical",
  "Pallidum"              = "Non-cortical",
  "NeoW_Frahm"            = "White matter",
  "Total_insula_volume_L" = "Cortical gray",
  "Nucleus_subthalamicus" = "Non-cortical",
  "Capsula_interna"       = "White matter",
  "Striatum"              = "Non-cortical",
  "ASG_Sousa"             = "Cortical gray",
  "NeoG_Frahm"            = "Cortical gray",
  "Mesencephalon"         = "Non-cortical",
  "Cerebellum"            = "Non-cortical",
  "Hippocampus"           = "Cortical gray"
)

# Defensive checks: every analyzed structure must have exactly one group in each scheme.
stopifnot(setequal(names(region_labels), names(region_groups_cortical2)))
stopifnot(setequal(names(region_labels), names(region_groups_cortical3)))

cortical_2_group_levels <- c("Cortical", "Non-cortical")
cortical_3_group_levels <- c("Cortical gray", "White matter", "Non-cortical")

# --- Robust helpers to convert raw variable names to display labels / anatomical groups ---
label_region <- function(x) {
  x_chr <- as.character(x)
  out <- unname(region_labels[x_chr])
  out[is.na(out)] <- x_chr[is.na(out)]
  out
}

group_region_cortical2 <- function(x) {
  x_chr <- as.character(x)
  out <- unname(region_groups_cortical2[x_chr])
  out[is.na(out)] <- "Unassigned"
  factor(out, levels = cortical_2_group_levels)
}

group_region_cortical3 <- function(x) {
  x_chr <- as.character(x)
  out <- unname(region_groups_cortical3[x_chr])
  out[is.na(out)] <- "Unassigned"
  factor(out, levels = cortical_3_group_levels)
}

# Backward-compatible alias used by the existing plotting code.
# The default RegionGroup is now the 2-way cortical/non-cortical split.
group_region <- group_region_cortical2

region_metadata <- tibble(
  Variable = target_cols,
  Structure = label_region(target_cols),
  CorticalGroup2 = group_region_cortical2(target_cols),
  CorticalGroup3 = group_region_cortical3(target_cols)
)

cortical_2_group_counts <- region_metadata %>%
  count(CorticalGroup2, name = "n_structures")

cortical_3_group_counts <- region_metadata %>%
  count(CorticalGroup3, name = "n_structures")

print(cortical_2_group_counts)
print(cortical_3_group_counts)

# --- Helper: transform VCV by lambda (off-diagonals scaled, diagonals fixed)
transform_vcv_lambda <- function(tree, lambda) {
  v <- vcv(tree, corr = TRUE)
  diag_vals <- diag(v)
  v <- v * lambda
  diag(v) <- diag_vals
  v
}

# --- Helper: prepare one structure (data + tree alignment + drop human)
prep_structure <- function(var_name,
                           data,
                           tr,
                           human = "Homo_sapiens",
                           predictor_col = "Preferred_brain_volume") {
  if (!var_name %in% names(data)) return(NULL)
  if (!predictor_col %in% names(data)) stop("predictor_col not found in data")

  d1 <- data[, c("Species", predictor_col, var_name)]
  d1 <- as.data.frame(na.omit(d1))

  # Part-whole correction: Rest-of-Brain (ROB) = Total - Part
  d1$Rest_of_Brain <- d1[[predictor_col]] - d1[[var_name]]
  d1 <- d1[d1$Rest_of_Brain > 0, ]
  rownames(d1) <- d1$Species

  if (!human %in% d1$Species) return(NULL)

  obs_val <- d1[d1$Species == human, var_name][1]
  xh_raw  <- d1[d1$Species == human, "Rest_of_Brain"][1]

  # Align to tree
  clean <- clean.data(d1, tr)
  tr1 <- clean$tree
  data2 <- clean$data
  if (!is.null(data2$Species)) rownames(data2) <- data2$Species

  # Fit excluding human (consistent with your original logic)
  tr2 <- drop.tip(tr1, human)
  data_model <- subset(data2, Species != human)

  fmla <- as.formula(paste0("log(", var_name, ") ~ log(Rest_of_Brain)"))

  list(
    var_name = var_name,
    human = human,
    d_full = d1,
    tr1 = tr1,
    tr2 = tr2,
    data2 = data2,
    data_model = data_model,
    obs_val = as.numeric(obs_val),
    xh_raw = as.numeric(xh_raw),
    fmla = fmla
  )
}

# --- Helper: safe GLS fit
fit_gls_safe <- function(fmla, cor_struct, data_model) {
  tryCatch(
    gls(
      fmla,
      correlation = cor_struct,
      data = data_model,
      control = glsControl(opt = "optim", msMaxIter = 1000, msTol = 1e-6)
    ),
    error = function(e) NULL
  )
}

######### INDIV SECTION PLOTS PT 1 START
# --- Helper: BM PGLS plot for each individual structure
plot_region_gls <- function(region, data, tr, human = "Homo_sapiens") {
  
  # Display label for plotting only. The raw region name is still used for modeling.
  pretty_region <- label_region(region)
  
  # Prep
  pp <- prep_structure(region, data = data, tr = tr, human = human)
  
  if (is.null(pp)) {
    warning("Skipping ", region, ": prep_structure() returned NULL")
    return(NULL)
  }
  
  # Phylogenetic GLS model: BM
  cor_struct <- corBrownian(1, form = ~Species, phy = pp$tr2)
  
  fit <- fit_gls_safe(pp$fmla, cor_struct, pp$data_model)
  
  if (is.null(fit)) {
    warning("Skipping ", region, ": fit_gls_safe() returned NULL")
    return(NULL)
  }
  
  # Training range excludes human
  x_min_train <- min(pp$data_model$Rest_of_Brain, na.rm = TRUE)
  x_max_train <- max(pp$data_model$Rest_of_Brain, na.rm = TRUE)
  
  # Plotting range includes human, so the line extends to the human x-value
  x_max_plot <- max(pp$d_full$Rest_of_Brain, na.rm = TRUE)
  
  # Log-spaced prediction grid for log-log plotting
  x_grid_raw <- exp(seq(
    from = log(x_min_train),
    to   = log(x_max_plot),
    length.out = 200
  ))
  
  # Model: log(y) = b0 + b1 * log(ROB) + error
  X <- cbind(1, log(x_grid_raw))
  colnames(X) <- names(coef(fit))
  
  beta <- coef(fit)
  Vb   <- vcov(fit)
  
  eta <- as.numeric(X %*% beta)
  se  <- sqrt(diag(X %*% Vb %*% t(X)))
  z   <- qnorm(0.975)
  
  pred_df <- tibble(
    Rest_of_Brain = x_grid_raw,
    fit = exp(eta),
    lwr = exp(eta - z * se),
    upr = exp(eta + z * se),
    extrapolated = Rest_of_Brain > x_max_train
  )
  
  plot_df <- pp$d_full %>%
    mutate(
      is_human = Species == pp$human
    )
  
  ggplot() +
    geom_ribbon(
      data = filter(pred_df, !extrapolated),
      aes(x = Rest_of_Brain, ymin = lwr, ymax = upr),
      alpha = 0.2
    ) +
    geom_ribbon(
      data = filter(pred_df, extrapolated),
      aes(x = Rest_of_Brain, ymin = lwr, ymax = upr),
      alpha = 0.1
    ) +
    geom_line(
      data = filter(pred_df, !extrapolated),
      aes(x = Rest_of_Brain, y = fit),
      linewidth = 1
    ) +
    geom_line(
      data = filter(pred_df, extrapolated),
      aes(x = Rest_of_Brain, y = fit),
      linewidth = 1,
      linetype = "dashed"
    ) +
    geom_point(
      data = filter(plot_df, !is_human),
      aes(x = Rest_of_Brain, y = .data[[region]]),
      size = 2,
      alpha = 0.8
    ) +
    geom_point(
      data = filter(plot_df, is_human),
      aes(x = Rest_of_Brain, y = .data[[region]]),
      size = 3,
      shape = 17
    ) +
    scale_x_log10() +
    scale_y_log10() +
    theme_bw() +
    labs(
      title = paste0("Phylogenetic GLS (BM): ", pretty_region, " vs Rest of Brain"),
      subtitle = "Solid = non-human range; dashed = extrapolation to human. Human shown as triangle; fit excludes human.",
      x = "Rest of Brain (Total brain − structure)",
      y = paste0(pretty_region, " volume")
    )
}
######### INDIV SECTION PLOTS PT 1 END


# --- Helper: extract lambda for each model spec
lambda_from_fit <- function(fit, model_type, fixed_lambda = NA_real_) {
  if (model_type == "PagelML") {
    as.numeric(coef(fit$modelStruct, unconstrained = FALSE))
  } else if (model_type == "Fixed") {
    as.numeric(fixed_lambda)
  } else if (model_type == "BM") {
    1.0
  } else {
    stop("Unknown model_type")
  }
}

# --- Helper: phylogenetic conditional mean correction mu and conditional variance ch
phylo_mu_ch <- function(v_full, v_reduced, human, X_vec) {
  cc <- which(rownames(v_full) == human)
  others <- rownames(v_reduced)
  Cih <- v_full[cc, others, drop = FALSE]

  inv_v <- solve(v_reduced)
  xbar <- mean(X_vec, na.rm = TRUE)

  mu <- Cih %*% inv_v %*% (X_vec - xbar)
  ch <- v_full[cc, cc] - Cih %*% inv_v %*% t(Cih)

  list(mu = as.numeric(mu), ch = as.numeric(ch))
}

# --- Helper: run a set of models for one structure
# mu_source controls what you used as X in the phylo correction:
#   "trait"     matches your Plot I/II function (X = log(trait))
#   "predictor" matches your Plot III loop as written (X = log(ROB))
run_models_for_structure <- function(pp,
                                     model_specs,
                                     want_ci = TRUE,
                                     mu_source = c("trait", "predictor")) {
  mu_source <- match.arg(mu_source)
  if (is.null(pp)) return(NULL)

  out <- list()

  for (spec in model_specs) {
    cor_struct <- spec$cor(pp$tr2)
    fit <- fit_gls_safe(pp$fmla, cor_struct, pp$data_model)
    if (is.null(fit)) next

    sigma <- fit$sigma
    lambda_est <- lambda_from_fit(fit, spec$type, spec$lambda)

    # Build vcv with lambda
    v_full <- transform_vcv_lambda(pp$tr1, lambda_est)
    v_red  <- transform_vcv_lambda(pp$tr2, lambda_est)

### CHECK THIS PART OF SCRIPT BELOW
    # Choose X for phylo correction (to preserve your original behavior)
    if (mu_source == "trait") {
      X <- log(pp$data_model[[pp$var_name]])
    } else {
      X <- log(pp$data_model$Rest_of_Brain)
    }

    mc <- phylo_mu_ch(v_full, v_red, pp$human, X)
    mu <- mc$mu
    ch <- mc$ch

    if (!is.finite(ch) || ch <= 0) ch <- 1  # preserve your “clamp” behavior

    xh <- log(pp$xh_raw)
    pred_log <- as.numeric(c(1, xh) %*% coef(fit) + mu)

    if (want_ci) {
      se <- sqrt(as.numeric(sigma^2 * ch))
      upper_log <- pred_log + qnorm(0.975) * se
      lower_log <- pred_log - qnorm(0.975) * se

      out[[spec$name]] <- data.frame(
        Variable  = pp$var_name,
        Model     = spec$name,
        Observed  = pp$obs_val,
        Lower     = exp(lower_log),
        Predicted = exp(pred_log),
        Upper     = exp(upper_log),
        lambda    = lambda_est,
        N         = nrow(pp$data_model),
        stringsAsFactors = FALSE
      )
    } else {
      out[[spec$name]] <- data.frame(
        Variable  = pp$var_name,
        Model     = spec$name,
        Observed  = pp$obs_val,
        Predicted = exp(pred_log),
        lambda    = lambda_est,
        N         = nrow(pp$data_model),
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(out) == 0) return(NULL)
  do.call(rbind, out)
}

## -------------------------
## MODEL SPECIFICATIONS
## -------------------------

specs_plot12 <- list(
  list(
    name = "Brownian (BM)",
    type = "BM",
    lambda = NA_real_,
    cor = function(tr2) corBrownian(1, form = ~Species, phy = tr2)
  ),
  list(
    name = "Pagel's lambda (ML)",
    type = "PagelML",
    lambda = NA_real_,
    cor = function(tr2) corPagel(1, form = ~Species, phy = tr2, fixed = FALSE)
  )
)

specs_plot3 <- list(
  list(
    name = "Independence (λ=0)",
    type = "Fixed",
    lambda = 0,
    cor = function(tr2) corPagel(0, form = ~Species, phy = tr2, fixed = TRUE)
  ),
  list(
    name = "Pagel's ML (Estimated)",
    type = "PagelML",
    lambda = NA_real_,
    cor = function(tr2) corPagel(1, form = ~Species, phy = tr2, fixed = FALSE)
  ),
  list(
    name = "Brownian (λ=1)",
    type = "BM",
    lambda = NA_real_,
    cor = function(tr2) corBrownian(1, form = ~Species, phy = tr2)
  )
)

## -------------------------
## PLOT I + II DATA (one pass)
## -------------------------

core_df <- purrr::map_dfr(target_cols, function(vn) {
  pp <- prep_structure(vn, data = data, tr = tr, human = "Homo_sapiens")
  run_models_for_structure(pp, specs_plot12, want_ci = TRUE, mu_source = "trait")
})

core_df$Variable <- factor(core_df$Variable, levels = target_cols)
core_df$VarLabel <- label_region(core_df$Variable)
core_df$CorticalGroup2 <- group_region_cortical2(core_df$Variable)
core_df$CorticalGroup3 <- group_region_cortical3(core_df$Variable)
core_df$RegionGroup <- core_df$CorticalGroup2

## -------------------------
## PLOT I (unchanged plotting code)
## -------------------------
p1 <- ggplot(core_df, aes(x = VarLabel)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, color = "#377eb8") +
  geom_point(aes(y = Predicted, shape = "Predicted"), color = "#377eb8", size = 2.5) +
  geom_point(aes(y = Observed, shape = "Observed"), color = "#e41a1c", size = 2.5) +
  scale_shape_manual(name = "Value", values = c("Predicted" = 16, "Observed" = 17)) +
  coord_flip() +
  facet_wrap(vars(RegionGroup, Model), scales = "free_x") +
  theme_bw() +
  labs(
    title = "Prediction: Brownian Motion vs Pagel's lambda (Part-Whole Corrected)",
    subtitle = "Predictor = Rest of Brain (Total - Structure)",
    y = "Volume / Value (Original Scale)",
    x = "Brain Structure"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 11),
    axis.text.y = element_text(size = 9)
  )

p1

## -------------------------
## PLOT II (normalization + plot unchanged)
## -------------------------
final_df_norm <- core_df %>%
  group_by(Variable) %>%
  mutate(
    local_min = min(c(Lower, Upper, Observed, Predicted)),
    local_max = max(c(Lower, Upper, Observed, Predicted)),
    Predicted_Sc = (Predicted - local_min) / (local_max - local_min),
    Observed_Sc  = (Observed  - local_min) / (local_max - local_min),
    Lower_Sc     = (Lower     - local_min) / (local_max - local_min),
    Upper_Sc     = (Upper     - local_min) / (local_max - local_min)
  ) %>%
  ungroup()


p2<-ggplot(final_df_norm, aes(x = VarLabel)) +
  geom_errorbar(aes(ymin = Lower_Sc, ymax = Upper_Sc), width = 0.2, color = "#377eb8") +
  geom_point(aes(y = Predicted_Sc, shape = "Predicted"), color = "#377eb8", size = 2.5) +
  geom_point(aes(y = Observed_Sc, shape = "Observed"), color = "#e41a1c", size = 2.5) +
  scale_shape_manual(name = "Value", values = c("Predicted" = 16, "Observed" = 17)) +
  coord_flip() +
  facet_wrap(vars(RegionGroup, Model), scales = "free_x") +
  theme_bw() +
  labs(
    title = "Standardized Prediction Error (0-1 Scale, Part-Whole Corrected)",
    subtitle = "0 = Min value for structure, 1 = Max value for structure (across CI and Obs). Predictor = Rest of Brain",
    y = "Standardized Position (0 to 1)",
    x = "Brain Structure"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 11),
    axis.text.y = element_text(size = 9)
  )

p2
  
## -------------------------
## PLOT III DATA (one pass, no CI)
## NOTE: mu_source = "predictor" preserves your Plot III as-written behavior.
## -------------------------
final_df <- purrr::map_dfr(target_cols, function(vn) {
  pp <- prep_structure(vn, data = data, tr = tr, human = "Homo_sapiens")
  run_models_for_structure(pp, specs_plot3, want_ci = FALSE, mu_source = "predictor")
})

final_df$PropDiff <- (final_df$Observed - final_df$Predicted) / final_df$Predicted
final_df$Direction <- ifelse(final_df$PropDiff > 0, "Larger than Predicted", "Smaller than Predicted")
final_df$Model <- factor(final_df$Model, levels = c("Independence (λ=0)", "Pagel's ML (Estimated)", "Brownian (λ=1)"))
final_df$VarLabel <- label_region(final_df$Variable)
final_df$CorticalGroup2 <- group_region_cortical2(final_df$Variable)
final_df$CorticalGroup3 <- group_region_cortical3(final_df$Variable)
final_df$RegionGroup <- final_df$CorticalGroup2

p3<-ggplot(final_df, aes(x = reorder(VarLabel, PropDiff), y = PropDiff, color = Direction)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_segment(aes(xend = reorder(VarLabel, PropDiff), yend = 0), size = 1.2) +
  geom_point(size = 4) +
  geom_text(aes(label = ifelse(Model == "Pagel's ML (Estimated)",
                               sprintf("λ=%.2f", lambda), "")),
            nudge_x = -0.4,
            size = 3,
            fontface = "italic",
            color = "black",
            show.legend = FALSE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_color_manual(values = c("Larger than Predicted" = "#e41a1c",
                                "Smaller than Predicted" = "#377eb8")) +
  coord_flip() +
  facet_wrap(vars(RegionGroup, Model), scales = "free_x") +
  theme_bw() +
  labs(
    title = "Human Brain Mosaicism (Corrected for Part-Whole)",
    subtitle = "Predictor = Rest of Brain (Total - Structure).",
    y = "Deviation from Prediction (%)",
    x = "Brain Structure",
    color = "Direction"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    axis.text.y = element_text(size = 9, face = "bold")
  )

p3

## -------------------------
## PLOT IV: lambda likelihood profiles (uses same prep_structure)
## -------------------------
lambda_seq <- seq(0, 1, by = 0.05)

profile_tmp <- purrr::map(target_cols, function(vn) {
  pp <- prep_structure(vn, data = data, tr = tr, human = "Homo_sapiens")
  if (is.null(pp)) return(list(scan = NULL, mle = NULL))

  # MLE (estimated lambda)
  fit_mle <- fit_gls_safe(
    pp$fmla,
    corPagel(1, form = ~Species, phy = pp$tr2, fixed = FALSE),
    pp$data_model
  )

  mle_row <- NULL
  if (!is.null(fit_mle)) {
    best_lambda <- as.numeric(coef(fit_mle$modelStruct, unconstrained = FALSE))
    best_loglik <- as.numeric(logLik(fit_mle))
    best_lambda_plot <- pmax(0, pmin(1, best_lambda))
    mle_row <- data.frame(
      Variable = vn,
      CorticalGroup2 = group_region_cortical2(vn),
      CorticalGroup3 = group_region_cortical3(vn),
      RegionGroup = group_region_cortical2(vn),
      lambda = best_lambda_plot,
      LogLik = best_loglik
    )
  }

  # Profile scan (fixed lambda)
  scan_df <- purrr::map_dfr(lambda_seq, function(val) {
    fit_scan <- fit_gls_safe(
      pp$fmla,
      corPagel(val, form = ~Species, phy = pp$tr2, fixed = TRUE),
      pp$data_model
    )
    if (is.null(fit_scan)) return(NULL)
    data.frame(
      Variable = vn,
      CorticalGroup2 = group_region_cortical2(vn),
      CorticalGroup3 = group_region_cortical3(vn),
      RegionGroup = group_region_cortical2(vn),
      lambda = val,
      LogLik = as.numeric(logLik(fit_scan))
    )
  })

  list(scan = scan_df, mle = mle_row)
})

df_profile <- dplyr::bind_rows(purrr::map(profile_tmp, "scan"))
df_mle     <- dplyr::bind_rows(purrr::map(profile_tmp, "mle"))

p4<-ggplot(df_profile, aes(x = lambda, y = LogLik)) +
  geom_line(color = "#377eb8", size = 1) +
  geom_point(data = df_mle, aes(x = lambda, y = LogLik), color = "red", size = 3) +
  geom_vline(data = df_mle, aes(xintercept = lambda), linetype = "dashed", color = "red", alpha = 0.5) +
  geom_text(data = df_mle, aes(label = sprintf("%.2f", lambda), x = 0.1, y = LogLik),
            color = "red", size = 3, hjust = 0, vjust = 1) +
  facet_wrap(
    vars(RegionGroup, Variable),
    scales = "free_y",
    labeller = labeller(Variable = function(x) label_region(x))
  ) +
  theme_bw() +
  labs(
    title = "lambda Likelihood Profiles (Part-Whole Corrected)",
    subtitle = "Predictor = Rest of Brain. Curve shows model fit at different lambda values. Red dot = MLE.",
    x = "lambda (Phylogenetic Signal)",
    y = "Log-Likelihood"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 9),
    axis.text.y = element_text(size = 7)
  )
  
p4  

# ============================================================
# Save plots
# ============================================================

ggsave("figs/s3/cortical/plot_prediction_BM_vs_pagel.png", p1, width = 10, height = 7, dpi = 300)
ggsave("figs/s3/cortical/plot_standardized_prediction_error.png", p2, width = 10, height = 7, dpi = 300)
ggsave("figs/s3/cortical/plot_human_brain_mosaicism.png", p3, width = 12, height = 8, dpi = 300)
ggsave("figs/s3/cortical/plot_pagel_lambda_profile.png", p4, width = 12, height = 9, dpi = 300)

ggsave("figs/s3/cortical/plot_prediction_BM_vs_pagel.pdf", p1, width = 10, height = 7)
ggsave("figs/s3/cortical/plot_standardized_prediction_error.pdf", p2, width = 10, height = 7)
ggsave("figs/s3/cortical/plot_human_brain_mosaicism.pdf", p3, width = 12, height = 8)
ggsave("figs/s3/cortical/plot_pagel_lambda_profile.pdf", p4, width = 12, height = 9)


######### INDIV SECTION PLOTS PT 2 START
## -------------------------
## PLOT V: individual BM PGLS plots by structure
## -------------------------

# Use target_cols to preserve your intended structure order.
# Raw names are used for modeling; pretty labels are used for plot names and axis labels.
plots_raw <- target_cols %>%
  set_names(label_region(.)) %>%
  map(
    ~ plot_region_gls(
      region = .x,
      data = data,
      tr = tr,
      human = "Homo_sapiens"
    )
  )

# Check failures, if any
failed_regions <- names(plots_raw)[map_lgl(plots_raw, is.null)]
failed_regions

# Keep only successful plots
plots <- compact(plots_raw)

# Confirm number of plots
length(plots)
names(plots)

# Print all plots one by one
walk(plots, print)

# Save all plots to a multipage PDF
pdf("figs/s3/cortical/phylo_gls_region_plots.pdf", width = 7, height = 5)

walk(plots, print)

dev.off()

# Save one presentation-ready JPG per structure
dir.create("figs/s3/cortical/phylo_gls_region_jpgs", showWarnings = FALSE)

iwalk(
  plots,
  function(p, region_label) {
    
    safe_region <- gsub("[^A-Za-z0-9]+", "_", region_label)
    safe_region <- gsub("_+$", "", safe_region)
    
    ggsave(
      filename = file.path(
        "figs/s3/cortical/phylo_gls_region_jpgs",
        paste0("phylo_gls_", safe_region, ".jpg")
      ),
      plot = p,
      width = 13.333,
      height = 7.5,
      units = "in",
      dpi = 300,
      device = "jpeg",
      bg = "white"
    )
  }
)
######### INDIV SECTION PLOTS PT 2 END


## -------------------------
## TABLE I: Predicted values + CIs (one pass, with CI) and rCMRGlc
## -------------------------
# --- 1) Build table with predicted values + CIs (from Plot I/II code, with some extra columns for later)

digits_lambda <- 3
digits_vals   <- 5   # Predicted, CIs, Observed

core_df_out <- core_df %>%
# 1) Column calculations
  mutate(
    Variable = factor(Variable, levels = target_cols),
    VarLabel = label_region(Variable),
    CorticalGroup2 = group_region_cortical2(Variable),
    CorticalGroup3 = group_region_cortical3(Variable),
    RegionGroup = CorticalGroup2,
    Diff.min = (Observed - Lower) / Observed,
    Diff.pre = (Observed - Predicted) / Observed,
    Diff.max = (Observed - Upper) / Observed
    ) %>%
# 2) Shape table + renaming
  transmute(
    CorticalGroup2,
    CorticalGroup3,
    RegionGroup,
    Structure       = VarLabel,
    Model,
    lambda,
    `95% CI min`    = Lower,
    Predicted,
    `95% CI max`    = Upper,
    Observed,
    Diff.min,
    Diff.pre,
    Diff.max,
    N
  ) %>%
# 3) Order rows by label order and model
  mutate(
    CorticalGroup2 = factor(CorticalGroup2, levels = cortical_2_group_levels),
    CorticalGroup3 = factor(CorticalGroup3, levels = cortical_3_group_levels),
    RegionGroup = factor(RegionGroup, levels = cortical_2_group_levels),
    Structure = factor(Structure, levels = label_region(target_cols))
  ) %>%
  arrange(CorticalGroup2, CorticalGroup3, Structure, Model) %>%
# 4) Limit significant figures & clean N
  mutate(
    lambda       = signif(lambda,       digits_lambda),
    `95% CI min` = signif(`95% CI min`, digits_vals),
    Predicted    = signif(Predicted,    digits_vals),
    `95% CI max` = signif(`95% CI max`, digits_vals),
    Observed     = signif(Observed,     digits_vals),
    N            = as.integer(N)
  )

# --- 2) Add rCMRGlc by matching Structure from Sup Table 1

# Keep only the relevant columns; ensure rCMRGlc is numeric and rounded to 1 decimal place
sup_clean <- heiss_stephan_tbl %>%
  select(
    `volume_term`,
    `rCMRGlc_mean_both_hemispheres`
  ) %>%
  rename(
    Structure = `volume_term`,
    rCMRGlc   = `rCMRGlc_mean_both_hemispheres`
  ) %>%
  mutate(
    rCMRGlc = round(as.numeric(rCMRGlc), 1)
  )

# Join rCMRGlc into the main table by Structure
core_with_rCMRGlc <- core_df_out %>%
  left_join(
    sup_clean %>% select(Structure, rCMRGlc),
    by = "Structure"
  )

# Save model + rCMRGlc table for downstream astrocyte plotting
dir.create("tables", showWarnings = FALSE, recursive = TRUE)

write_csv(
  core_with_rCMRGlc,
  "checks/s3/cortical/core_with_rCMRGlc_predicted_volumes.csv"
)

# --- 3) Print-ready table (BM vs ML grouped columns)

metric_cols <- c(
  "lambda", "95% CI min", "Predicted", "95% CI max",
  "Observed", "Diff.min", "Diff.pre", "Diff.max"
)

# model-specific metrics only
model_wide <- core_with_rCMRGlc %>%
  mutate(
    ModelGroup = case_when(
      grepl("brownian|\\bbm\\b", tolower(Model)) ~ "Brownian (BM)",
      grepl("pagel|lambda|\\bml\\b", tolower(Model)) ~ "Pagel's lambda (ML)",
      TRUE ~ as.character(Model)
    ),
    ModelGroup = factor(ModelGroup, levels = c("Brownian (BM)", "Pagel's lambda (ML)"))
  ) %>%
  select(CorticalGroup2, CorticalGroup3, RegionGroup, Structure, ModelGroup, all_of(metric_cols)) %>%
  pivot_wider(
    id_cols = c(CorticalGroup2, CorticalGroup3, RegionGroup, Structure),
    names_from  = ModelGroup,
    values_from = all_of(metric_cols),
    names_glue  = "{ModelGroup}__{.value}"
  )

# shared columns (same across BM/ML for each structure)
shared_cols <- core_with_rCMRGlc %>%
  group_by(CorticalGroup2, CorticalGroup3, RegionGroup, Structure) %>%
  summarise(
    N = dplyr::first(N),
    rCMRGlc = dplyr::first(rCMRGlc),
    .groups = "drop"
  )

# enforce final column order
bm_cols <- paste0("Brownian (BM)__", metric_cols)
ml_cols <- paste0("Pagel's lambda (ML)__", metric_cols)

tbl_body <- model_wide %>%
  left_join(shared_cols, by = c("CorticalGroup2", "CorticalGroup3", "RegionGroup", "Structure")) %>%
  arrange(CorticalGroup2, CorticalGroup3, Structure) %>%
  select(CorticalGroup2, CorticalGroup3, RegionGroup, Structure, all_of(bm_cols), all_of(ml_cols), N, rCMRGlc)

# from your existing tbl_body (Structure + BM block + ML block + N + rCMRGlc)

tbl_body2 <- tbl_body %>%
  select(-any_of(c("Model", "ModelGroup")))   # removes model column if present

# row 1: show each model name only once (first column of each block)
header_top <- c(
  "CorticalGroup2",
  "CorticalGroup3",
  "RegionGroup",
  "Structure",
  "Brownian (BM)", rep("", length(metric_cols) - 1),
  "Pagel's lambda (ML)", rep("", length(metric_cols) - 1),
  "N", "rCMRGlc" 
)

# row 2: subcolumn labels
header_sub <- c("", "", "", "", metric_cols, metric_cols, "", "")

print_ready_df <- as.data.frame(
  rbind(
    header_top,
    header_sub,
    as.matrix(tbl_body2)
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

names(print_ready_df) <- rep("", ncol(print_ready_df))
# Save the table with data about predicted values, CIs, and rCMRGlc to Excel
write_xlsx(print_ready_df, "tables/s3/cortical/Table 4 Predicted volumes for Brownian and Pagel's lambda models.xlsx", col_names = FALSE)

## -------------------------
## TO REPORT IN TEXT OR SUP: Correlation analysis between rCMRGlc and prediction error (Diff.pre)
## Analyses are run separately for both grouping schemes:
##   1) Cortical vs non-cortical
##   2) Cortical gray vs white matter vs non-cortical
## -------------------------

diff_cols <- c(
  "Brownian (BM)__Diff.min",
  "Brownian (BM)__Diff.pre",
  "Brownian (BM)__Diff.max",
  "Pagel's lambda (ML)__Diff.min",
  "Pagel's lambda (ML)__Diff.pre",
  "Pagel's lambda (ML)__Diff.max"
)

num <- function(x) readr::parse_number(as.character(x))

safe_cor_test <- function(x, y, method = c("pearson", "spearman")) {
  method <- match.arg(method)
  ok <- complete.cases(x, y)
  n_ok <- sum(ok)

  # Need at least 3 complete observations and non-zero variance in both vectors.
  if (n_ok < 3 || sd(x[ok]) == 0 || sd(y[ok]) == 0) {
    return(tibble(
      n = n_ok,
      estimate = NA_real_,
      p_value = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_
    ))
  }

  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = method))

  tibble(
    n = n_ok,
    estimate = unname(ct$estimate),
    p_value = ct$p.value,
    ci_low = if (!is.null(ct$conf.int)) ct$conf.int[1] else NA_real_,
    ci_high = if (!is.null(ct$conf.int)) ct$conf.int[2] else NA_real_
  )
}

# ============================================================
# rCMRGlc grouped analysis helper
# ============================================================

run_rCMRGlc_group_analyses <- function(core_with_rCMRGlc,
                                       tbl_body,
                                       group_col,
                                       group_levels,
                                       output_suffix,
                                       group_label) {
  group_col <- as.character(group_col)

  # -------------------------
  # Correlations by group
  # -------------------------
  cor_df <- tbl_body %>%
    transmute(
      Group = factor(as.character(.data[[group_col]]), levels = group_levels),
      Structure,
      rCMRGlc = parse_number(as.character(rCMRGlc)),
      across(all_of(diff_cols), ~ as.numeric(as.character(.x)))
    )

  cor_pearson <- cor_df %>%
    group_by(Group) %>%
    group_modify(~ map_dfr(diff_cols, function(col) {
      safe_cor_test(.x[[col]], .x$rCMRGlc, method = "pearson") %>%
        mutate(Diff_column = col, .before = 1)
    })) %>%
    ungroup() %>%
    rename(r = estimate)

  cor_spearman <- cor_df %>%
    group_by(Group) %>%
    group_modify(~ map_dfr(diff_cols, function(col) {
      safe_cor_test(.x[[col]], .x$rCMRGlc, method = "spearman") %>%
        mutate(Diff_column = col, .before = 1)
    })) %>%
    ungroup() %>%
    rename(rho = estimate)

  write_csv(
    cor_pearson,
    paste0("checks/s3/cortical/correlation_rCMRGlc_prediction_error_by_", output_suffix, "_pearson.csv")
  )
  write_csv(
    cor_spearman,
    paste0("checks/s3/correlation_rCMRGlc_prediction_error_by_", output_suffix, "_spearman.csv")
  )

  # -------------------------
  # Data for rCMRGlc plots
  # -------------------------
  df_bm <- core_with_rCMRGlc %>%
    filter(Model == "Brownian (BM)") %>%
    transmute(
      Group = factor(as.character(.data[[group_col]]), levels = group_levels),
      Structure,
      rCMRGlc  = num(rCMRGlc),
      Diff.pre = num(Diff.pre),
      Diff.min = num(Diff.min),
      Diff.max = num(Diff.max),
      excluded = Structure == "Neocortex white"
    ) %>%
    drop_na(rCMRGlc, Diff.pre)

  df_fit <- df_bm %>%
    filter(!excluded)

  df_neo <- df_bm %>%
    filter(excluded)

  xg <- seq(10, 40, length.out = 300)

  x_lab <- "rCMRGlc (µmol/100 g/min)"
  y_lab <- "Difference from prediction (BM)"

  theme_paper <- theme_bw(base_size = 16) +
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

  add_points_and_labels <- list(
    geom_point(
      data = df_fit,
      aes(x = rCMRGlc, y = Diff.pre),
      inherit.aes = FALSE,
      pch = 16,
      size = 2.2
    ),
    geom_text(
      data = df_fit,
      aes(x = rCMRGlc, y = Diff.pre, label = Structure),
      inherit.aes = FALSE,
      hjust = -0.08,
      size = 2.8,
      check_overlap = TRUE
    ),
    geom_point(
      data = df_neo,
      aes(x = rCMRGlc, y = Diff.pre),
      inherit.aes = FALSE,
      shape = 24,
      fill = "firebrick2",
      color = "firebrick4",
      size = 4.5
    ),
    geom_text(
      data = df_neo,
      aes(x = rCMRGlc, y = Diff.pre, label = "Neocortex white (excluded)"),
      inherit.aes = FALSE,
      hjust = -0.08,
      size = 3.1,
      color = "firebrick4"
    )
  )

  # ============================================================
  # Plot: LOESS with Diff.min / Diff.max error bars, fit separately by group
  # ============================================================

  p_loess_bars <- ggplot() +
    geom_errorbar(
      data = df_bm,
      aes(x = rCMRGlc, ymin = Diff.min, ymax = Diff.max),
      inherit.aes = FALSE,
      width = 0,
      na.rm = TRUE
    ) +
    geom_smooth(
      data = df_fit,
      aes(x = rCMRGlc, y = Diff.pre),
      method = "loess",
      se = FALSE,
      span = 0.75,
      color = "steelblue4",
      linewidth = 1.2
    ) +
    add_points_and_labels +
    facet_wrap(~Group) +
    coord_cartesian(
      xlim = c(10, 40),
      ylim = c(-3.5, 1.5),
      clip = "on"
    ) +
    labs(
      title = "LOESS fit with error bars",
      subtitle = paste0("Fits are stratified by ", group_label, "; Neocortex white is plotted but excluded from fits"),
      x = x_lab,
      y = y_lab
    ) +
    theme_paper

  # ============================================================
  # Plot: Polynomial comparison, fit separately by group
  # Degree is capped within each group to avoid overfitting/singular models.
  # ============================================================

  fit_poly_by_group <- function(df_group, max_degree = 5) {
    group_name <- as.character(df_group$Group[1])
    n_obs <- nrow(df_group)
    n_unique_x <- n_distinct(df_group$rCMRGlc)
    max_allowed <- min(max_degree, n_obs - 2, n_unique_x - 1)

    if (!is.finite(max_allowed) || max_allowed < 1) {
      return(list(curves = NULL, stats = NULL))
    }

    degrees <- seq_len(max_allowed)

    fit_poly <- function(d) {
      form <- if (d == 1) {
        Diff.pre ~ rCMRGlc
      } else {
        as.formula(sprintf("Diff.pre ~ poly(rCMRGlc, %d, raw = TRUE)", d))
      }

      lm(form, data = df_group)
    }

    fits <- map(degrees, fit_poly)
    names(fits) <- as.character(degrees)

    adj_r2 <- map_dbl(fits, ~ summary(.x)$adj.r.squared)
    best_degree <- degrees[which.max(adj_r2)]
    best_fit <- fits[[as.character(best_degree)]]

    best_sum <- summary(best_fit)
    fstat <- best_sum$fstatistic
    model_p <- if (!is.null(fstat)) pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE) else NA_real_

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

    poly_stats <- paste(
      paste0("Best degree: ", best_degree),
      paste0("Adj. R²: ", sprintf("%.4f", adj_r2[which(degrees == best_degree)])),
      paste0("Model p: ", format.pval(model_p, digits = 3, eps = 1e-4)),
      make_eqn(coef(best_fit)),
      sep = "\n"
    )

    curves <- map_dfr(degrees, function(d) {
      tibble(
        Group = group_name,
        degree = d,
        is_best = d == best_degree,
        degree_label = paste0(
          "Degree ", d,
          "  Adj. R² = ", sprintf("%.3f", adj_r2[which(degrees == d)])
        ),
        rCMRGlc = xg,
        Diff.pre = predict(
          fits[[as.character(d)]],
          newdata = data.frame(rCMRGlc = xg)
        )
      )
    })

    stats <- tibble(
      Group = group_name,
      label = poly_stats,
      x = 10.6,
      y = 1.35
    )

    list(curves = curves, stats = stats)
  }

  poly_results <- df_fit %>%
    group_split(Group) %>%
    map(fit_poly_by_group)

  poly_curves <- bind_rows(map(poly_results, "curves"))
  poly_stats_df <- bind_rows(map(poly_results, "stats"))

  p_poly_compare <- ggplot() +
    geom_line(
      data = poly_curves,
      aes(x = rCMRGlc, y = Diff.pre, color = degree_label),
      linewidth = 0.9
    ) +
    geom_line(
      data = filter(poly_curves, is_best),
      aes(x = rCMRGlc, y = Diff.pre),
      linewidth = 1.2,
      linetype = "dashed",
      color = "black"
    ) +
    add_points_and_labels +
    geom_text(
      data = poly_stats_df,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 3.2,
      color = "navy"
    ) +
    facet_wrap(~Group) +
    coord_cartesian(
      xlim = c(10, 40),
      ylim = c(-3.5, 1.5),
      clip = "on"
    ) +
    labs(
      title = paste0("Polynomial fits by ", group_label),
      subtitle = "Models exclude Neocortex white; maximum degree is capped within each group",
      x = x_lab,
      y = y_lab,
      color = NULL
    ) +
    theme_paper +
    theme(legend.position = "right")

  # ============================================================
  # Plot: Quadratic fit with 95% confidence interval, fit separately by group
  # ============================================================

  fit_quadratic_by_group <- function(df_group) {
    group_name <- as.character(df_group$Group[1])

    if (nrow(df_group) < 4 || n_distinct(df_group$rCMRGlc) < 3) {
      return(list(pred = NULL, stats = NULL))
    }

    quad_fit <- lm(Diff.pre ~ poly(rCMRGlc, 2, raw = TRUE), data = df_group)

    quad_sum <- summary(quad_fit)
    quad_adj_r2 <- quad_sum$adj.r.squared
    quad_fstat <- quad_sum$fstatistic
    quad_p <- if (!is.null(quad_fstat)) pf(quad_fstat[1], quad_fstat[2], quad_fstat[3], lower.tail = FALSE) else NA_real_

    b <- coef(quad_fit)

    quad_eqn <- paste0(
      "y = ", sprintf("%.4f", b[1]),
      ifelse(b[2] >= 0, " + ", " - "), sprintf("%.4f", abs(b[2])), "x",
      ifelse(b[3] >= 0, " + ", " - "), sprintf("%.4f", abs(b[3])), "x²"
    )

    quad_stats <- paste0(
      "Adj. R² = ", sprintf("%.3f", quad_adj_r2),
      "\nModel p = ", format.pval(quad_p, digits = 2, eps = 1e-4),
      "\n", quad_eqn
    )

    quad_pred <- data.frame(rCMRGlc = xg)

    quad_pred <- cbind(
      quad_pred,
      as.data.frame(
        predict(
          quad_fit,
          newdata = quad_pred,
          interval = "confidence"
        )
      )
    ) %>%
      as_tibble() %>%
      mutate(Group = group_name, .before = 1)

    quad_stats_df <- tibble(
      Group = group_name,
      label = quad_stats,
      x = 10.6,
      y = 1.35
    )

    list(pred = quad_pred, stats = quad_stats_df)
  }

  quad_results <- df_fit %>%
    group_split(Group) %>%
    map(fit_quadratic_by_group)

  quad_pred <- bind_rows(map(quad_results, "pred"))
  quad_stats_df <- bind_rows(map(quad_results, "stats"))

  p_quad_ci <- ggplot() +
    geom_ribbon(
      data = quad_pred,
      aes(x = rCMRGlc, ymin = lwr, ymax = upr),
      inherit.aes = FALSE,
      fill = "steelblue",
      alpha = 0.2
    ) +
    geom_line(
      data = quad_pred,
      aes(x = rCMRGlc, y = fit),
      inherit.aes = FALSE,
      color = "steelblue4",
      linewidth = 1.2
    ) +
    add_points_and_labels +
    geom_text(
      data = quad_stats_df,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 3.2,
      color = "navy"
    ) +
    facet_wrap(~Group) +
    coord_cartesian(
      xlim = c(10, 40),
      ylim = c(-3.5, 1.5),
      clip = "on"
    ) +
    labs(
      title = paste0("Quadratic fits by ", group_label),
      subtitle = "Fits exclude Neocortex white; plotted points include it",
      x = x_lab,
      y = y_lab
    ) +
    theme_paper

  # ============================================================
  # Save grouped outputs
  # ============================================================

  ggsave(paste0("figs/s3/cortical/plot_loess_errorbars_BM_rCMRGlc_by_", output_suffix, ".png"), p_loess_bars, width = 9, height = 5, dpi = 300)
  ggsave(paste0("figs/s3/cortical/plot_polynomial_comparison_BM_rCMRGlc_by_", output_suffix, ".png"), p_poly_compare, width = 11, height = 5.5, dpi = 300)
  ggsave(paste0("figs/s3/cortical/plot_quadratic_CI_BM_rCMRGlc_by_", output_suffix, ".png"), p_quad_ci, width = 9, height = 5, dpi = 300)

  ggsave(paste0("figs/s3/cortical/plot_loess_errorbars_BM_rCMRGlc_by_", output_suffix, ".pdf"), p_loess_bars, width = 9, height = 5)
  ggsave(paste0("figs/s3/cortical/plot_polynomial_comparison_BM_rCMRGlc_by_", output_suffix, ".pdf"), p_poly_compare, width = 11, height = 5.5)
  ggsave(paste0("figs/s3/cortical/plot_quadratic_CI_BM_rCMRGlc_by_", output_suffix, ".pdf"), p_quad_ci, width = 9, height = 5)

  invisible(list(
    cor_pearson = cor_pearson,
    cor_spearman = cor_spearman,
    p_loess_bars = p_loess_bars,
    p_poly_compare = p_poly_compare,
    p_quad_ci = p_quad_ci
  ))
}

# Run both requested grouping schemes.
rCMRGlc_results_cortical2 <- run_rCMRGlc_group_analyses(
  core_with_rCMRGlc = core_with_rCMRGlc,
  tbl_body = tbl_body,
  group_col = "CorticalGroup2",
  group_levels = cortical_2_group_levels,
  output_suffix = "cortical2_group",
  group_label = "cortical vs non-cortical group"
)

rCMRGlc_results_cortical3 <- run_rCMRGlc_group_analyses(
  core_with_rCMRGlc = core_with_rCMRGlc,
  tbl_body = tbl_body,
  group_col = "CorticalGroup3",
  group_levels = cortical_3_group_levels,
  output_suffix = "cortical3_group",
  group_label = "cortical gray / white matter / non-cortical group"
)
