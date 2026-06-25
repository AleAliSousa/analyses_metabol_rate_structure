  # ============================================================
# Study 3 — Predicted vs observed brain-region volumes in Homo sapiens
# Phylogenetic GLS (part-whole corrected) + relationship to rCMRGlc
#
# Tidied / compacted version of s3_predicValuesPGLS_01062026.R (2026-06-04).
# Modelling logic is UNCHANGED; this version reorganises, de-duplicates,
# removes dead code, routes all figure saves through one helper, and wires
# the Phase 1 fixes in as CONFIG toggles so previous results stay reproducible.
#
# This variant ADDS a "Remainder" structure (= Preferred_brain_volume minus the
# sum of the modelled regions) that is processed alongside the others. See the
# REMAINDER section below for the exact constituent regions and exclusions.
#
# Phase 1 fixes (see metadata/PHASE1_missing_data_strategy.md):
#   1.1  Common-species handling is explicit and configurable (CONFIG$restrict_species).
#        Strict listwise across all 13 regions leaves only 4 species, so the
#        DEFAULT remains max-data per-region fitting; see the diagnostic script
#        s3_0_missingness_clade_diagnostic_04062026.R to choose a common set,
#        and the MI script s3_1_phylo_multiple_imputation_04062026.R for the
#        balanced-sample sensitivity analysis.
#   1.2  Neocortex grey is de-overlapped:
#           Neocortex grey (remaining) = NeoG_Frahm - V1(ASG_Sousa) - Insula(grey)
#        V1 is present for all 24 neocortex species (free to subtract); insula
#        is missing for 15 of them and is only ~0.5-1.1% of neocortex grey, so
#        CONFIG$neocortex_subtract_insula lets you keep n=24 (V1-only) or honour
#        the full definition at n=9. See memo for the recommended route (impute
#        the ~1% insula term to keep n=24 AND remove both sub-regions).
# ============================================================

setwd("~/Library/CloudStorage/Dropbox/COLLABORATIVE/Do expensive brain regions increase less in humans/analyses_metabol_rate_structure")

library(here)
library(ape)
library(nlme)        # Mac note: ML optimizer differs from PC; opt = "optim" set below.
library(tidyverse)
library(dispRity)
library(scales)
library(readxl)
library(writexl)

options(scipen = 999)

# ============================================================
# CONFIG — edit here, not in the body
# ============================================================
CONFIG <- list(
  human         = "Homo_sapiens",
  predictor_col = "Preferred_brain_volume",
  
  # Phase 1.2 — neocortex de-overlap.
  # DEFAULT OFF: this script is the faithful TIDY of 03062026 and reproduces its
  # numbers. The "replace neocortex grey with the remainder" analysis (which also
  # needs the upstream rCMRGlc reweighted) lives in the separate v2 pipeline.
  deoverlap_neocortex      = FALSE,  # TRUE => NeoG - V1 (- insula); FALSE => raw NeoG_Frahm
  neocortex_subtract_insula = TRUE,  # TRUE => NeoG - V1 - insula (n=9, full definition)
  # FALSE => NeoG - V1 only (n=24, insula ~1% left in)
  
  # Phase 1.1 — restrict to a common species set (NULL = max-data per region)
  restrict_species = NULL,           # e.g. a character vector of Species to force a common set
  
  # I/O
  dir_figs   = "figs/s3/all",
  dir_checks = "checks/s3/all",
  dir_tables = "tables/s3/all"
)

# Create output directories up front (previously scattered / partly missing)
for (d in c(CONFIG$dir_figs,
            file.path(CONFIG$dir_figs, "phylo_gls_region_jpgs"),
            CONFIG$dir_checks,
            CONFIG$dir_tables)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# ============================================================
# LOAD
# ============================================================
heiss_stephan_tbl <- read.csv("data_intermediate/Heiss_Stephan_data.csv")
tr               <- read.tree("data_raw/species.nwk")
Stephan_primates <- read.csv("data_raw/Stephan_primates.csv")

# Drop legacy empty/index columns if still present
data_clean <- Stephan_primates[, !names(Stephan_primates) %in% c("X", "order")]

# Drop blank/junk rows: the source CSV ends with an all-empty row (no Species).
# (There are 59 real species; the 60th parsed row is that phantom blank line.)
data_clean <- data_clean[!is.na(data_clean$Species) & trimws(data_clean$Species) != "", ]

# Unified "preferred brain volume": first non-NA of the three candidate columns
data <- data_clean %>%
  mutate(
    Preferred_brain_volume = coalesce(Brain_volume, Brainvol, Total_brain_net_volume)
  )

# Sanity check: every named species should get a predictor. Warn, don't halt —
# any species still missing all three volume columns simply drops out of each fit.
n_missing_pred <- sum(is.na(data$Preferred_brain_volume))
if (n_missing_pred > 0) {
  warning(n_missing_pred, " species lack Brain_volume/Brainvol/Total_brain_net_volume ",
          "and will be excluded from fits: ",
          paste(data$Species[is.na(data$Preferred_brain_volume)], collapse = ", "))
}

# ============================================================
# PHASE 1.2 — de-overlap neocortex grey
#   V1 = ASG_Sousa (area striata grey); Insula = Total_insula_volume_L (grey).
#   NOTE: all 24 NeoG_Frahm species have V1; only 9 also have insula.
# ============================================================
neo_region <- "NeoG_Frahm"   # default modelled neocortex column

if (isTRUE(CONFIG$deoverlap_neocortex)) {
  insula_term <- if (isTRUE(CONFIG$neocortex_subtract_insula)) {
    coalesce(data$Total_insula_volume_L, NA_real_)   # NA where missing -> remainder NA (n=9)
  } else {
    0                                                # V1-only remainder (n=24)
  }
  data$Neocortex_grey_remaining <- data$NeoG_Frahm - data$ASG_Sousa - insula_term
  # Guard: a valid remainder must be positive
  bad <- which(!is.na(data$Neocortex_grey_remaining) & data$Neocortex_grey_remaining <= 0)
  if (length(bad)) {
    warning("Non-positive neocortex remainder set to NA for: ",
            paste(data$Species[bad], collapse = ", "))
    data$Neocortex_grey_remaining[bad] <- NA_real_
  }
  neo_region <- "Neocortex_grey_remaining"
  message(sprintf("Phase 1.2: neocortex de-overlapped (%s); usable n = %d",
                  if (CONFIG$neocortex_subtract_insula) "minus V1 + insula" else "minus V1 only",
                  sum(!is.na(data$Neocortex_grey_remaining))))
}

# ============================================================
# REMAINDER — brain volume not accounted for by the modelled regions.
#   Remainder = Preferred_brain_volume - sum(constituent regions below)
#   Uses the same neocortex-grey column as the analysis (neo_region), so it
#   honours the de-overlap toggle. By design this EXCLUDES Area striata grey
#   (ASG_Sousa) and Insular cortex grey (Total_insula_volume_L).
#   na.rm = FALSE on purpose: a species missing ANY constituent region gets
#   Remainder = NA — we do NOT silently treat missing sub-volumes as zero.
# ============================================================
remainder_regions <- c(
  "LGN_Sousa",             # Corpus geniculatum laterale
  "Amygdala",              # Amygdala
  "Pallidum",              # Pallidum
  "NeoW_Frahm",            # Neocortex white
  "Nucleus_subthalamicus", # Nucleus subthalamicus Luysi
  "Capsula_interna",       # Capsula interna
  "Striatum",              # Striatum
  neo_region,              # Neocortex grey (raw or de-overlapped per CONFIG)
  "Mesencephalon",         # Mesencephalon
  "Cerebellum",            # Cerebellum
  "Hippocampus"            # Hippocampus
)

missing_remainder_cols <- setdiff(remainder_regions, names(data))
if (length(missing_remainder_cols)) {
  stop("Remainder: missing expected columns: ",
       paste(missing_remainder_cols, collapse = ", "))
}

data$Remainder <- data[[CONFIG$predictor_col]] -
  rowSums(data[, remainder_regions], na.rm = FALSE)

# Guard: a valid remainder must be positive
bad_rem <- which(!is.na(data$Remainder) & data$Remainder <= 0)
if (length(bad_rem)) {
  warning("Non-positive Remainder set to NA for: ",
          paste(data$Species[bad_rem], collapse = ", "))
  data$Remainder[bad_rem] <- NA_real_
}

message(sprintf("Remainder computed; usable n = %d (of %d species require all %d regions)",
                sum(!is.na(data$Remainder)), nrow(data), length(remainder_regions)))

# ============================================================
# REGIONS — raw column -> display label, in plotting/table order.
#   The neocortex entry uses neo_region but keeps the label "Neocortex grey"
#   so the rCMRGlc join (volume_term = "Neocortex grey") still matches.
# ============================================================
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
  "__NEOCORTEX_GREY__"    = "Neocortex grey",   # raw column set below (de-overlap toggle)
  "Mesencephalon"         = "Mesencephalon",
  "Cerebellum"            = "Cerebellum",
  "Hippocampus"           = "Hippocampus",
  "Remainder"             = "Remainder"   # brain minus the modelled regions (computed above)
)
# Point the neocortex slot at the right raw column while preserving row order
names(region_labels)[names(region_labels) == "__NEOCORTEX_GREY__"] <- neo_region

target_cols <- names(region_labels)

# Robust raw-name -> label (safe with factors)
label_region <- function(x) {
  x_chr <- as.character(x)
  out <- unname(region_labels[x_chr])
  out[is.na(out)] <- x_chr[is.na(out)]
  out
}

# ============================================================
# HELPERS
# ============================================================

# Transform a correlation VCV by Pagel's lambda (scale off-diagonals, keep diagonal)
transform_vcv_lambda <- function(tree, lambda) {
  v <- vcv(tree, corr = TRUE)
  d <- diag(v)
  v <- v * lambda
  diag(v) <- d
  v
}

# Species with complete data across a set of region columns (+ predictor)
get_common_species <- function(data, region_cols, predictor_col = CONFIG$predictor_col) {
  cc <- complete.cases(data[, c(predictor_col, region_cols), drop = FALSE])
  as.character(data$Species[cc])
}

# Prepare one structure: select, drop NA, part-whole correct, align to tree, split human.
#   restrict_species = optional common species set (Phase 1.1). NULL keeps max data.
prep_structure <- function(var_name, data, tr,
                           human = CONFIG$human,
                           predictor_col = CONFIG$predictor_col,
                           restrict_species = CONFIG$restrict_species) {
  if (!var_name %in% names(data)) return(NULL)
  if (!predictor_col %in% names(data)) stop("predictor_col not found in data")
  
  d1 <- data[, c("Species", predictor_col, var_name)]
  d1 <- as.data.frame(na.omit(d1))
  
  if (!is.null(restrict_species)) {
    d1 <- d1[d1$Species %in% restrict_species, , drop = FALSE]
  }
  
  # Part-whole correction: Rest-of-Brain (ROB) = Total - Part
  d1$Rest_of_Brain <- d1[[predictor_col]] - d1[[var_name]]
  d1 <- d1[d1$Rest_of_Brain > 0, ]
  rownames(d1) <- d1$Species
  
  if (!human %in% d1$Species) return(NULL)
  
  obs_val <- d1[d1$Species == human, var_name][1]
  xh_raw  <- d1[d1$Species == human, "Rest_of_Brain"][1]
  
  # Align to tree
  clean <- clean.data(d1, tr)
  tr1   <- clean$tree
  data2 <- clean$data
  if (!is.null(data2$Species)) rownames(data2) <- data2$Species
  
  # Fit excludes human (original logic)
  tr2        <- drop.tip(tr1, human)
  data_model <- subset(data2, Species != human)
  
  fmla <- as.formula(paste0("log(", var_name, ") ~ log(Rest_of_Brain)"))
  
  list(
    var_name = var_name, human = human,
    d_full = d1, tr1 = tr1, tr2 = tr2,
    data2 = data2, data_model = data_model,
    obs_val = as.numeric(obs_val), xh_raw = as.numeric(xh_raw),
    fmla = fmla
  )
}

# Safe GLS fit (NULL on failure)
fit_gls_safe <- function(fmla, cor_struct, data_model) {
  tryCatch(
    gls(fmla, correlation = cor_struct, data = data_model,
        control = glsControl(opt = "optim", msMaxIter = 1000, msTol = 1e-6)),
    error = function(e) NULL
  )
}

# Extract lambda for a given model spec
lambda_from_fit <- function(fit, model_type, fixed_lambda = NA_real_) {
  switch(model_type,
         PagelML = as.numeric(coef(fit$modelStruct, unconstrained = FALSE)),
         Fixed   = as.numeric(fixed_lambda),
         BM      = 1.0,
         stop("Unknown model_type")
  )
}

# Phylogenetic conditional mean correction (mu) and conditional variance (ch)
phylo_mu_ch <- function(v_full, v_reduced, human, X_vec) {
  cc     <- which(rownames(v_full) == human)
  others <- rownames(v_reduced)
  Cih    <- v_full[cc, others, drop = FALSE]
  
  inv_v <- solve(v_reduced)
  xbar  <- mean(X_vec, na.rm = TRUE)
  
  mu <- Cih %*% inv_v %*% (X_vec - xbar)
  ch <- v_full[cc, cc] - Cih %*% inv_v %*% t(Cih)
  list(mu = as.numeric(mu), ch = as.numeric(ch))
}

# Run a set of models for one structure.
#   mu_source: "trait" (X = log trait; Plot I/II) or "predictor" (X = log ROB; Plot III)
run_models_for_structure <- function(pp, model_specs, want_ci = TRUE,
                                     mu_source = c("trait", "predictor")) {
  mu_source <- match.arg(mu_source)
  if (is.null(pp)) return(NULL)
  
  out <- list()
  for (spec in model_specs) {
    cor_struct <- spec$cor(pp$tr2)
    fit <- fit_gls_safe(pp$fmla, cor_struct, pp$data_model)
    if (is.null(fit)) next
    
    sigma      <- fit$sigma
    lambda_est <- lambda_from_fit(fit, spec$type, spec$lambda)
    
    v_full <- transform_vcv_lambda(pp$tr1, lambda_est)
    v_red  <- transform_vcv_lambda(pp$tr2, lambda_est)
    
    X  <- if (mu_source == "trait") log(pp$data_model[[pp$var_name]]) else log(pp$data_model$Rest_of_Brain)
    mc <- phylo_mu_ch(v_full, v_red, pp$human, X)
    mu <- mc$mu
    ch <- mc$ch
    if (!is.finite(ch) || ch <= 0) ch <- 1   # preserve original clamp
    
    xh       <- log(pp$xh_raw)
    pred_log <- as.numeric(c(1, xh) %*% coef(fit) + mu)
    
    base <- data.frame(
      Variable = pp$var_name, Model = spec$name, Observed = pp$obs_val,
      stringsAsFactors = FALSE
    )
    if (want_ci) {
      se <- sqrt(as.numeric(sigma^2 * ch))
      out[[spec$name]] <- cbind(base, data.frame(
        Lower = exp(pred_log - qnorm(0.975) * se),
        Predicted = exp(pred_log),
        Upper = exp(pred_log + qnorm(0.975) * se),
        lambda = lambda_est, N = nrow(pp$data_model)
      ))
    } else {
      out[[spec$name]] <- cbind(base, data.frame(
        Predicted = exp(pred_log), lambda = lambda_est, N = nrow(pp$data_model)
      ))
    }
  }
  if (length(out) == 0) return(NULL)
  do.call(rbind, out)
}

# Save a ggplot to png + pdf (+ optional jpg) in one call
save_plot <- function(plot, name, width = 7, height = 5, dpi = 300, jpg = FALSE,
                      dir = CONFIG$dir_figs) {
  ggsave(file.path(dir, paste0(name, ".png")), plot, width = width, height = height, dpi = dpi)
  ggsave(file.path(dir, paste0(name, ".pdf")), plot, width = width, height = height)
  if (jpg) {
    ggsave(file.path(dir, paste0(name, ".jpg")), plot, width = width, height = height,
           dpi = dpi, device = "jpeg", bg = "white")
  }
  invisible(plot)
}

# ============================================================
# MODEL SPECIFICATIONS (defined once, reused everywhere)
# ============================================================
spec_independence <- list(name = "Independence (λ=0)", type = "Fixed", lambda = 0,
                          cor = function(tr2) corPagel(0, form = ~Species, phy = tr2, fixed = TRUE))
spec_pagelML <- list(name = "Pagel's ML (Estimated)", type = "PagelML", lambda = NA_real_,
                     cor = function(tr2) corPagel(1, form = ~Species, phy = tr2, fixed = FALSE))
spec_BM <- list(name = "Brownian (λ=1)", type = "BM", lambda = NA_real_,
                cor = function(tr2) corBrownian(1, form = ~Species, phy = tr2))

# Plot I/II uses BM + ML but labels ML "Pagel's lambda (ML)"
spec_pagelML_plot12 <- modifyList(spec_pagelML, list(name = "Pagel's lambda (ML)"))
spec_BM_plot12      <- modifyList(spec_BM,      list(name = "Brownian (BM)"))

specs_plot12  <- list(spec_BM_plot12, spec_pagelML_plot12)
specs_plot3   <- list(spec_independence, spec_pagelML, spec_BM)
specs_compare <- specs_plot3   # 12-way robustness reuses the Plot III specs

# ============================================================
# PLOT I + II DATA (one pass, with CI; mu_source = "trait")
# ============================================================
# add a new column: mu
core_df <- map_dfr(target_cols, function(vn)
  run_models_for_structure(prep_structure(vn, data, tr), specs_plot12, TRUE, "trait") %>%
    mutate(mu = {
      pp <- prep_structure(vn, data, tr)
      if (is.null(pp)) return(NA_real_)
      v_full <- transform_vcv_lambda(pp$tr1, ifelse(Model == "Brownian (BM)", 1,
                                                    ifelse(grepl("Pagel", Model), lambda, 0)))
      v_red  <- transform_vcv_lambda(pp$tr2, ifelse(Model == "Brownian (BM)", 1,
                                                    ifelse(grepl("Pagel", Model), lambda, 0)))
      X <- log(pp$data_model[[vn]])
      phylo_mu_ch(v_full, v_red, pp$human, X)$mu
    })
)

# core_df <- map_dfr(target_cols, function(vn) {
#   run_models_for_structure(prep_structure(vn, data, tr), specs_plot12,
#                            want_ci = TRUE, mu_source = "trait")
# })

core_df$Variable <- factor(core_df$Variable, levels = target_cols)
core_df$VarLabel <- label_region(core_df$Variable)

# ---- PLOT I: predicted vs observed, BM vs Pagel ----
p1 <- ggplot(core_df, aes(x = VarLabel)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, color = "#377eb8") +
  geom_point(aes(y = Predicted, shape = "Predicted"), color = "#377eb8", size = 2.5) +
  geom_point(aes(y = Observed, shape = "Observed"), color = "#e41a1c", size = 2.5) +
  scale_shape_manual(name = "Value", values = c("Predicted" = 16, "Observed" = 17)) +
  coord_flip() +
  facet_wrap(~Model, scales = "free_x") +
  theme_bw() +
  labs(title = "Prediction: Brownian Motion vs Pagel's lambda (Part-Whole Corrected)",
       subtitle = "Predictor = Rest of Brain (Total - Structure)",
       y = "Volume / Value (Original Scale)", x = "Brain Structure") +
  theme(strip.text = element_text(face = "bold", size = 11),
        axis.text.y = element_text(size = 9))
p1

# ---- PLOT II: per-structure 0-1 normalized prediction error ----
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

p2 <- ggplot(final_df_norm, aes(x = VarLabel)) +
  geom_errorbar(aes(ymin = Lower_Sc, ymax = Upper_Sc), width = 0.2, color = "#377eb8") +
  geom_point(aes(y = Predicted_Sc, shape = "Predicted"), color = "#377eb8", size = 2.5) +
  geom_point(aes(y = Observed_Sc, shape = "Observed"), color = "#e41a1c", size = 2.5) +
  scale_shape_manual(name = "Value", values = c("Predicted" = 16, "Observed" = 17)) +
  coord_flip() +
  facet_wrap(~Model) +
  theme_bw() +
  labs(title = "Standardized Prediction Error (0-1 Scale, Part-Whole Corrected)",
       subtitle = "0 = min, 1 = max value per structure (across CI and Obs). Predictor = Rest of Brain",
       y = "Standardized Position (0 to 1)", x = "Brain Structure") +
  theme(strip.text = element_text(face = "bold", size = 11),
        axis.text.y = element_text(size = 9))
p2

# ============================================================
# PLOT III DATA (one pass, no CI; mu_source = "predictor")
# ============================================================
final_df <- map_dfr(target_cols, function(vn) {
  
  pp <- prep_structure(vn, data, tr)
  
  run_models_for_structure(pp, specs_plot3,
                           want_ci = FALSE, mu_source = "predictor") %>%
    mutate(mu = purrr::map2_dbl(Model, lambda, function(m, lam) {
      
      if (is.null(pp)) return(NA_real_)
      
      lam_use <- if (m == "Brownian (λ=1)") {
        1
      } else if (m == "Pagel's ML (Estimated)") {
        lam
      } else {
        0
      }
      
      v_full <- transform_vcv_lambda(pp$tr1, lam_use)
      v_red  <- transform_vcv_lambda(pp$tr2, lam_use)
      
      X <- log(pp$data_model$Rest_of_Brain)
      
      phylo_mu_ch(v_full, v_red, pp$human, X)$mu
    }))
})

check_mu <- final_df %>%
  filter(Model == "Independence (λ=0)") %>%
  select(Variable, mu)

print(check_mu)




# final_df <- map_dfr(target_cols, function(vn) {
#   run_models_for_structure(prep_structure(vn, data, tr), specs_plot3,
#                            want_ci = FALSE, mu_source = "predictor")
# })
final_df$PropDiff  <- (final_df$Observed - final_df$Predicted) / final_df$Predicted
final_df$Direction <- ifelse(final_df$PropDiff > 0, "Larger than Predicted", "Smaller than Predicted")
final_df$Model     <- factor(final_df$Model,
                             levels = c("Independence (λ=0)", "Pagel's ML (Estimated)", "Brownian (λ=1)"))
final_df$VarLabel  <- label_region(final_df$Variable)

p3 <- ggplot(final_df, aes(x = reorder(VarLabel, PropDiff), y = PropDiff, color = Direction)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_segment(aes(xend = reorder(VarLabel, PropDiff), yend = 0), size = 1.2) +
  geom_point(size = 4) +
  geom_text(aes(label = ifelse(Model == "Pagel's ML (Estimated)", sprintf("λ=%.2f", lambda), "")),
            nudge_x = -0.4, size = 3, fontface = "italic", color = "black", show.legend = FALSE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_color_manual(values = c("Larger than Predicted" = "#e41a1c",
                                "Smaller than Predicted" = "#377eb8")) +
  coord_flip() +
  facet_wrap(~Model) +
  theme_bw() +
  labs(title = "Human Brain Mosaicism (Corrected for Part-Whole)",
       subtitle = "Predictor = Rest of Brain (Total - Structure).",
       y = "Deviation from Prediction (%)", x = "Brain Structure", color = "Direction") +
  theme(strip.text = element_text(face = "bold", size = 10),
        legend.position = "bottom",
        axis.text.y = element_text(size = 9, face = "bold"))
p3

# ============================================================
# PLOT IV — lambda likelihood profiles
# ============================================================
lambda_seq <- seq(0, 1, by = 0.05)

profile_data <- map_dfr(target_cols, function(vn) {
  pp <- prep_structure(vn, data, tr)
  if (is.null(pp)) return(NULL)
  
  fit_mle <- fit_gls_safe(pp$fmla, corPagel(1, form = ~Species, phy = pp$tr2, fixed = FALSE), pp$data_model)
  mle_row <- NULL
  if (!is.null(fit_mle)) {
    best_lambda <- as.numeric(coef(fit_mle$modelStruct, unconstrained = FALSE))
    mle_row <- data.frame(Variable = vn,
                          lambda = pmax(0, pmin(1, best_lambda)),
                          LogLik = as.numeric(logLik(fit_mle)))
  }
  scan_df <- map_dfr(lambda_seq, function(val) {
    fit_scan <- fit_gls_safe(pp$fmla, corPagel(val, form = ~Species, phy = pp$tr2, fixed = TRUE), pp$data_model)
    if (is.null(fit_scan)) return(NULL)
    data.frame(Variable = vn, lambda = val, LogLik = as.numeric(logLik(fit_scan)))
  })
  list(scan = scan_df, mle = mle_row)
})
df_profile <- profile_data$scan
df_mle     <- profile_data$mle

p4 <- ggplot(df_profile, aes(x = lambda, y = LogLik)) +
  geom_line(color = "#377eb8", size = 1) +
  geom_point(data = df_mle, aes(x = lambda, y = LogLik), color = "red", size = 3) +
  geom_vline(data = df_mle, aes(xintercept = lambda), linetype = "dashed", color = "red", alpha = 0.5) +
  geom_text(data = df_mle, aes(label = sprintf("%.2f", lambda), x = 0.1, y = LogLik),
            color = "red", size = 3, hjust = 0, vjust = 1) +
  facet_wrap(~Variable, scales = "free_y",
             labeller = labeller(Variable = function(x) label_region(x))) +
  theme_bw() +
  labs(title = "lambda Likelihood Profiles (Part-Whole Corrected)",
       subtitle = "Predictor = Rest of Brain. Curve = fit at fixed lambda; red dot = MLE.",
       x = "lambda (Phylogenetic Signal)", y = "Log-Likelihood") +
  theme(strip.text = element_text(face = "bold", size = 9),
        axis.text.y = element_text(size = 7))
p4

# ---- Save Plots I-IV ----
save_plot(p1, "plot_prediction_BM_vs_pagel")
save_plot(p2, "plot_standardized_prediction_error")
save_plot(p3, "plot_human_brain_mosaicism")
save_plot(p4, "plot_pagel_lambda_profile")

# ============================================================
# PLOT V — individual BM PGLS plots per structure
# ============================================================
plot_region_gls <- function(region, data, tr, human = CONFIG$human) {
  pretty_region <- label_region(region)
  pp <- prep_structure(region, data, tr, human = human)
  if (is.null(pp)) { warning("Skipping ", region, ": prep_structure() returned NULL"); return(NULL) }
  
  fit <- fit_gls_safe(pp$fmla, corBrownian(1, form = ~Species, phy = pp$tr2), pp$data_model)
  if (is.null(fit)) { warning("Skipping ", region, ": fit_gls_safe() returned NULL"); return(NULL) }
  
  x_min_train <- min(pp$data_model$Rest_of_Brain, na.rm = TRUE)
  x_max_train <- max(pp$data_model$Rest_of_Brain, na.rm = TRUE)
  x_max_plot  <- max(pp$d_full$Rest_of_Brain,   na.rm = TRUE)
  
  x_grid_raw <- exp(seq(log(x_min_train), log(x_max_plot), length.out = 200))
  X <- cbind(1, log(x_grid_raw)); colnames(X) <- names(coef(fit))
  
  eta <- as.numeric(X %*% coef(fit))
  se  <- sqrt(diag(X %*% vcov(fit) %*% t(X)))
  z   <- qnorm(0.975)
  
  pred_df <- tibble(Rest_of_Brain = x_grid_raw, fit = exp(eta),
                    lwr = exp(eta - z * se), upr = exp(eta + z * se),
                    extrapolated = Rest_of_Brain > x_max_train)
  plot_df <- pp$d_full %>% mutate(is_human = Species == pp$human)
  
  ggplot() +
    geom_ribbon(data = filter(pred_df, !extrapolated), aes(Rest_of_Brain, ymin = lwr, ymax = upr), alpha = 0.2) +
    geom_ribbon(data = filter(pred_df,  extrapolated), aes(Rest_of_Brain, ymin = lwr, ymax = upr), alpha = 0.1) +
    geom_line(data = filter(pred_df, !extrapolated), aes(Rest_of_Brain, fit), linewidth = 1) +
    geom_line(data = filter(pred_df,  extrapolated), aes(Rest_of_Brain, fit), linewidth = 1, linetype = "dashed") +
    geom_point(data = filter(plot_df, !is_human), aes(Rest_of_Brain, .data[[region]]), size = 2, alpha = 0.8) +
    geom_point(data = filter(plot_df,  is_human), aes(Rest_of_Brain, .data[[region]]), size = 3, shape = 17) +
    scale_x_log10() + scale_y_log10() + theme_bw() +
    labs(title = paste0("Phylogenetic GLS (BM): ", pretty_region, " vs Rest of Brain"),
         subtitle = "Solid = non-human range; dashed = extrapolation. Human = triangle; fit excludes human.",
         x = "Rest of Brain (Total brain - structure)", y = paste0(pretty_region, " volume"))
}

plots <- target_cols %>%
  set_names(label_region(.)) %>%
  map(~ plot_region_gls(.x, data, tr)) %>%
  compact()

failed_regions <- setdiff(label_region(target_cols), names(plots))
if (length(failed_regions)) message("Region plots skipped: ", paste(failed_regions, collapse = ", "))

walk(plots, print)

# Multipage PDF
pdf(file.path(CONFIG$dir_figs, "phylo_gls_region_plots.pdf"), width = 7, height = 5)
walk(plots, print)
dev.off()

# One JPG per structure
iwalk(plots, function(p, region_label) {
  safe_region <- gsub("_+$", "", gsub("[^A-Za-z0-9]+", "_", region_label))
  ggsave(file.path(CONFIG$dir_figs, "phylo_gls_region_jpgs", paste0("phylo_gls_", safe_region, ".jpg")),
         plot = p, width = 13.333, height = 7.5, units = "in", dpi = 300, device = "jpeg", bg = "white")
})

# ============================================================
# TABLE I — predicted values + CIs, joined to rCMRGlc
# ============================================================
digits_lambda <- 3
digits_vals   <- 5

core_df_out <- core_df %>%
  mutate(
    Variable = factor(Variable, levels = target_cols),
    VarLabel = label_region(Variable),
    Diff.min = (Observed - Lower)     / Observed,
    Diff.pre = (Observed - Predicted) / Observed,
    Diff.max = (Observed - Upper)     / Observed
  ) %>%
  transmute(
    Structure = VarLabel, Model, lambda,
    `95% CI min` = Lower, Predicted, `95% CI max` = Upper, Observed,
    Diff.min, Diff.pre, Diff.max, N
  ) %>%
  mutate(Structure = factor(Structure, levels = label_region(target_cols))) %>%
  arrange(Structure, Model) %>%
  mutate(
    lambda       = signif(lambda,       digits_lambda),
    `95% CI min` = signif(`95% CI min`, digits_vals),
    Predicted    = signif(Predicted,    digits_vals),
    `95% CI max` = signif(`95% CI max`, digits_vals),
    Observed     = signif(Observed,     digits_vals),
    N            = as.integer(N)
  )

# rCMRGlc by Structure (volume_term) from Supplementary Table 1
sup_clean <- heiss_stephan_tbl %>%
  select(volume_term, rCMRGlc_mean_both_hemispheres) %>%
  rename(Structure = volume_term, rCMRGlc = rCMRGlc_mean_both_hemispheres) %>%
  mutate(rCMRGlc = round(as.numeric(rCMRGlc), 1))

core_with_rCMRGlc <- core_df_out %>%
  left_join(sup_clean %>% select(Structure, rCMRGlc), by = "Structure")

write_csv(core_with_rCMRGlc, file.path(CONFIG$dir_checks, "core_with_rCMRGlc_predicted_volumes.csv"))

# ---- Print-ready table (BM vs ML grouped columns) ----
metric_cols <- c("lambda", "95% CI min", "Predicted", "95% CI max",
                 "Observed", "Diff.min", "Diff.pre", "Diff.max")

model_wide <- core_with_rCMRGlc %>%
  mutate(
    ModelGroup = case_when(
      grepl("brownian|\\bbm\\b", tolower(Model)) ~ "Brownian (BM)",
      grepl("pagel|lambda|\\bml\\b", tolower(Model)) ~ "Pagel's lambda (ML)",
      TRUE ~ as.character(Model)
    ),
    ModelGroup = factor(ModelGroup, levels = c("Brownian (BM)", "Pagel's lambda (ML)"))
  ) %>%
  select(Structure, ModelGroup, all_of(metric_cols)) %>%
  pivot_wider(names_from = ModelGroup, values_from = all_of(metric_cols),
              names_glue = "{ModelGroup}__{.value}")

shared_cols <- core_with_rCMRGlc %>%
  group_by(Structure) %>%
  summarise(N = dplyr::first(N), rCMRGlc = dplyr::first(rCMRGlc), .groups = "drop")

bm_cols <- paste0("Brownian (BM)__", metric_cols)
ml_cols <- paste0("Pagel's lambda (ML)__", metric_cols)

tbl_body <- model_wide %>%
  left_join(shared_cols, by = "Structure") %>%
  arrange(Structure) %>%
  select(Structure, all_of(bm_cols), all_of(ml_cols), N, rCMRGlc)

tbl_body2 <- tbl_body %>% select(-any_of(c("Model", "ModelGroup")))

header_top <- c("Structure",
                "Brownian (BM)",       rep("", length(metric_cols) - 1),
                "Pagel's lambda (ML)", rep("", length(metric_cols) - 1),
                "N", "rCMRGlc")
header_sub <- c("", metric_cols, metric_cols, "", "")

print_ready_df <- as.data.frame(
  rbind(header_top, header_sub, as.matrix(tbl_body2)),
  stringsAsFactors = FALSE, check.names = FALSE
)
names(print_ready_df) <- rep("", ncol(print_ready_df))
write_xlsx(print_ready_df,
           file.path(CONFIG$dir_tables, "Table 4 Predicted volumes for Brownian and Pagel's lambda models.xlsx"),
           col_names = FALSE)

# ============================================================
# REPORT — correlation between rCMRGlc and prediction error across structures
# ============================================================
diff_cols <- c(paste0("Brownian (BM)__",       c("Diff.min", "Diff.pre", "Diff.max")),
               paste0("Pagel's lambda (ML)__", c("Diff.min", "Diff.pre", "Diff.max")))

cor_df <- tbl_body %>%
  transmute(Structure,
            rCMRGlc = parse_number(as.character(rCMRGlc)),
            across(all_of(diff_cols), ~ as.numeric(as.character(.x))))

cor_by_method <- function(method) {
  est_name <- if (method == "pearson") "r" else "rho"
  map_dfr(diff_cols, function(col) {
    ok <- complete.cases(cor_df[[col]], cor_df$rCMRGlc)
    ct <- cor.test(cor_df[[col]][ok], cor_df$rCMRGlc[ok], method = method)
    tibble(Diff_column = col, n = sum(ok),
           !!est_name := unname(ct$estimate), p_value = ct$p.value,
           ci_low = if (!is.null(ct$conf.int)) ct$conf.int[1] else NA_real_,
           ci_high = if (!is.null(ct$conf.int)) ct$conf.int[2] else NA_real_)
  })
}
cor_pearson  <- cor_by_method("pearson")
cor_spearman <- cor_by_method("spearman")

# ============================================================
# rCMRGlc vs prediction-error plots (fit excludes Neocortex white)
# ============================================================
num <- function(x) readr::parse_number(as.character(x))

df_bm <- core_with_rCMRGlc %>%
  filter(Model == "Brownian (BM)") %>%
  transmute(Structure, rCMRGlc = num(rCMRGlc),
            Diff.pre = num(Diff.pre), Diff.min = num(Diff.min), Diff.max = num(Diff.max),
            excluded = Structure == "Neocortex white") %>%
  drop_na(rCMRGlc, Diff.pre)
df_fit <- filter(df_bm, !excluded)
df_neo <- filter(df_bm,  excluded)

x_lab <- "rCMRGlc (µmol/100 g/min)"
y_lab <- "Difference from prediction (BM)"

theme_paper <- theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18), axis.text = element_text(size = 15),
        plot.title = element_text(size = 18, face = "bold"), plot.subtitle = element_text(size = 15),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        plot.margin = margin(8, 8, 8, 8))

add_points_and_labels <- list(
  geom_point(data = df_fit, aes(rCMRGlc, Diff.pre), inherit.aes = FALSE, pch = 16, size = 2.2),
  geom_text(data = df_fit, aes(rCMRGlc, Diff.pre, label = Structure), inherit.aes = FALSE,
            hjust = -0.08, size = 2.8, check_overlap = TRUE),
  geom_point(data = df_neo, aes(rCMRGlc, Diff.pre), inherit.aes = FALSE,
             shape = 24, fill = "firebrick2", color = "firebrick4", size = 4.5),
  geom_text(data = df_neo, aes(rCMRGlc, Diff.pre, label = "Neocortex white (excluded)"),
            inherit.aes = FALSE, hjust = -0.08, size = 3.1, color = "firebrick4")
)

# ---- LOESS with Diff.min/Diff.max error bars ----
p_loess_bars <- ggplot() +
  geom_errorbar(data = df_bm, aes(rCMRGlc, ymin = Diff.min, ymax = Diff.max),
                inherit.aes = FALSE, width = 0, na.rm = TRUE) +
  geom_smooth(data = df_fit, aes(rCMRGlc, Diff.pre), method = "loess", se = FALSE,
              span = 0.75, color = "steelblue4", linewidth = 1.2) +
  add_points_and_labels +
  coord_cartesian(xlim = c(10, 40), ylim = c(-3.5, 1.5), clip = "on") +
  labs(title = "LOESS fit with error bars",
       subtitle = "Fit excludes Neocortex white; plotted points include it", x = x_lab, y = y_lab) +
  theme_paper
p_loess_bars

# ---- Polynomial comparison, degrees 1-5 ----
xg <- seq(10, 40, length.out = 300)
fit_poly <- function(d) {
  # NB: keep the braces — a bare formula as an if-consequent before `else`
  # is a parse error in R (the `~` and `else` collide on precedence).
  form <- if (d == 1) {
    Diff.pre ~ rCMRGlc
  } else {
    as.formula(sprintf("Diff.pre ~ poly(rCMRGlc, %d, raw = TRUE)", d))
  }
  lm(form, data = df_fit)
}
fits        <- map(1:5, fit_poly)
adj_r2      <- map_dbl(fits, ~ summary(.x)$adj.r.squared)
best_degree <- which.max(adj_r2)
best_fit    <- fits[[best_degree]]
best_sum    <- summary(best_fit)
fstat       <- best_sum$fstatistic
model_p     <- pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)

fmt <- function(x) formatC(x, digits = 4, format = "f")
make_eqn <- function(coefs) {
  terms <- map_chr(seq_along(coefs), function(i) {
    if (i == 1) return(fmt(coefs[i]))
    pow <- i - 1
    paste0(ifelse(coefs[i] >= 0, " + ", " - "), fmt(abs(coefs[i])), "·x",
           ifelse(pow > 1, paste0("^", pow), ""))
  })
  paste0("y = ", paste0(terms, collapse = ""))
}
poly_stats <- paste(paste0("Best degree: ", best_degree),
                    paste0("Adj. R²: ", sprintf("%.4f", adj_r2[best_degree])),
                    paste0("Model p: ", format.pval(model_p, digits = 3, eps = 1e-4)),
                    make_eqn(coef(best_fit)), sep = "\n")
poly_curves <- map_dfr(1:5, function(d) {
  tibble(degree = d,
         degree_label = paste0("Degree ", d, "  Adj. R² = ", sprintf("%.3f", adj_r2[d])),
         rCMRGlc = xg, Diff.pre = predict(fits[[d]], newdata = data.frame(rCMRGlc = xg)))
})

p_poly_compare <- ggplot() +
  geom_line(data = poly_curves, aes(rCMRGlc, Diff.pre, color = degree_label), linewidth = 0.9) +
  geom_line(data = filter(poly_curves, degree == best_degree), aes(rCMRGlc, Diff.pre),
            linewidth = 1.2, linetype = "dashed", color = "black") +
  add_points_and_labels +
  annotate("text", x = 10.6, y = 1.35, label = poly_stats, hjust = 0, vjust = 1, size = 3.2, color = "navy") +
  coord_cartesian(xlim = c(10, 40), ylim = c(-3.5, 1.5), clip = "on") +
  labs(title = "Polynomial fits, degree 1 to 5",
       subtitle = paste0("Models exclude Neocortex white; best fit: degree ", best_degree),
       x = x_lab, y = y_lab, color = NULL) +
  theme_paper + theme(legend.position = "right")
p_poly_compare

# ---- Quadratic fit with 95% CI ----
quad_fit <- lm(Diff.pre ~ poly(rCMRGlc, 2, raw = TRUE), data = df_fit)
quad_sum <- summary(quad_fit)
quad_fstat <- quad_sum$fstatistic
quad_p <- pf(quad_fstat[1], quad_fstat[2], quad_fstat[3], lower.tail = FALSE)
b <- coef(quad_fit)
quad_eqn <- paste0("y = ", sprintf("%.4f", b[1]),
                   ifelse(b[2] >= 0, " + ", " - "), sprintf("%.4f", abs(b[2])), "x",
                   ifelse(b[3] >= 0, " + ", " - "), sprintf("%.4f", abs(b[3])), "x²")
quad_stats <- paste0("Adj. R² = ", sprintf("%.3f", quad_sum$adj.r.squared),
                     "\nModel p = ", format.pval(quad_p, digits = 2, eps = 1e-4), "\n", quad_eqn)
quad_pred <- data.frame(rCMRGlc = xg)
quad_pred <- cbind(quad_pred, as.data.frame(predict(quad_fit, newdata = quad_pred, interval = "confidence")))

p_quad_ci <- ggplot() +
  geom_ribbon(data = quad_pred, aes(rCMRGlc, ymin = lwr, ymax = upr), inherit.aes = FALSE,
              fill = "steelblue", alpha = 0.2) +
  geom_line(data = quad_pred, aes(rCMRGlc, fit), inherit.aes = FALSE, color = "steelblue4", linewidth = 1.2) +
  add_points_and_labels +
  annotate("text", x = 10.6, y = 1.35, label = quad_stats, hjust = 0, vjust = 1, size = 3.2, color = "navy") +
  coord_cartesian(xlim = c(10, 40), ylim = c(-3.5, 1.5), clip = "on") +
  labs(x = x_lab, y = y_lab) + theme_paper
p_quad_ci

save_plot(p_loess_bars,   "plot_loess_errorbars_BM_rCMRGlc")
save_plot(p_poly_compare, "plot_polynomial_comparison_BM_rCMRGlc", width = 9, height = 5.5)
save_plot(p_quad_ci,      "plot_quadratic_CI_BM_rCMRGlc")

# ============================================================
# PLOT VII + TABLE — 12-way quadratic robustness
#   4 structure configs x 3 PGLS evolution models.
#   mu_source = "trait" so the (Exclude Neocortex white, Brownian) cell
#   reproduces the published quad_fit above.
# ============================================================
structure_configs <- list(
  "All structures"            = character(0),
  "Exclude Neocortex white"   = "Neocortex white",
  "Exclude Capsula interna"   = "Capsula interna",
  "Exclude both white matter" = c("Neocortex white", "Capsula interna")
)
evolution_levels <- vapply(specs_compare, function(s) s$name, character(1))
structure_levels <- names(structure_configs)

# PGLS per (evolution model x structure); structure exclusion only affects the lm() below
diff_pre_long <- map_dfr(specs_compare, function(spec) {
  per_struct <- map_dfr(target_cols, function(vn) {
    run_models_for_structure(prep_structure(vn, data, tr), list(spec),
                             want_ci = FALSE, mu_source = "trait")
  })
  if (is.null(per_struct) || nrow(per_struct) == 0) return(NULL)
  per_struct %>%
    transmute(evolution_model = spec$name, Structure = label_region(Variable),
              Observed, Predicted, Diff.pre = (Observed - Predicted) / Observed)
})

diff_pre_with_rCMRGlc <- diff_pre_long %>%
  left_join(sup_clean %>% select(Structure, rCMRGlc), by = "Structure") %>%
  mutate(rCMRGlc = as.numeric(rCMRGlc)) %>%
  drop_na(rCMRGlc, Diff.pre)

# One quadratic per cell -> 12-row stats table
quad_grid <- map_dfr(structure_levels, function(sc_name) {
  excl <- structure_configs[[sc_name]]
  map_dfr(evolution_levels, function(ev_name) {
    sub    <- diff_pre_with_rCMRGlc %>% filter(evolution_model == ev_name) %>%
      mutate(excluded = Structure %in% excl)
    fit_df <- filter(sub, !excluded)
    n_used <- nrow(fit_df); n_excl <- sum(sub$excluded)
    
    if (n_used < 4) {
      return(tibble(structure_config = sc_name, evolution_model = ev_name,
                    n_used = n_used, n_excluded = n_excl,
                    intercept = NA_real_, b_linear = NA_real_, b_quadratic = NA_real_,
                    p_intercept = NA_real_, p_linear = NA_real_, p_quadratic = NA_real_,
                    adj_R2 = NA_real_, AIC = NA_real_, BIC = NA_real_, model_p_value = NA_real_))
    }
    fit   <- lm(Diff.pre ~ poly(rCMRGlc, 2, raw = TRUE), data = fit_df)
    summ  <- summary(fit); coefs <- coef(summ); fstat <- summ$fstatistic
    model_p <- if (!is.null(fstat)) unname(pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)) else NA_real_
    tibble(structure_config = sc_name, evolution_model = ev_name,
           n_used = n_used, n_excluded = n_excl,
           intercept = unname(coefs[1, "Estimate"]), b_linear = unname(coefs[2, "Estimate"]),
           b_quadratic = unname(coefs[3, "Estimate"]),
           p_intercept = unname(coefs[1, "Pr(>|t|)"]), p_linear = unname(coefs[2, "Pr(>|t|)"]),
           p_quadratic = unname(coefs[3, "Pr(>|t|)"]),
           adj_R2 = summ$adj.r.squared, AIC = AIC(fit), BIC = BIC(fit), model_p_value = model_p)
  })
})

# Prediction grid (curve + 95% CI) per cell
xg_compare <- seq(10, 40, length.out = 300)
quad_pred_grid <- map_dfr(structure_levels, function(sc_name) {
  excl <- structure_configs[[sc_name]]
  map_dfr(evolution_levels, function(ev_name) {
    fit_df <- diff_pre_with_rCMRGlc %>% filter(evolution_model == ev_name) %>%
      mutate(excluded = Structure %in% excl) %>% filter(!excluded)
    if (nrow(fit_df) < 4) return(NULL)
    fit  <- lm(Diff.pre ~ poly(rCMRGlc, 2, raw = TRUE), data = fit_df)
    pred <- as.data.frame(predict(fit, newdata = data.frame(rCMRGlc = xg_compare), interval = "confidence"))
    tibble(structure_config = sc_name, evolution_model = ev_name,
           rCMRGlc = xg_compare, fit = pred$fit, lwr = pred$lwr, upr = pred$upr)
  })
})

# Per-cell point data
points_grid <- map_dfr(structure_levels, function(sc_name) {
  diff_pre_with_rCMRGlc %>%
    mutate(structure_config = sc_name, excluded = Structure %in% structure_configs[[sc_name]])
})

# Lock facet ordering
points_grid$structure_config    <- factor(points_grid$structure_config,    levels = structure_levels)
points_grid$evolution_model     <- factor(points_grid$evolution_model,     levels = evolution_levels)
quad_pred_grid$structure_config <- factor(quad_pred_grid$structure_config, levels = structure_levels)
quad_pred_grid$evolution_model  <- factor(quad_pred_grid$evolution_model,  levels = evolution_levels)
quad_grid$structure_config      <- factor(quad_grid$structure_config,      levels = structure_levels)
quad_grid$evolution_model       <- factor(quad_grid$evolution_model,       levels = evolution_levels)

p_quad_compare_12way <- ggplot() +
  geom_ribbon(data = quad_pred_grid, aes(rCMRGlc, ymin = lwr, ymax = upr),
              inherit.aes = FALSE, fill = "steelblue", alpha = 0.2) +
  geom_line(data = quad_pred_grid, aes(rCMRGlc, fit), inherit.aes = FALSE,
            color = "steelblue4", linewidth = 1) +
  geom_point(data = filter(points_grid, !excluded), aes(rCMRGlc, Diff.pre),
             inherit.aes = FALSE, pch = 16, size = 1.8) +
  geom_point(data = filter(points_grid, excluded), aes(rCMRGlc, Diff.pre),
             inherit.aes = FALSE, shape = 24, fill = "firebrick2", color = "firebrick4", size = 3) +
  facet_grid(structure_config ~ evolution_model) +
  coord_cartesian(xlim = c(10, 40), ylim = c(-3.5, 1.5)) +
  labs(title = "Quadratic robustness: 4 structure configurations x 3 PGLS evolution models",
       subtitle = "Fit + 95% CI on included subset; excluded structures as red triangles",
       x = x_lab, y = y_lab) +
  theme_paper +
  theme(strip.text = element_text(face = "bold", size = 10),
        panel.spacing = unit(0.8, "lines"),
        plot.title = element_text(size = 14, face = "bold"), plot.subtitle = element_text(size = 11))
p_quad_compare_12way

save_plot(p_quad_compare_12way, "plot_quadratic_comparison_12way", width = 14, height = 12, jpg = TRUE)
write_xlsx(quad_grid %>% arrange(structure_config, evolution_model),
           file.path(CONFIG$dir_tables, "Table_S_quadratic_comparison_12way.xlsx"))

You
Sat 6/13/2026 9:03 AM
Prompt: I want to make an updated script based on s3_predicValuesPGLS_13062026_mucheck.R (Attached) It should have an additional structure processed with the others, that is called "Remainder". The volume of remainder should be calculated. It is equal to Preferred_brain_volume