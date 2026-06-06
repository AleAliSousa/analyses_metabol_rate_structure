# ============================================================
# Study 3 — Phase 1.1 SENSITIVITY: phylogenetic multiple imputation (MI)
#
# Question this answers: do the human prediction errors (Diff.pre) survive
# putting every region on the SAME, balanced species set, instead of each
# region using whatever subset happens to have data?
#
# Design:
#   * Engine: Rphylopars (multivariate Brownian model on the tree). It imputes
#     missing trait values by borrowing strength across (a) phylogeny and
#     (b) correlations with the other traits — the right tool for continuous
#     traits on a tree, unlike plain `mice` which assumes independent rows.
#   * Traits imputed on the LOG scale (the analysis scale).
#   * Homo sapiens is observed for all regions, so it is NEVER imputed — we only
#     impute non-human species to balance the comparison sample. The human value
#     being tested is untouched.
#   * Multiple imputation: draw m datasets from the per-tip predictive
#     distribution, refit the BM PGLS on each, and POOL with Rubin's rules so the
#     between-imputation variance (your "keep all the variance") is retained.
#
# CAVEATS (report these):
#   * Amygdala (~70% missing) and the Frahm neocortex pair (~60%) are imputed
#     heavily; their imputed slopes are model-driven. Imputation cannot create
#     information — treat the sparsest regions cautiously. `frac_imputed` per
#     region is reported so heavily-imputed cells can be discounted.
#   * Draws use a per-tip normal approximation (independent across tips/traits).
#     This is adequate for a sensitivity analysis; for the conditional
#     covariance among missing tips, see the note at the bottom.
#
# NOTE: This script could not be executed in the build sandbox (no R there).
#       Run locally and sanity-check `frac_imputed` and `pooled` before use.
# ============================================================

setwd("~/Library/CloudStorage/Dropbox/COLLABORATIVE/Do expensive brain regions increase less in humans/analyses_metabol_rate_structure")

library(tidyverse)
library(ape)
library(nlme)
library(Rphylopars)   # install.packages("Rphylopars")

set.seed(1)

# ------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------
MI <- list(
  human         = "Homo_sapiens",
  predictor_col = "Preferred_brain_volume",
  m             = 25,            # number of imputations
  model         = "BM",         # Rphylopars evolutionary model for imputation
  out_dir       = "checks/s3/phase1"
)
dir.create(MI$out_dir, showWarnings = FALSE, recursive = TRUE)

# Regions to put on a common footing (raw column -> label)
region_labels <- c(
  LGN_Sousa = "Corpus geniculatum laterale", Amygdala = "Amygdala", Pallidum = "Pallidum",
  NeoW_Frahm = "Neocortex white", Total_insula_volume_L = "Insular cortex (grey)",
  Nucleus_subthalamicus = "Nucleus subthalamicus Luysi", Capsula_interna = "Capsula interna",
  Striatum = "Striatum", ASG_Sousa = "Area striata grey", NeoG_Frahm = "Neocortex grey",
  Mesencephalon = "Mesencephalon", Cerebellum = "Cerebellum", Hippocampus = "Hippocampus"
)
target_cols <- names(region_labels)

# ------------------------------------------------------------
# DATA + TREE
# ------------------------------------------------------------
tr <- read.tree("data_raw/species.nwk")
Stephan <- read.csv("data_raw/Stephan_primates.csv") %>%
  mutate(Preferred_brain_volume = coalesce(Brain_volume, Brainvol, Total_brain_net_volume))

# Keep species present in both data and tree
keep <- intersect(Stephan$Species, tr$tip.label)
tr   <- drop.tip(tr, setdiff(tr$tip.label, keep))
dat  <- Stephan %>% filter(Species %in% keep)

# Log-scale trait matrix: regions + predictor as auxiliary (aids imputation)
log_cols <- c(target_cols, MI$predictor_col)
trait_log <- dat %>%
  transmute(species = Species,
            across(all_of(log_cols), ~ log(.x), .names = "{.col}")) %>%
  as.data.frame()

# Record which entries were originally missing (so we only impute those, and
# never touch observed values such as all of Homo sapiens)
missing_mask <- is.na(trait_log[, log_cols, drop = FALSE])
rownames(missing_mask) <- trait_log$species

frac_imputed <- tibble(
  raw_column = log_cols,
  Structure  = ifelse(log_cols %in% names(region_labels), region_labels[log_cols], log_cols),
  n_missing  = colSums(missing_mask),
  frac_imputed = round(colSums(missing_mask) / nrow(missing_mask), 3)
)
print(frac_imputed)

# ------------------------------------------------------------
# Fit the phylogenetic imputation model once
# ------------------------------------------------------------
pe <- phylopars(trait_data = trait_log, tree = tr, model = MI$model)

# Reconstructed means + variances for tips (missing entries get BLUPs).
# Guard: confirm the API returns species-named rows (versions differ); fail loudly if not.
stopifnot(all(trait_log$species %in% rownames(pe$anc_recon)),
          all(log_cols %in% colnames(pe$anc_recon)))
recon_mean <- pe$anc_recon[trait_log$species, log_cols, drop = FALSE]
recon_var  <- pe$anc_var [trait_log$species, log_cols, drop = FALSE]
recon_var[recon_var < 0] <- 0   # numerical guard

# ------------------------------------------------------------
# Compact BM PGLS: predict the human value for one region on one dataset.
# Mirrors the main script (BM, mu_source = "trait"); see s3_predicValuesPGLS_04062026.R.
# ------------------------------------------------------------
pgls_predict_human <- function(d, var_name, tree, human, predictor_col = MI$predictor_col) {
  d <- d[stats::complete.cases(d[, c(predictor_col, var_name)]), ]
  d$Rest_of_Brain <- d[[predictor_col]] - d[[var_name]]
  d <- d[d$Rest_of_Brain > 0, ]
  rownames(d) <- d$Species
  if (!human %in% d$Species) return(NULL)

  obs <- d[d$Species == human, var_name]
  xh  <- log(d[d$Species == human, "Rest_of_Brain"])

  tr_full <- drop.tip(tree, setdiff(tree$tip.label, d$Species))
  tr_red  <- drop.tip(tr_full, human)
  d_model <- d[d$Species != human, ]

  fmla <- as.formula(paste0("log(", var_name, ") ~ log(Rest_of_Brain)"))
  fit <- tryCatch(
    gls(fmla, correlation = corBrownian(1, form = ~Species, phy = tr_red), data = d_model,
        control = glsControl(opt = "optim", msMaxIter = 1000, msTol = 1e-6)),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)

  # Phylogenetic conditional-mean correction (mu), X = log(trait), as in main script
  v_full <- vcv(tr_full, corr = TRUE)
  v_red  <- vcv(tr_red,  corr = TRUE)
  cc     <- which(rownames(v_full) == human)
  others <- rownames(v_red)
  Cih    <- v_full[cc, others, drop = FALSE]
  X      <- log(d_model[[var_name]])
  mu     <- as.numeric(Cih %*% solve(v_red) %*% (X - mean(X)))

  pred_log <- as.numeric(c(1, xh) %*% coef(fit) + mu)
  tibble(Variable = var_name, slope = unname(coef(fit)[2]),
         Observed = obs, Predicted = exp(pred_log),
         Diff.pre = (obs - exp(pred_log)) / obs)
}

# ------------------------------------------------------------
# Generate m imputations, refit per region, collect estimates
# ------------------------------------------------------------
imp_results <- map_dfr(seq_len(MI$m), function(mi) {
  # Draw missing log-values; keep observed values exactly
  drawn <- trait_log
  for (col in log_cols) {
    miss <- which(missing_mask[, col])
    if (length(miss)) {
      drawn[miss, col] <- rnorm(length(miss), recon_mean[miss, col], sqrt(recon_var[miss, col]))
    }
  }
  # Back-transform to original scale; rebuild a data frame for PGLS
  d_imp <- dat
  for (col in log_cols) d_imp[[col]] <- exp(drawn[[col]])

  map_dfr(target_cols, function(vn) {
    res <- pgls_predict_human(d_imp, vn, tr, MI$human)
    if (!is.null(res)) mutate(res, imputation = mi)
    else NULL
  })
})

# ------------------------------------------------------------
# Pool with Rubin's rules (per region)
#   pooled mean = mean over imputations
#   total var   = within + (1 + 1/m) * between   (between-imputation variance is
#   the part that preserves "all the original variance" you wanted to keep)
# Each imputation here yields a point estimate, so we pool the point estimates
# and report the between-imputation SD as the imputation-driven uncertainty.
# ------------------------------------------------------------
pooled <- imp_results %>%
  group_by(Variable) %>%
  summarise(
    Structure = region_labels[dplyr::first(Variable)],
    slope_pooled    = mean(slope),    slope_between_sd    = sd(slope),
    Diff.pre_pooled = mean(Diff.pre), Diff.pre_between_sd = sd(Diff.pre),
    Diff.pre_total_se = sqrt((1 + 1/n()) * var(Diff.pre)),
    n_imputations = n(),
    .groups = "drop"
  ) %>%
  left_join(frac_imputed %>% select(raw_column, frac_imputed),
            by = c("Variable" = "raw_column")) %>%
  arrange(match(Variable, target_cols))

write_csv(imp_results, file.path(MI$out_dir, "s3_MI_per_imputation_estimates.csv"))
write_csv(pooled,      file.path(MI$out_dir, "s3_MI_pooled_Diff_pre.csv"))
write_csv(frac_imputed,file.path(MI$out_dir, "s3_MI_fraction_imputed.csv"))

message("Phylogenetic MI complete: ", MI$m, " imputations. Pooled results in ", MI$out_dir)
print(pooled)

# ------------------------------------------------------------
# Interpretation aid: compare pooled (balanced-sample) Diff.pre to the main
# max-data Diff.pre. If the sign and rough magnitude hold for the well-covered
# regions, the cross-region pattern is not an artefact of differential sampling.
# Heavily-imputed regions (high frac_imputed) should be flagged, not trusted.
# ------------------------------------------------------------

# --- OPTIONAL, more rigorous draw (conditional covariance among missing tips) ---
# The loop above draws each missing tip independently. To draw jointly (so
# imputed values respect the phylogenetic covariance among missing relatives),
# replace the per-column rnorm() with a multivariate draw per trait using the
# conditional covariance from the BM model:
#   Sigma_mm - Sigma_mo %*% solve(Sigma_oo) %*% Sigma_om   (per trait, scaled by
#   the fitted rate). This is heavier to code and rarely changes a sensitivity
#   conclusion, but is the technically complete version.
