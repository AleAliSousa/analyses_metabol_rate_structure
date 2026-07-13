# ============================================================
# Study 3 — Predicted vs observed brain-region volumes in Homo sapiens
# MERGED-DATA VARIANT (data-source / source-filter / clade switchable)
#
# Adapted from s3_predicValuesPGLS_16062026_..._PATCHED.R.
# The ENTIRE modelling body (HELPERS onward) is byte-identical to the
# original; ONLY the CONFIG + data-loading header below has changed.
#
# What this variant adds (all via CONFIG, no edits to the body):
#   * data_source   = "stephan" | "merged"
#       "stephan"  -> original Stephan_primates.csv behaviour (reproduces prior results)
#       "merged"   -> the compiled volumes_wide.csv / volumes_long.csv merge
#   * source_filter = "all" | "histological"   (merged only)
#       "histological" keeps only rows whose contributing Teams are ALL in
#       CONFIG$histological_teams (i.e. purely histological/stereological values),
#       EXCLUDING MRI (RillingInsel). See keep_human_always below.
#   * keep_human_always = TRUE
#       The human's merged cells are team-blends (e.g. total brain =
#       "Bush; RillingInsel; Stephan_collection; Zilles"). A strict histological
#       filter would delete the human from every region and Study 3 REQUIRES the
#       human, so the human row is retained regardless of the source filter.
#   * clade_restrict = "none" | "anthropoids"
#       "anthropoids" keeps only Platyrrhini + Cercopithecoidea + Hominoidea
#       (includes Homo). Use this for the frontal-lobe / vermis-type structures
#       that are only measured in anthropoids, and for all-primate vs anthropoid
#       human-prediction comparisons.
#
# Downstream code is unchanged because the merged region columns are RENAMED to
# the original Stephan raw column names (LGN_Sousa, NeoG_Frahm, ...), so
# region_labels / target_cols / the rCMRGlc join all match exactly as before.
#
# The three sensitivity configurations are driven by s3_run_and_compare_configs.R,
# which sets .s3_config_override and source()s this file once per configuration.
# ============================================================

setwd(local({ d <- normalizePath(getwd()); while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d); d }))  # repo root (portable; replaces hardcoded path -- see R/project_root.R)

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
CONFIG_default <- list(
  human         = "Homo_sapiens",
  predictor_col = "Preferred_brain_volume",

  # Phase 1.1 — restrict to a common species set (NULL = max-data per region)
  restrict_species = NULL,

  # ---- NEW switches ----
  data_source   = "merged",          # "stephan" | "merged"
  source_filter = "all",             # "all" | "histological"  (merged only)
  keep_human_always = TRUE,          # retain human even if its cells blend non-histological teams
  clade_restrict = "anthropoids",           # "none" | "anthropoids"
  restrict_to_tree = TRUE,           # keep only species on species.nwk (= primates); see note at the filter

  # ---- Extra regions beyond the core 14 (default: none) ----
  # Each entry: list(raw=<internal name>, merged_col=<volumes_wide column>,
  #                  label=<display label = rCMRGlc volume_term>, rCMRGlc=<value>)
  # When non-empty they are appended to merged_region_cols + region_labels, and a
  # matching rCMRGlc row is added so the body's volume_term join finds them.
  extra_regions = list(),
  extra_tag     = "ext",             # run_tag suffix used when extra_regions is non-empty

  # Teams counted as histological/stereological (Stephan-comparable).
  # RillingInsel (MRI) and Ashwell (marsupial series) are intentionally excluded.
  # NOTE: Bush/Sherwood/Stimpson/Barks are treated as histological/stereological here;
  #       edit this vector if you classify any of them differently.
  histological_teams = c("Stephan_collection", "Zilles", "Sherwood",
                         "Stimpson", "Barks", "Bush"),

  # Species-name reconciliation (merged spelling -> tree/Stephan spelling).
  # Pongo_sp. is a single merged row mapped to one tree tip; adjust if desired.
  species_crosswalk = c(
    "Lagothrix_lagothricha" = "Lagothrix_lagotricha",
    "Pongo_sp."             = "Pongo_pygmaeus"
  ),

  # I/O paths
  path_stephan     = "data_raw/Stephan_primates.csv",
  path_merged_wide = "data_intermediate/volumes_wide.csv",
  path_merged_long = "data_intermediate/volumes_long.csv",
  path_tree        = "data_raw/species.nwk",
  path_heiss       = "data_intermediate/Heiss_Stephan_data.csv"
)

# Allow s3_run_and_compare_configs.R to override selected fields before source().
if (exists(".s3_config_override") && is.list(.s3_config_override)) {
  CONFIG <- modifyList(CONFIG_default, .s3_config_override)
} else {
  CONFIG <- CONFIG_default
}

# Output run tag (keeps each configuration's figures/tables separate)
CONFIG$run_tag <- if (CONFIG$data_source == "stephan") "stephan" else
                    paste0("merged_", CONFIG$source_filter)
if (CONFIG$clade_restrict == "anthropoids")
  CONFIG$run_tag <- paste0(CONFIG$run_tag, "_anthro")
if (length(CONFIG$extra_regions) > 0)
  CONFIG$run_tag <- paste0(CONFIG$run_tag, "_", CONFIG$extra_tag)

CONFIG$dir_figs   <- file.path("figs",   "s3", CONFIG$run_tag)
CONFIG$dir_checks <- file.path("checks", "s3", CONFIG$run_tag)
CONFIG$dir_tables <- file.path("tables", "s3", CONFIG$run_tag)

message(sprintf("[Study 3] data_source=%s | source_filter=%s | clade=%s -> run_tag=%s",
                CONFIG$data_source, CONFIG$source_filter, CONFIG$clade_restrict, CONFIG$run_tag))

# Create output directories up front
for (d in c(CONFIG$dir_figs,
            file.path(CONFIG$dir_figs, "phylo_gls_region_jpgs"),
            CONFIG$dir_checks,
            CONFIG$dir_tables)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# ============================================================
# TAXONOMY — genus -> primate grade + the master primate genus list.
#   These genus vectors are the SINGLE source of truth for "is this a primate?"
#   They cover every primate genus present in volumes_wide.csv (verified against
#   all 192 merged genera), including 8 off-tree primates that grade_of() used to
#   miss (Callicebus, Carlito, Cebuella, Lemur, Mandrillus, Prolemur, Sapajus,
#   Semnopithecus). Primate-adjacent NON-primates are deliberately excluded:
#   Cynocephalus (colugo/Dermoptera), Tupaia & Urogale (treeshrews/Scandentia).
# ============================================================
STREPSIRRHINI    <- c("Avahi","Cheirogaleus","Daubentonia","Eulemur","Galago","Galagoides",
                      "Indri","Lemur","Lepilemur","Loris","Microcebus","Nycticebus","Otolemur",
                      "Perodicticus","Prolemur","Propithecus","Varecia")
TARSIIFORMES     <- c("Tarsius","Carlito")
PLATYRRHINI      <- c("Alouatta","Aotus","Ateles","Callicebus","Callimico","Callithrix",
                      "Cebuella","Cebus","Lagothrix","Leontopithecus","Pithecia",
                      "Plecturocebus","Saguinus","Saimiri","Sapajus")
CERCOPITHECOIDEA <- c("Cercocebus","Cercopithecus","Colobus","Erythrocebus","Lophocebus",
                      "Macaca","Mandrillus","Miopithecus","Nasalis","Papio","Piliocolobus",
                      "Pygathrix","Semnopithecus","Trachypithecus")
HOMINOIDEA       <- c("Gorilla","Homo","Hylobates","Pan","Pongo","Symphalangus")

PRIMATE_GENERA <- unique(c(STREPSIRRHINI, TARSIIFORMES, PLATYRRHINI,
                           CERCOPITHECOIDEA, HOMINOIDEA))

grade_of <- function(species) {
  genus <- sub("_.*$", "", species)
  dplyr::case_when(
    genus %in% TARSIIFORMES     ~ "Tarsiiformes",
    genus %in% STREPSIRRHINI    ~ "Strepsirrhini",
    genus %in% PLATYRRHINI      ~ "Platyrrhini",
    genus %in% CERCOPITHECOIDEA ~ "Cercopithecoidea",
    genus %in% HOMINOIDEA       ~ "Hominoidea",
    TRUE                        ~ "Unassigned"
  )
}
ANTHROPOID_GRADES <- c("Platyrrhini", "Cercopithecoidea", "Hominoidea")

# TRUE iff the species' genus is a primate genus (the authoritative primate test).
is_primate <- function(species) sub("_.*$", "", species) %in% PRIMATE_GENERA

# ============================================================
# REGION COLUMN MAP (merged column -> original Stephan raw name)
#   region_labels below is IDENTICAL to the original, so downstream code
#   (target_cols, label_region, rCMRGlc join) is unchanged. For merged data we
#   simply rename the merged volume columns to these raw names.
# ============================================================
merged_region_cols <- c(
  LGN_Sousa                 = "Corpus_geniculatum_laterale_Vol.mm3",
  Amygdala                  = "Amygdala_Vol.mm3",
  Pallidum                  = "Pallidum_Vol.mm3",
  NeoW_Frahm                = "Neocortex_white_matter_Vol.mm3",
  Total_insula_volume_L     = "Insula_left_Vol.mm3",   # LEFT insula: matches Stephan Total_insula_volume_L 1:1 (Insula_Vol.mm3 is bilateral = 2x)
  Nucleus_subthalamicus     = "Nucleus_subthalamicus_Vol.mm3",
  Capsula_interna           = "Capsula_interna_Vol.mm3",
  Striatum                  = "Striatum_Vol.mm3",
  ASG_Sousa                 = "Area_striata_grey_matter_Vol.mm3",
  NeoG_Frahm                = "Neocortex_grey_matter_Vol.mm3",
  Mesencephalon             = "Mesencephalon_Vol.mm3",
  Cerebellum                = "Cerebellum_Vol.mm3",
  Hippocampus               = "Hippocampus_Vol.mm3",
  Lateral_cerebellar_nuclei = "Lateral_cerebellar_nuclei_Vol.mm3"
)
merged_predictor_col_raw <- "Total_brain_net_volume_Vol.mm3"

# Append any extra regions' merged columns (so build_merged_data() pulls them)
for (er in CONFIG$extra_regions) merged_region_cols[er$raw] <- er$merged_col

# ============================================================
# LOAD
# ============================================================
heiss_stephan_tbl <- read.csv(CONFIG$path_heiss)
tr                <- read.tree(CONFIG$path_tree)

# Add rCMRGlc rows for extra regions so the body's volume_term join finds them.
# (Only volume_term + rCMRGlc_mean_both_hemispheres are used by sup_clean; other
#  columns are set NA.)
if (length(CONFIG$extra_regions) > 0) {
  .tmpl <- heiss_stephan_tbl[1, , drop = FALSE]; .tmpl[] <- NA
  for (er in CONFIG$extra_regions) {
    .r <- .tmpl
    .r$volume_term                    <- er$label
    .r$rCMRGlc_mean_both_hemispheres  <- er$rCMRGlc
    heiss_stephan_tbl <- rbind(heiss_stephan_tbl, .r)
  }
  rm(.tmpl)
}

apply_crosswalk <- function(sp) {
  xw <- CONFIG$species_crosswalk
  hit <- sp %in% names(xw)
  sp[hit] <- unname(xw[sp[hit]])
  sp
}

# ---- Stephan loader (reproduces original behaviour) ----
build_stephan_data <- function() {
  Stephan_primates <- read.csv(CONFIG$path_stephan)
  dc <- Stephan_primates[, !names(Stephan_primates) %in% c("X", "order")]
  dc <- dc[!is.na(dc$Species) & trimws(dc$Species) != "", ]
  dc$Species <- gsub(" ", "_", dc$Species)
  dc <- dc[is_primate(dc$Species), , drop = FALSE]   # PRIMATE CROP at source
  dc %>% mutate(
    Preferred_brain_volume = coalesce(Brain_volume, Brainvol, Total_brain_net_volume)
  )
}

# ---- Merged loader (all sources, or histological-only via the long file) ----
build_merged_data <- function() {
  if (CONFIG$source_filter == "histological") {
    long <- read.csv(CONFIG$path_merged_long, check.names = FALSE, stringsAsFactors = FALSE)
    long$Species <- apply_crosswalk(gsub(" ", "_", long$Species))
    # keep rows whose contributing teams are ALL histological
    is_hist <- vapply(long$Teams, function(ts) {
      parts <- trimws(unlist(strsplit(as.character(ts), "[;,]")))
      length(parts) > 0 && all(parts %in% CONFIG$histological_teams)
    }, logical(1))
    keep <- is_hist
    if (isTRUE(CONFIG$keep_human_always)) keep <- keep | (long$Species == CONFIG$human)
    long <- long[keep, , drop = FALSE]
    wide <- long %>%
      dplyr::select(Species, Variable, Value) %>%
      dplyr::distinct(Species, Variable, .keep_all = TRUE) %>%
      tidyr::pivot_wider(names_from = Variable, values_from = Value)
    wide <- as.data.frame(wide, check.names = FALSE)
  } else {
    wide <- read.csv(CONFIG$path_merged_wide, check.names = FALSE, stringsAsFactors = FALSE)
    wide$Species <- apply_crosswalk(gsub(" ", "_", wide$Species))
  }

  # PRIMATE CROP at source — drop every non-primate row before building `data`.
  wide <- wide[is_primate(wide$Species), , drop = FALSE]

  out <- data.frame(Species = wide$Species, stringsAsFactors = FALSE)
  num <- function(z) suppressWarnings(as.numeric(as.character(z)))
  out$Preferred_brain_volume <- if (merged_predictor_col_raw %in% names(wide))
    num(wide[[merged_predictor_col_raw]]) else NA_real_
  for (raw in names(merged_region_cols)) {
    mc <- merged_region_cols[[raw]]
    out[[raw]] <- if (mc %in% names(wide)) num(wide[[mc]]) else NA_real_
  }
  out
}

data <- if (CONFIG$data_source == "stephan") build_stephan_data() else build_merged_data()

# ============================================================
# SPECIES CROP — done HERE, at the first step, before regions/labels/plots are
# defined and before `data` is used anywhere. Everything downstream (prep_structure,
# d_full, every point on every plot) is built from this cropped `data`, so nothing
# non-primate can reach a figure. Three sequential crops, each reported:
#   (1) PRIMATE crop  — genus must be a primate genus (already enforced in the
#       loaders; re-applied here as a guard).
#   (2) TREE crop     — keep species on species.nwk (the analysable primate set).
#   (3) ANTHROPOID crop (optional) — Platyrrhini + Cercopithecoidea + Hominoidea,
#       human always retained.
# A hard stopifnot() at the end aborts the run if any non-primate survived.
# ============================================================
message(sprintf("[crop] loaded %d species from %s.", nrow(data), CONFIG$data_source))

# (1) PRIMATE crop
.nonprimate <- unique(data$Species[!is_primate(data$Species)])
if (length(.nonprimate) > 0) {
  message(sprintf("[crop] (1) primate crop: dropping %d non-primate species: %s%s",
                  length(.nonprimate), paste(head(.nonprimate, 12), collapse = ", "),
                  if (length(.nonprimate) > 12) ", ..." else ""))
  data <- data[is_primate(data$Species), , drop = FALSE]
} else {
  message("[crop] (1) primate crop: no non-primates present.")
}

# (2) TREE crop
if (isTRUE(CONFIG$restrict_to_tree)) {
  .offtree <- setdiff(data$Species, tr$tip.label)
  data <- data[data$Species %in% tr$tip.label, , drop = FALSE]
  message(sprintf("[crop] (2) tree crop: kept %d species on species.nwk; dropped %d off-tree primate(s)%s.",
                  nrow(data), length(.offtree),
                  if (length(.offtree)) paste0(": ", paste(head(.offtree, 12), collapse = ", ")) else ""))
}

# (3) ANTHROPOID crop (optional; human always retained)
if (CONFIG$clade_restrict == "anthropoids") {
  keep_row <- grade_of(data$Species) %in% ANTHROPOID_GRADES | data$Species == CONFIG$human
  data <- data[keep_row, , drop = FALSE]
  message(sprintf("[crop] (3) anthropoid crop: kept %d species (Platyrrhini+Cercopithecoidea+Hominoidea + human).",
                  nrow(data)))
}

# ---- CONFIRM: hard stop if anything non-primate slipped through ----
.bad <- unique(data$Species[!is_primate(data$Species)])
if (length(.bad) > 0) stop("Non-primate species survived the crop: ", paste(.bad, collapse = ", "))
if (isTRUE(CONFIG$restrict_to_tree) && !all(data$Species %in% tr$tip.label))
  stop("Species off the tree survived the tree crop.")
.grade_tab <- table(grade_of(data$Species))
message(sprintf("[crop] CONFIRMED primates only: %d species | %s",
                nrow(data),
                paste(names(.grade_tab), as.integer(.grade_tab), sep = ":", collapse = ", ")))

# Sanity check on the predictor
n_missing_pred <- sum(is.na(data$Preferred_brain_volume))
if (n_missing_pred > 0) {
  warning(n_missing_pred, " species lack a predictor (Preferred_brain_volume) ",
          "and will drop out of each fit.")
}

# ============================================================
# REGIONS — raw column -> display label, in plotting/table order.
#   IDENTICAL to the original so the rCMRGlc join (volume_term) still matches.
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
  "NeoG_Frahm"    = "Neocortex grey",
  "Mesencephalon"         = "Mesencephalon",
  "Cerebellum"            = "Cerebellum",
  "Hippocampus"           = "Hippocampus",
  "Lateral_cerebellar_nuclei" = "Nucleus dentatus cerebelli"
)

# Append any extra regions (label used here must equal the volume_term added to
# heiss_stephan_tbl above, so the rCMRGlc join matches).
for (er in CONFIG$extra_regions) region_labels[er$raw] <- er$label

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
#   mu_source: "trait" (X = log trait; Plot I/II) or "predictor" (X = log ROB).
#   add_phylo_mu = TRUE returns the conditional phylogenetic prediction:
#     log(predicted y) = beta0 + beta1 * log(human ROB) + mu.
#     This is the prediction used by the mosaicism plots below.
#   add_phylo_mu = FALSE returns the ordinary fitted GLS line prediction only:
#     log(predicted y) = beta0 + beta1 * log(human ROB).
run_models_for_structure <- function(pp, model_specs, want_ci = TRUE,
                                      mu_source = c("trait", "predictor"),
                                      add_phylo_mu = TRUE) {
  mu_source <- match.arg(mu_source)
  if (is.null(pp)) return(NULL)

  out <- list()
  for (spec in model_specs) {
    cor_struct <- spec$cor(pp$tr2)
    fit <- fit_gls_safe(pp$fmla, cor_struct, pp$data_model)
    if (is.null(fit)) next

    sigma      <- fit$sigma
    lambda_est <- lambda_from_fit(fit, spec$type, spec$lambda)

    mu <- 0
    ch <- 1
    if (isTRUE(add_phylo_mu)) {
      v_full <- transform_vcv_lambda(pp$tr1, lambda_est)
      v_red  <- transform_vcv_lambda(pp$tr2, lambda_est)

      X  <- if (mu_source == "trait") log(pp$data_model[[pp$var_name]]) else log(pp$data_model$Rest_of_Brain)
      mc <- phylo_mu_ch(v_full, v_red, pp$human, X)
      mu <- mc$mu
      ch <- mc$ch
      if (!is.finite(ch) || ch <= 0) ch <- 1   # preserve original clamp
    }

    xh            <- log(pp$xh_raw)
    fitted_log    <- as.numeric(c(1, xh) %*% coef(fit))
    pred_log      <- fitted_log + mu

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
        Predicted_log = pred_log,
        FittedLine_log = fitted_log,
        PhyloMu_log = mu,
        lambda = lambda_est, N = nrow(pp$data_model)
      ))
    } else {
      out[[spec$name]] <- cbind(base, data.frame(
        Predicted = exp(pred_log),
        Predicted_log = pred_log,
        FittedLine_log = fitted_log,
        PhyloMu_log = mu,
        lambda = lambda_est, N = nrow(pp$data_model)
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
core_df <- map_dfr(target_cols, function(vn) {
  run_models_for_structure(prep_structure(vn, data, tr), specs_plot12,
                           want_ci = TRUE, mu_source = "trait")
})
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
# PLOT III DATA (one pass, no CI; conditional phylogenetic prediction)
# ============================================================
# These mosaicism plots use the held-out-human prediction equation with the
# phylogenetic conditional mean correction:
#   log(predicted y) = beta0 + beta1 * log(human Rest_of_Brain) + mu.
# The BM PGLS region plots below are drawn with the same +mu shift, so the
# human point's above/below-line status matches the BM mosaicism sign.
# Models are fit with natural logs via log(y) ~ log(Rest_of_Brain); plot axes
# are displayed as log10 values.
final_df <- map_dfr(target_cols, function(vn) {
  run_models_for_structure(prep_structure(vn, data, tr), specs_plot3,
                           want_ci = FALSE, mu_source = "predictor",
                           add_phylo_mu = TRUE)
}) %>%
  mutate(
    PropDiff        = (Observed - Predicted) / Predicted,
    Observed_log10  = log10(Observed),
    Predicted_log10 = log10(Predicted),
    Log10Diff       = Observed_log10 - Predicted_log10,
    Direction       = ifelse(Log10Diff > 0, "Larger than Predicted", "Smaller than Predicted"),
    Model           = factor(Model,
                             levels = c("Independence (λ=0)", "Pagel's ML (Estimated)", "Brownian (λ=1)")),
    VarLabel        = label_region(Variable)
  )

plot_mosaicism <- function(df, title, subtitle) {
  df <- df %>% mutate(Model = fct_drop(Model))
  lambda_label_df <- df %>% filter(as.character(Model) == "Pagel's ML (Estimated)")

  ggplot(df, aes(x = reorder(VarLabel, PropDiff), y = PropDiff, color = Direction)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    geom_segment(aes(xend = reorder(VarLabel, PropDiff), yend = 0), size = 1.2) +
    geom_point(size = 4) +
    geom_text(data = lambda_label_df,
              aes(label = sprintf("λ=%.2f", lambda)),
              nudge_x = -0.4, size = 3, fontface = "italic",
              color = "black", show.legend = FALSE) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_color_manual(values = c("Larger than Predicted" = "#e41a1c",
                                  "Smaller than Predicted" = "#377eb8")) +
    coord_flip() +
    facet_wrap(~Model) +
    theme_bw() +
    labs(title = title, subtitle = subtitle,
         y = "Deviation from conditional prediction (%)", x = "Brain Structure", color = "Direction") +
    theme(strip.text = element_text(face = "bold", size = 10),
          legend.position = "bottom",
          axis.text.y = element_text(size = 9, face = "bold"))
}

p3 <- plot_mosaicism(
  final_df,
  title = "Human Brain Mosaicism (Corrected for Part-Whole)",
  subtitle = "Predictor = Rest of Brain (Total - Structure). Predictions include the human-specific phylogenetic conditional mean correction (mu)."
)
p3

# Requested additional versions of the Human Brain Mosaicism plot
p3_pagel_bm <- plot_mosaicism(
  final_df %>% filter(as.character(Model) %in% c("Pagel's ML (Estimated)", "Brownian (λ=1)")),
  title = "Human Brain Mosaicism (Corrected for Part-Whole): Pagel's ML and Brownian",
  subtitle = "Models shown = Pagel's ML and Brownian. Predictions include the human-specific phylogenetic conditional mean correction (mu)."
)
p3_pagel_bm

p3_bm_only <- plot_mosaicism(
  final_df %>% filter(as.character(Model) == "Brownian (λ=1)"),
  title = "Human Brain Mosaicism (Corrected for Part-Whole): Brownian Only",
  subtitle = "Model shown = Brownian. Prediction includes the same BM human-specific phylogenetic mu used in the BM PGLS plots."
)
p3_bm_only

# Combined version of the three requested mosaicism views, kept as one figure.
mosaicism_versions_df <- bind_rows(
  final_df %>% mutate(MosaicismVersion = "All models"),
  final_df %>%
    filter(as.character(Model) %in% c("Pagel's ML (Estimated)", "Brownian (λ=1)")) %>%
    mutate(MosaicismVersion = "Pagel ML + Brownian"),
  final_df %>%
    filter(as.character(Model) == "Brownian (λ=1)") %>%
    mutate(MosaicismVersion = "Brownian only")
) %>%
  mutate(
    MosaicismVersion = factor(MosaicismVersion,
                              levels = c("All models", "Pagel ML + Brownian", "Brownian only")),
    MosaicismPanel = factor(
      paste(MosaicismVersion, Model, sep = " | "),
      levels = c(
        "All models | Independence (λ=0)",
        "All models | Pagel's ML (Estimated)",
        "All models | Brownian (λ=1)",
        "Pagel ML + Brownian | Pagel's ML (Estimated)",
        "Pagel ML + Brownian | Brownian (λ=1)",
        "Brownian only | Brownian (λ=1)"
      )
    )
  )

# p3_requested_versions <- ggplot(
#   mosaicism_versions_df,
#   aes(x = reorder(VarLabel, PropDiff), y = PropDiff, color = Direction)
# ) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
#   geom_segment(aes(xend = reorder(VarLabel, PropDiff), yend = 0), size = 1.0) +
#   geom_point(size = 3) +
#   geom_text(data = filter(mosaicism_versions_df, as.character(Model) == "Pagel's ML (Estimated)"),
#             aes(label = sprintf("λ=%.2f", lambda)),
#             nudge_x = -0.4, size = 2.6, fontface = "italic",
#             color = "black", show.legend = FALSE) +
#   scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
#   scale_color_manual(values = c("Larger than Predicted" = "#e41a1c",
#                                 "Smaller than Predicted" = "#377eb8")) +
#   coord_flip() +
#   facet_wrap(~MosaicismPanel, ncol = 3) +
#   theme_bw() +
#   labs(title = "Human Brain Mosaicism (Corrected for Part-Whole): All Requested Views",
#        subtitle = "Includes the original all-model view, Pagel's ML + Brownian, and Brownian-only view; predictions include the human-specific phylogenetic conditional mean correction (mu).",
#        y = "Deviation from conditional prediction (%)", x = "Brain Structure", color = "Direction") +
#   theme(strip.text = element_text(face = "bold", size = 8),
#         legend.position = "bottom",
#         axis.text.y = element_text(size = 7, face = "bold"))
# p3_requested_versions

# Diagnostic table: for BM, this sign is now the same sign you see by comparing
# the human triangle with the +mu conditional prediction line in the BM PGLS plots.
bm_mosaicism_mu_alignment_check <- final_df %>%
  filter(as.character(Model) == "Brownian (λ=1)") %>%
  transmute(
    Structure = VarLabel,
    Observed,
    Predicted_from_BM_mu_adjusted_equation = Predicted,
    Observed_log10,
    Predicted_log10,
    Log10Diff,
    PropDiff,
    Direction,
    FittedLine_log,
    PhyloMu_log,
    Predicted_log,
    lambda,
    N
  )
write.csv(bm_mosaicism_mu_alignment_check,
          file.path(CONFIG$dir_checks, "bm_mosaicism_mu_adjusted_alignment_check.csv"),
          row.names = FALSE)

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

# ---- Save Plots I-IV plus requested mosaicism variants ----
save_plot(p1, "plot_prediction_BM_vs_pagel")
save_plot(p2, "plot_standardized_prediction_error")
save_plot(p3, "plot_human_brain_mosaicism")
save_plot(p3_pagel_bm, "plot_human_brain_mosaicism_pagelML_BM")
save_plot(p3_bm_only, "plot_human_brain_mosaicism_BM_only")
#save_plot(p3_requested_versions, "plot_human_brain_mosaicism_all_requested_views", width = 14, height = 8)
save_plot(p4, "plot_pagel_lambda_profile")

# ============================================================
# PLOT V - individual BM PGLS plots per structure + combined faceted plot
#          These plots use the same BM + human-specific mu equation as the
#          Brownian panel in the Human Brain Mosaicism plot.
# ============================================================
build_region_gls_plot_data <- function(region, data, tr, human = CONFIG$human) {
  pretty_region <- label_region(region)
  pp <- prep_structure(region, data, tr, human = human)
  if (is.null(pp)) { warning("Skipping ", region, ": prep_structure() returned NULL"); return(NULL) }

  fit <- fit_gls_safe(pp$fmla, corBrownian(1, form = ~Species, phy = pp$tr2), pp$data_model)
  if (is.null(fit)) { warning("Skipping ", region, ": fit_gls_safe() returned NULL"); return(NULL) }

  # Brownian model uses lambda = 1.  This is the same human-specific mu term
  # used by the Brownian rows in the Human Brain Mosaicism plot.
  lambda_est <- 1.0
  v_full <- transform_vcv_lambda(pp$tr1, lambda_est)
  v_red  <- transform_vcv_lambda(pp$tr2, lambda_est)
  mc <- phylo_mu_ch(v_full, v_red, pp$human, log(pp$data_model$Rest_of_Brain))
  mu <- mc$mu
  ch <- mc$ch
  if (!is.finite(ch) || ch <= 0) ch <- 1

  x_min_train <- min(pp$data_model$Rest_of_Brain, na.rm = TRUE)
  x_max_train <- max(pp$data_model$Rest_of_Brain, na.rm = TRUE)
  x_max_plot  <- max(pp$d_full$Rest_of_Brain,   na.rm = TRUE)

  x_grid_raw <- exp(seq(log(x_min_train), log(x_max_plot), length.out = 200))
  X <- cbind(1, log(x_grid_raw)); colnames(X) <- names(coef(fit))

  eta_fitted <- as.numeric(X %*% coef(fit))
  eta        <- eta_fitted + mu
  se         <- sqrt(diag(X %*% vcov(fit) %*% t(X)))
  z          <- qnorm(0.975)

  # Equation of the plotted BM + human-mu line on the log10 axes.
  # Original fitted model is: ln(region) = b0 + b1*ln(ROB).
  # The plotted line adds mu, then divides by ln(10):
  # log10(region) = (b0 + mu)/ln(10) + b1*log10(ROB).
  gls_coef <- coef(fit)
  line_intercept_log10 <- unname((gls_coef[1] + mu) / log(10))
  line_slope_log10 <- unname(gls_coef[2])

  fmt_eq_num <- function(z, digits = 3) {
    formatC(z, format = "f", digits = digits)
  }

  line_eqn <- paste0(
    "log10(", pretty_region, " volume) = ",
    fmt_eq_num(line_intercept_log10),
    ifelse(line_slope_log10 >= 0, " + ", " - "),
    fmt_eq_num(abs(line_slope_log10)),
    " * log10(Rest of Brain)"
  )

  pred_df <- tibble(Rest_of_Brain = x_grid_raw,
                    fit = exp(eta),
                    fit_unadjusted = exp(eta_fitted),
                    lwr = exp(eta - z * se),
                    upr = exp(eta + z * se),
                    extrapolated = Rest_of_Brain > x_max_train) %>%
    mutate(Variable = region,
           VarLabel = pretty_region,
           PhyloMu_log = mu,
           x_log10 = log10(Rest_of_Brain),
           fit_log10 = eta / log(10),
           fit_unadjusted_log10 = eta_fitted / log(10),
           lwr_log10 = (eta - z * se) / log(10),
           upr_log10 = (eta + z * se) / log(10))

  human_eta_fitted <- as.numeric(c(1, log(pp$xh_raw)) %*% coef(fit))
  human_eta        <- human_eta_fitted + mu
  human_pred_df <- tibble(
    Variable = region,
    VarLabel = pretty_region,
    Rest_of_Brain = pp$xh_raw,
    x_log10 = log10(pp$xh_raw),
    pred_log10 = human_eta / log(10),
    fitted_line_log10 = human_eta_fitted / log(10),
    PhyloMu_log = mu
  )

  plot_df <- pp$d_full %>%
    filter(Rest_of_Brain > 0, .data[[region]] > 0) %>%
    mutate(Variable = region,
           VarLabel = pretty_region,
           is_human = Species == pp$human,
           x_log10 = log10(Rest_of_Brain),
           y_log10 = log10(.data[[region]]))

  list(
    region = region,
    pretty_region = pretty_region,
    pred_df = pred_df,
    plot_df = plot_df,
    human_pred_df = human_pred_df,
    PhyloMu_log = mu,
    line_intercept_log10 = line_intercept_log10,
    line_slope_log10 = line_slope_log10,
    line_eqn = line_eqn
  )
}

make_region_gls_plot <- function(plot_data) {
  pred_df <- plot_data$pred_df
  plot_df <- plot_data$plot_df
  human_pred_df <- plot_data$human_pred_df
  pretty_region <- plot_data$pretty_region

  eqn_label <- plot_data$line_eqn
  if (is.null(eqn_label) || length(eqn_label) == 0 || is.na(eqn_label)) {
    eqn_label <- paste0("log10(", pretty_region, " volume) = fitted BM + human-mu line")
  }

  ggplot() +
    geom_ribbon(
      data = filter(pred_df, !extrapolated),
      aes(x = x_log10, ymin = lwr_log10, ymax = upr_log10),
      alpha = 0.2
    ) +
    geom_ribbon(
      data = filter(pred_df, extrapolated),
      aes(x = x_log10, ymin = lwr_log10, ymax = upr_log10),
      alpha = 0.1
    ) +
    geom_line(
      data = filter(pred_df, !extrapolated),
      aes(x = x_log10, y = fit_log10),
      linewidth = 1
    ) +
    geom_line(
      data = filter(pred_df, extrapolated),
      aes(x = x_log10, y = fit_log10),
      linewidth = 1,
      linetype = "dashed"
    ) +
    geom_point(
      data = filter(plot_df, !is_human),
      aes(x = x_log10, y = y_log10),
      size = 2,
      alpha = 0.8
    ) +
    geom_point(
      data = filter(plot_df, is_human),
      aes(x = x_log10, y = y_log10),
      size = 3,
      shape = 17
    ) +
    # Uncomment this if you want to show the predicted human value as an x.
    # geom_point(
    #   data = human_pred_df,
    #   aes(x = x_log10, y = pred_log10),
    #   size = 3,
    #   shape = 4,
    #   stroke = 1.1
    # ) +
    theme_bw() +
    labs(
      title = paste0(
        "Phylogenetic GLS (BM + human mu): ",
        pretty_region,
        " vs Rest of Brain"
      ),
      subtitle = paste0(
        eqn_label,
        "\nHuman = triangle; fit excludes human; dashed line = extrapolated segment."
      ),
      x = "log10(Rest of Brain [Total brain - structure])",
      y = paste0("log10(", pretty_region, " volume)")
    )
}

plot_region_gls <- function(region, data, tr, human = CONFIG$human) {
  plot_data <- build_region_gls_plot_data(region, data, tr, human = human)
  if (is.null(plot_data)) return(NULL)
  make_region_gls_plot(plot_data)
}

gls_plot_data <- target_cols %>%
  set_names(label_region(.)) %>%
  map(~ build_region_gls_plot_data(.x, data, tr)) %>%
  compact()

plots <- gls_plot_data %>% map(make_region_gls_plot)

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

# Requested single figure containing all BM PGLS region plots together.
combined_pred_df <- map_dfr(gls_plot_data, "pred_df")
combined_plot_df <- map_dfr(gls_plot_data, "plot_df")
combined_human_pred_df <- map_dfr(gls_plot_data, "human_pred_df")

p5_combined_gls <- ggplot() +
  geom_ribbon(data = filter(combined_pred_df, !extrapolated),
              aes(x = x_log10, ymin = lwr_log10, ymax = upr_log10), alpha = 0.2) +
  geom_ribbon(data = filter(combined_pred_df,  extrapolated),
              aes(x = x_log10, ymin = lwr_log10, ymax = upr_log10), alpha = 0.1) +
  geom_line(data = filter(combined_pred_df, !extrapolated),
            aes(x = x_log10, y = fit_log10, group = VarLabel), linewidth = 0.8) +
  geom_line(data = filter(combined_pred_df,  extrapolated),
            aes(x = x_log10, y = fit_log10, group = VarLabel), linewidth = 0.8, linetype = "dashed") +
  geom_point(data = filter(combined_plot_df, !is_human), aes(x = x_log10, y = y_log10), size = 1.4, alpha = 0.8) +
  geom_point(data = filter(combined_plot_df,  is_human), aes(x = x_log10, y = y_log10), size = 2.2, shape = 17) +
  #geom_point(data = combined_human_pred_df, aes(x = x_log10, y = pred_log10), size = 2.0, shape = 4, stroke = 0.9) +
  facet_wrap(~VarLabel, scales = "free", ncol = 4) +
  theme_bw() +
  labs(title = "Phylogenetic GLS (BM + human mu): All Brain Structures vs Rest of Brain",
       subtitle = "nonhuman primates = points; human = triangle",
       x = "log10(Rest of Brain [Total brain - structure])",
       y = "log10(Structure volume)") +
  theme(strip.text = element_text(face = "bold", size = 8),
        axis.text = element_text(size = 7))
p5_combined_gls

save_plot(p5_combined_gls, "phylo_gls_region_plots_combined", width = 14, height = 11, jpg = TRUE)
save_plot(p5_combined_gls, "phylo_gls_region_plots_combined_mu_adjusted", width = 14, height = 11, jpg = TRUE)

# Extra diagnostic: at each human x-value, the BM + mu line used in the PGLS plot
# should exactly match the Brownian prediction used in the mosaicism plot.
bm_gls_mu_plot_alignment_check <- combined_human_pred_df %>%
  transmute(Structure = VarLabel,
            PGLS_plot_predicted_log10_at_human_x = pred_log10,
            PGLS_plot_fitted_line_log10_at_human_x = fitted_line_log10,
            PhyloMu_log) %>%
  left_join(
    bm_mosaicism_mu_alignment_check %>%
      transmute(Structure,
                Mosaicism_predicted_log10 = Predicted_log10,
                Mosaicism_observed_log10 = Observed_log10,
                Mosaicism_log10_diff = Log10Diff,
                Mosaicism_direction = Direction),
    by = "Structure"
  ) %>%
  mutate(abs_log10_difference_between_plot_and_mosaicism = abs(PGLS_plot_predicted_log10_at_human_x - Mosaicism_predicted_log10))

write.csv(bm_gls_mu_plot_alignment_check,
          file.path(CONFIG$dir_checks, "bm_gls_plot_line_matches_mosaicism_mu_prediction_check.csv"),
          row.names = FALSE)

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

# Structures included as plotted points but excluded from rCMRGlc fits/correlations.
always_excluded_from_fit <- c("Neocortex grey", "Cerebellum")
always_excluded_label <- paste(always_excluded_from_fit, collapse = " and ")

cor_df <- tbl_body %>%
  filter(!Structure %in% always_excluded_from_fit) %>%
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
# rCMRGlc vs prediction-error plots (fit excludes Neocortex grey + Cerebellum)
# ============================================================
num <- function(x) readr::parse_number(as.character(x))

df_bm <- core_with_rCMRGlc %>%
  filter(Model == "Brownian (BM)") %>%
  transmute(Structure, rCMRGlc = num(rCMRGlc),
            Diff.pre = num(Diff.pre), Diff.min = num(Diff.min), Diff.max = num(Diff.max),
            excluded = Structure %in% always_excluded_from_fit) %>%
  drop_na(rCMRGlc, Diff.pre)
df_fit <- filter(df_bm, !excluded)
df_excluded <- filter(df_bm,  excluded)

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
  geom_point(data = df_excluded, aes(rCMRGlc, Diff.pre), inherit.aes = FALSE,
             shape = 22, fill = "deepskyblue1", color = "blue4", size = 4.5, stroke = 0.9),
  geom_text(data = df_excluded, aes(rCMRGlc, Diff.pre, label = paste0(Structure, " (excluded)")),
            inherit.aes = FALSE, hjust = -0.08, size = 3.1, color = "blue4")
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
       subtitle = paste0("Fit excludes ", always_excluded_label, "; plotted points include them"),
       x = x_lab, y = y_lab) +
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
       subtitle = paste0("Models exclude ", always_excluded_label, "; best fit: degree ", best_degree),
       x = x_lab, y = y_lab, color = NULL) +
  theme_paper + theme(legend.position = "right")
p_poly_compare

# ---- Shared Spearman correlation + optional linear regression objects ----
# Both analyses use df_fit, so Neocortex grey and Cerebellum are excluded
# from the fitted statistics but still shown as excluded points.

fmt_p <- function(p, digits = 2, eps = 1e-4) {
  format.pval(p, digits = digits, eps = eps)
}

# ---- Optional linear regression, kept only for descriptive linear-fit panel ----
linear_fit <- lm(Diff.pre ~ rCMRGlc, data = df_fit)
linear_sum <- summary(linear_fit)
linear_fstat <- linear_sum$fstatistic
linear_p <- pf(linear_fstat[1], linear_fstat[2], linear_fstat[3], lower.tail = FALSE)

linear_coef <- coef(linear_fit)
linear_eqn <- paste0(
  "y = ", sprintf("%.4f", linear_coef[1]),
  ifelse(linear_coef[2] >= 0, " + ", " - "),
  sprintf("%.4f", abs(linear_coef[2])), "x"
)

linear_pred <- data.frame(rCMRGlc = xg)
linear_pred <- cbind(
  linear_pred,
  as.data.frame(predict(linear_fit, newdata = linear_pred, interval = "confidence"))
)

# ---- Spearman correlation ----
# exact = FALSE avoids warnings/problems when there are ties
cor_test <- cor.test(
  df_fit$rCMRGlc,
  df_fit$Diff.pre,
  method = "spearman",
  exact = FALSE
)

cor_stats <- paste0(
  "Spearman \u03c1 = ", sprintf("%.3f", unname(cor_test$estimate)),
  "\np = ", fmt_p(cor_test$p.value),
  "\nn = ", sum(complete.cases(df_fit$rCMRGlc, df_fit$Diff.pre))
)

linear_stats <- paste0(
  "Adj. R\u00b2 = ", sprintf("%.3f", linear_sum$adj.r.squared),
  "\np = ", fmt_p(linear_p)
)

# ---- p_scatter: scatterplot for Spearman correlation; no fitted line ----
p_scatter <- ggplot() +
  add_points_and_labels +
  annotate(
    "text",
    x = 10.6, y = 1.35,
    label = cor_stats,
    hjust = 0, vjust = 1,
    size = 3.2,
    color = "navy"
  ) +
  coord_cartesian(xlim = c(10, 40), ylim = c(-3.5, 1.5), clip = "on") +
  labs(
    title = "Scatterplot: rCMRGlc and BM prediction error",
    subtitle = paste0("Spearman correlation excludes ", always_excluded_label),
    x = x_lab,
    y = y_lab
  ) +
  theme_paper

p_scatter

# ---- p_linear: optional linear fit with 95% CI ribbon and model p ----
# Keep this only if you want to show a descriptive linear model.
# Note: this is separate from the Spearman rank correlation above.
p_linear <- ggplot() +
  geom_ribbon(
    data = linear_pred,
    aes(rCMRGlc, ymin = lwr, ymax = upr),
    inherit.aes = FALSE,
    fill = "steelblue",
    alpha = 0.2
  ) +
  geom_line(
    data = linear_pred,
    aes(rCMRGlc, fit),
    inherit.aes = FALSE,
    color = "steelblue4",
    linewidth = 1.2
  ) +
  add_points_and_labels +
  annotate(
    "text",
    x = 10.6, y = 1.35,
    label = linear_stats,
    hjust = 0, vjust = 1,
    size = 3.2,
    color = "navy"
  ) +
  coord_cartesian(xlim = c(10, 40), ylim = c(-3.5, 1.5), clip = "on") +
  labs(
    title = "Linear regression: rCMRGlc and BM prediction error",
    subtitle = paste0(
      "Linear model excludes ",
      always_excluded_label,
      "; plotted points include them"
    ),
    x = x_lab,
    y = y_lab
  ) +
  theme_paper

p_linear

save_plot(p_loess_bars,   "plot_loess_errorbars_BM_rCMRGlc")
save_plot(p_poly_compare, "plot_polynomial_comparison_BM_rCMRGlc", width = 9, height = 5.5)
save_plot(p_scatter,      "plot_scatter_spearman_BM_rCMRGlc")
save_plot(p_linear,       "plot_linear_CI_BM_rCMRGlc")


# ============================================================
# PLOT VII + TABLE — 12-way cubic / degree-3 robustness
#   4 structure configs x 3 PGLS evolution models.
#   Neocortex grey and Cerebellum are always excluded from the rCMRGlc regressions and
#   are plotted afterwards as distinct bright-blue squares.
# ============================================================

# always_excluded_from_fit is defined above as c("Neocortex grey", "Cerebellum")
structure_configs <- list(
  "All structures"            = character(0),
  "Exclude Neocortex white"   = "Neocortex white",
  "Exclude Capsula interna"   = "Capsula interna",
  "Exclude both white matter" = c("Neocortex white", "Capsula interna")
)

evolution_levels <- vapply(specs_compare, function(s) s$name, character(1))
structure_levels <- names(structure_configs)

# PGLS per evolution model x structure;
# structure exclusion only affects the lm() below.

diff_pre_long <- map_dfr(specs_compare, function(spec) {
  per_struct <- map_dfr(target_cols, function(vn) {
    run_models_for_structure(
      prep_structure(vn, data, tr),
      list(spec),
      want_ci = FALSE,
      mu_source = "trait"
    )
  })
  
  if (is.null(per_struct) || nrow(per_struct) == 0) return(NULL)
  per_struct %>%
    transmute(
      evolution_model = spec$name,
      Structure = label_region(Variable),
      Observed,
      Predicted,
      Diff.pre = (Observed - Predicted) / Observed
    )
})

diff_pre_with_rCMRGlc <- diff_pre_long %>%
  left_join(
    sup_clean %>% select(Structure, rCMRGlc),
    by = "Structure"
  ) %>%
  
  mutate(rCMRGlc = as.numeric(rCMRGlc)) %>%
  drop_na(rCMRGlc, Diff.pre)

# One cubic / degree-3 model per cell -> 12-row stats table

cubic_grid <- map_dfr(structure_levels, function(sc_name) {
  excl <- unique(c(always_excluded_from_fit, structure_configs[[sc_name]]))
  map_dfr(evolution_levels, function(ev_name) {
    sub <- diff_pre_with_rCMRGlc %>%
      filter(evolution_model == ev_name) %>%
      mutate(excluded = Structure %in% excl)
    fit_df <- filter(sub, !excluded)
    n_used <- nrow(fit_df)
    n_excl <- sum(sub$excluded)
    
    # Cubic has 4 coefficients: intercept, x, x^2, x^3.
    # Need at least 5 observations to have residual df for p-values.
    
    if (n_used < 5) {
      return(tibble(
        structure_config = sc_name,
        evolution_model = ev_name,
        n_used = n_used,
        n_excluded = n_excl,
        intercept = NA_real_,
        b_linear = NA_real_,
        b_quadratic = NA_real_,
        b_cubic = NA_real_,
        p_intercept = NA_real_,
        p_linear = NA_real_,
        p_quadratic = NA_real_,
        p_cubic = NA_real_,
        adj_R2 = NA_real_,
        AIC = NA_real_,
        BIC = NA_real_,
        model_p_value = NA_real_
      ))
    }
    
    fit <- lm(
      Diff.pre ~ rCMRGlc + I(rCMRGlc^2) + I(rCMRGlc^3),
      data = fit_df
    )
    
    summ <- summary(fit)
    coefs <- coef(summ)
    fstat <- summ$fstatistic
    model_p <- if (!is.null(fstat)) {
      unname(pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE))
    } else {
      NA_real_
    }
    
    tibble(
      structure_config = sc_name,
      evolution_model = ev_name,
      n_used = n_used,
      n_excluded = n_excl,
      intercept = unname(coefs["(Intercept)", "Estimate"]),
      b_linear = unname(coefs["rCMRGlc", "Estimate"]),
      b_quadratic = unname(coefs["I(rCMRGlc^2)", "Estimate"]),
      b_cubic = unname(coefs["I(rCMRGlc^3)", "Estimate"]),
      p_intercept = unname(coefs["(Intercept)", "Pr(>|t|)"]),
      p_linear = unname(coefs["rCMRGlc", "Pr(>|t|)"]),
      p_quadratic = unname(coefs["I(rCMRGlc^2)", "Pr(>|t|)"]),
      p_cubic = unname(coefs["I(rCMRGlc^3)", "Pr(>|t|)"]),
      adj_R2 = summ$adj.r.squared,
      AIC = AIC(fit),
      BIC = BIC(fit),
      model_p_value = model_p
    )
  })
})
# ============================================================
# PLOT VII + TABLE — 12-way linear robustness
#   4 structure configs x 3 PGLS evolution models.
#   One linear model per cell.
# ============================================================

# One linear model per cell -> 12-row stats table
linear_grid <- map_dfr(structure_levels, function(sc_name) {
  excl <- unique(c(always_excluded_from_fit, structure_configs[[sc_name]]))
  
  map_dfr(evolution_levels, function(ev_name) {
    sub <- diff_pre_with_rCMRGlc %>%
      filter(evolution_model == ev_name) %>%
      mutate(excluded = Structure %in% excl)
    
    fit_df <- filter(sub, !excluded)
    
    n_used <- nrow(fit_df)
    n_excl <- sum(sub$excluded)
    
    # Linear model has 2 coefficients: intercept and slope.
    # Need at least 3 observations for residual df and p-values.
    if (n_used < 3) {
      return(tibble(
        structure_config = sc_name,
        evolution_model = ev_name,
        n_used = n_used,
        n_excluded = n_excl,
        intercept = NA_real_,
        b_linear = NA_real_,
        p_intercept = NA_real_,
        p_linear = NA_real_,
        adj_R2 = NA_real_,
        R2 = NA_real_,
        AIC = NA_real_,
        BIC = NA_real_,
        model_p_value = NA_real_
      ))
    }
    
    fit <- lm(Diff.pre ~ rCMRGlc, data = fit_df)
    
    summ <- summary(fit)
    coefs <- coef(summ)
    fstat <- summ$fstatistic
    
    model_p <- if (!is.null(fstat)) {
      unname(pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE))
    } else {
      NA_real_
    }
    
    tibble(
      structure_config = sc_name,
      evolution_model = ev_name,
      n_used = n_used,
      n_excluded = n_excl,
      
      intercept = unname(coefs["(Intercept)", "Estimate"]),
      b_linear = unname(coefs["rCMRGlc", "Estimate"]),
      
      p_intercept = unname(coefs["(Intercept)", "Pr(>|t|)"]),
      p_linear = unname(coefs["rCMRGlc", "Pr(>|t|)"]),
      
      adj_R2 = summ$adj.r.squared,
      R2 = summ$r.squared,
      AIC = AIC(fit),
      BIC = BIC(fit),
      model_p_value = model_p
    )
  })
})

# Prediction grid: linear fit + 95% CI per cell
xg_compare <- seq(10, 40, length.out = 300)

linear_pred_grid <- map_dfr(structure_levels, function(sc_name) {
  excl <- unique(c(always_excluded_from_fit, structure_configs[[sc_name]]))
  
  map_dfr(evolution_levels, function(ev_name) {
    fit_df <- diff_pre_with_rCMRGlc %>%
      filter(evolution_model == ev_name) %>%
      mutate(excluded = Structure %in% excl) %>%
      filter(!excluded)
    
    if (nrow(fit_df) < 3) return(NULL)
    
    fit <- lm(Diff.pre ~ rCMRGlc, data = fit_df)
    
    pred <- as.data.frame(
      predict(
        fit,
        newdata = data.frame(rCMRGlc = xg_compare),
        interval = "confidence"
      )
    )
    
    tibble(
      structure_config = sc_name,
      evolution_model = ev_name,
      rCMRGlc = xg_compare,
      fit = pred$fit,
      lwr = pred$lwr,
      upr = pred$upr
    )
  })
})

# Per-cell point data.
# Neocortex grey and Cerebellum are always excluded from the regression fits,
# but plotted afterwards as their own point class. White-matter exclusions are panel-specific.
points_grid <- map_dfr(structure_levels, function(sc_name) {
  panel_excl <- structure_configs[[sc_name]]
  fit_excl <- unique(c(always_excluded_from_fit, panel_excl))

  diff_pre_with_rCMRGlc %>%
    mutate(
      structure_config = sc_name,
      always_excluded_structure = Structure %in% always_excluded_from_fit,
      excluded_white_matter = Structure %in% panel_excl,
      excluded = Structure %in% fit_excl
    )
})

# Lock facet ordering
points_grid$structure_config      <- factor(points_grid$structure_config,      levels = structure_levels)
points_grid$evolution_model       <- factor(points_grid$evolution_model,       levels = evolution_levels)
linear_pred_grid$structure_config <- factor(linear_pred_grid$structure_config, levels = structure_levels)
linear_pred_grid$evolution_model  <- factor(linear_pred_grid$evolution_model,  levels = evolution_levels)
linear_grid$structure_config      <- factor(linear_grid$structure_config,      levels = structure_levels)
linear_grid$evolution_model       <- factor(linear_grid$evolution_model,       levels = evolution_levels)

p_linear_compare_12way <- ggplot() +
  geom_ribbon(
    data = linear_pred_grid,
    aes(rCMRGlc, ymin = lwr, ymax = upr),
    inherit.aes = FALSE,
    fill = "steelblue",
    alpha = 0.2
  ) +
  geom_line(
    data = linear_pred_grid,
    aes(rCMRGlc, fit),
    inherit.aes = FALSE,
    color = "steelblue4",
    linewidth = 1
 ) +
  # Points used in the regression fit
  geom_point(
    data = filter(points_grid, !excluded),
    aes(rCMRGlc, Diff.pre),
    inherit.aes = FALSE,
    shape = 16,
    size = 1.8
  ) +
  # Panel-specific excluded white-matter structures
  geom_point(
    data = filter(points_grid, excluded_white_matter),
    aes(rCMRGlc, Diff.pre),
    inherit.aes = FALSE,
    shape = 24,
    fill = "firebrick2",
    color = "firebrick4",
    size = 3
  ) +
  # Neocortex grey and Cerebellum: always excluded from regression and overlaid afterwards
  geom_point(
    data = filter(points_grid, always_excluded_structure),
    aes(rCMRGlc, Diff.pre),
    inherit.aes = FALSE,
    shape = 22,
    fill = "deepskyblue1",
    color = "blue4",
    size = 3.4,
    stroke = 0.9
  ) +
  facet_grid(structure_config ~ evolution_model) +
  coord_cartesian(xlim = c(10, 40), ylim = c(-3.5, 1.5)) +
  labs(
    title = "Linear robustness: 4 structure configurations x 3 PGLS evolution models",
    subtitle = "Linear fit + 95% CI on included subset; white-matter exclusions = red triangles; Neocortex grey and Cerebellum = bright blue squares",
    x = x_lab,
    y = y_lab
  ) +
  theme_paper +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    panel.spacing = unit(0.8, "lines"),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11)
  )

p_linear_compare_12way

save_plot(
  p_linear_compare_12way,
  "plot_linear_comparison_12way",
  width = 14,
  height = 12,
  jpg = TRUE
)

# Full stats table for the 12 linear models
write_xlsx(
  linear_grid %>% arrange(structure_config, evolution_model),
  file.path(CONFIG$dir_tables, "Table_S_linear_comparison_12way.xlsx")
)

# P-value-focused table for the 12 linear models
linear_pvals_12way <- linear_grid %>%
  arrange(structure_config, evolution_model) %>%
  select(
    structure_config,
    evolution_model,
    n_used,
    n_excluded,
    p_intercept,
    p_linear,
    model_p_value,
    R2,
    adj_R2,
    AIC,
    BIC
  )

write_xlsx(
  linear_pvals_12way,
  file.path(CONFIG$dir_tables, "Table_S_linear_comparison_12way_pvalues.xlsx")
)