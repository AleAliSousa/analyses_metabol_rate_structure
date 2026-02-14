setwd("~/Library/CloudStorage/Dropbox/COLLABORATIVE/Do expensive brain regions increase less in humans/1A analyses metabol rates regions/R/original with humans _18122025")

### Loop over structures and produce a combined plot
# Refactor of your Striatum-only script into a loop over all structures

library(ape)
library(nlme)
library(tidyverse)
library(dispRity)
library(scales)

# --- Load tree and data
tr <- read.tree("species.nwk")
Stephan_primates <- read.csv("Stephan_primates.csv")

# --- Clean columns
data <- subset(Stephan_primates, select = -c(X, order))

# --- Structures list (as in your script)
structures <- c(
  "LGN_Sousa",
  "Amygdala",
  "Pallidum",
  "NeoW_Frahm",
  "Medulla_oblongata",
  "Nucleus_subthalamicus",
  "Capsula_interna",
  "Striatum",
  "Diencephalon",
  "ASG_Sousa",
  "NeoG_Frahm",
  "Mesencephalon",
  "Cerebellum",
  "Hippocampus",
  "Total_insula_volume_L"
)

# --- Sanity check: structures present?
missing_cols <- setdiff(structures, colnames(data))
if (length(missing_cols) > 0) {
  stop(paste("These structures are missing from data:", paste(missing_cols, collapse = ", ")))
}

# --- Required columns present?
req_cols <- c("Species", "Brain_weight")
if (any(!req_cols %in% colnames(data))) {
  stop("Missing required columns: Species and/or Brain_weight")
}

# --- Utility: one-structure computation (same logic you used for Striatum)
predict_one_structure <- function(structure_name, data, tr, focal_species = "Homo_sapiens") {
  
  # --- Pull required columns and coerce numeric
  d1 <- data[, c("Species", "Brain_weight", structure_name)]
  d1$Brain_weight        <- suppressWarnings(as.numeric(d1$Brain_weight))
  d1[[structure_name]]   <- suppressWarnings(as.numeric(d1[[structure_name]]))
  
  # --- Keep only finite, positive values (prevents log(0), NA, Inf)
  keep <- is.finite(d1$Brain_weight) & is.finite(d1[[structure_name]]) &
    d1$Brain_weight > 0 & d1[[structure_name]] > 0
  d1 <- d1[keep, , drop = FALSE]
  
  if (nrow(d1) < 5) return(NULL)
  if (!(focal_species %in% d1$Species)) return(NULL)
  
  # --- Match/prune tree and data
  # Put taxa names in rownames for clean.data()
  rownames(d1) <- d1$Species
  d1$Species <- NULL
  
  clean <- dispRity::clean.data(d1, tr)
  tr1   <- clean$tree
  
  data2 <- as.data.frame(clean$data)
  data2$Species <- rownames(data2)
  rownames(data2) <- NULL
  
  if (!(focal_species %in% data2$Species)) return(NULL)
  
  # --- Focal observed and predictor from CLEANED data
  focal_row <- data2[data2$Species == focal_species, , drop = FALSE]
  if (nrow(focal_row) != 1) return(NULL)
  
  observed_unlog <- as.numeric(focal_row[[structure_name]])
  bw_focal       <- as.numeric(focal_row$Brain_weight)
  if (!is.finite(observed_unlog) || observed_unlog <= 0) return(NULL)
  if (!is.finite(bw_focal) || bw_focal <= 0) return(NULL)
  
  xh <- log(bw_focal)
  
  # --- Drop focal species for model fit
  tr2 <- ape::drop.tip(tr1, focal_species)
  data2_nofocal <- data2[data2$Species != focal_species, , drop = FALSE]
  if (nrow(data2_nofocal) < 3) return(NULL)
  
  # --- Ensure data order matches tree order (critical for vcv/solve)
  data2_nofocal <- data2_nofocal[match(tr2$tip.label, data2_nofocal$Species), , drop = FALSE]
  if (anyNA(data2_nofocal$Species)) return(NULL)
  
  # --- Precompute logs for stability and transparency
  data2_nofocal$log_bw <- log(as.numeric(data2_nofocal$Brain_weight))
  data2_nofocal$log_y  <- log(as.numeric(data2_nofocal[[structure_name]]))
  
  keep2 <- is.finite(data2_nofocal$log_bw) & is.finite(data2_nofocal$log_y)
  data2_nofocal <- data2_nofocal[keep2, , drop = FALSE]
  if (nrow(data2_nofocal) < 3) return(NULL)
  
  # If any rows were dropped, prune tree accordingly and re-match order
  if (length(setdiff(tr2$tip.label, data2_nofocal$Species)) > 0) {
    tr2 <- ape::drop.tip(tr2, setdiff(tr2$tip.label, data2_nofocal$Species))
    data2_nofocal <- data2_nofocal[match(tr2$tip.label, data2_nofocal$Species), , drop = FALSE]
  }
  
  data2_nofocal$Species <- factor(data2_nofocal$Species, levels = tr2$tip.label)
  
  # --- Fit Pagel-lambda GLS with proper tryCatch
  fit <- tryCatch(
    nlme::gls(
      log_y ~ log_bw,
      correlation = ape::corPagel(1, phy = tr2, form = ~Species, fixed = FALSE),
      data = data2_nofocal,
      method = "ML"
    ),
    error = function(e) {
      message("Pagel fit failed for ", structure_name, ": ", e$message)
      return(NULL)
    }
  )
  if (is.null(fit)) return(NULL)
  
  lambda_hat <- as.numeric(coef(fit$modelStruct$corStruct, unconstrained = FALSE))
  sigma      <- fit$sigma
  
  # --- Brownian correlation matrices
  R_full  <- ape::vcv(tr1, corr = TRUE)  # includes focal
  R_other <- ape::vcv(tr2, corr = TRUE)  # excludes focal
  
  # Apply Pagel's lambda to off-diagonals
  off_full <- lower.tri(R_full) | upper.tri(R_full)
  R_full[off_full] <- lambda_hat * R_full[off_full]
  diag(R_full) <- 1
  
  off_other <- lower.tri(R_other) | upper.tri(R_other)
  R_other[off_other] <- lambda_hat * R_other[off_other]
  diag(R_other) <- 1
  
  # --- Conditional mean/variance for focal species
  cc <- which(rownames(R_full) == focal_species)
  if (length(cc) != 1) return(NULL)
  
  Cih <- R_full[cc, -cc, drop = FALSE]
  
  y_other <- data2_nofocal$log_y
  xbar    <- mean(y_other, na.rm = TRUE)
  
  mu <- Cih %*% solve(R_other, (y_other - xbar))
  ch <- R_full[cc, cc] - Cih %*% solve(R_other, t(Cih))
  
  predicted_log <- as.numeric(c(1, xh) %*% coef(fit) + mu)
  se_log        <- sqrt(as.numeric(sigma^2 * ch))
  
  upper_log <- predicted_log + qnorm(0.975) * se_log
  lower_log <- predicted_log - qnorm(0.975) * se_log
  
  data.frame(
    Structure = structure_name,
    Predicted = exp(predicted_log),
    Observed  = observed_unlog,
    Lower     = exp(lower_log),
    Upper     = exp(upper_log),
    Lambda    = lambda_hat,
    N_used    = nrow(data2_nofocal),
    stringsAsFactors = FALSE
  )
}

# --- Run loop
results_list <- lapply(structures, predict_one_structure, data = data, tr = tr, focal_species = "Homo_sapiens")
results <- do.call(rbind, results_list)

if (is.null(results) || nrow(results) == 0) {
  stop("No results produced. Check that Homo_sapiens is present and structures have data after NA omission.")
}

# Optional: write results
write.csv(results, "predicted_all_structures.csv", row.names = FALSE)

# Human-readable labels to match your figure
label_map <- c(
  "Total_insula_volume_L" = "Insular cortex grey",
  "Hippocampus"           = "Hippocampus",
  "Cerebellum"            = "Cerebellum",
  "Mesencephalon"         = "Mesencephalon",
  "NeoG_Frahm"            = "Neocortex grey",
  "ASG_Sousa"             = "Area striata grey",
  "Striatum"              = "Striatum",
  "Capsula_interna"       = "Capsula interna",
  "Nucleus_subthalamicus" = "Nucleus subthalamicus Luysi",
  "NeoW_Frahm"            = "Neocortex white",
  "Pallidum"              = "Pallidum",
  "Amygdala"              = "Amygdala",
  "LGN_Sousa"             = "Corpus geniculatum laterale",
  "Medulla_oblongata"     = "Medulla oblongata",
  "Diencephalon"          = "Diencephalon"
)

results_plot <- results %>%
  mutate(
    Structure_label = label_map[Structure],
    Structure_label = factor(
      Structure_label,
      levels = Structure_label[order(Predicted)]
    )
  )
# --- One combined plot (all structures together)
p_pred_obs <- ggplot(results_plot,
            aes(x = Structure_label, y = Predicted)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                width = 0.15, colour = "black") +
  geom_point(size = 2, shape = 4, colour = "black") +
  geom_point(aes(y = Observed), size = 2, shape = 1, colour = "red") +
  scale_y_log10(
    labels = label_comma(accuracy = 1)
  ) +
  xlab("Structure") +
  ylab("Observed (red) and predicted (black) volumes") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p_pred_obs)
ggsave("predicted_vs_observed_all_structures.png", plot = p_pred_obs, width = 10, height = 5, dpi = 300)

# Proportional difference (Predicted - Observed) / Observed
results2 <- results %>%
  mutate(
    prop_diff = (Predicted - Observed) / Observed,
    Structure_label = unname(label_map[Structure]),
    Structure_label = ifelse(is.na(Structure_label), Structure, Structure_label)
  ) %>%
  arrange(prop_diff) %>%
  mutate(Structure_label = factor(Structure_label, levels = rev(Structure_label)))

## Plot (ggplot2), styled like your example
p_propdiff <- ggplot(results2, aes(x = prop_diff, y = Structure_label)) +
  geom_col(width = 0.8, fill = "#69b3a2", color = "black", linewidth = 0.4) +
  geom_vline(xintercept = 0, linewidth = 0.6) +
  scale_x_continuous(breaks = seq(-1, 1, 0.5)) +
  coord_cartesian(xlim = c(-1, 1)) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank()
  )

print(p_propdiff)
ggsave("Fig_proportional_difference_all_structures.png", p_propdiff, width = 8, height = 6, dpi = 300)


# Positive values mean the model overpredicts the observed structure size.
# Negative values mean the model underpredicts it.
