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
  
  # Subset and drop NA
  d1 <- subset(data, select = c("Species", "Brain_weight", structure_name))
  d1 <- as.data.frame(na.omit(d1))
  rownames(d1) <- d1$Species
  
  # Must have focal species in this structure data
  if (!(focal_species %in% d1$Species)) {
    return(NULL)
  }
  
  # Match tree and data
  clean <- clean.data(d1, tr)
  tr1 <- clean$tree
  data2 <- clean$data
  
  # Must still contain focal species after cleaning
  if (!(focal_species %in% rownames(data2))) {
    return(NULL)
  }
  
  # Observed values (original units)
  observed_unlog <- as.numeric(data[data$Species == focal_species, structure_name])
  
  # Predictor for focal species (log brain weight)
  xh <- log(as.numeric(data[data$Species == focal_species, "Brain_weight"]))
  
  # Drop focal species from tree and data for prediction
  tr2 <- drop.tip(tr1, focal_species)
  data2_nofocal <- subset(data2, Species != focal_species)
  
  # If too few points, skip
  if (nrow(data2_nofocal) < 3) {
    return(NULL)
  }
  
  form <- as.formula(paste0("log(", structure_name, ") ~ log(Brain_weight)"))

  bad <- lapply(structures, function(s) {
    d <- data[, c("Species", "Brain_weight", s)]
    d <- na.omit(d)
    any(d$Brain_weight <= 0 | d[[s]] <= 0)
  })
  structures[unlist(bad)]
  
  # Fit GLS with Brownian lambda = 1
  fit <- gls(
    form,
    correlation = corBrownian(1, form = ~Species, phy = tr2),
    data = data2_nofocal
  )

  # Fit GLS with Pagel's lambda
  # fit <- gls(
  #   form,
  #   correlation = corPagel(1, phy = tr2, form = ~Species, fixed = FALSE),
  #   data = data2_nofocal,
  #   method = "ML"
  # )
  
  # # Extract lambda
  # lambda_hat <- coef(fit$modelStruct$corStruct, unconstrained = FALSE) # added for Pagel's lambda
  
  sigma <- fit$sigma
  
  # Tree covariance (Brownian, correlation scale as in your script)
  v.full <- vcv(tr1, corr = TRUE)
  cc <- which(rownames(v.full) == focal_species)
  
  # If focal species is not in v.full (should not happen if checks passed), skip
  if (length(cc) != 1) {
    return(NULL)
  }
  
  Cih <- v.full[cc, -cc, drop = FALSE]
  
  # Response vector for other species (log scale)
  y_other <- log(as.numeric(data2_nofocal[[structure_name]]))
  xbar <- mean(y_other, na.rm = TRUE)
  
  # Conditional mean + conditional variance (as in your code)
  V_other <- vcv(tr2, corr = TRUE)
  mu <- Cih %*% solve(V_other) %*% (y_other - xbar)
  ch <- v.full[cc, cc] - Cih %*% solve(V_other) %*% t(Cih)
  
  # Linear predictor + bias
  predicted <- c(1, xh) %*% coef(fit) + mu
  std.error <- sqrt(as.numeric(sigma^2 * ch))
  
  # 95% CI (log scale), then back-transform
  upper <- predicted + qnorm(0.975) * std.error
  lower <- predicted - qnorm(0.975) * std.error
  
  out <- data.frame(
    Structure = structure_name,
    Predicted = as.numeric(exp(predicted)),
    Observed  = observed_unlog,
    Lower     = as.numeric(exp(lower)),
    Upper     = as.numeric(exp(upper)),
#    Lambda    = lambda_hat, # added for Pagel's lambda
    N_used    = nrow(data2_nofocal),
    stringsAsFactors = FALSE
  )
  
  return(out)
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
