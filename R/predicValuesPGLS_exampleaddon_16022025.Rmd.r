## -------------------------
## ONE EXAMPLE REGION: Hippocampus
## -------------------------

region <- "Hippocampus"

# Prep (includes part-whole correction ROB = Total - Part, and drops human for fitting)
pp <- prep_structure(region, data = data, tr = tr, human = "Homo_sapiens")
stopifnot(!is.null(pp))

# Choose ONE phylo model for the regression line (BM is a good default)
cor_struct <- corBrownian(1, form = ~Species, phy = pp$tr2)

fit <- fit_gls_safe(pp$fmla, cor_struct, pp$data_model)
stopifnot(!is.null(fit))

# Build a prediction grid across Rest_of_Brain (ROB) observed range (non-human)
x_grid_raw <- seq(
  from = min(pp$data_model$Rest_of_Brain, na.rm = TRUE),
  to   = max(pp$data_model$Rest_of_Brain, na.rm = TRUE),
  length.out = 200
)

# Model is: log(y) = b0 + b1 * log(ROB) + error
X <- cbind(1, log(x_grid_raw))
colnames(X) <- names(coef(fit))  # usually "(Intercept)" and "log(Rest_of_Brain)"

beta <- coef(fit)
Vb   <- vcov(fit)                # covariance of fixed effects

eta  <- as.numeric(X %*% beta)                          # mean on log scale
se   <- sqrt(diag(X %*% Vb %*% t(X)))                  # SE of mean on log scale
z    <- qnorm(0.975)

pred_df <- tibble::tibble(
  Rest_of_Brain = x_grid_raw,
  fit  = exp(eta),
  lwr  = exp(eta - z * se),
  upr  = exp(eta + z * se)
)

# Data for plotting (include human for points; show human separately)
plot_df <- pp$d_full %>%
  dplyr::mutate(
    is_human = (Species == pp$human)
  )

# Plot on log-log axes (matches the fitted relationship)
ggplot() +
  geom_ribbon(
    data = pred_df,
    aes(x = Rest_of_Brain, ymin = lwr, ymax = upr),
    alpha = 0.2
  ) +
  geom_line(
    data = pred_df,
    aes(x = Rest_of_Brain, y = fit),
    linewidth = 1
  ) +
  geom_point(
    data = dplyr::filter(plot_df, !is_human),
    aes(x = Rest_of_Brain, y = .data[[region]]),
    size = 2,
    alpha = 0.8
  ) +
  geom_point(
    data = dplyr::filter(plot_df, is_human),
    aes(x = Rest_of_Brain, y = .data[[region]]),
    size = 3,
    shape = 17
  ) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  labs(
    title = paste0("Phylogenetic GLS (BM): ", region, " vs Rest of Brain (ROB)"),
    subtitle = "Line = fitted mean; band = 95% CI for mean. Human shown as triangle; fit excludes human.",
    x = "Rest of Brain (Total brain − structure)",
    y = paste0(region, " volume")
  )
