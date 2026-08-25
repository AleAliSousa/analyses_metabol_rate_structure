library(ape)
library(nlme)
library(dispRity)

# ============================================================
# PHYLOGENETIC REGRESSION FOR ONE BRAIN REGION
# ============================================================

fit_region_pgls <- function(
    region,
    data,
    tree,
    human = "Homo_sapiens",
    predictor_col = "Preferred_brain_volume"
) {
  
  # ----------------------------------------------------------
  # 1. Keep only the species, total brain volume, and the
  #    selected brain-region volume.
  # ----------------------------------------------------------
  model_data <- data[, c("Species", predictor_col, region)]
  
  # Remove rows missing either the predictor or region volume.
  model_data <- na.omit(model_data)
  
  
  # ----------------------------------------------------------
  # 2. Part-whole correction
  #
  # Instead of predicting a region from total brain volume,
  # predict it from the remainder of the brain:
  #
  # Rest of Brain = Total brain volume - Region volume
  #
  # This avoids correlating a structure directly with a total
  # that already contains that structure.
  # ----------------------------------------------------------
  model_data$Rest_of_Brain <-
    model_data[[predictor_col]] - model_data[[region]]
  
  # Logarithms require strictly positive values.
  model_data <- model_data[
    model_data$Rest_of_Brain > 0 &
      model_data[[region]] > 0,
  ]
  
  rownames(model_data) <- model_data$Species
  
  
  # ----------------------------------------------------------
  # 3. Align the data and phylogenetic tree
  #
  # clean.data() removes unmatched species and returns a tree
  # and dataset containing the same taxa.
  # ----------------------------------------------------------
  aligned <- clean.data(model_data, tree)
  
  model_tree <- aligned$tree
  model_data <- aligned$data
  
  rownames(model_data) <- model_data$Species
  
  
  # ----------------------------------------------------------
  # 4. Hold Homo sapiens out of the regression
  #
  # The PGLS model is fitted to nonhuman species. The fitted
  # model can then be used to predict the human region volume.
  # ----------------------------------------------------------
  if (!human %in% model_data$Species) {
    stop("Human is not present after tree-data alignment.")
  }
  
  human_data <- model_data[
    model_data$Species == human,
    ,
    drop = FALSE
  ]
  
  model_data <- model_data[
    model_data$Species != human,
    ,
    drop = FALSE
  ]
  
  model_tree <- drop.tip(model_tree, human)
  
  
  # ----------------------------------------------------------
  # 5. Define the phylogenetic regression
  #
  # Both variables are natural-log transformed.
  # ----------------------------------------------------------
  model_formula <- as.formula(
    paste0(
      "log(", region, ") ~ log(Rest_of_Brain)"
    )
  )
  
  
  # ----------------------------------------------------------
  # 6. Define the phylogenetic correlation structure
  #
  # fixed = FALSE means Pagel's lambda is estimated by
  # maximum likelihood rather than fixed in advance.
  # ----------------------------------------------------------
  phylo_correlation <- corPagel(
    value = 1,
    form = ~Species,
    phy = model_tree,
    fixed = FALSE
  )
  
  
  # ----------------------------------------------------------
  # 7. Fit the phylogenetic generalized least-squares model
  # ----------------------------------------------------------
  fit <- gls(
    model = model_formula,
    correlation = phylo_correlation,
    data = model_data,
    control = glsControl(
      opt = "optim",
      msMaxIter = 1000,
      msTol = 1e-6
    )
  )
  
  
  # ----------------------------------------------------------
  # 8. Extract the fitted Pagel lambda
  # ----------------------------------------------------------
  lambda <- as.numeric(
    coef(
      fit$modelStruct,
      unconstrained = FALSE
    )
  )
  
  
  # ----------------------------------------------------------
  # 9. Ordinary model-line prediction for the held-out human
  #
  # This is the prediction from the fitted PGLS regression
  # coefficients. It does not yet include the additional
  # human-specific conditional phylogenetic mu correction used
  # elsewhere in your full script.
  # ----------------------------------------------------------
  human_log_prediction <- as.numeric(
    c(1, log(human_data$Rest_of_Brain)) %*% coef(fit)
  )
  
  human_prediction <- exp(human_log_prediction)
  
  
  # ----------------------------------------------------------
  # 10. Return the model and key results
  # ----------------------------------------------------------
  list(
    region = region,
    fit = fit,
    coefficients = coef(fit),
    lambda = lambda,
    n_nonhuman_species = nrow(model_data),
    human_observed = human_data[[region]],
    human_predicted = human_prediction,
    human_rest_of_brain = human_data$Rest_of_Brain,
    model_data = model_data,
    model_tree = model_tree
  )
}