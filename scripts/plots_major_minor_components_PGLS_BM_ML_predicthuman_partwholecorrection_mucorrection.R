# ------------------------------------------------------------
# plot_scaling_partwhole_mu_FINAL.R
# ------------------------------------------------------------

setwd("~/Library/CloudStorage/Dropbox/COLLABORATIVE/Do expensive brain regions increase less in humans/analyses_metabol_rate_structure")

library(ape)
library(caper)
library(dplyr)

# ---- Load data ----
df <- read.csv("data_raw/Stephan_primates.csv", stringsAsFactors = FALSE)
df$Species <- gsub(" ", "_", df$Species)
df$Total_brain <- df$Total_brain_net_volume

tree <- read.tree("data_raw/species.nwk")

# ============================================================
# HELPERS
# ============================================================

transform_vcv_lambda <- function(tree, lambda) {
  v <- vcv(tree, corr=TRUE)
  d <- diag(v)
  v <- v * lambda
  diag(v) <- d
  v
}

compute_mu <- function(tree_full, tree_reduced, human, X, lambda) {
  
  v_full <- transform_vcv_lambda(tree_full, lambda)
  v_red  <- transform_vcv_lambda(tree_reduced, lambda)
  
  # enforce shared species
  common <- intersect(names(X), rownames(v_red))
  X <- X[common]
  v_red <- v_red[common, common]
  
  if (length(common) < 3) return(0)
  
  cc <- which(rownames(v_full) == human)
  
  Cih <- v_full[cc, common, drop=FALSE]
  inv_v <- solve(v_red)
  
  xbar <- mean(X)
  
  mu <- Cih %*% inv_v %*% (X - xbar)
  
  as.numeric(mu)
}

# ============================================================
# MODEL FITTING
# ============================================================

fit_model <- function(d, model) {
  
  if (model == "OLS") {
    return(list(fit = lm(y ~ x, data=d), lambda = 0))
  }
  
  comp_data <- comparative.data(tree, d,
                                names.col="Species",
                                vcv=TRUE,
                                warn.dropped=FALSE)
  
  if (model == "PGLS_BM") {
    fit <- pgls(y ~ x, data=comp_data, lambda=1)
    return(list(fit=fit, lambda=1))
  }
  
  fit <- pgls(y ~ x, data=comp_data)
  return(list(fit=fit, lambda=fit$param["lambda"]))
}

# ============================================================
# CORE COMPONENT FUNCTION
# ============================================================

run_component <- function(comp, model) {
  
  d_all <- data.frame(
    Species=df$Species,
    comp=df[[comp]],
    total=df$Total_brain
  )
  
  d_all$ROB <- d_all$total - d_all$comp
  d_all$x <- log10(d_all$ROB)
  d_all$y <- log10(d_all$comp)
  
  d_all <- d_all[complete.cases(d_all), ]
  d_all <- d_all[d_all$ROB > 0, ]
  
  # keep only species in tree
  d_all <- d_all[d_all$Species %in% tree$tip.label, ]
  
  # split
  d <- subset(d_all, Species != "Homo_sapiens")
  human <- subset(d_all, Species == "Homo_sapiens")
  
  if (nrow(d) < 5 || nrow(human) != 1) return(NULL)
  
  # prune tree to match data
  tree_red <- drop.tip(tree, setdiff(tree$tip.label, d$Species))
  
  # fit
  model_out <- fit_model(d[,c("Species","x","y")], model)
  fit <- model_out$fit
  lambda <- model_out$lambda
  
  # build X aligned to tree
  X <- d$y
  names(X) <- d$Species
  X <- X[tree_red$tip.label]
  
  if (any(is.na(X))) return(NULL)
  
  mu <- compute_mu(tree, tree_red, "Homo_sapiens", X, lambda)
  
  list(d=d, human=human, fit=fit, mu=mu)
}

# ============================================================
# PLOTTING
# ============================================================

plot_scaling <- function(components, cols, filename, model, labels_map=NULL) {
  
  dev.new()
  
  # axis limits
  x_vals <- c()
  y_vals <- c()
  
  for (comp in components) {
    
    tmp <- data.frame(
      comp=df[[comp]],
      total=df$Total_brain
    )
    
    tmp$ROB <- tmp$total - tmp$comp
    tmp$x <- log10(tmp$ROB)
    tmp$y <- log10(tmp$comp)
    
    tmp <- tmp[complete.cases(tmp), ]
    tmp <- tmp[tmp$ROB > 0, ]
    
    x_vals <- c(x_vals, tmp$x)
    y_vals <- c(y_vals, tmp$y)
  }
  
  x_vals <- x_vals[is.finite(x_vals)]
  y_vals <- y_vals[is.finite(y_vals)]
  
  plot(NULL,
       xlim=range(x_vals),
       ylim=range(y_vals),
       xlab="log10 Rest of Brain",
       ylab="log10 component volume")
  
  slopes <- c()
  
  for (i in seq_along(components)) {
    
    comp <- components[i]
    res <- run_component(comp, model)
    if (is.null(res)) next
    
    d <- res$d
    human <- res$human
    fit <- res$fit
    mu <- res$mu
    
    beta <- coef(fit)
    slopes[i] <- beta[2]
    
    newx <- seq(min(d$x), max(d$x), length.out=100)
    preds <- beta[1] + beta[2]*newx + mu
    
    lines(newx, preds, col=cols[i], lwd=2)
    points(d$x, d$y, col=cols[i], pch=16, cex=0.8)
    
    points(human$x, human$y,
           pch=1, col="black", cex=1.6, lwd=2)
  }
  
  nice <- if (!is.null(labels_map)) labels_map[components] else components
  
  legend("topleft",
         legend=paste0(nice," (",signif(slopes,2),")"),
         col=cols, pch=16, lwd=2, cex=0.7, bty="n",
         title=model)
  
  dev.copy(png, filename, width=900, height=700)
  dev.off()
}

# ============================================================
# RESIDUALS
# ============================================================

get_residuals <- function(components, model, labels_map=NULL) {
  
  out <- data.frame()
  
  for (comp in components) {
    
    res <- run_component(comp, model)
    if (is.null(res)) next
    
    human <- res$human
    fit <- res$fit
    mu <- res$mu
    
    beta <- coef(fit)
    
    pred <- beta[1] + beta[2]*human$x + mu
    obs <- human$y
    
    name <- if (!is.null(labels_map)) labels_map[comp] else comp
    
    out <- rbind(out, data.frame(
      Region=name,
      Model=model,
      Residual=obs - pred
    ))
  }
  
  out
}

# ============================================================
# RUN
# ============================================================

components1 <- c("Medulla_oblongata","Cerebellum","Mesencephalon",
                 "Diencephalon","Telencephalon")

region_labels <- c(
  "LGN_Sousa"="Corpus geniculatum laterale",
  "Amygdala"="Amygdala",
  "Pallidum"="Pallidum",
  "NeoW_Frahm"="Neocortex white",
  "Total_insula_volume_L"="Insular cortex (grey)",
  "Nucleus_subthalamicus"="Nucleus subthalamicus Luysi",
  "Capsula_interna"="Capsula interna",
  "Striatum"="Striatum",
  "ASG_Sousa"="Area striata grey",
  "NeoG_Frahm"="Neocortex grey",
  "Mesencephalon"="Mesencephalon",
  "Cerebellum"="Cerebellum",
  "Hippocampus"="Hippocampus"
)

components2 <- names(region_labels)

cols1 <- c("steelblue","darkorange","forestgreen","purple","firebrick")
cols2 <- rainbow(length(components2))

models <- c("OLS","PGLS_ML","PGLS_BM")

for (m in models) {
  
  plot_scaling(components1, cols1,
               paste0("figs/traits/components_major_partwhole_mu_",m,".png"),
               model=m)
  
  plot_scaling(components2, cols2,
               paste0("figs/traits/components_regions_partwhole_mu_",m,".png"),
               model=m,
               labels_map=region_labels)
}

# ============================================================
# OUTPUT
# ============================================================

res_all <- do.call(rbind, lapply(models, function(m) {
  rbind(
    get_residuals(components1, m),
    get_residuals(components2, m, region_labels)
  )
}))

write.csv(res_all,
          "figs/traits/human_residuals_partwhole_mu.csv",
          row.names=FALSE)
