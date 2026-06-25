# ------------------------------------------------------------
# plot_scaling_with_residuals.R
# ------------------------------------------------------------

setwd("~/Library/CloudStorage/Dropbox/COLLABORATIVE/Do expensive brain regions increase less in humans/analyses_metabol_rate_structure")

library(ape)
library(caper)
library(dplyr)

# ---- Load data ----
df <- read.csv("data_raw/Stephan_primates.csv", stringsAsFactors = FALSE)

# Clean species names
df$Species <- gsub(" ", "_", df$Species)

# ---- Preferred brain volume ----
df <- df %>%
  mutate(
    Preferred_brain_volume = coalesce(Brain_volume, Brainvol, Total_brain_net_volume),
    logPreferred = log10(Preferred_brain_volume),
    logTotal = log10(Total_brain_net_volume)
  )

# ---- Load phylogeny ----
tree <- read.tree("data_raw/species.nwk")

# ============================================================
# SETS
# ============================================================

components1 <- c(
  "Medulla_oblongata",
  "Cerebellum",
  "Mesencephalon",
  "Diencephalon",
  "Telencephalon"
)

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

components2 <- names(region_labels)

cols1 <- c("steelblue", "darkorange", "forestgreen", "purple", "firebrick")
cols2 <- rainbow(length(components2))

# ============================================================
# OUTPUT
# ============================================================

dir.create("figs/traits", recursive = TRUE, showWarnings = FALSE)

# ============================================================
# CORE FUNCTION (now flexible predictor)
# ============================================================

plot_scaling <- function(components, cols, filename,
                         predictor = "logTotal",
                         use_pgls = FALSE,
                         labels_map = NULL) {
  
  dev.new()
  
  x_all <- df[[predictor]]
  y_all <- log10(as.matrix(df[, components]))
  
  plot(NULL,
       xlim = range(x_all, na.rm = TRUE),
       ylim = range(y_all, na.rm = TRUE),
       xlab = expression(log[10]*" brain volume (mm"^3*")"),
       ylab = expression(log[10]*" component volume (mm"^3*")"))
  
  slopes <- c()
  
  for (i in seq_along(components)) {
    
    comp <- components[i]
    
    y <- log10(df[[comp]])
    x <- df[[predictor]]
    
    d <- data.frame(Species=df$Species, x=x, y=y)
    d <- d[complete.cases(d), ]
    
    if (nrow(d) < 5) next
    
    if (use_pgls) {
      
      comp_data <- comparative.data(tree, d,
                                    names.col="Species",
                                    vcv=TRUE,
                                    warn.dropped=FALSE)
      
      fit <- pgls(y ~ x, data = comp_data)
      slope <- coef(fit)[2]
      
      newx <- seq(min(d$x), max(d$x), length.out=100)
      preds <- predict(fit, newdata=data.frame(x=newx))
      
      lines(newx, preds, col=cols[i], lwd=2)
      
    } else {
      
      fit <- lm(y ~ x, data = d)
      slope <- coef(fit)[2]
      
      abline(fit, col=cols[i], lwd=2)
    }
    
    slopes[i] <- slope
    
    points(d$x, d$y, col=cols[i], pch=16, cex=0.8)
  }
  
  if (!is.null(labels_map)) {
    nice <- labels_map[components]
  } else {
    nice <- components
  }
  
  legend_labels <- paste0(nice, " (", signif(slopes,2), ")")
  
  legend("topleft",
         legend=legend_labels,
         col=cols,
         pch=16,
         lwd=2,
         cex=0.7,
         bty="n",
         title=ifelse(use_pgls,"PGLS slope","OLS slope"))
  
  dev.copy(png, filename, width=900, height=700)
  dev.off()
}

# ============================================================
# HUMAN RESIDUALS (also uses predictor)
# ============================================================

get_human_residuals <- function(components,
                                predictor = "logTotal",
                                labels_map=NULL) {
  
  results <- data.frame()
  
  for (comp in components) {
    
    y <- log10(df[[comp]])
    x <- df[[predictor]]
    
    d <- data.frame(Species=df$Species, x=x, y=y)
    d <- d[complete.cases(d), ]
    
    if (!"Homo_sapiens" %in% d$Species) next
    
    # OLS
    fit_lm <- lm(y ~ x, data=d)
    pred_lm <- predict(fit_lm, newdata=subset(d, Species=="Homo_sapiens"))
    obs <- subset(d, Species=="Homo_sapiens")$y
    resid_lm <- obs - pred_lm
    
    # PGLS
    comp_data <- comparative.data(tree, d,
                                  names.col="Species",
                                  vcv=TRUE,
                                  warn.dropped=FALSE)
    
    fit_pgls <- pgls(y ~ x, data=comp_data)
    
    pred_pgls <- predict(fit_pgls,
                         newdata=data.frame(x=subset(d, Species=="Homo_sapiens")$x))
    
    resid_pgls <- obs - pred_pgls
    
    name <- if (!is.null(labels_map)) labels_map[comp] else comp
    
    results <- rbind(results, data.frame(
      Region = name,
      Residual_OLS = resid_lm,
      Residual_PGLS = resid_pgls
    ))
  }
  
  return(results)
}

# ============================================================
# RUN PLOTS
# ============================================================

# ✅ Set 1 → strict total brain
plot_scaling(components1, cols1,
             "figs/traits/components_major_OLS.png",
             predictor = "logTotal",
             use_pgls = FALSE)

plot_scaling(components1, cols1,
             "figs/traits/components_major_PGLS.png",
             predictor = "logTotal",
             use_pgls = TRUE)

# ✅ Set 2 → preferred brain volume
plot_scaling(components2, cols2,
             "figs/traits/components_regions_OLS.png",
             predictor = "logPreferred",
             use_pgls = FALSE,
             labels_map = region_labels)

plot_scaling(components2, cols2,
             "figs/traits/components_regions_PGLS.png",
             predictor = "logPreferred",
             use_pgls = TRUE,
             labels_map = region_labels)

# ============================================================
# HUMAN RESIDUALS
# ============================================================

res1 <- get_human_residuals(components1, predictor="logTotal")

res2 <- get_human_residuals(components2,
                            predictor="logPreferred",
                            labels_map=region_labels)

res_all <- rbind(res1, res2)

write.csv(res_all,
          "figs/traits/human_residuals.csv",
          row.names = FALSE)
