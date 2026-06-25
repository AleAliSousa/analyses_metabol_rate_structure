# ------------------------------------------------------------
# plot_scaling_full_models_nohumans.R
# ------------------------------------------------------------

setwd("~/Library/CloudStorage/Dropbox/COLLABORATIVE/Do expensive brain regions increase less in humans/analyses_metabol_rate_structure")

library(ape)
library(caper)
library(dplyr)

# ---- Load data ----
df <- read.csv("data_raw/Stephan_primates.csv", stringsAsFactors = FALSE)
df$Species <- gsub(" ", "_", df$Species)

# ---- Brain size variables ----
df <- df %>%
  mutate(
    Preferred_brain_volume = coalesce(Brain_volume, Brainvol, Total_brain_net_volume),
    logTotal = log10(Total_brain_net_volume),
    logPreferred = log10(Preferred_brain_volume)
  )

# ---- Tree ----
tree <- read.tree("data_raw/species.nwk")

# ============================================================
# SETS
# ============================================================

components1 <- c(
  "Medulla_oblongata","Cerebellum","Mesencephalon",
  "Diencephalon","Telencephalon"
)

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

dir.create("figs/traits", recursive=TRUE, showWarnings=FALSE)

# ============================================================
# MODEL FUNCTION (NO HUMANS EVER)
# ============================================================

fit_model <- function(d, model) {
  
  if (model == "OLS") {
    return(lm(y ~ x, data=d))
  }
  
  comp_data <- comparative.data(tree, d,
                                names.col="Species",
                                vcv=TRUE,
                                warn.dropped=FALSE)
  
  if (model == "PGLS_BM") {
    return(pgls(y ~ x, data=comp_data, lambda=1))
  } else {
    return(pgls(y ~ x, data=comp_data))
  }
}

# ============================================================
# PLOTTING FUNCTION
# ============================================================

plot_scaling <- function(components, cols, filename,
                         predictor, model,
                         labels_map=NULL) {
  
  dev.new()
  
  x_all <- df[[predictor]]
  y_all <- log10(as.matrix(df[,components]))
  
  plot(NULL,
       xlim=range(x_all,na.rm=TRUE),
       ylim=range(y_all,na.rm=TRUE),
       xlab="log10 brain volume",
       ylab="log10 component volume")
  
  slopes <- c()
  
  for (i in seq_along(components)) {
    
    comp <- components[i]
    
    # ---- full dataset ----
    d_all <- data.frame(
      Species=df$Species,
      x=df[[predictor]],
      y=log10(df[[comp]])
    )
    
    d_all <- d_all[complete.cases(d_all), ]
    
    # ---- split ----
    d <- subset(d_all, Species != "Homo_sapiens")
    human <- subset(d_all, Species == "Homo_sapiens")
    
    if (nrow(d) < 5) next
    
    # ---- fit WITHOUT humans ----
    fit <- fit_model(d, model)
    
    slope <- coef(fit)[2]
    slopes[i] <- slope
    
    # ---- line ----
    newx <- seq(min(d$x), max(d$x), length.out=100)
    preds <- predict(fit, newdata=data.frame(x=newx))
    
    lines(newx, preds, col=cols[i], lwd=2)
    
    # ---- non-human points ----
    points(d$x, d$y, col=cols[i], pch=16, cex=0.8)
    
    # ---- human point LAST ----
    if (nrow(human) == 1) {
      points(human$x, human$y,
             pch=1, col="black", cex=1.6, lwd=2)
    }
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
# HUMAN RESIDUALS (STRICTLY OUT-OF-SAMPLE)
# ============================================================

get_residuals <- function(components, predictor, model, labels_map=NULL) {
  
  out <- data.frame()
  
  for (comp in components) {
    
    d_all <- data.frame(
      Species=df$Species,
      x=df[[predictor]],
      y=log10(df[[comp]])
    )
    
    d_all <- d_all[complete.cases(d_all), ]
    
    if (!"Homo_sapiens" %in% d_all$Species) next
    
    d <- subset(d_all, Species != "Homo_sapiens")
    human <- subset(d_all, Species == "Homo_sapiens")
    
    fit <- fit_model(d, model)
    
    obs <- human$y
    pred <- predict(fit, newdata=data.frame(x=human$x))
    
    name <- if (!is.null(labels_map)) labels_map[comp] else comp
    
    out <- rbind(out, data.frame(
      Region=name,
      Model=model,
      Residual=obs - pred
    ))
  }
  
  return(out)
}

# ============================================================
# RUN PLOTS
# ============================================================

models <- c("OLS","PGLS_ML","PGLS_BM")

for (m in models) {
  
  plot_scaling(components1, cols1,
               paste0("figs/traits/components_major_",m,".png"),
               predictor="logTotal",
               model=m)
  
  plot_scaling(components2, cols2,
               paste0("figs/traits/components_regions_",m,".png"),
               predictor="logPreferred",
               model=m,
               labels_map=region_labels)
}

# ============================================================
# RESIDUALS TABLE
# ============================================================

res_all <- data.frame()

for (m in models) {
  
  res_all <- rbind(res_all,
                   get_residuals(components1, "logTotal", m),
                   get_residuals(components2, "logPreferred", m, region_labels))
}

write.csv(res_all,
          "figs/traits/human_residuals_all_models.csv",
          row.names=FALSE)
