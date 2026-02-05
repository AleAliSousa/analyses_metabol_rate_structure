library(nlme)
library(ape)
library(ggplot2)
library(dplyr)

# --- CONFIGURATION ---
target_cols <- c(
  "LGN_Sousa", "Amygdala", "Pallidum", "NeoW_Frahm", "Medulla_oblongata",
  "Nucleus_subthalamicus", "Capsula_interna", "Striatum", "Diencephalon",
  "ASG_Sousa", "NeoG_Frahm", "Mesencephalon", "Cerebellum", "Hippocampus",
  "Total_insula_volume_L"
)

# --- PLOTTING FUNCTION ---
plot_regressions <- function(model_type = "lambda_ml") {
  
  all_points <- list()
  all_lines <- list()
  
  print(paste("Generating plots for model:", model_type))
  
  for (var_name in target_cols) {
    if (!var_name %in% names(data)) next
    
    # 1. PREPARE DATA
    d1 <- subset(data, select=c("Species", "Brain_weight", var_name))
    d1 <- as.data.frame(na.omit(d1))
    
    # Rest of Brain Calculation
    # Note: Use [[var_name]] to access the column dynamically
    d1$Rest_of_Brain <- d1$Brain_weight - d1[[var_name]]
    d1 <- d1[d1$Rest_of_Brain > 0, ] 
    
    rownames(d1) <- d1$Species
    
    clean <- clean.data(d1, tr) 
    tr1 <- clean$tree
    data_all <- clean$data
    if(!is.null(data_all$Species)) rownames(data_all) <- data_all$Species
    
    tr2 <- drop.tip(tr1, "Homo_sapiens")
    data_model <- subset(data_all, Species != "Homo_sapiens")
    
    fmla <- as.formula(paste("log(", var_name, ") ~ log(Rest_of_Brain)"))
    
    # 2. FIT MODEL
    fit <- NULL
    if (model_type == "lambda_0") {
      fit <- tryCatch({ gls(fmla, correlation = corPagel(0, form=~Species, phy=tr2, fixed=TRUE), data = data_model) }, error = function(e) NULL)
    } else if (model_type == "lambda_ml") {
      fit <- tryCatch({ gls(fmla, correlation = corPagel(1, form=~Species, phy=tr2, fixed=FALSE), data = data_model) }, error = function(e) NULL)
    } else if (model_type == "lambda_1") {
      fit <- tryCatch({ gls(fmla, correlation = corBrownian(1, form=~Species, phy=tr2), data = data_model) }, error = function(e) NULL)
    }
    
    if (is.null(fit)) next
    
    # 3. STANDARDIZE DATA FOR PLOTTING (The Fix)
    # We create a new dataframe with FIXED column names
    df_temp <- data.frame(
      Species = data_all$Species,
      Structure_Name = var_name,          # Store the name as a value, not a column header
      Log_X = log(data_all$Rest_of_Brain),
      Log_Y = log(data_all[[var_name]]),  # Dynamically grab the value
      Is_Human = ifelse(data_all$Species == "Homo_sapiens", "Human", "Primate"),
      stringsAsFactors = FALSE
    )
    
    all_points[[var_name]] <- df_temp
    
    # 4. STORE REGRESSION LINE INFO
    coefs <- coef(fit)
    
    if (model_type == "lambda_ml") {
      lambda_val <- as.numeric(coef(fit$modelStruct, unconstrained = FALSE))
      # Clamp for display if needed
      lambda_val <- max(0, min(1, lambda_val)) 
      label_txt <- sprintf("λ=%.2f", lambda_val)
    } else if (model_type == "lambda_0") {
      label_txt <- "λ=0"
    } else {
      label_txt <- "λ=1"
    }
    
    all_lines[[var_name]] <- data.frame(
      Structure_Name = var_name,
      Intercept = coefs[1],
      Slope = coefs[2],
      Label = label_txt
    )
  }
  # 5. COMBINE
  df_plot_points <- do.call(rbind, all_points)
  df_plot_lines <- do.call(rbind, all_lines)
  
  # 6. GGPLOT
  p <- ggplot(df_plot_points, aes(x = Log_X, y = Log_Y)) +
    
    geom_abline(data = df_plot_lines, aes(intercept = Intercept, slope = Slope), 
                color = "#377eb8", size = 0.8) +
    
    geom_point(data = subset(df_plot_points, Is_Human == "Primate"), 
               color = "gray60", alpha = 0.6, size = 1.5) +
    
    geom_point(data = subset(df_plot_points, Is_Human == "Human"), 
               color = "#e41a1c", shape = 17, size = 3) +
    
    geom_text(data = df_plot_lines, aes(label = Label), 
              x = -Inf, y = Inf, hjust = -0.2, vjust = 1.5, 
              size = 3, fontface = "italic") +
    
    # Facet by the standardized column "Structure_Name"
    facet_wrap(~Structure_Name, scales = "free") +
    
    theme_bw() +
    labs(
      title = paste("Regression: Log(Structure) vs Log(Rest of Brain) |", toupper(model_type)),
      subtitle = "Blue Line = PGLS Regression (excluding Human). Red Triangle = Human Observed.",
      x = "Log(Rest of Brain Weight)",
      y = "Log(Structure Volume)"
    ) +
    theme(strip.text = element_text(face="bold"))
  
  print(p)
}