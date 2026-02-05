setwd("~/Library/CloudStorage/Dropbox/COLLABORATIVE/Do expensive brain regions increase less in humans/1A analyses metabol rates regions/R/original with humans _18122025")

library(ape)
library(nlme)  
# Mac note: ML optimizer is different in Mac than PC so extra code was added
library(tidyverse)
library(dispRity)
library(scales)
library(readxl) 
library(writexl)   # for writing Excel .xlsx files
# --- SETUP ALL ---
# --- Load tree and data
tr <- read.tree("species.nwk")
Stephan_primates <- read.csv("Stephan_primates.csv")

# --- Clean columns
data_clean <- subset(Stephan_primates, select = -c(X, order))

#  --- Create a unified “preferred brain volume” column
data <- data_clean %>%
  mutate(
    Preferred_brain_volume = coalesce(
      Brain_volume,              # 1st choice
      Brainvol,               # 2nd choice
      Total_brain_net_volume  # 3rd choice
    )
  )

#  --- Sanity check: did everyone get a value?
sum(is.na(data$Preferred_brain_volume))

# --- Parameters
options(scipen=999)

# --- Structures list
target_cols <- c(
  #"Diencephalon", "Medulla_oblongata",
  "LGN_Sousa", "Amygdala", "Pallidum", "NeoW_Frahm", "Total_insula_volume_L",
  "Nucleus_subthalamicus", "Capsula_interna", "Striatum", 
  "ASG_Sousa", "NeoG_Frahm", "Mesencephalon", "Cerebellum", "Hippocampus"
)

# Human-readable labels to match your figure
label_map <- c(
  #  "Diencephalon"          = "Diencephalon",
  #  "Medulla_oblongata"     = "Medulla oblongata",
  "Total_insula_volume_L" = "Insular cortex (grey)",
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
  "LGN_Sousa"             = "Corpus geniculatum laterale"
)

# --- Robust label helper to convert raw variable names → pretty names (safe with factors) ---
pretty_var <- function(x) {
  x_chr <- as.character(x)
  out <- unname(label_map[x_chr])
  out[is.na(out)] <- x_chr[is.na(out)]
  out
}
## Analysis A: Brain-weight predictor (BM vs Pagel ML)
# --- Helper function to transform VCV by Lambda
# (Multiplies off-diagonals by lambda, keeps diagonal variances constant)
transform_vcv_lambda <- function(tree, lambda) {
  v <- vcv(tree, corr = TRUE)
  diag_vals <- diag(v) # Save diagonal (variances)
  v <- v * lambda      # Shrink covariances
  diag(v) <- diag_vals # Restore diagonal
  return(v)
}

# --- LOOP for Plots I, II, IV ---
predict_human_one_var <- function(var_name,
                                  predictor_col = "Preferred_brain_volume",
                                  data,
                                  tr,
                                  human = "Homo_sapiens") {
  # Safety
  if (!var_name %in% names(data)) return(NULL)
  if (!predictor_col %in% names(data)) stop("predictor_col not found in data")
  
  # Prepare minimal data
  d1 <- data[, c("Species", predictor_col, var_name)]
  d1 <- as.data.frame(na.omit(d1))
  
  # --- CALCULATE REST OF BRAIN ---
  # ROB = Total - Part
  d1$Rest_of_Brain <- d1[[predictor_col]] - d1[[var_name]]
  
  # Safety check: if structure > brain (data error), ROB will be negative/zero, log will fail
  d1 <- d1[d1$Rest_of_Brain > 0, ]
  
  rownames(d1) <- d1$Species
  
  if (!human %in% d1$Species) return(NULL)
  
  # Observed human + predictor (now Rest_of_Brain)
  obs_val <- d1[d1$Species == human, var_name][1]
  xh_raw  <- d1[d1$Species == human, "Rest_of_Brain"][1]
  
  # Clean + align with tree (dispRity::clean.data)
  clean <- clean.data(d1, tr)
  tr1 <- clean$tree
  data2 <- clean$data
  if (!is.null(data2$Species)) rownames(data2) <- data2$Species
  
  # Fit without human
  tr2 <- drop.tip(tr1, human)
  data_model <- subset(data2, Species != human)
  N <- nrow(data_model)
  
  # --- NEW FORMULA: Uses Rest_of_Brain ---
  fmla <- as.formula(paste0("log(", var_name, ") ~ log(Rest_of_Brain)"))
  
  models_to_run <- list(
    list(name="Brownian (BM)", cor=corBrownian(1, form=~Species, phy=tr2)),
    list(name="Pagel's Lambda (ML)", cor=corPagel(1, form=~Species, phy=tr2, fixed=FALSE))
  )
  
  out <- list()
  
  for (mod in models_to_run) {
    
    # Mac note
    # tryCatch tries to run the code even when there is a problem, so it conceals the problem because you don't get a warning
    fit <- tryCatch(
      gls(fmla, 
          correlation = mod$cor, 
          data = data_model, 
          control = glsControl(opt = "optim", msMaxIter = 1000)), 
      error = function(e) NULL
    )    
    if (is.null(fit)) next
    
    sigma <- fit$sigma
    
    # lambda + VCVs
    if (mod$name == "Pagel's Lambda (ML)") {
      lambda_est <- as.numeric(coef(fit$modelStruct, unconstrained = FALSE))
      v.full    <- transform_vcv_lambda(tr1, lambda_est)
      v.reduced <- transform_vcv_lambda(tr2, lambda_est)
    } else {
      lambda_est <- 1.0
      v.full    <- vcv(tr1, corr = TRUE)
      v.reduced <- vcv(tr2, corr = TRUE)
    }
    
    # Phylogenetic conditional mean correction
    cc <- which(rownames(v.full) == human)
    others_names <- rownames(v.reduced)
    Cih <- v.full[cc, others_names, drop = FALSE]
    
    # IMPORTANT: X here is the RESPONSE vector (log trait), as in your code
    X <- log(data_model[[var_name]])
    xbar <- mean(X, na.rm = TRUE)
    
    inv_v_reduced <- solve(v.reduced)
    
    mu <- Cih %*% inv_v_reduced %*% (X - xbar)
    ch <- v.full[cc, cc] - Cih %*% inv_v_reduced %*% t(Cih)
    
    # Safety check: ch should be positive (variance can't be negative)
    # Small negative values are numerical errors - clamp to small positive

    if (ch <= 0) {
      ch <- 1  # Tweak CI based in lambda. If estimation fails, multiply CI by 1, so it doesn't change CI based on lambda 
    }
    # for example, we know observed is smaller than predicted but it cannot be negative variance
    
    xh <- log(xh_raw)
    predicted_log <- as.numeric(c(1, xh) %*% coef(fit) + mu)
    std.error <- sqrt(as.numeric(sigma^2 * ch))
    
    upper_log <- predicted_log + qnorm(0.975) * std.error
    lower_log <- predicted_log - qnorm(0.975) * std.error    
    out[[mod$name]] <- data.frame(
      Variable = var_name,
      Model = mod$name,
      Observed = as.numeric(obs_val),
      Lower = exp(lower_log),
      Predicted = exp(predicted_log),
      Upper = exp(upper_log),
      Lambda = lambda_est,
      stringsAsFactors = FALSE,
      N = nrow(data_model)
    )
  }
  
  if (length(out) == 0) return(NULL)
  do.call(rbind, out)
}

# Compute once for all structures (predictor = Rest_of_Brain)
core_df <- do.call(rbind, lapply(target_cols, predict_human_one_var,
                                 predictor_col = "Preferred_brain_volume",
                                 data = data,
                                 tr = tr))
core_df$Variable <- factor(core_df$Variable, levels = target_cols)
core_df$VarLabel <- pretty_var(core_df$Variable) #relabel

#### PLOT I: Prediction: Brownian Motion vs Pagel's Lambda (raw scale)

# --- PLOT ---
ggplot(core_df, aes(x = VarLabel)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, color = "#377eb8") +
  geom_point(aes(y = Predicted, shape = "Predicted"), color = "#377eb8", size = 2.5) +
  geom_point(aes(y = Observed, shape = "Observed"), color = "#e41a1c", size = 2.5) +
  scale_shape_manual(name = "Value", values = c("Predicted" = 16, "Observed" = 17)) +
  coord_flip() +
  facet_wrap(~Model, scales = "free_x") + # TWO PANELS comparison
  theme_bw() +
  labs(
    title = "Prediction: Brownian Motion vs Pagel's Lambda (Part-Whole Corrected)",
    subtitle = "Predictor = Rest of Brain (Total - Structure)",
    y = "Volume / Value (Original Scale)",
    x = "Brain Structure"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 11),
    axis.text.y = element_text(size = 9)
  )

#### PLOT II: Standard Prediction Error (0-1 Scale) (min–max scaling, for a given structure in dataset 0 = min , 1 = max)

# --- NORMALIZATION (0 to 1) ---
# Group by Variable to find the min/max range for that specific brain structure
# across BOTH models so the comparison remains valid.

final_df_norm <- core_df %>%
  group_by(Variable) %>%
  mutate(
    # Find global min and max for this variable (including CIs and Observed)
    local_min = min(c(Lower, Upper, Observed, Predicted)),
    local_max = max(c(Lower, Upper, Observed, Predicted)),
    
    # Apply Min-Max Scaling
    Predicted_Sc = (Predicted - local_min) / (local_max - local_min),
    Observed_Sc = (Observed - local_min) / (local_max - local_min),
    Lower_Sc = (Lower - local_min) / (local_max - local_min),
    Upper_Sc = (Upper - local_min) / (local_max - local_min)
  ) %>%
  ungroup()
final_df_norm$VarLabel <- pretty_var(final_df_norm$Variable)
# --- PLOT ---
# We use the _Sc (Scaled) columns now

ggplot(final_df_norm, aes(x = VarLabel)) +
  # Prediction Interval
  geom_errorbar(aes(ymin = Lower_Sc, ymax = Upper_Sc), width = 0.2, color = "#377eb8") +
  # Predicted Point
  geom_point(aes(y = Predicted_Sc, shape = "Predicted"), color = "#377eb8", size = 2.5) +
  # Observed Point
  geom_point(aes(y = Observed_Sc, shape = "Observed"), color = "#e41a1c", size = 2.5) +
  scale_shape_manual(name = "Value", values = c("Predicted" = 16, "Observed" = 17)) +
  coord_flip() +
  facet_wrap(~Model) + 
  theme_bw() +
  labs(
    title = "Standardized Prediction Error (0-1 Scale, Part-Whole Corrected)",
    subtitle = "0 = Min value for structure, 1 = Max value for structure (across CI and Obs). Predictor = Rest of Brain",
    y = "Standardized Position (0 to 1)",
    x = "Brain Structure"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 11),
    axis.text.y = element_text(size = 9)
  )

#### PLOT IV: Lambda Likelihood Profiles

# --- SETUP ---
profile_data <- list()
mle_data <- list()

# Define the sequence of Lambdas to test (Scan from 0 to 1)
lambda_seq <- seq(0, 1, by = 0.05)

# --- MAIN LOOP ---
for (var_name in target_cols) {
  
  if (!var_name %in% names(data)) next
  
  # A. Prepare Data
  d1 <- subset(data, select=c("Species", "Preferred_brain_volume", var_name))
  d1 <- as.data.frame(na.omit(d1))
  
  # --- CALCULATE REST OF BRAIN ---
  d1$Rest_of_Brain <- d1$Preferred_brain_volume - d1[[var_name]]
  
  # Safety check: if structure > brain (data error), ROB will be negative/zero, log will fail
  d1 <- d1[d1$Rest_of_Brain > 0, ]
  
  rownames(d1) <- d1$Species
  
  clean <- clean.data(d1, tr) 
  tr1 <- clean$tree
  data2 <- clean$data
  if(!is.null(data2$Species)) rownames(data2) <- data2$Species
  
  # Drop human for the fitting (consistent with your prediction model)
  tr2 <- drop.tip(tr1, "Homo_sapiens")
  data_model <- subset(data2, Species != "Homo_sapiens")
  
  # --- NEW FORMULA: Uses Rest_of_Brain ---
  fmla <- as.formula(paste("log(", var_name, ") ~ log(Rest_of_Brain)"))
  
  # B. Run Maximum Likelihood Estimate (Best point)
  fit_mle <- tryCatch({
    gls(fmla, correlation = corPagel(1, form=~Species, phy=tr2, fixed=FALSE), data = data_model, control = glsControl(opt = "optim", msMaxIter = 1000, msTol = 1e-6))
  }, error = function(e) NULL)
  
  if (!is.null(fit_mle)) {
    best_lambda <- as.numeric(coef(fit_mle$modelStruct, unconstrained = FALSE))
    best_loglik <- as.numeric(logLik(fit_mle))
    
    # Clamp visual to 0-1 range if estimation goes slightly out
    best_lambda_plot <- pmax(0, pmin(1, best_lambda))
    
    mle_data[[var_name]] <- data.frame(
      Variable = var_name,
      Lambda = best_lambda_plot,
      LogLik = best_loglik
    )
  }
  
  # C. Run Profile Scan (The Curve)
  for (val in lambda_seq) {
    fit_scan <- tryCatch({
      gls(fmla, correlation = corPagel(val, form=~Species, phy=tr2, fixed=TRUE), data = data_model, control = glsControl(opt = "optim", msMaxIter = 1000, msTol = 1e-6))
    }, error = function(e) NULL)
    
    if (!is.null(fit_scan)) {
      profile_data[[paste(var_name, val)]] <- data.frame(
        Variable = var_name,
        Lambda = val,
        LogLik = as.numeric(logLik(fit_scan))
      )
    }
  }
}

# --- COMBINE DATA ---
df_profile <- do.call(rbind, profile_data)
df_mle <- do.call(rbind, mle_data)

# --- PLOT ---
ggplot(df_profile, aes(x = Lambda, y = LogLik)) +
  
  # The Likelihood Curve
  geom_line(color = "#377eb8", size = 1) +
  
  # The MLE Point (Red Dot)
  geom_point(data = df_mle, aes(x = Lambda, y = LogLik), color = "red", size = 3) +
  
  # Vertical Line at MLE
  geom_vline(data = df_mle, aes(xintercept = Lambda), linetype = "dashed", color = "red", alpha=0.5) +
  
  # Text Label for the Best Lambda
  geom_text(data = df_mle, aes(label = sprintf("%.2f", Lambda), x = 0.1, y = LogLik), 
            color = "red", size = 3, hjust = 0, vjust = 1) +
  
  # Facet Wrap creates the panel
  facet_wrap(
    ~Variable, 
    scales = "free_y",
    labeller = labeller(Variable = function(x) pretty_var(x))
  ) +
  theme_bw() +
  labs(
    title = "Lambda Likelihood Profiles (Part-Whole Corrected)",
    subtitle = "Predictor = Rest of Brain. Curve shows model fit at different Lambda values. Red dot = MLE.",
    x = "Lambda (Phylogenetic Signal)",
    y = "Log-Likelihood"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 9),
    axis.text.y = element_text(size = 7)
  )

#### PLOT V: Human Brain Mosaicism (Corrected for Part-Whole) (ANOTHER New SECTION )

# --- SETUP ---
#### Analysis B: Part–whole–corrected predictions (Rest_of_Brain)
results_rob <- list()   # <-- DEFINE HERE

# --- MAIN LOOP ---
for (var_name in target_cols) {
  if (!var_name %in% names(data)) next
  
  # A. PREPARE DATA
  d1 <- subset(data, select=c("Species", "Preferred_brain_volume", var_name))
  d1 <- as.data.frame(na.omit(d1))
  
  # --- NEW: CALCULATE REST OF BRAIN ---
  # ROB = Total - Part
  # Ensure units match before running this!
  d1$Rest_of_Brain <- d1$Preferred_brain_volume - d1[[var_name]]
  
  # Safety check: if structure > brain (data error), ROB will be negative/zero, log will fail
  d1 <- d1[d1$Rest_of_Brain > 0, ]
  
  rownames(d1) <- d1$Species
  
  # B. EXTRACT HUMAN VALUES
  if ("Homo_sapiens" %in% d1$Species) {
    obs_val <- d1[d1$Species == "Homo_sapiens", var_name]
    
    # The predictor for humans is now Rest-of-Brain, not Total Brain
    xh_raw <- d1[d1$Species == "Homo_sapiens", "Rest_of_Brain"] 
  } else { next }
  
  clean <- clean.data(d1, tr) 
  tr1 <- clean$tree
  data2 <- clean$data
  if(!is.null(data2$Species)) rownames(data2) <- data2$Species
  
  tr2 <- drop.tip(tr1, "Homo_sapiens")
  data_model <- subset(data2, Species != "Homo_sapiens")
  
  # --- NEW FORMULA: Uses Rest_of_Brain ---
  fmla <- as.formula(paste("log(", var_name, ") ~ log(Rest_of_Brain)"))
  
  # C. MODELS (3-Panel Setup)
  models_to_run <- list(
    list(name="Independence (λ=0)", cor=corPagel(0, form=~Species, phy=tr2, fixed=TRUE)),
    list(name="Pagel's ML (Estimated)", cor=corPagel(1, form=~Species, phy=tr2, fixed=FALSE)),
    list(name="Brownian (λ=1)", cor=corBrownian(1, form=~Species, phy=tr2))
  )
  
  for (mod in models_to_run) {
    fit <- tryCatch({ gls(fmla, correlation = mod$cor, data = data_model, control = glsControl(opt = "optim", msMaxIter = 1000, msTol = 1e-6)) }, 
                    error = function(e) return(NULL))
    if(is.null(fit)) next
    
    sigma <- fit$sigma
    
    if (mod$name == "Pagel's ML (Estimated)") {
      lambda_est <- as.numeric(coef(fit$modelStruct, unconstrained = FALSE))
    } else if (mod$name == "Independence (λ=0)") {
      lambda_est <- 0
    } else {
      lambda_est <- 1.0 
    }
    
    v.full <- transform_vcv_lambda(tr1, lambda_est)
    v.reduced <- transform_vcv_lambda(tr2, lambda_est)
    
    cc <- which(rownames(v.full) == "Homo_sapiens")
    others_names <- rownames(v.reduced)
    Cih <- v.full[cc, others_names] 
    
    # Predictor vector X is now Log(Rest_of_Brain)
    X <- log(data_model$Rest_of_Brain)
    xbar <- mean(X, na.rm=TRUE)
    
    inv_v_reduced <- solve(v.reduced)
    mu <- Cih %*% inv_v_reduced %*% (X - xbar)
    
    # xh is also Log(Rest_of_Brain) for human
    xh <- log(xh_raw)
    predicted_log <- c(1, xh) %*% coef(fit) + mu 
    
    results_rob[[paste(var_name, mod$name)]] <- data.frame(
      Variable = var_name,
      Model = mod$name,
      Observed = obs_val,
      Predicted = exp(predicted_log),
      Lambda = lambda_est
    )
  }
}

final_df <- do.call(rbind, results_rob)
final_df$PropDiff <- (final_df$Observed - final_df$Predicted) / final_df$Predicted
final_df$Direction <- ifelse(final_df$PropDiff > 0, "Larger than Predicted", "Smaller than Predicted")
final_df$Model <- factor(final_df$Model, levels = c("Independence (λ=0)", "Pagel's ML (Estimated)", "Brownian (λ=1)"))
final_df$VarLabel <- pretty_var(final_df$Variable)
# --- PLOT ---
ggplot(final_df, aes(x = reorder(VarLabel, PropDiff), y = PropDiff, color = Direction)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_segment(aes(xend = reorder(VarLabel, PropDiff), yend = 0), size = 1.2) +
  geom_point(size = 4) +
  geom_text(aes(label = ifelse(Model == "Pagel's ML (Estimated)", 
                               sprintf("λ=%.2f", Lambda), 
                               "")), 
            nudge_x = -0.4, 
            size = 3, 
            fontface = "italic",
            color = "black", 
            show.legend = FALSE) +
  
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_color_manual(values = c("Larger than Predicted" = "#e41a1c", 
                                "Smaller than Predicted" = "#377eb8")) +
  
  coord_flip() +
  facet_wrap(~Model) +
  theme_bw() +
  labs(
    title = "Human Brain Mosaicism (Corrected for Part-Whole)",
    subtitle = "Predictor = Rest of Brain (Total - Structure).",
    y = "Deviation from Prediction (%)",
    x = "Brain Structure",
    color = "Direction"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    axis.text.y = element_text(size = 9, face="bold")
  )

#### TABLE 

# --- 1) Build core_df_out in one pipeline ---------------------------

digits_lambda <- 3
digits_vals   <- 5   # Predicted, CIs, Observed

core_df_out <- core_df %>%
  # 1) Column calculations
  mutate(
    Variable = factor(Variable, levels = target_cols),
    VarLabel = pretty_var(Variable),
    Diff.min = (Observed - Lower) / Observed,
    Diff.pre = (Observed - Predicted) / Observed,
    Diff.max = (Observed - Upper) / Observed
    ) %>%
  # 2) Shape table + renaming
  transmute(
    Structure       = VarLabel,
    Model,
    Lambda,
    `95% CI min`    = Lower,
    Predicted,
    `95% CI max`    = Upper,
    Observed,
    Diff.min,
    Diff.pre,
    Diff.max,
    N
  ) %>%
  # 3) Order rows by label order and model
  mutate(
    Structure = factor(Structure, levels = pretty_var(target_cols))
  ) %>%
  arrange(Structure, Model) %>%
  # 4) Limit significant figures & clean N
  mutate(
    Lambda       = signif(Lambda,       digits_lambda),
    `95% CI min` = signif(`95% CI min`, digits_vals),
    Predicted    = signif(Predicted,    digits_vals),
    `95% CI max` = signif(`95% CI max`, digits_vals),
    Observed     = signif(Observed,     digits_vals),
    N            = as.integer(N)
  )

# --- 2) Add rCMRGIc by matching Structure from Sup Table 1

# Read Supplementary Table 1
sup_tbl_raw <- read_excel(
  "/Users/crossmodal/Library/CloudStorage/Dropbox/COLLABORATIVE/Do expensive brain regions increase less in humans/MS Is human brain organization economical/Supplementary Figs and Tables/Sup Table 1 Regional cerebral metabolic rates and data sources using in Study 1.xlsx",
  sheet     = 1,
  col_names = FALSE
)

# Use row 3 as column names
colnames(sup_tbl_raw) <- as.character(unlist(sup_tbl_raw[3, ]))

# Drop rows 1–4 and keep only the relevant columns
sup_clean <- sup_tbl_raw[-c(1:4), ] %>%
  select(
    `Volume term`,
    `rCMRGIc (µmol/100 g/min.)`
  ) %>%
  rename(
    Structure = `Volume term`,
    rCMRGIc   = `rCMRGIc (µmol/100 g/min.)`
  )

# Join rCMRGIc into the main table by Structure
core_with_rCMRGIc <- core_df_out %>%
  left_join(
    sup_clean %>% select(Structure, rCMRGIc),
    by = "Structure"
  )

write.csv(core_with_rCMRGIc,"core_df_out_with_rCMRGIc.csv",row.names = FALSE)
write_xlsx(core_with_rCMRGIc,"core_df_out_with_rCMRGIc.xlsx")

