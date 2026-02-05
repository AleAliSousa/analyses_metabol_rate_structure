# this script was adapted from the script saved on 6 Nov 2020 named "original with humans/predicValuesPGLS.R" 

library(ape)
library(MASS)
library(nlme)
library(dispRity)
library(ggplot2)

# --- Load tree and data
tr=read.tree("species.nwk")
Stephan_primates <- read.csv ("Stephan_primates.csv")

# --- Clean columns
colnames(Stephan_primates) # Check column names first
data <- subset(Stephan_primates, select = -c(X, order))

# --- Make list of structures
structures <- list(
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

setdiff(structures, colnames(data)) #check for those columns in data

# --- Prepare one dependent variable (Striatum) for illustration
d1 <- subset(data, select=c("Species", "Brain_weight", "Striatum")) #The last var - the dependent - should be changed here, in the pgsl and X
d1 <- as.data.frame(na.omit(d1))
rownames(d1) <- d1$Species

# --- Ensure tree and data match
clean <- clean.data(d1, tr) #same species in data and tree
tr1 <- clean$tree
data2 <- clean$data
setdiff(d1$Species, rownames(data2)) #check which species were dropped

# --- Check the observed values for Homo_sapiens
data[data$Species == "Homo_sapiens", c("Brain_weight", "Striatum")]

# --- Remove Homo_sapiens from data and tree (we will predict it)
tr2 <- drop.tip(tr1, "Homo_sapiens")
data2 <- subset(data2, Species!="Homo_sapiens") # remove Homo sapiens
which(rownames(data2) == "Homo_sapiens") # check that Homo sapiens was removed

## Do gls analysis ###

# --- Predictor value for Homo 
xh <- log(data[data$Species == "Homo_sapiens", "Brain_weight"]) # x for human

# --- Fit GLS with Pagel's lambda estimated
fit <- gls(log(Striatum) ~ log(Brain_weight), correlation=corBrownian (1, form=~Species, phy=tr2), data=data2)

# --- Extract residual SD and estimated Pagel's lambda
sigma <- fit$sigma  #get the res std err
lambda <-  fit$modelStruct$corStruct[1]
v.full <- vcv(tr1, corr=TRUE)
cc <- which(rownames(v.full) == "Homo_sapiens")

Cih <- v.full[cc, -cc] # Covariance of species with humans

# --- Response vector for "other species" on log scale, aligned to tree tip order
X <- log(data2$Striatum) ## covariate data for the "other species"
is.numeric(X) # Check if numeric
anyNA(X) # Check if there are any NA values
xbar <- mean(X, na.rm=TRUE)

# --- Bias (conditional mean) and conditional variance with λ
mu <- Cih %*% solve(vcv(tr2, corr=TRUE)) %*% (X-xbar) # conditional mean

ch <- v.full[cc,cc] - Cih %*% solve(vcv(tr2, corr=TRUE)) %*% Cih # conditional variance

# --- Linear predictor (log scale) + bias
predicted <- c(1, xh) %*% coef(fit) + mu ## mu is the BIAS
std.error <- sqrt(sigma^2*ch)
# --- 95% CI (log scale)
upper <- predicted + qnorm(0.975)*std.error
lower <- predicted - qnorm(0.975) * std.error
c(exp(lower), exp(predicted), exp(upper))


## NB didn't account for lambda.


# --- Compare predicted vs observed on original scale
exp(predicted)  # back-transform predicted
data[data$Species == "Homo_sapiens", "Striatum"]

# --- Back-transform to original Striatum scale and compare with observed
predicted_unlog <- as.numeric(exp(predicted))
lower_unlog     <- as.numeric(exp(lower))
upper_unlog     <- as.numeric(exp(upper))
observed_unlog  <- as.numeric(data[data$Species == "Homo_sapiens", "Striatum"])

cat("\nPredicted (orig units):", predicted_unlog,
    "\nObserved  (orig units):", observed_unlog,
    "\n95% CI (orig units):   [", lower_unlog, ", ", upper_unlog, "]\n")

# --- Minimal plot (one point with error bar)
plot_df <- data.frame(
  Structure = "Striatum",
  Predicted = predicted_unlog,
  Observed  = observed_unlog,
  Lower     = lower_unlog,
  Upper     = upper_unlog
)


#Plot with original units

ggplot(plot_df, aes(x = Structure, y = Predicted)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.15, colour = "black") +
  geom_point(size = 2, shape = 4, colour = "black") +
  geom_point(aes(y = Observed), size = 2, shape = 1, colour = "red") +
  xlab("Structure") +
  ylab("Observed and Predicted size with 95% CI (original units)") +
  theme_bw()

### ### 
###  these values do not match the output files such as "original with humans/predicted.csv"
###  ### 


# 
# ##PLOT all structures on one figure
# ggplot(data, aes(x = reorder(Area, n), y=predicted)) + 
#   geom_errorbar(aes(ymin=min, ymax=max), colour="black", width=.1) +
#   geom_point(size=2, shape=4, colour = "black") +
#   geom_point(aes(y = observed), size = 2, shape = 1, colour = "red") +
#   xlab("Area") +
#   ylab("Observed and Predicted size with 95% CI (log)") +
#   expand_limits(y=0) +                        # Expand y range
#   scale_y_continuous(breaks=0:20*5) +     
#   scale_x_discrete(breaks=0:10*1) + # Set tick every 4
#   #scale_y_continuous(breaks = NULL) +
#   #scale_x_discrete(breaks = NULL) +
#   theme_bw() +
#   theme(legend.justification=c(1,0),
#         legend.position=c(1,0),
#         text = element_text(size = 13),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text.x = element_text(angle = 90, hjust = 1))
# 
# barplot(predicted_forplotting$diff.pre, names=predicted_forplotting$Area, col="#69b3a2",horiz=T , las=2)