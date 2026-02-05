library(ape)
library(MASS)
library(nlme)
require(dispRity)


tr=read.tree("species.nwk")
data <- read.csv ("Stephan_primates.csv")
sapply(data, class)
###Subsetting

d1 <- subset(data, select=c("Species", "Brain_weight", "Striatum")) #The last var - the dependent - should be changed here, in the pgsl and X
d1 <- as.data.frame(na.omit(d1))
rownames(d1) <- d1$Species
clean <- clean.data(d1, tr)
tr1 <- clean$tree
data2 <- clean$data
#which(rownames(v.full) == "Homo_sapiens")
#data2 <- data2[-c(44), ]
data2 <- subset(data2, Species!="Homo_sapiens")
which(rownames(data2) == "Homo_sapiens")
rownames(data2) <- data2@Species

## Will predict "t1" therefore remove it from the full matrix
tr2 <- drop.tip(tr1, "Homo_sapiens")

#x <- as.data.frame(t(mvrnorm(2, rep(0, 50), Sigma=vcv(tr2, corr=TRUE)))) ## use reduced tree
#x$t <- rownames(x)

## Do gls analysis ###

fit <- gls(log(Striatum) ~ log(Brain_weight), correlation=corBrownian (1, form=~Species, phy=tr2), data=data2)

sigma <- fit$sigma  #get the res std err
lambda <-  fit$modelStruct$corStruct[1]
v.full <- vcv(tr1, corr=TRUE)
cc <- which(rownames(v.full) == "Homo_sapiens")

Cih <- v.full[cc, -cc]

xh <- 14.10069 ## Say our "human" had a trait value of 0.1

X <- log(data2$Striatum) ## covariate data for the "other species"
#X <- na.omit (X) 
#X <- as.numeric(X)
xbar <- mean(X, na.rm=TRUE)
mu <- Cih %*% solve(vcv(tr2, corr=TRUE)) %*% (X-xbar)

ch <- v.full[cc,cc] - Cih %*% solve(vcv(tr2, corr=TRUE)) %*% Cih

## Therefore, the predicted value for V1 is:

predicted <- c(1, xh) %*% coef(fit) + mu ## mu is the BIAS
std.error <- sqrt(sigma^2*ch)

upper <- predicted + qnorm(0.975)*std.error
lower <- predicted - qnorm(0.975) * std.error

c(exp(lower), exp(predicted), exp(upper))

## NB didn't account for lambda.


##PLOT

ggplot(data, aes(x = reorder(Area, n), y=predicted)) + 
  geom_errorbar(aes(ymin=min, ymax=max), colour="black", width=.1) +
  geom_point(size=2, shape=4, colour = "black") +
  geom_point(aes(y = observed), size = 2, shape = 1, colour = "red") +
  xlab("Area") +
  ylab("Observed and Predicted size with 95% CI (log)") +
  expand_limits(y=0) +                        # Expand y range
  scale_y_continuous(breaks=0:20*5) +     
  scale_x_discrete(breaks=0:10*1) + # Set tick every 4
  #scale_y_continuous(breaks = NULL) +
  #scale_x_discrete(breaks = NULL) +
  theme_bw() +
  theme(legend.justification=c(1,0),
        legend.position=c(1,0),
        text = element_text(size = 13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))
  
 

barplot(predicted_forplotting$diff.pre, names=predicted_forplotting$Area, col="#69b3a2",horiz=T , las=2)