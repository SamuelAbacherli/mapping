
# ------------------------------------------------------------------------------
# SET UP
# ------------------------------------------------------------------------------

rm(list = ls())

library(AER)
library(spdep)
library(spatGeoDa)
library(spatialreg)
library(psych)
library(MASS)
library(rowr)

# Load data
data(airbnb)
airbnb <- subset(airbnb, !is.na(airbnb$price_pp))  # remove rows where price_pp is missing

set.seed(1)

# ------------------------------------------------------------------------------
# CALCULATE STANDARD ERRORS MANUALLY
# ------------------------------------------------------------------------------

# Convert the polygon data to a listw object with queen contiguity and row-standardised
nb <- poly2nb(airbnb, queen = TRUE)
W <- nb2listw(nb, style = "W")

# Compute a SAR model based on the airbnb data and its contiguity matrix W
lm1 <- price_pp ~ unemployed + poverty + without_hs
sar1 <- lagsarlm(lm1, airbnb, W, tol.solve = 1.0e-30)

# Estimate effects manually
rho <- as.numeric(coefficients(sar1)["rho"])
beta.hat <- coefficients(sar1)[-1]
N <- nrow(airbnb)
I <- diag(N)
I_rhoWinv <- solve(I - rho * listw2mat(W))

# Build matrix of partial derivatives S for the beta variable unemployed
S_W.unemployed <- I_rhoWinv %*% (I * beta.hat["unemployed"])
S_W.poverty <- I_rhoWinv %*% (I * beta.hat["poverty"])
S_W.without_hs <- I_rhoWinv %*% (I * beta.hat["without_hs"])

# Direct effect is the sum of the diagonal divded by N
Direct.unemployed <- tr(S_W.unemployed) / N
Direct.poverty <- tr(S_W.poverty) / N
Direct.without_hs <- tr(S_W.without_hs) / N

# Total effect is the sum of all derivatives
Total.unemployed <- sum(S_W.unemployed) / N
Total.poverty <- sum(S_W.poverty) / N
Total.without_hs <- sum(S_W.without_hs) / N

# Indirect effects is their difference between total and direct effect
Indirect.unemployed <- Total.unemployed - Direct.unemployed
Indirect.poverty <- Total.poverty - Direct.poverty
Indirect.without_hs <- Total.without_hs - Direct.without_hs

effects <- cbind("Direct" = c(Direct.unemployed, Direct.poverty, Direct.without_hs),
                "Indirect" = c(Indirect.unemployed, Indirect.poverty, Indirect.without_hs),
                "Total" = c(Total.unemployed, Total.poverty, Total.without_hs))
rownames(effects) <- c("Unemployed", "Poverty", "Without_HS") ; effects

# Compare with automatically calculated effects
impacts <- impacts(sar1, listw = W, R = 100)
summary(impacts, zstats = T, short = T)

# Draw parameters from a multivariate normal distribution and recalulate impacts
R <- mvrnorm(n = 100, mu = c("rho" = rho, beta.hat), Sigma = sar1$resvar[-1,-1])
M <- matrix(rep(0, 100*9), nrow = 100)
colnames(M) <- c("Direct Unemployed", "Indirect Unemployed", "Total Unemployed",
                 "Direct Poverty", "Indirect Poverty", "Total Poverty",
                 "Direct Without_HS", "Indirect Without_HS", "Total Without_HS")

# Calculating the impacts for each random sample of coefficients
for (r in 1:nrow(R)) {
  rho <- as.numeric(R[r, 1])
  beta.hat <- as.numeric(R[r, 3:5])
  N <- nrow(airbnb)
  I <- diag(N)
  I_rhoWinv <- solve(I - rho * listw2mat(W))
  
  # Build matrix of partial derivatives S for the beta variable unemployed
  S_W.unemployed <- I_rhoWinv %*% (I * beta.hat[1])
  S_W.poverty <- I_rhoWinv %*% (I * beta.hat[2])
  S_W.without_hs <- I_rhoWinv %*% (I * beta.hat[3])
  
  # Direct effect is the sum of the diagonal divded by N
  M[r, "Direct Unemployed"] <- tr(S_W.unemployed) / N
  M[r, "Direct Poverty"] <- tr(S_W.poverty) / N
  M[r, "Direct Without_HS"] <- tr(S_W.without_hs) / N
  
  # Total effect is the sum of all derivatives
  M[r, "Total Unemployed"] <- sum(S_W.unemployed) / N
  M[r, "Total Poverty"] <- sum(S_W.poverty) / N
  M[r, "Total Without_HS"] <- sum(S_W.without_hs) / N
  
  # Indirect effects is their difference between total and direct effect
  M[r, "Indirect Unemployed"] <- M[r, "Total Unemployed"] - M[r, "Direct Unemployed"]
  M[r, "Indirect Poverty"] <- M[r, "Total Poverty"] - M[r, "Direct Poverty"]
  M[r, "Indirect Without_HS"] <- M[r, "Total Without_HS"] - M[r, "Direct Without_HS"]
}

# Calculating the standard erros of the impacts
se <- matrix(c(sd(M[,1]), sd(M[,2]), sd(M[,3]), sd(M[,4]), sd(M[,5]), sd(M[,6]), 
               sd(M[,7]), sd(M[,8]), sd(M[,9])), nrow = 3, byrow = T)
colnames(se) <- c("Direct", "Indirect", "Total")
rownames(se) <- c("unemployed", "poverty", "without_hs")
se

# ------------------------------------------------------------------------------
# ESTIMATE A S2SLS BY HAND
# ------------------------------------------------------------------------------

# Define the dependent variable
y <- c(airbnb$price_pp)

# Create instruments using lag.listw()
X <- cbind(airbnb$unemployed, airbnb$poverty, airbnb$without_hs)
WX <- lag.listw(W, X)
W2X <- lag.listw(W, WX)

# Use ivreg() from the AER package
round(ivreg(formula = y ~ listw2mat(W) %*% y + X | X + WX + W2X)$coeff, 4)

# Use stsls() from the spatialreg package
round(stsls(formula = y ~ X, listw = W)$coeff, 4)


# ------------------------------------------------------------------------------
# ESTIMATING SPATIAL MODELS
# ------------------------------------------------------------------------------

# OLS
ols <- lm(lm1, data = airbnb)

# SAR
sar_ml <- lagsarlm(lm1, airbnb, W, tol.solve = 1.0e-30)

# SLX
slx <- lmSLX(lm1, airbnb, W)

# SEM
sem_ml <- errorsarlm(lm1, airbnb, W, tol.solve = 1.0e-30)
sem_gm <- GMerrorsar(lm1, airbnb, W)

# SAC
sac_ml <- sacsarlm(lm1, airbnb, W, tol.solve = 1.0e-30)
sac_gmm <- gstsls(lm1, airbnb, W) ; names(sac_gmm$coefficients)[1] <- "rho"

# SDM
sdm <- lagsarlm(lm1, airbnb, W, tol.solve = 1.0e-30, Durbin = T)

# SDEM
sdem_ml <- errorsarlm(lm1, airbnb, W, tol.solve = 1.0e-30, Durbin = T)

# GNS
gns_ml <- sacsarlm(lm1, airbnb, W, tol.solve = 1.0e-30, Durbin = T)

# Merger of all results the hard way to compare results
x <- merge(data.frame("GNS" = coef(gns_ml)),
           data.frame("SDEM" = coef(sdem_ml)),
           by = "row.names", all.x = T)
rownames(x) <- x[,1] ; x <- x[, -1]

x <- merge(x, data.frame("SAC_GMM" = coef(sac_gmm)),
           by = "row.names", all.x = T)
rownames(x) <- x[,1] ; x <- x[, -1]

x <- merge(x, data.frame("SAC_ML" = coef(sac_ml)),
           by = "row.names", all.x = T)
rownames(x) <- x[,1] ; x <- x[, -1]

x <- merge(x, data.frame("SEM_GM" = coef(sem_gm)),
           by = "row.names", all.x = T)
rownames(x) <- x[,1] ; x <- x[, -1]

x <- merge(x, data.frame("SEM_ML" = coef(sem_ml)),
           by = "row.names", all.x = T)
rownames(x) <- x[,1] ; x <- x[, -1]

x <- merge(x, data.frame("SLX" = coef(slx)),
           by = "row.names", all.x = T)
rownames(x) <- x[,1] ; x <- x[, -1]

x <- merge(x, data.frame("SAR" = coef(sar_ml)),
           by = "row.names", all.x = T)
rownames(x) <- x[,1] ; x <- x[, -1]

x <- merge(x, data.frame("OLS" = coef(ols)),
           by = "row.names", all.x = T)
rownames(x) <- x[,1] ; x <- x[, -1]

x  # the results



