#- Estimating OLS
library(MASS)
library(tidyverse)
library(lmtest)
library(imager)
data(Boston)

## Task 1
# Assign variables
y <- Boston$medv
X <- Boston %>% 
  dplyr::select(crim, indus, age) %>% # select dependent variables
  as.matrix()  # convert from data.frame to matrix
X <- cbind(intercept = 1, X)  # ad the intercept

# Linear model without weights
m.lm <- lm(medv ~ crim + indus + age, data = Boston)
bptest(m.lm)  # heteroskedasticity at the 10% significance level
dwtest(m.lm)  # there is autocorrelation in the model

# Plot of residuals versus fitted values
plot(m.lm$residuals ~ m.lm$fitted.values, pch = 19, cex = 0.7, col = rgb(0,0,0,0.5),
     main = "Heteroskedasticity 1", xlab = "Fitted Values", ylab = "Residuals")

# Creating weight matrix with zn using O = P'P
P <- diag(sqrt(Boston$zn / sum(Boston$zn)))

# Adjust the data with the weight matrix
y.star <- P %*% y
X.star <- P %*% X

# Use ordinary OLS on the transformed data
m.est <- solve(t(X.star) %*% X.star) %*% t(X.star) %*% y.star %>% 
  as.vector()
names(m.est) <- c("(Intercept)", "crim", "indus", "age") ; m.est

# Check the heteroskedasticity
d.star <- as.data.frame(cbind(y.star, X.star))
names(d.star)[1] <- "medv"  # name the y variable
m.weighted <- lm(medv ~ crim + indus + age, data = d.star)
bptest(m.weighted)  # the weight matrix did not get rid of the heteroskedasticity, especially as 73.5% of the observations were adjusted to 0
dwtest(m.weighted)  # the weighting also didn't get rid of autocorrelation

# Beware that 73.5% of all observations lie at (0,0)
plot(m.weighted$residuals ~ m.weighted$fitted.values, pch = 19, cex = 0.7, 
     col = rgb(0,0,0,0.5), main = "Heteroskedasticity 2", 
     xlab = "Fitted Values", ylab = "Residuals")


## Task 2
# Use ordinary ML on the transformed data
logL <- function(par, XX = X.star, yy = y.star) {
  sigma <- par[length(par)]
  beta <- par[-length(par)]
  sum(-0.5*log(sigma) - 0.5*log(2*pi) - 
        (1/(2*sigma)) * (yy - XX %*% matrix(beta))^2)
}

# Optimization function using the logL function
r <- optim(rep(0.1, 5), logL, method = "BFGS",
           control = list(maxit = 100000, abstol = 10^-6, fnscale = -1))  # use negative fnscale to maximize instead of minimize
names(r$par) <-  c("(Intercept)", "crim", "indus", "age", "sigma") ; r$par[1:4]


## Task 3
ref.point <- load.image("data/Hermannskogel.png")
plot(ref.point, xaxt = "n", axes = F, main = "Reference Point of MGI Austria")
