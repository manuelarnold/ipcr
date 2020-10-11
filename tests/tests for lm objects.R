# Set Up -----
rm(list = ls())
library(ipcr)



# Generate data -----
x1 <- rnorm(1000, mean = -1, sd = sqrt(2))
x2 <- rnorm(1000, mean = -2, sd = sqrt(0.5))

z1 <- rep(c(0, 1), each = 500)
z2 <- rnorm(1000)
covariates <- data.frame(z1 = z1, z2 = z2)

y <- c(3 + 0.4 * x1[1:500], 3 + 0.7 * x1[501:1000])
y <- y + 0.1 * z2 * x2 + rnorm(1000)



# Fit models ----
m1 <- lm(y ~ x1)
m2 <- lm(y ~ x1 + x2)



# Test get_scores ----
scores_m1 <- get_scores(m1)
head(scores_m1)
colMeans(scores_m1)

scores_m2 <- get_scores(m2)
head(scores_m2)
colMeans(scores_m2)



# Test get_ipcs ----
# m1
ipcs_m1 <- get_ipcs(m1)
head(ipcs_m1)

colMeans(ipcs_m1)
coef(m1)

sqrt(1 / nrow(ipcs_m1)) * apply(ipcs_m1, MARGIN = 2, FUN = sd)
summary(m1)$coefficients[ ,2]
sqrt(diag(sandwich::sandwich(m1)))

# m2
ipcs_m2 <- get_ipcs(m2)
head(ipcs_m2)

colMeans(ipcs_m2)
coef(m2)

sqrt(1 / nrow(ipcs_m2)) * apply(ipcs_m2, MARGIN = 2, FUN = sd)
summary(m2)$coefficients[ ,2]
sqrt(diag(sandwich::sandwich(m2)))



# Test ipcr without covariates ----
#m1
ipcr_m1 <- ipcr(m1)
ipcr_m1
head(ipcr_m1$IPCs)
head(ipcs_m1)

#m2
ipcr_m2 <- ipcr(m2)
ipcr_m2
head(ipcr_m2$IPCs)
head(ipcs_m2)



# IPC regression
ipcr_m11 <- ipcr(m1, covariates = z1)
summary(ipcr_m11)
ipcr_m12 <- ipcr(m1, covariates = covariates)
summary(ipcr_m12)
ipcr_m13 <- ipcr(m1, covariates = covariates)
summary(ipcr_m13)
