#set.seed(3)

# Friends and stability

# 1 0.8
# 2 0.7
# 3 0.6
# 4 0.5
# 5 0.4
# 6 0.3

library(lavaan)

n <- 100

m_gen <- '
d2 ~ %s*d1
d3 ~ %s*d2
d4 ~ %s*d3
d5 ~ %s*d4

d1 ~~ 1*d1
d2 ~~ 0.6975*d2
d3 ~~ 0.6975*d3
d4 ~~ 0.6975*d4
d5 ~~ 0.6975*d5

d1 ~ 0*1
d2 ~ 0*1
d3 ~ 0*1
d4 ~ 0*1
d5 ~ 0*1
'

d <- rbind(
  lavaan::simulateData(model = sprintf(m_gen, 0.8, 0.8, 0.8, 0.8), sample.nobs = n),
  lavaan::simulateData(model = sprintf(m_gen, 0.7, 0.7, 0.7, 0.7), sample.nobs = n),
  lavaan::simulateData(model = sprintf(m_gen, 0.6, 0.6, 0.6, 0.6), sample.nobs = n),
  lavaan::simulateData(model = sprintf(m_gen, 0.5, 0.5, 0.5, 0.5), sample.nobs = n),
  lavaan::simulateData(model = sprintf(m_gen, 0.4, 0.4, 0.4, 0.4), sample.nobs = n),
  lavaan::simulateData(model = sprintf(m_gen, 0.3, 0.3, 0.3, 0.3), sample.nobs = n)
)

d <- scale(x = d, center = TRUE, scale = FALSE)

m <- '
d2 ~ beta*d1
d3 ~ beta*d2
d4 ~ beta*d3
d5 ~ beta*d4

d1 ~~ phi*d1
d2 ~~ psi*d2
d3 ~~ psi*d3
d4 ~~ psi*d4
d5 ~~ psi*d5

d1 ~ 0*1
d2 ~ 0*1
d3 ~ 0*1
d4 ~ 0*1
d5 ~ 0*1
'
library(OpenMx)

m <- mxModel(manifestVars = c("d1", "d2", "d3", "d4", "d5"),
             type = "RAM",
             mxData(observed = d, type = "raw"),
             mxPath(from = c("d1", "d2", "d3", "d4"), to = c("d2", "d3", "d4", "d5"),
                    connect = "single", arrows = 1, free = TRUE, values = 0.5,
                    labels = "beta"),
             mxPath(from = "d1", arrows = 2, free = TRUE, values = 1, labels = "phi"),
             mxPath(from = c("d2", "d3", "d4", "d5"), connect = "single", arrows = 2,
                    free = TRUE, values = 0.6975, labels = "psi"),
             mxPath(from = "one", to = c("d1", "d2", "d3", "d4", "d5"), arrows = 1,
                    free = FALSE, values = 0))

fit <- mxTryHard(m)
summary(fit)

Covariates <- data.frame(gender = rep(c(0, 1), times = 3*n),
                         contacts = rep(c(2,3,4,5,6,7), each = 100),
                         age = rnorm(n = 6*n, mean = 40, sd = 10))



res <- ipcr(fit, covariates = Covariates)
AIC(res)
coef(res)

coef(res, parameter = "phi")

coef(res)
summary(res)

predict(res)

sdf <- cv.glmnet(x = as.matrix(Covariates), y = res$IPCs$beta)
a <- coef(sdf, "lambda.min")
a@x
