rm(list = ls())

#### Load ipcr function ####
library(ipcr)


#### Generate Data ####
x <- rnorm(1000, mean = -1, sd = sqrt(2))
y <- c(3 + 0.4 * x[1:500], 3 + 0.7 * x[501:1000])
y <- y + rnorm(1000)
cov <- c(rep(c(0, 1), each = 500))
noise <- rep(c(0,1), times = 500)

data_model <- data.frame(x, y)
covariates <- data.frame(cov, noise)

conv = 0.0001
max_it = 50
linear = TRUE
IPC <- list("info" = list(name = deparse(substitute(x)),
                          class = class(x)))

#### OpenMx ####
# 2 variables mean cov
m2mc_mx <- OpenMx::mxModel(manifestVars = c("x", "y"),
                           type = "RAM",
                           OpenMx::mxData(observed = data_model, type = "raw"),
                           OpenMx::mxPath(from = "x", to = "y", arrows = 1, free = TRUE,
                                  values = 0.5, labels = "beta"),
                           OpenMx::mxPath(from = c("x", "y"), connect = "single", arrows = 2,
                                  free = TRUE, values = 1, labels = c("var_x", "var_y")),
                           OpenMx::mxPath(from = "one", to = c("x", "y"), arrows = 1, free = TRUE,
                                  values = 0, labels = c("m_x", "m_y")))
x_mx <- OpenMx::mxTryHard(model = m2mc_mx)

m2mc_lav <- '
y ~ x
y ~~ y
x ~~ x
x ~ 1
y ~ 1
'
x_lav <- lavaan(model = m2mc_lav, data = data_model)


res1_mx <- ipcr3(fit = x_mx, covariates = covariates, iterated = TRUE, linear = TRUE)
res2_mx <- ipcr3(fit = x_mx, covariates = covariates, iterated = TRUE, linear = FALSE)

res1_lav <- ipcr3(fit = x_lav, covariates = covariates, iterated = TRUE)
summary(res1_lav)

sapply(X = res1$regression[["beta"]], FUN = print)

summary(res1, method = "standard")


# 2 variables cov
m2c <- mxModel(manifestVars = c("x", "y"),
                type = "RAM",
                mxData(observed = data_model, type = "raw"),
                mxPath(from = "x", to = "y", arrows = 1, free = TRUE,
                       values = 0.5, labels = "beta"),
                mxPath(from = c("x", "y"), connect = "single", arrows = 2,
                       free = TRUE, values = 1, labels = c("var_x", "var_y")),
                mxPath(from = "one", to = c("x", "y"), arrows = 1, free = FALSE,
                       values = 0, labels = c("m_x", "m_y")))
x <- mxTryHard(model = m2c)

res1 <- ipcr3(fit = x, covariates = covariates, iterated = TRUE, linear = TRUE)
res2 <- ipcr3(fit = x, covariates = covariates, iterated = TRUE, linear = FALSE)

fitted(res1, parameter = "beta")
fitted.values(res1)


predict(res2)
sapply(X = res1$regression, FUN = predict, newdata = data.frame(cov = 3, noise = 3))

predict(res1, method = "standard", parameter = "beta")

predict(res1, method = "standard", parameter = "beta", newdata = data.frame(cov = 3, noise = 3))
nobs(res1)

# 1 variables mean cov
m1mc <- mxModel(manifestVars = "x",
               type = "RAM",
               mxData(observed = data_model, type = "raw"),
               mxPath(from = "x", connect = "single", arrows = 2,
                      free = TRUE, values = 1, labels = "var_x"),
               mxPath(from = "one", to = "x", arrows = 1, free = TRUE,
                      values = 0, labels = "m_x"))
x <- mxTryHard(model = m1mc)

res1 <- ipcr3(fit = x, covariates = covariates, iterated = TRUE, linear = TRUE)
res2 <- ipcr3(fit = x, covariates = covariates, iterated = TRUE, linear = FALSE)


# 1 variables cov
m1c <- mxModel(manifestVars = "x",
                type = "RAM",
                mxData(observed = data_model, type = "raw"),
                mxPath(from = "x", connect = "single", arrows = 2,
                       free = TRUE, values = 1, labels = "var_x"),
                mxPath(from = "one", to = "x", arrows = 1, free = FALSE,
                       values = 0, labels = "m_x"))
x <- mxTryHard(model = m1c)

res1 <- ipcr3(fit = x, covariates = covariates, iterated = TRUE, linear = TRUE)
res2 <- ipcr3(fit = x, covariates = covariates, iterated = TRUE, linear = FALSE)


conv = 0.0001
max_it = 50
linear = FALSE
IPC <- list("info" = list(name = deparse(substitute(x)),
                          class = class(x)))


res <- iterated_ipcr.MxRAMModel3(x, IPC, covariates, conv, max_it, linear = FALSE)

res_ipcr3(fit = x, covariates = covariates, iterated = TRUE, linear = FALSE)

res_ipcr2 <- ipcr2(fit = x, covariates = covariates, iterated = TRUE)
res_ipcr3 <- ipcr3(fit = x, covariates = covariates, iterated = TRUE, linear = TRUE)

initial_ipc <- ipcr3(fit = x, covariates = covariates)
initial_ipc$regression

system.time(ipcr2(fit = x, covariates = covariates, iterated = TRUE))
system.time(ipcr3(fit = x, covariates = covariates, iterated = TRUE, linear = TRUE))

res_ipcr2
res_ipcr3$iterated_regression$beta


