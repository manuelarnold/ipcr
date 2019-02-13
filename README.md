# Introduction

ipcr is an R package that makes it easy for researchers to study heterogeneity in the parameter estimates of a structural equation model (SEM)  with [individual parameter contribution (IPC) regression](https://psyarxiv.com/sbyux/). IPC regression allows regressing SEM parameter estimates on additional covariates.
ipcr requires that the SEM is fitted with R package [OpenMx](https://openmx.ssri.psu.edu/). As of now, only single-group RAM-type OpenMx models are supported. Moreover, the data used to estimate the SEM and the additional covariates need to be complete without any missing values.
This package is still under development. Please report any bugs.

# Installation
Use the following commands in R to install ipcr from GitHub:

```{r, eval=FALSE}
if (!require(devtools)) {install.packages("devtools")}
devtools::install_github("manuelarnold/ipcr")
```

# Example
``` r
# Specify an OpenMx model
m <- mxModel(model = "CFA",
             manifestVars = c("x1", "x2", "x3"),
             latentVars = "f",
             type = "RAM",
             mxData(observed = ipcr_data, type = "raw"),
             mxPath(from = "f", to = c("x1", "x2", "x3"), arrows = 1,
                    free = c(FALSE, TRUE, TRUE), values = 1,
                    labels = c("l1", "l2", "l3")),
             mxPath(from = "f", arrows = 2, free = TRUE, values = 0.75,
                    labels = "var_f"),
             mxPath(from = c("x1", "x2", "x3"), arrows = 2, free = TRUE,
                    values = 0.25, labels = c("e1", "e2", "e3")),
             mxPath(from = "one", to = c("x1", "x2", "x3"), arrows = 1,
                    free = FALSE, values = 0))

# Fit the model
fit <- mxTryHard(model = m)

# Investigate the parameter estimates with IPC regression
IPC_reg <- ipcr(fit = fit, covariates = ipcr_covariates, iterated = TRUE)

# Get an overview about parameter differences
plot(IPC_reg)

# Show IPC regression output
summary(object = IPC_reg)
```
