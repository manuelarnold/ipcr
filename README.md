# Introduction

ipcr is an R package for predicting and explaining individual differences in model parameters with [individual parameter contribution (IPC) regression](https://doi.org/10.1080/10705511.2019.1667240). IPC regression allows regressing model parameters on covariates. IPC regression can be used as an alternative to other methods such as random-effects models or multi-group models.

ipcr was mainly written for structural equation models (SEMs) estimated with the [lavaan](https://lavaan.ugent.be/) or [OpenMx](https://openmx.ssri.psu.edu/) package. However, icpr can also be used to investigate models fitted with R's lm and glm function.

As of now, the data used to estimate the model and covariates need to be complete without any missing values.

This package is still under development. Please report any bugs.

# Installation
Use the following commands in R to install ipcr from GitHub:

```{r, eval=FALSE}
if (!require(devtools)) {install.packages("devtools")}
devtools::install_github("manuelarnold/ipcr")
```

# Example
``` r
# Structural equation model example using the lavaan package

## Load Holzinger and Swineford (1939) data provided by the lavaan package
HS_data <- lavaan::HolzingerSwineford1939

## Remove observations with missing values
HS_data <- HS_data[stats::complete.cases(HS_data), ]

## lavaan model syntac for a single group model
m <- 'visual =~ x1 + x2 + x3
      textual =~ x4 + x5 + x6
      speed =~ x7 + x8 + x9'

## Fit the model
fit <- lavaan::cfa(model = m, data = HS_data)

## Prepare a data.frame with covariates
covariates <- HS_data[, c("sex", "ageyr", "agemo", "school", "grade")]

## Regress parameters on covariates with the ipcr function
res <- ipcr(fit = fit, covariates = covariates)

## Plot heatmap with the correlation between parameters and predictors
plot(res)

## Show results (standard IPC regression)
summary(res)

## IPC regression with LASSO regularization
res_reg <- ipcr(fit = fit, covariates = covariates, regularization = TRUE)

## Show results (regularized standard IPC regression)
summary(res_reg)
```
