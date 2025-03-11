# ipcr: Individual Parameter Contribution Regression

## Introduction

`ipcr` is an R package for predicting and explaining individual differences in model parameters using [individual parameter contribution regression (IPCR)](https://doi.org/10.1080/10705511.2019.1667240). IPCR allows estimated model parameters to be regressed on predictors, enabling the study of parameter heterogeneity. It serves as an alternative to methods such as random-effects models and multi-group models.

`ipcr` is primarily designed for **structural equation models (SEMs)** estimated with [`lavaan`](https://lavaan.ugent.be/) or [`OpenMx`](https://openmx.ssri.psu.edu/). However, it can also be applied to models fitted using base R’s [`lm()`](https://rdrr.io/r/stats/lm.html) and [`glm()`](https://rdrr.io/r/stats/glm.html) functions, as well as [`lmer()`](https://cran.r-project.org/package=lme4) from the `lme4` package.

This package is still under development. If you encounter any bugs, please report them in the [GitHub Issues](https://github.com/manuelarnold/ipcr/issues) section.

**Note:** This is the development branch, which reworks `ipcr` version 0.2 as described in [individual parameter contribution regression (IPCR)](https://doi.org/10.3390/psych3030027}  and introduces new features. Currently, iterated IPCR (`ipcr_it`) and regularized IPCR (`ipcr_reg`) are not functional.

---

## Installation

You can install `ipcr` from GitHub using the `remotes` package:

```r
install.packages("remotes") # Install remotes if not already installed
remotes::install_github("manuelarnold/ipcr")
```


## Example

```r
# Load Holzinger and Swineford (1939) dataset from the lavaan package
HS_data <- lavaan::HolzingerSwineford1939

# Remove observations with missing values to ensure complete data
HS_data <- HS_data[stats::complete.cases(HS_data), ]

# Define a confirmatory factor analysis (CFA) model using lavaan syntax
# This model specifies three latent factors: visual, textual, and speed
m <- 'visual  =~ x1 + x2 + x3
      textual =~ x4 + x5 + x6
      speed   =~ x7 + x8 + x9'

# Fit the CFA model
fit <- lavaan::cfa(model = m, data = HS_data)

# Select predictors for parameter heterogeneity analysis
predictors <- HS_data[, c("sex", "ageyr", "agemo", "school", "grade")]

# Perform Individual Parameter Contribution Regression (IPCR)
res <- ipcr(fit = fit, predictors = predictors)

# Plot a heatmap showing correlations between parameters and predictors
#plot(res) # not functional

# Display a summary of the IPC regression results
summary(res)

# Perform IPC regression with LASSO regularization
#res_reg <- ipcr(fit = fit, predictors = predictors, regularization = TRUE) # not functional

# Display results of regularized IPC regression
#summary(res_reg) # not functional
```
