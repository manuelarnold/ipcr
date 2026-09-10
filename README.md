# ipcr: Individual Parameter Contribution Regression

## Introduction

**ipcr** is an R package for predicting and explaining individual differences in model parameters using [individual parameter contribution regression (IPCR)](https://doi.org/10.1080/10705511.2019.1667240). IPCR allows estimated model parameters to be regressed on predictors, enabling the study of parameter heterogeneity. It serves as an alternative to methods such as random-effects models and multi-group models.

The primary **ipcr** backends are structural equation models estimated with [`OpenMx`](https://openmx.ssri.psu.edu/) and [`lavaan`](https://lavaan.ugent.be/). Additional officially supported backends are base R's [`lm()`](https://rdrr.io/r/stats/lm.html) and [`glm()`](https://rdrr.io/r/stats/glm.html) functions and `lmer()`/`glmer()` from the [`lme4`](https://cran.r-project.org/package=lme4) package. Mixed-model IPCs currently use one row per level of a single grouping factor; `lmer()` models must use `REML = FALSE`, and `glmer()` support is limited to binomial and Poisson models. Other fitted-model classes are attempted through compatible `sandwich` methods with a warning but are not officially supported.

lavaan support includes continuous ML, GLS, ULS, and WLS models, multigroup models, FIML, and simple linear equality constraints. Ordered indicators currently require full WLS because lavaan does not provide compatible casewise scores for the default WLSMV estimator. Unsupported estimators and nonlinear or inequality constraints are rejected.

OpenMx support requires a meaningful row-wise likelihood based on raw data. RAM and non-RAM models use numerical OpenMx row gradients by default. Continuous raw-data RAM models can optionally use faster analytical scores; these support missing manifest observations, definition variables in the A, S, and M matrices, and equality restrictions represented by shared labels. For standard `ipcr()`, a request for analytical scores falls back to numerical scores with a warning when the RAM model contains unsupported features such as explicit `mxConstraint()` objects or algebras. Analytical `ipcr_it()` instead rejects unsupported RAM features; users can explicitly select `analytic = FALSE` for the much slower numerical iteration.

The predictor argument is consistently named `predictors`. It accepts numeric, logical, and factor predictors; character columns must first be converted to factors or numeric indicators. Predictor rows are matched positionally to IPC units. If the first-stage model omits observations, `predictors` must contain only the retained IPC units, in the order used by the fitted model. Rows with missing predictor values are omitted from the second-stage regression, while the complete IPC and score matrices remain available in the result. In an iterated analysis, IPCs and scores for these excluded rows retain their standard-IPCR values and are not recalculated.

Second-stage summaries include an overall F-test and R-squared value for every IPC response. When multiple IPC responses are present and the multivariate test is estimable, `car::Manova()` additionally provides Type-II tests using Pillai's trace, so each predictor term is tested conditionally on the remaining terms rather than sequentially. If the residual IPC response matrix is rank deficient or the multivariate calculation otherwise fails, the MANCOVA is unavailable but the separate IPC regressions and their F-tests remain available.

The returned `IPCs` and `scores` components are numeric matrices. The result also retains the full normalized predictor data in `predictors` and a logical `complete_rows` vector identifying the rows used by the second-stage analysis. The `get_ipcs()` and `get_scores()` convenience functions return the same matrix representation.

This package is research software under active development. If you encounter any bugs, please report them in the [GitHub Issues](https://github.com/manuelarnold/ipcr/issues) section.

Version 0.4 substantially reworks **ipcr** version 0.2 as described in [individual parameter contribution regression (IPCR)](https://doi.org/10.3390/psych3030027) and introduces new features. Iterated IPCR (`ipcr_it()`) is available for supported lavaan and raw-data OpenMx models. OpenMx RAM models can use fast analytical scores; general OpenMx models use a much slower numerical fallback. Regularized IPCR (`ipcr_reg()`) applies a cross-validated multivariate group elastic net to the IPCs from a standard or iterated IPCR result.

---

## Installation

You can install `ipcr` from GitHub using the `remotes` package:

```r
# Install remotes package if not already installed
install.packages("remotes")

# Install ipcr package
remotes::install_github("manuelarnold/ipcr")

# Load ipcr package
library(ipcr)
```

---

## Supported backend examples

### OpenMx

OpenMx is a primary backend. Numerical row scores are the general default;
this simple raw-data RAM model also supports the faster analytical scores.

```r
set.seed(1)
n <- 60L
openmx_data <- data.frame(
  y = stats::rnorm(n),
  group = rep(c(0, 1), each = n / 2)
)

openmx_model <- OpenMx::mxModel(
  "one_variable",
  type = "RAM",
  manifestVars = "y",
  OpenMx::mxPath(
    from = "y", arrows = 2, free = TRUE, values = 1,
    labels = "variance", lbound = 0.0001
  ),
  OpenMx::mxPath(
    from = "one", to = "y", arrows = 1, free = TRUE,
    values = 0, labels = "mean"
  ),
  OpenMx::mxData(openmx_data["y"], type = "raw")
)
openmx_fit <- OpenMx::mxTryHard(openmx_model, extraTries = 3)
openmx_ipcr <- ipcr(
  openmx_fit,
  predictors = openmx_data["group"],
  analytic = TRUE
)
summary(openmx_ipcr)
```

### lavaan

```r
lavaan_data <- lavaan::HolzingerSwineford1939
lavaan_fit <- lavaan::cfa(
  "visual =~ x1 + x2 + x3",
  data = lavaan_data
)
lavaan_predictors <- data.frame(sex = lavaan_data$sex - 1)

lavaan_ipcr <- ipcr(lavaan_fit, predictors = lavaan_predictors)
summary(lavaan_ipcr)
```

### lm

```r
set.seed(2)
n <- 100L
group <- rep(c(0, 1), each = n / 2)
x <- stats::rnorm(n)
y <- 1 + 0.5 * x * group + stats::rnorm(n)

lm_fit <- stats::lm(y ~ x)
lm_ipcr <- ipcr(lm_fit, predictors = data.frame(group = group))
summary(lm_ipcr)
```

## Iterated IPCR

`ipcr_it()` tries each supplied learning rate in order and returns the first
converged analysis. If no attempt converges, it signals a hard error and does
not return coefficients or a partial IPCR result. The stopping rule is based
on coefficient changes relative to their standard-IPCR standard errors. A
separate learning-rate-adjusted likelihood diagnostic, `C_l / learning_rate`,
produces a warning when it exceeds 0.001 at convergence; `print()` and
`summary()` repeat that warning because this diagnostic is not part of the
stopping rule.

```r
lavaan_data <- lavaan::HolzingerSwineford1939
lavaan_fit <- lavaan::cfa(
  "visual =~ x1 + x2 + x3",
  data = lavaan_data
)
lavaan_ipcr_it <- ipcr_it(
  lavaan_fit,
  predictors = data.frame(sex = lavaan_data$sex - 1)
)
summary(lavaan_ipcr_it)
```

## Regularized IPCR

`ipcr_reg()` takes an existing standard or iterated IPCR result and applies a
cross-validated multivariate group elastic net to its stored IPCs. Predictors
are selected jointly across the first-stage model parameters. Multi-level
factors use the same dummy-variable contrasts as the original IPC regression.

```r
set.seed(3)
regularization_predictors <- data.frame(
  group = group,
  age = stats::rnorm(n),
  noise = stats::rnorm(n)
)
lm_ipcr_many <- ipcr(lm_fit, predictors = regularization_predictors)
lm_ipcr_regularized <- ipcr_reg(lm_ipcr_many)
summary(lm_ipcr_regularized)
```
