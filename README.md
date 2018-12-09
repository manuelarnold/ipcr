# ipcr

Identify and estimate parameter heterogeneity in structural equation 
models (SEM) fitted with OpenMx using individual parameter contribution (IPC)
regression. 

This package is still under development. Pleaser report any bugs.

# Installation
If you want to install ipcr from GitHub, use the following commands in R:

```{r, eval=FALSE}
if (!require(devtools)) {install.packages("devtools")}
devtools::install_github("manuelarnold/ipcr")
```

# Example
``` r
# Specify an OpenMx model
Model <- mxModel(manifestVars = c("x1", "x2", "x3", "y1", "y2", "y3"),
                latentVars = c("xi", "eta"),
                type = "RAM",
                mxData(observed = ipcr_data, type = "raw"),
                mxPath(from = "xi", to = c("x1", "x2", "x3"),
                       connect = "single", arrows = 1,
                       free = c(FALSE, TRUE, TRUE), values = 1,
                       labels = c("l11", "l22", "l33")),
                mxPath(from = "eta", to = c("y1", "y2", "y3"),
                connect = "single", arrows = 1, free = c(FALSE, TRUE, TRUE),
                values = 1, labels = c("l44", "l55", "l66")),
                mxPath(from = "xi", to = "eta", arrows = 1, free = TRUE,
                values = 0.5, labels = "beta"),
                mxPath(from = c("xi", "eta", "x1", "x2", "x3", "y1", "y2", "y3"),
                       connect = "single", arrows = 2, free = TRUE,
                       values = c(0.75, 0.555, rep(0.25, 6)),
                       labels = c("phi11", "phi22", "e11", "e22", "e33", "e44", "e55", "e66")),
                mxPath(from = "one", to = c("x1", "x2", "x3", "y1", "y2", "y3"),
                       connect = "single", arrows = 1, free = TRUE, values = 0,
                        labels = c("m_x1", "m_x2", "m_x3", "m_y1", "m_y2", "m_y3"))
)

# Fit the model
Fit <- mxTryHard(model = Model)

# Run the model
IPC_reg <- ipcr(fit = Fit, formula = ~ z1 + z2, covariates = ipcr_data, iterate = TRUE)

# Display the iterated individual parameter contribution regression estimates of the
regression coefficient beta
summary(object = IPC_reg, parameter = "beta", method = "iterated"
```

# Future plans
* Interactions between covariates
* Polynomials of covariates
* Support for missing data
