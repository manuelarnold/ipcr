# Checks for MxRAMModel ----

make_openmx_lgcm_population_model <- function(residual_variance) {
  manifest_variables <- paste0("y", 1:4)

  OpenMx::mxModel(
    "PopulationModel",
    type = "RAM",
    manifestVars = manifest_variables,
    latentVars = c("intercept", "slope"),
    OpenMx::mxPath(
      from = "intercept",
      to = manifest_variables,
      arrows = 1,
      free = FALSE,
      values = 1
    ),
    OpenMx::mxPath(
      from = "slope",
      to = manifest_variables,
      arrows = 1,
      free = FALSE,
      values = 0:3
    ),
    OpenMx::mxPath(
      from = manifest_variables,
      arrows = 2,
      free = FALSE,
      values = residual_variance
    ),
    OpenMx::mxPath(
      from = "intercept",
      arrows = 2,
      free = FALSE,
      values = 1
    ),
    OpenMx::mxPath(
      from = "slope",
      arrows = 2,
      free = FALSE,
      values = 0.1
    ),
    OpenMx::mxPath(
      from = "intercept",
      to = "slope",
      arrows = 2,
      free = FALSE,
      values = 0.2
    ),
    OpenMx::mxPath(
      from = "one",
      to = c("intercept", "slope"),
      arrows = 1,
      free = FALSE,
      values = c(1, 0.1)
    )
  )
}

make_openmx_lgcm_sample_model <- function(data) {
  manifest_variables <- paste0("y", 1:4)

  OpenMx::mxModel(
    "SampleModel",
    type = "RAM",
    manifestVars = manifest_variables,
    latentVars = c("intercept", "slope"),
    OpenMx::mxPath(
      from = "intercept",
      to = manifest_variables,
      arrows = 1,
      free = FALSE,
      values = 1
    ),
    OpenMx::mxPath(
      from = "slope",
      to = manifest_variables,
      arrows = 1,
      free = FALSE,
      values = 0:3
    ),
    OpenMx::mxPath(
      from = manifest_variables,
      arrows = 2,
      free = TRUE,
      values = 0.4,
      labels = "residual",
      lbound = 0.0001
    ),
    OpenMx::mxPath(
      from = "intercept",
      arrows = 2,
      free = TRUE,
      values = 1,
      labels = "VarI",
      lbound = 0.0001
    ),
    OpenMx::mxPath(
      from = "slope",
      arrows = 2,
      free = TRUE,
      values = 0.1,
      labels = "VarS",
      lbound = 0.0001
    ),
    OpenMx::mxPath(
      from = "intercept",
      to = "slope",
      arrows = 2,
      free = TRUE,
      values = 0.2,
      labels = "CovIS"
    ),
    OpenMx::mxPath(
      from = "one",
      to = c("intercept", "slope"),
      arrows = 1,
      free = TRUE,
      values = c(1, 0.1),
      labels = c("MeanI", "MeanS")
    ),
    OpenMx::mxData(observed = data, type = "raw")
  )
}

set.seed(23484)

n_per_group <- 100L
population_group_1 <- make_openmx_lgcm_population_model(0.3)
population_group_2 <- make_openmx_lgcm_population_model(0.5)
data_group_1 <- OpenMx::mxGenerateData(
  population_group_1,
  nrows = n_per_group
)
data_group_2 <- OpenMx::mxGenerateData(
  population_group_2,
  nrows = n_per_group
)
colnames(data_group_1) <- paste0("y", 1:4)
colnames(data_group_2) <- paste0("y", 1:4)
data <- rbind(data_group_1, data_group_2)
predictors <- data.frame(
  group = rep(c(0, 1), each = n_per_group)
)

fit <- OpenMx::mxTryHard(
  make_openmx_lgcm_sample_model(data)
)

# IPCR based on numerically computed scores
result_numeric <- ipcr(x = fit, predictors = predictors, analytic = FALSE)

# IPCR based on analytically computed scores
result_analytic <- ipcr(x = fit, predictors = predictors, analytic = TRUE)

# Compare computation time of numerically and analytically computed IPCs
run_time_numeric <- system.time(ipcr(x = fit, predictors = predictors,
                                         analytic = FALSE))
run_time_analytical <- system.time(ipcr(x = fit, predictors = predictors,
                                         analytic = TRUE))
c(numeric = run_time_numeric["elapsed"],
  analytical = run_time_analytical["elapsed"])
# Analytically computed IPCs should be much faster

fitted_parameters <- OpenMx::omxGetParameters(fit)
class(result_numeric) # expect: ipcr (numeric IPCs)
class(result_analytic) # expect: ipcr (analytic IPCs)
class(fit) # "MxRAMModel"

# check dimensions of IPC matrix
dim(result_numeric$IPCs) # dimensions IPC matrix (numeric IPCs)
dim(result_analytic$IPCs) # dimensions IPC matrix (analytic IPCs)
c(2L * n_per_group, length(fitted_parameters)) # target dimensions

# check names of IPC matrix
colnames(result_numeric$IPCs) # names of IPC matrix (numeric IPCs)
colnames(result_analytic$IPCs) # names of IPC matrix (analytic IPCs)
names(fitted_parameters) # target names

# check if the means of the IPCs are aproximately identical to the parameters
colMeans(result_numeric$IPCs) # means of the IPCs (numeric IPCs)
colMeans(result_analytic$IPCs) # means of the IPCs (analytic IPCs)
fitted_parameters # target values

# Check if the sums of the scores are almost zero
colSums(result_numeric$scores) # means of the scores (numeric scores)
colSums(result_analytic$scores) # means of the scores (analytic scores)

# check if the first six IPCs are the same
head(result_numeric$IPCs) # numeric IPCs
head(result_analytic$IPCs) # analytic IPCs
