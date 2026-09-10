fit_openmx_test_model <- function(model) {
  fit <- NULL
  invisible(utils::capture.output(
    fit <- OpenMx::mxTryHard(
      model,
      extraTries = 3,
      bestInitsOutput = FALSE,
      silent = TRUE,
      verbose = 0
    ),
    type = "output"
  ))
  fit
}

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

test_that("ipcr fits an OpenMx latent growth curve model", {
  skip_if_not_installed("OpenMx")
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

  fit <- fit_openmx_test_model(make_openmx_lgcm_sample_model(data))
  result <- suppressMessages(
    ipcr(x = fit, predictors = predictors)
  )

  fitted_parameters <- OpenMx::omxGetParameters(fit)
  expect_s3_class(result, "ipcr")
  expect_s4_class(fit, "MxRAMModel")
  expect_identical(
    dim(result$IPCs),
    c(2L * n_per_group, length(fitted_parameters))
  )
  expect_identical(colnames(result$IPCs), names(fitted_parameters))
  expect_equal(
    colMeans(result$IPCs),
    fitted_parameters,
    tolerance = 1e-4
  )
  expect_identical(stats::nobs(fit), 2L * n_per_group)
})

make_openmx_nonram_normal_model <- function(data) {
  OpenMx::mxModel(
    "NonRAMNormalModel",
    OpenMx::mxMatrix(
      type = "Symm",
      nrow = 1,
      ncol = 1,
      free = TRUE,
      values = 1,
      labels = "variance",
      lbound = 0.0001,
      name = "expected_covariance"
    ),
    OpenMx::mxMatrix(
      type = "Full",
      nrow = 1,
      ncol = 1,
      free = TRUE,
      values = 0,
      labels = "mean",
      name = "expected_means"
    ),
    OpenMx::mxData(observed = data, type = "raw"),
    OpenMx::mxExpectationNormal(
      covariance = "expected_covariance",
      means = "expected_means",
      dimnames = "y"
    ),
    OpenMx::mxFitFunctionML()
  )
}

test_that("a non-RAM MxModel uses numerical row scores", {
  skip_if_not_installed("OpenMx")
  set.seed(817)
  n <- 80L
  data <- data.frame(y = stats::rnorm(n))
  fit <- fit_openmx_test_model(make_openmx_nonram_normal_model(data))
  predictors <- data.frame(group = rep(c(0, 1), each = n / 2L))

  result <- ipcr(fit, predictors)
  fitted_parameters <- OpenMx::omxGetParameters(fit)

  expect_s4_class(fit, "MxModel")
  expect_false(inherits(fit, "MxRAMModel"))
  expect_identical(dim(result$IPCs), c(n, length(fitted_parameters)))
  expect_identical(stats::nobs(fit), n)
  expect_equal(
    colMeans(result$IPCs),
    fitted_parameters,
    tolerance = 1e-5
  )
  expect_warning(
    analytical_request <- ipcr(fit, predictors, analytic = TRUE),
    "using numerical OpenMx row scores"
  )
  expect_equal(analytical_request$IPCs, result$IPCs, tolerance = 1e-10)
  expect_output(print(result), "Model class")
})

test_that("analytical and numerical OpenMx RAM scores agree", {
  skip_if_not_installed("OpenMx")
  set.seed(818)
  population <- make_openmx_lgcm_population_model(0.4)
  data <- OpenMx::mxGenerateData(population, nrows = 80L)
  colnames(data) <- paste0("y", 1:4)
  fit <- fit_openmx_test_model(make_openmx_lgcm_sample_model(data))
  predictors <- data.frame(group = rep(c(0, 1), each = 40L))

  numerical_scores <- get_scores(fit, analytic = FALSE)
  analytical_scores <- get_scores(fit, analytic = TRUE)
  result <- ipcr(fit, predictors, analytic = TRUE)

  expect_equal(analytical_scores, numerical_scores, tolerance = 1e-7)
  expect_lt(max(abs(colSums(analytical_scores))), 5e-4)
  expect_identical(dim(result$IPCs), c(80L, 6L))
  expect_true(all(is.finite(result$IPCs)))
  expect_equal(result$scores, as.matrix(analytical_scores), tolerance = 1e-12)
})

test_that("analytical OpenMx RAM scores handle missing observations", {
  skip_if_not_installed("OpenMx")
  set.seed(819)
  population <- make_openmx_lgcm_population_model(0.4)
  data <- OpenMx::mxGenerateData(population, nrows = 80L)
  colnames(data) <- paste0("y", 1:4)
  data[seq(2L, 80L, by = 5L), "y1"] <- NA_real_
  data[seq(3L, 80L, by = 7L), "y3"] <- NA_real_
  fit <- fit_openmx_test_model(make_openmx_lgcm_sample_model(data))

  numerical_scores <- get_scores(fit, analytic = FALSE)
  analytical_scores <- get_scores(fit, analytic = TRUE)
  result <- ipcr(
    fit,
    predictors = data.frame(group = rep(c(0, 1), each = 40L)),
    analytic = TRUE
  )

  expect_identical(dim(analytical_scores), c(80L, 6L))
  expect_equal(analytical_scores, numerical_scores, tolerance = 1e-7)
  expect_lt(max(abs(colSums(analytical_scores))), 5e-4)
  expect_identical(dim(result$IPCs), c(80L, 6L))
  expect_equal(
    colMeans(result$IPCs),
    OpenMx::omxGetParameters(fit),
    tolerance = 1e-4
  )
})

make_openmx_definition_model <- function(data) {
  OpenMx::mxModel(
    "DefinitionModel",
    type = "RAM",
    manifestVars = c("y1", "y2"),
    OpenMx::mxPath(
      from = "y1",
      arrows = 2,
      free = TRUE,
      values = 0.7,
      labels = "var1"
    ),
    OpenMx::mxPath(
      from = "y2",
      arrows = 2,
      free = FALSE,
      labels = "data.resvar"
    ),
    OpenMx::mxPath(
      from = "y1",
      to = "y2",
      arrows = 1,
      free = FALSE,
      labels = "data.beta"
    ),
    OpenMx::mxPath(
      from = "one",
      to = "y1",
      arrows = 1,
      free = FALSE,
      labels = "data.mu"
    ),
    OpenMx::mxPath(
      from = "one",
      to = "y2",
      arrows = 1,
      free = TRUE,
      values = 0.2,
      labels = "mean2"
    ),
    OpenMx::mxData(observed = data, type = "raw")
  )
}

test_that("analytical scores handle definition variables and missingness", {
  skip_if_not_installed("OpenMx")
  set.seed(820)
  n <- 60L
  definition_mean <- rep(c(-0.4, 0.2), each = n / 2L)
  definition_path <- rep(c(0.2, 0.6), length.out = n)
  definition_variance <- rep(c(0.7, 1.1), length.out = n)
  y1 <- definition_mean + stats::rnorm(n, sd = 0.8)
  y2 <- 0.3 + definition_path * y1 +
    stats::rnorm(n, sd = sqrt(definition_variance))
  data <- data.frame(
    y1 = y1,
    y2 = y2,
    mu = definition_mean,
    beta = definition_path,
    resvar = definition_variance
  )
  data[seq(2L, n, by = 6L), "y1"] <- NA_real_
  data[seq(4L, n, by = 7L), "y2"] <- NA_real_
  fit <- fit_openmx_test_model(make_openmx_definition_model(data))

  numerical_scores <- get_scores(fit, analytic = FALSE)
  analytical_scores <- get_scores(fit, analytic = TRUE)

  expect_identical(dim(analytical_scores), c(n, 2L))
  expect_equal(analytical_scores, numerical_scores, tolerance = 1e-7)
  expect_lt(max(abs(colSums(analytical_scores))), 5e-4)
})

test_that("analytical scores handle unlabelled free RAM parameters", {
  skip_if_not_installed("OpenMx")
  set.seed(821)
  data <- data.frame(y = stats::rnorm(50L))
  model <- OpenMx::mxModel(
    "Unlabelled",
    type = "RAM",
    manifestVars = "y",
    OpenMx::mxPath(from = "y", arrows = 2, free = TRUE, values = 1),
    OpenMx::mxPath(from = "one", to = "y", arrows = 1,
                   free = TRUE, values = 0),
    OpenMx::mxData(observed = data, type = "raw")
  )
  fit <- fit_openmx_test_model(model)

  numerical_scores <- get_scores(fit, analytic = FALSE)
  analytical_scores <- get_scores(fit, analytic = TRUE)

  expect_identical(colnames(analytical_scores), names(fit$output$estimate))
  expect_equal(analytical_scores, numerical_scores, tolerance = 1e-7)
})

test_that("explicit OpenMx constraints use numerical scores", {
  skip_if_not_installed("OpenMx")
  data <- data.frame(y = stats::rnorm(20L))
  model <- OpenMx::mxModel(
    "Constrained",
    type = "RAM",
    manifestVars = "y",
    OpenMx::mxPath(from = "y", arrows = 2, free = TRUE,
                   values = 1, labels = "variance"),
    OpenMx::mxPath(from = "one", to = "y", arrows = 1,
                   free = TRUE, values = 0, labels = "mean"),
    OpenMx::mxConstraint(M[1, 1] == 0, name = "zero_mean"),
    OpenMx::mxData(observed = data, type = "raw")
  )

  expect_warning(
    analytic <- ipcr:::normalize_analytic_ipcr(model, TRUE),
    "mxConstraint"
  )
  expect_false(analytic)
})

test_that("OpenMx algebras use numerical scores", {
  skip_if_not_installed("OpenMx")
  data <- data.frame(y = stats::rnorm(20L))
  model <- OpenMx::mxModel(
    "WithAlgebra",
    type = "RAM",
    manifestVars = "y",
    OpenMx::mxPath(from = "y", arrows = 2, free = TRUE,
                   values = 1, labels = "variance"),
    OpenMx::mxPath(from = "one", to = "y", arrows = 1,
                   free = TRUE, values = 0, labels = "mean"),
    OpenMx::mxAlgebra(A + A, name = "unused_algebra"),
    OpenMx::mxData(observed = data, type = "raw")
  )

  expect_warning(
    analytic <- ipcr:::normalize_analytic_ipcr(model, TRUE),
    "algebras"
  )
  expect_false(analytic)
})

make_openmx_iterated_location_scale_model <- function(data) {
  OpenMx::mxModel(
    "IteratedLocationScale",
    type = "RAM",
    manifestVars = "y",
    OpenMx::mxPath(
      from = "y",
      arrows = 2,
      free = TRUE,
      values = 1,
      labels = "variance",
      lbound = 0.0001
    ),
    OpenMx::mxPath(
      from = "one",
      to = "y",
      arrows = 1,
      free = TRUE,
      values = 0,
      labels = "mean"
    ),
    OpenMx::mxData(observed = data, type = "raw")
  )
}

make_openmx_fixed_shared_variance_model <- function(data) {
  manifest_variables <- c("y1", "y2", "y3")
  OpenMx::mxModel(
    "FixedSharedVariance",
    type = "RAM",
    manifestVars = manifest_variables,
    OpenMx::mxPath(
      from = "y1", arrows = 2, free = FALSE, values = 1
    ),
    OpenMx::mxPath(
      from = c("y2", "y3"), arrows = 2, free = TRUE,
      values = 0.8, labels = "shared_variance", lbound = 0.0001
    ),
    OpenMx::mxPath(
      from = "one", to = manifest_variables, arrows = 1,
      free = TRUE, values = 0,
      labels = paste0("mean", seq_along(manifest_variables))
    ),
    OpenMx::mxData(observed = data, type = "raw")
  )
}

test_that("cached analytical RAM components preserve iterated calculations", {
  skip_if_not_installed("OpenMx")
  set.seed(900)
  data <- data.frame(y = stats::rnorm(30L))
  fit <- fit_openmx_test_model(
    make_openmx_iterated_location_scale_model(data)
  )
  parameter_estimates <- coef_ipcr(fit)
  parameter_values <- matrix(
    parameter_estimates,
    nrow = nrow(data),
    ncol = length(parameter_estimates),
    byrow = TRUE
  )
  rows <- seq_len(nrow(data))
  context <- analytical_iterated_ram_context_ipcr(fit)

  uncached <- analytical_iterated_ram_components_ipcr(
    x = fit,
    rows = rows,
    parameter_values = parameter_values,
    learning_rate = 0.5
  )
  cached <- analytical_iterated_ram_components_ipcr(
    x = fit,
    rows = rows,
    parameter_values = parameter_values,
    learning_rate = 0.5,
    context = context
  )
  likelihood <- analytical_iterated_ram_components_ipcr(
    x = fit,
    rows = rows,
    parameter_values = parameter_values,
    learning_rate = 0.5,
    context = context,
    likelihood_only = TRUE
  )

  expect_equal(cached, uncached, tolerance = 1e-12)
  expect_identical(names(likelihood), "log_likelihood")
  expect_equal(
    likelihood$log_likelihood,
    cached$log_likelihood,
    tolerance = 1e-12
  )
})

test_that("analytical RAM components use the log-variance chain rule", {
  skip_if_not_installed("OpenMx")
  set.seed(909)
  data <- data.frame(y = stats::rnorm(40L, sd = 0.8))
  fit <- fit_openmx_test_model(
    make_openmx_iterated_location_scale_model(data)
  )
  rows <- seq_len(nrow(data))
  natural_estimates <- coef_ipcr(fit)
  natural_values <- matrix(
    natural_estimates,
    nrow = nrow(data),
    ncol = length(natural_estimates),
    byrow = TRUE
  )
  original <- analytical_iterated_ram_components_ipcr(
    x = fit,
    rows = rows,
    parameter_values = natural_values,
    learning_rate = 0.5
  )

  spec <- ram_variance_parameterization_ipcr(fit, "log")
  log_values <- matrix(
    spec$parameter_estimates,
    nrow = nrow(data),
    ncol = length(spec$parameter_estimates),
    byrow = TRUE
  )
  transformed <- analytical_iterated_ram_components_ipcr(
    x = fit,
    rows = rows,
    parameter_values = log_values,
    learning_rate = 0.5,
    context = analytical_iterated_ram_context_ipcr(
      fit,
      parameterization = spec
    )
  )

  variance <- unname(natural_estimates["variance"])
  expect_equal(
    transformed$scores[, "log(variance)"],
    original$scores[, "variance"] * variance,
    tolerance = 1e-10
  )
  expect_equal(
    transformed$scores[, "mean"],
    original$scores[, "mean"],
    tolerance = 1e-10
  )
  expect_equal(
    transformed$IPCs[, "log(variance)"] - log(variance),
    (original$IPCs[, "variance"] - variance) / variance,
    tolerance = 1e-10
  )
  expect_equal(
    transformed$log_likelihood,
    original$log_likelihood,
    tolerance = 1e-10
  )
})

test_that("ipcr_it converges for an analytical MxRAMModel", {
  skip_if_not_installed("OpenMx")
  set.seed(901)
  n <- 100L
  group <- rep(c(0, 1), each = n / 2L)
  y <- 0.7 * group + stats::rnorm(
    n,
    sd = ifelse(group == 0, 0.8, 1)
  )
  fit <- fit_openmx_test_model(
    make_openmx_iterated_location_scale_model(data.frame(y = y))
  )
  named_ram_fit <- fit

  result <- ipcr_it(
    named_ram_fit,
    predictors = data.frame(group = group),
    conv = 1e-7,
    learning_rate = c(1, 0.5),
    max_it = 50,
    iteration_info = TRUE
  )
  group_means <- tapply(y, group, mean)
  group_variances <- tapply(
    y,
    group,
    function(values) mean((values - mean(values))^2)
  )

  expect_s3_class(result, "ipcr_it")
  expect_s3_class(result, "ipcr")
  expect_identical(result$info$name, "named_ram_fit")
  expect_true(result$convergence$converged)
  expect_lte(result$convergence$C_SE, 1e-7)
  expect_true(is.finite(result$convergence$C_l))
  expect_true(all(c("C_SE", "C_l") %in% names(result$iteration_history)))
  expect_equal(
    result$convergence$C_SE,
    utils::tail(result$iteration_history$C_SE, 1L)
  )
  expect_equal(
    result$convergence$C_l,
    utils::tail(result$iteration_history$C_l, 1L)
  )
  expect_equal(
    result$convergence$log_likelihood,
    utils::tail(result$iteration_history$log_likelihood, 1L)
  )
  if (NROW(result$iteration_history) > 1L) {
    expect_equal(
      result$iteration_history$C_l[-1L],
      abs(diff(result$iteration_history$log_likelihood)) /
        sum(result$complete_rows)
    )
  }
  expect_identical(result$convergence$attempts$learning_rate, 1)
  expect_equal(
    stats::coef(result)["group", "mean"],
    unname(diff(group_means)),
    tolerance = 1e-6
  )
  expect_equal(
    stats::coef(result)["group", "variance"],
    unname(diff(group_variances)),
    tolerance = 1e-6
  )
  expect_identical(dim(result$IPCs), c(n, 2L))
  expect_identical(dim(result$scores), c(n, 2L))
  expect_identical(tail(result$iteration_history$status, 1L), "converged")
  expect_output(print(result), "Iterated IPCR status: Converged")
})

test_that("ipcr_it reports log and natural-scale variance effects", {
  skip_if_not_installed("OpenMx")
  set.seed(910)
  n <- 120L
  group <- rep(c(0, 1), each = n / 2L)
  y <- 0.4 * group + stats::rnorm(
    n,
    sd = ifelse(group == 0, sqrt(0.6), 1)
  )
  fit <- fit_openmx_test_model(
    make_openmx_iterated_location_scale_model(data.frame(y = y))
  )

  result <- ipcr_it(
    fit,
    predictors = data.frame(group = group),
    conv = 1e-7,
    learning_rate = 1,
    max_it = 50,
    iteration_info = TRUE,
    variance_parameterization = "log"
  )
  group_variances <- tapply(
    y,
    group,
    function(values) mean((values - mean(values))^2)
  )
  expected_log_effect <- unname(diff(log(group_variances)))
  expected_variance_change <- unname(diff(group_variances))
  effects <- result$output$natural_scale_effects

  expect_true(result$convergence$converged)
  expect_identical(
    result$info$parameters,
    c("log(variance)", "mean")
  )
  expect_identical(
    result$parameterization$log_variance_parameters,
    "variance"
  )
  expect_equal(
    stats::coef(result)["group", "log(variance)"],
    expected_log_effect,
    tolerance = 1e-6
  )
  expect_identical(colnames(result$natural_scale_IPCs), c("variance", "mean"))
  expect_true(all(result$natural_scale_IPCs[, "variance"] > 0))
  expect_equal(effects$Variance_Ratio, exp(expected_log_effect), tolerance = 1e-6)
  expect_equal(
    effects$Reference_Variance,
    unname(group_variances[1L]),
    tolerance = 1e-6
  )
  expect_equal(
    effects$Variance_Change,
    expected_variance_change,
    tolerance = 1e-6
  )
  expect_output(
    summary(result, what = "ipcr"),
    "Back-Transformed Variance Effects"
  )
})

test_that("log variance handles fixed and shared S-diagonal elements", {
  skip_if_not_installed("OpenMx")
  set.seed(911)
  n <- 80L
  data <- data.frame(
    y1 = stats::rnorm(n, sd = 1),
    y2 = stats::rnorm(n, sd = sqrt(0.7)),
    y3 = stats::rnorm(n, sd = sqrt(0.7))
  )
  fit <- fit_openmx_test_model(
    make_openmx_fixed_shared_variance_model(data)
  )
  result <- ipcr_it(
    fit,
    predictors = data.frame(group = rep(c(0, 1), each = n / 2L)),
    conv = 1e6,
    learning_rate = 1,
    max_it = 1,
    variance_parameterization = "log"
  )

  expect_true(result$convergence$converged)
  expect_lte(result$convergence$C_SE, 1e6)
  expect_identical(
    result$parameterization$log_variance_parameters,
    "shared_variance"
  )
  expect_true("log(shared_variance)" %in% colnames(result$IPCs))
  expect_false(any(grepl("y1", colnames(result$IPCs), fixed = TRUE)))
  expect_true(all(result$natural_scale_IPCs[, "shared_variance"] > 0))
})

test_that("log variance is limited to analytical MxRAMModel iteration", {
  skip_if_not_installed("OpenMx")
  set.seed(912)
  data <- data.frame(y = stats::rnorm(20L))
  ram_fit <- fit_openmx_test_model(
    make_openmx_iterated_location_scale_model(data)
  )
  nonram_fit <- fit_openmx_test_model(
    make_openmx_nonram_normal_model(data)
  )

  expect_error(
    ipcr_it(
      ram_fit,
      predictors = data.frame(group = rep(c(0, 1), each = 10L)),
      analytic = FALSE,
      variance_parameterization = "log"
    ),
    "only with analytic = TRUE"
  )
  expect_error(
    ipcr_it(
      nonram_fit,
      predictors = data.frame(group = rep(c(0, 1), each = 10L)),
      variance_parameterization = "log"
    ),
    "available only for analytical OpenMx MxRAMModel"
  )
})

test_that("ipcr_it updates shared-label MxRAMModel parameters analytically", {
  skip_if_not_installed("OpenMx")
  set.seed(904)
  n_per_group <- 40L
  data_group_1 <- OpenMx::mxGenerateData(
    make_openmx_lgcm_population_model(0.35),
    nrows = n_per_group
  )
  data_group_2 <- OpenMx::mxGenerateData(
    make_openmx_lgcm_population_model(0.45),
    nrows = n_per_group
  )
  data <- rbind(data_group_1, data_group_2)
  colnames(data) <- paste0("y", 1:4)
  fit <- fit_openmx_test_model(make_openmx_lgcm_sample_model(data))

  result <- ipcr_it(
    fit,
    predictors = data.frame(
      group = rep(c(0, 1), each = n_per_group)
    ),
    conv = 1e6,
    learning_rate = c(1, 0.5),
    max_it = 2
  )
  log_result <- ipcr_it(
    fit,
    predictors = data.frame(
      group = rep(c(0, 1), each = n_per_group)
    ),
    conv = 1e-5,
    learning_rate = c(1, 0.5),
    max_it = 50,
    variance_parameterization = "log"
  )

  expect_true(result$convergence$converged)
  expect_identical(result$convergence$attempts$learning_rate, 1)
  expect_identical(colnames(result$IPCs), names(coef_ipcr(fit)))
  expect_true(all(is.finite(result$IPCs)))
  expect_true("residual" %in% colnames(result$IPCs))
  expect_true(log_result$convergence$converged)
  expect_identical(
    log_result$parameterization$log_variance_parameters,
    c("residual", "VarI", "VarS")
  )
  expect_true(all(c(
    "log(residual)", "log(VarI)", "log(VarS)"
  ) %in% colnames(log_result$IPCs)))
})

test_that("ipcr_it handles partially missing OpenMx raw observations", {
  skip_if_not_installed("OpenMx")
  set.seed(905)
  n <- 80L
  data <- OpenMx::mxGenerateData(
    make_openmx_lgcm_population_model(0.4),
    nrows = n
  )
  colnames(data) <- paste0("y", 1:4)
  data[seq_len(8L), "y1"] <- NA_real_
  data[9:16, "y3"] <- NA_real_
  fit <- fit_openmx_test_model(make_openmx_lgcm_sample_model(data))

  expect_warning(
    result <- ipcr_it(
      fit,
      predictors = data.frame(group = rep(c(0, 1), each = n / 2L)),
      conv = 1e6,
      learning_rate = c(1, 0.5),
      max_it = 2
    ),
    "learning-rate-adjusted likelihood change",
    fixed = TRUE
  )

  expect_true(result$convergence$converged)
  expect_false(result$convergence$likelihood_stable)
  expect_identical(dim(result$IPCs), c(n, 6L))
  expect_true(all(is.finite(result$IPCs)))
})

test_that("ipcr_it retains analytical RAM definition-variable handling", {
  skip_if_not_installed("OpenMx")
  set.seed(906)
  n <- 60L
  definition_mean <- rep(c(-0.3, 0.3), each = n / 2L)
  definition_path <- rep(c(0.2, 0.5), length.out = n)
  definition_variance <- rep(c(0.8, 1.1), length.out = n)
  y1 <- definition_mean + stats::rnorm(n, sd = 0.8)
  y2 <- 0.2 + definition_path * y1 + stats::rnorm(
    n,
    sd = sqrt(definition_variance)
  )
  data <- data.frame(
    y1 = y1,
    y2 = y2,
    mu = definition_mean,
    beta = definition_path,
    resvar = definition_variance
  )
  fit <- fit_openmx_test_model(make_openmx_definition_model(data))

  result <- ipcr_it(
    fit,
    predictors = data.frame(group = rep(c(0, 1), each = n / 2L)),
    conv = 1e6,
    learning_rate = c(1, 0.5),
    max_it = 2
  )

  expect_true(result$convergence$converged)
  expect_identical(dim(result$IPCs), c(n, 2L))
  expect_true(all(is.finite(result$scores)))
})

test_that("ipcr_it tries every rate before signaling nonconvergence", {
  skip_if_not_installed("OpenMx")
  set.seed(902)
  n <- 80L
  predictor <- seq(-1, 1, length.out = n)
  y <- 0.4 + 0.8 * predictor + stats::rnorm(
    n,
    sd = exp(0.2 + 0.2 * predictor)
  )
  fit <- fit_openmx_test_model(
    make_openmx_iterated_location_scale_model(data.frame(y = y))
  )

  nonconvergence <- tryCatch(
    ipcr_it(
      fit,
      predictors = data.frame(predictor = predictor),
      conv = 1e-15,
      learning_rate = c(1, 0.5),
      max_it = 1,
      iteration_info = TRUE
    ),
    ipcr_it_nonconvergence = function(error) error
  )

  expect_s3_class(nonconvergence, "ipcr_it_nonconvergence")
  expect_s3_class(nonconvergence, "error")
  expect_match(
    conditionMessage(nonconvergence),
    "did not converge for any supplied learning rate",
    fixed = TRUE
  )
  expect_s3_class(nonconvergence$attempts, "data.frame")
  expect_identical(
    names(nonconvergence$attempts),
    c("learning_rate", "iterations", "C_SE", "C_l", "status", "reason")
  )
  expect_identical(
    nonconvergence$attempts$learning_rate,
    c(1, 0.5)
  )
  expect_true(all(
    nonconvergence$attempts$status == "maximum iterations reached"
  ))
  expect_true(all(is.finite(nonconvergence$attempts$C_SE)))
  expect_true(all(is.finite(nonconvergence$attempts$C_l)))
  expect_s3_class(nonconvergence$iteration_history, "data.frame")
  expect_identical(
    names(nonconvergence$iteration_history),
    c(
      "learning_rate", "iteration", "C_SE", "C_l", "log_likelihood",
      "status"
    )
  )
  expect_identical(NROW(nonconvergence$iteration_history), 2L)
  expect_false(inherits(nonconvergence, "ipcr_it"))
  expect_false(any(
    c(
      "IPCs", "scores", "mlm", "output", "mancova", "coefficients",
      "tests", "result"
    ) %in% names(nonconvergence)
  ))

  default_nonconvergence <- tryCatch(
    ipcr_it(
      fit,
      predictors = data.frame(predictor = predictor),
      conv = 1e-15,
      max_it = 1
    ),
    ipcr_it_nonconvergence = function(error) error
  )

  expect_s3_class(default_nonconvergence, "ipcr_it_nonconvergence")
  expect_identical(
    default_nonconvergence$attempts$learning_rate,
    c(1, 0.5, 0.1)
  )
  expect_null(default_nonconvergence$iteration_history)

  malformed_backend <- tryCatch(
    ipcr:::run_iterated_ipcr_backend(
      fit = fit,
      predictors = data.frame(predictor = predictor),
      analytic = TRUE,
      conv = 0.01,
      learning_rate = c(0.8, 0.4),
      max_it = 1,
      iteration_info = FALSE,
      component_evaluator = function(...) {
        list(log_likelihood = NA_real_)
      },
      parameterization = ipcr:::ram_variance_parameterization_ipcr(
        fit,
        "original"
      )
    ),
    ipcr_it_nonconvergence = function(error) error
  )
  expect_s3_class(malformed_backend, "ipcr_it_nonconvergence")
  expect_identical(
    malformed_backend$attempts$learning_rate,
    c(0.8, 0.4)
  )
  expect_true(all(malformed_backend$attempts$status == "numerical failure"))
  expect_true(all(grepl(
    "invalid log-likelihood",
    malformed_backend$attempts$reason,
    fixed = TRUE
  )))
})

test_that("numerical iterated IPCR has a minimal MxRAMModel smoke test", {
  skip_if_not_installed("OpenMx")
  set.seed(907)
  n <- 12L
  group <- rep(c(0, 1), each = n / 2L)
  data <- data.frame(y = 0.4 * group + stats::rnorm(n))
  fit <- fit_openmx_test_model(
    make_openmx_iterated_location_scale_model(data)
  )

  expect_warning(
    result <- ipcr_it(
      fit,
      predictors = data.frame(group = group),
      analytic = FALSE,
      conv = 1e6,
      learning_rate = 1,
      max_it = 1
    ),
    "very long run time"
  )

  expect_s3_class(result, "ipcr_it")
  expect_true(result$convergence$converged)
  expect_false(result$info$analytic)
  expect_identical(dim(result$IPCs), c(n, 2L))
})

test_that("a general MxModel has a minimal numerical iteration smoke test", {
  skip_if_not_installed("OpenMx")
  set.seed(908)
  n <- 12L
  group <- rep(c(0, 1), each = n / 2L)
  data <- data.frame(y = 0.4 * group + stats::rnorm(n))
  fit <- fit_openmx_test_model(make_openmx_nonram_normal_model(data))
  named_general_fit <- fit

  expect_warning(
    result <- ipcr_it(
      named_general_fit,
      predictors = data.frame(group = group),
      analytic = TRUE,
      conv = 1e6,
      learning_rate = 1,
      max_it = 1,
      iteration_info = TRUE
    ),
    "Analytical scores are unavailable.*very long run time"
  )

  expect_s3_class(result, "ipcr_it")
  expect_identical(result$info$name, "named_general_fit")
  expect_true(result$convergence$converged)
  expect_false(result$info$analytic)
  expect_identical(dim(result$IPCs), c(n, 2L))
  expect_identical(dim(result$scores), c(n, 2L))
  expect_true(is.finite(result$convergence$log_likelihood))
})

# Numerical iteration tests intentionally use only 12 observations, one binary
# predictor, one learning rate, and one iteration. Do not add continuous IPC
# predictors to automatic numerical OpenMx tests: they can create a distinct
# predicted parameter vector for every observation and make the checks take a
# very long time. Use
# tests/manual/Checks/check-iterated-ipcr-backends.R instead.
test_that("ipcr_it rejects unsupported model classes", {
  expect_error(
    ipcr_it(stats::lm(mpg ~ wt, data = mtcars), predictors = mtcars$cyl),
    "currently implemented only for lavaan and OpenMx MxModel"
  )
})
