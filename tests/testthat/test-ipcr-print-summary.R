make_print_summary_example <- function(model = c("lm", "glm")) {
  model <- match.arg(model)
  set.seed(819)
  n <- 80L
  data <- data.frame(
    x = stats::rnorm(n),
    predictor = rep(c(0, 1), each = n / 2L)
  )

  if (model == "glm") {
    probability <- stats::plogis(-0.3 + 0.6 * data$x)
    data$y <- stats::rbinom(n, size = 1, prob = probability)
    fit <- stats::glm(y ~ x, data = data, family = stats::binomial())
  } else {
    data$y <- 1 + data$x + stats::rnorm(n)
    fit <- stats::lm(y ~ x, data = data)
  }

  ipcr(fit, data["predictor"])
}

collapse_output <- function(expr) {
  paste(utils::capture.output(force(expr)), collapse = "\n")
}

test_that("print.ipcr handles models with multiple classes", {
  result <- make_print_summary_example("glm")
  returned <- NULL
  output <- utils::capture.output(
    returned <- withVisible(print(result))
  )

  expect_match(paste(output, collapse = "\n"), "Model class: glm, lm", fixed = TRUE)
  expect_false(any(grepl("Linear OpenMx Model", output, fixed = TRUE)))
  expect_identical(returned$value, result)
  expect_false(returned$visible)
})

test_that("ipcr_it methods use the default learning-rate schedule", {
  expected <- c(1, 0.5, 0.1)
  expect_identical(eval(formals(ipcr_it)$learning_rate), expected)
  methods <- c("default", "MxRAMModel", "MxModel", "lavaan")
  for (method in methods) {
    method_formals <- formals(getS3method("ipcr_it", method))
    expect_identical(eval(method_formals$learning_rate), expected)
    expect_false("convergence_rule" %in% names(method_formals))
  }
  expect_false("convergence_rule" %in% names(formals(ipcr_it)))
})

test_that("C_SE is invariant to coefficient units and learning-rate damping", {
  current <- matrix(
    0,
    nrow = 2L,
    ncol = 2L,
    dimnames = list(c("(Intercept)", "x"), c("theta_1", "theta_2"))
  )
  updated <- matrix(c(0.02, 0.04, 0.03, 0.05), nrow = 2L)
  standard_errors <- matrix(c(1, 4, 2, 5), nrow = 2L)

  reference <- ipcr:::iterated_coefficient_change_ipcr(
    current,
    updated,
    standard_errors,
    learning_rate = 1
  )

  rescaled_updated <- updated
  rescaled_updated[, 2L] <- 100 * rescaled_updated[, 2L]
  rescaled_standard_errors <- standard_errors
  rescaled_standard_errors[, 2L] <-
    100 * rescaled_standard_errors[, 2L]
  rescaled <- ipcr:::iterated_coefficient_change_ipcr(
    current,
    rescaled_updated,
    rescaled_standard_errors,
    learning_rate = 1
  )
  damped <- ipcr:::iterated_coefficient_change_ipcr(
    current,
    0.5 * updated,
    standard_errors,
    learning_rate = 0.5
  )
  expect_equal(reference, rescaled)
  expect_equal(reference, damped)
})

test_that("C_l is an absolute per-observation likelihood change", {
  expect_equal(
    ipcr:::log_likelihood_change_per_observation_ipcr(
      current = -100,
      previous = -102,
      n = 50
    ),
    0.04
  )
  expect_true(is.na(
    ipcr:::log_likelihood_change_per_observation_ipcr(
      current = -100,
      previous = NA_real_,
      n = 50
    )
  ))
})

test_that("likelihood stability is adjusted for the learning rate", {
  stable <- ipcr:::likelihood_stability_ipcr(
    C_l = 5e-4,
    learning_rate = 0.5
  )
  unstable <- ipcr:::likelihood_stability_ipcr(
    C_l = 0.5801,
    learning_rate = 0.1
  )

  expect_named(
    stable,
    c("C_l_rate_adjusted", "threshold", "stable"),
    ignore.order = FALSE
  )
  expect_equal(stable$C_l_rate_adjusted, 1e-3)
  expect_equal(stable$threshold, 1e-3)
  expect_true(stable$stable)
  expect_equal(unstable$C_l_rate_adjusted, 5.801)
  expect_false(unstable$stable)
})

test_that("a converged backend result warns when likelihood is unstable", {
  set.seed(828)
  n <- 60L
  data <- data.frame(
    x = stats::rnorm(n),
    group = rep(c(0, 1), each = n / 2L)
  )
  data$y <- 1 + 0.5 * data$x + stats::rnorm(n)
  fit <- stats::lm(y ~ x, data = data)
  predictors <- data["group"]
  standard <- ipcr(fit, predictors)
  parameter_names <- colnames(standard$IPCs)
  parameter_estimates <- colMeans(standard$IPCs)
  parameterization <- list(
    type = "original",
    natural_names = parameter_names,
    parameter_names = parameter_names,
    transformed = rep(FALSE, length(parameter_names)),
    natural_estimates = parameter_estimates,
    parameter_estimates = parameter_estimates
  )
  component_evaluator <- function(
      rows, parameter_values, learning_rate, likelihood_only = FALSE) {
    if (isTRUE(likelihood_only)) {
      return(list(log_likelihood = n * 0.02))
    }
    list(
      IPCs = standard$IPCs[rows, , drop = FALSE],
      scores = standard$scores[rows, , drop = FALSE],
      log_likelihood = 0
    )
  }

  expect_warning(
    result <- ipcr:::run_iterated_ipcr_backend(
      fit = fit,
      predictors = predictors,
      analytic = TRUE,
      conv = 1e6,
      learning_rate = 0.1,
      max_it = 1L,
      iteration_info = FALSE,
      component_evaluator = component_evaluator,
      parameterization = parameterization
    ),
    "learning-rate-adjusted likelihood change",
    fixed = TRUE
  )

  expect_equal(result$convergence$C_l, 0.02)
  expect_equal(result$convergence$C_l_rate_adjusted, 0.2)
  expect_equal(result$convergence$likelihood_warning_threshold, 1e-3)
  expect_false(result$convergence$likelihood_stable)
})

test_that("print and summary flag likelihood instability", {
  result <- make_print_summary_example("lm")
  result$info$variance_parameterization <- "original"
  result$convergence <- list(
    status = "Converged after 4 iterations with learning rate 0.1.",
    C_l_rate_adjusted = 5.801,
    likelihood_warning_threshold = 1e-3,
    likelihood_stable = FALSE
  )
  class(result) <- c("ipcr_it", "ipcr")

  print_output <- collapse_output(print(result))
  summary_output <- collapse_output(summary(result, what = "ipcr"))

  expect_match(print_output, "Likelihood stability warning", fixed = TRUE)
  expect_match(summary_output, "Likelihood stability warning", fixed = TRUE)
  expect_match(print_output, "C_l / learning_rate = 5.801", fixed = TRUE)
  expect_match(summary_output, "C_l / learning_rate = 5.801", fixed = TRUE)

  result$convergence$likelihood_stable <- TRUE
  expect_false(grepl(
    "Likelihood stability warning",
    collapse_output(print(result)),
    fixed = TRUE
  ))
  expect_false(grepl(
    "Likelihood stability warning",
    collapse_output(summary(result, what = "ipcr")),
    fixed = TRUE
  ))
})

test_that("summary.ipcr honors every output mode", {
  result <- make_print_summary_example("lm")

  everything <- collapse_output(summary(result, what = "everything"))
  marginal <- collapse_output(summary(result, what = "marginal"))
  coefficients <- collapse_output(summary(result, what = "ipcr"))

  expect_match(everything, "Marginal Linear Effects", fixed = TRUE)
  expect_match(everything, "MANCOVA", fixed = TRUE)
  expect_match(everything, "F-Tests", fixed = TRUE)
  expect_match(
    everything,
    "Individual Parameter Contribution Regression Coefficients",
    fixed = TRUE
  )

  expect_match(marginal, "Marginal Linear Effects", fixed = TRUE)
  expect_match(marginal, "MANCOVA", fixed = TRUE)
  expect_match(marginal, "F-Tests", fixed = TRUE)
  expect_false(grepl(
    "Individual Parameter Contribution Regression Coefficients",
    marginal,
    fixed = TRUE
  ))

  expect_false(grepl("Marginal Linear Effects", coefficients, fixed = TRUE))
  expect_false(grepl("MANCOVA", coefficients, fixed = TRUE))
  expect_false(grepl("F-Tests", coefficients, fixed = TRUE))
  expect_match(
    coefficients,
    "Individual Parameter Contribution Regression Coefficients",
    fixed = TRUE
  )

  returned <- NULL
  utils::capture.output(returned <- withVisible(summary(result)))
  expect_null(returned$value)
  expect_false(returned$visible)
})

test_that("summary.ipcr validates what and digits", {
  result <- make_print_summary_example("lm")

  expect_error(summary(result, what = "unknown"), "'what' must be one of")
  expect_error(summary(result, what = c("ipcr", "marginal")), "'what'")
  expect_error(summary(result, what = NA_character_), "'what'")
  expect_error(summary(result, digits = -1), "'digits'")
  expect_error(summary(result, digits = 1.5), "'digits'")
  expect_error(summary(result, digits = Inf), "'digits'")
  expect_error(summary(result, digits = c(2, 3)), "'digits'")
  expect_no_error(utils::capture.output(summary(result, digits = 0)))
})

test_that("summary.ipcr formats tables by column name", {
  result <- make_print_summary_example("lm")
  result$output$F_tests <- result$output$F_tests[, c(
    "p_Value", "DF2", "Parameter", "F_Value", "DF1", "Rsquared"
  )]
  result$output$ipcr <- result$output$ipcr[, c(
    "p_Value", "Predictor", "Std_Error", "Parameter", "t_Value", "Estimate"
  )]

  expect_no_error(utils::capture.output(summary(result, what = "everything")))
})
