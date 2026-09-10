make_core_lm_example <- function(n = 40L) {
  set.seed(812)
  data <- data.frame(
    x = stats::rnorm(n),
    predictor = rep(c(0, 1), length.out = n)
  )
  data$y <- 1 + data$x + stats::rnorm(n)

  list(
    fit = stats::lm(y ~ x, data = data),
    predictors = data["predictor"]
  )
}

test_that("ipcr accessors delegate to the second-stage regression", {
  example <- make_core_lm_example()
  result <- ipcr(example$fit, example$predictors)
  new_predictors <- data.frame(predictor = c(0, 1))

  expect_equal(stats::coef(result), stats::coef(result$mlm))
  expect_equal(stats::coefficients(result), stats::coefficients(result$mlm))
  expect_equal(stats::fitted(result), stats::fitted(result$mlm))
  expect_identical(stats::nobs(result), stats::nobs(result$mlm))
  expect_equal(stats::predict(result), stats::predict(result$mlm))
  expect_equal(
    stats::predict(result, newdata = new_predictors),
    stats::predict(result$mlm, newdata = new_predictors)
  )
  expect_equal(stats::residuals(result), stats::residuals(result$mlm))
  expect_equal(stats::sigma(result), stats::sigma(result$mlm))
  expect_equal(stats::confint(result), stats::confint(result$mlm))
  expect_equal(stats::effects(result), stats::effects(result$mlm))
  expect_equal(stats::vcov(result), stats::vcov(result$mlm))

  expect_true(all(is.finite(result$output$F_tests$Rsquared)))
  expect_true(all(is.finite(result$output$F_tests$F_Value)))
  expect_true(all(is.finite(result$output$F_tests$p_Value)))
  expect_s3_class(result$mancova, "Anova.mlm")
  expect_identical(result$mancova$type, "II")
  expect_identical(result$mancova$test, "Pillai")
  expect_output(summary(result), "Type II MANOVA Tests", fixed = TRUE)
})

test_that("ipcr supports a model with one fitted parameter", {
  set.seed(813)
  data <- data.frame(
    x = seq_len(30L),
    predictor = rep(c(0, 1), each = 15L)
  )
  data$y <- 2 * data$x + stats::rnorm(nrow(data))
  fit <- stats::nls(y ~ slope * x, data = data, start = list(slope = 1))

  expect_warning(
    result <- ipcr(x = fit, predictors = data["predictor"]),
    "not officially supported"
  )

  expect_s3_class(result, "ipcr")
  expect_identical(dim(result$IPCs), c(nrow(data), 1L))
  expect_identical(colnames(result$IPCs), "slope")
  expect_null(result$mancova)
  expect_null(result$output$mancova)
  expect_identical(nrow(result$output$ipcr), 2L)
  expect_identical(nrow(result$output$F_tests), 1L)
  marginal_output <- paste(
    utils::capture.output(summary(result, what = "marginal")),
    collapse = "\n"
  )
  coefficient_output <- paste(
    utils::capture.output(summary(result, what = "ipcr")),
    collapse = "\n"
  )
  expect_match(marginal_output, "Not available for this result", fixed = TRUE)
  expect_false(grepl(
    "Individual Parameter Contribution Regression Coefficients",
    marginal_output,
    fixed = TRUE
  ))
  expect_match(
    coefficient_output,
    "Individual Parameter Contribution Regression Coefficients",
    fixed = TRUE
  )
})

test_that("experimental model subclasses warn but are attempted", {
  example <- make_core_lm_example()
  class(example$fit) <- c("experimental_lm", class(example$fit))

  expect_warning(
    result <- ipcr(example$fit, example$predictors),
    "not officially supported"
  )
  expect_s3_class(result, "ipcr")
})

test_that("the three public IPCR functions consistently use x", {
  expect_identical(names(formals(ipcr))[1L], "x")
  expect_identical(names(formals(ipcr_it))[1L], "x")
  expect_identical(names(formals(ipcr_reg))[1L], "x")
})

test_that("rank-deficient MANOVA does not prevent an ipcr result", {
  data <- data.frame(
    x = c(-1, 0, 1, 2),
    y = c(0, 1, 1, 4),
    predictor = c(0, 0, 1, 1)
  )
  fit <- stats::lm(y ~ x, data = data)

  expect_warning(
    result <- ipcr(x = fit, predictors = data["predictor"]),
    "residual IPC matrix has rank"
  )

  expect_s3_class(result, "ipcr")
  expect_s3_class(result$mancova, "Anova.mlm")
  expect_null(result$output$mancova)
  expect_output(summary(result), "Not available for this result", fixed = TRUE)
})

test_that("unavailable overall F-statistics are reported clearly", {
  example <- make_core_lm_example()
  warnings <- character()

  result <- withCallingHandlers(
    ipcr(
      example$fit,
      data.frame(constant_predictor = rep(1, nrow(example$predictors)))
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_true(any(grepl("overall F-statistic is unavailable", warnings)))
  expect_true(all(is.na(result$output$F_tests$F_Value)))
  expect_true(all(is.na(result$output$F_tests$p_Value)))
  expect_no_error(utils::capture.output(
    summary(result, what = "marginal")
  ))
})

test_that("IPC component validation rejects incompatible inputs", {
  parameters <- c(alpha = 1, beta = 2)
  valid_scores <- matrix(
    c(1, 0, 0, 1, -1, -1),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(NULL, names(parameters))
  )
  valid_bread <- diag(2)
  dimnames(valid_bread) <- list(names(parameters), names(parameters))

  expect_error(
    ipcr:::validate_ipcr_components(
      parameters,
      valid_scores[, 1, drop = FALSE],
      valid_bread,
      n = 3L
    ),
    "score matrix must have 2 columns"
  )
  expect_error(
    ipcr:::validate_ipcr_components(
      parameters,
      valid_scores,
      matrix(1, nrow = 2, ncol = 1),
      n = 3L
    ),
    "bread matrix must be a 2 by 2 matrix"
  )
  expect_error(
    ipcr:::validate_ipcr_components(
      parameters,
      valid_scores,
      matrix(c(1, 0, 0, 0), nrow = 2),
      n = 3L
    ),
    "bread matrix is singular"
  )

  invalid_names <- valid_scores
  colnames(invalid_names) <- rev(names(parameters))
  expect_error(
    ipcr:::validate_ipcr_components(
      parameters,
      invalid_names,
      valid_bread,
      n = 3L
    ),
    "Score column names do not match"
  )

  nonfinite_scores <- valid_scores
  nonfinite_scores[1, 1] <- NA_real_
  expect_error(
    ipcr:::validate_ipcr_components(
      parameters,
      nonfinite_scores,
      valid_bread,
      n = 3L
    ),
    "score matrix contains non-finite values"
  )
})

test_that("confirmed model nonconvergence is rejected", {
  nonconverged_glm <- structure(
    list(converged = FALSE),
    class = c("glm", "lm")
  )

  expect_error(
    ipcr:::check_model_convergence_ipcr(nonconverged_glm),
    "did not converge"
  )
})
