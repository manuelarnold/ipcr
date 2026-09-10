make_regularized_ipcr_example <- function(
    n = 120L, parameters = 2L, include_factor = FALSE) {
  set.seed(1201)
  predictors <- data.frame(
    signal_1 = stats::rnorm(n),
    signal_2 = stats::rnorm(n),
    noise = stats::rnorm(n)
  )
  if (include_factor) {
    predictors$site <- factor(
      rep(c("A", "B", "C"), length.out = n),
      levels = c("A", "B", "C")
    )
  }

  response_1 <- 2.5 * predictors$signal_1 +
    0.15 * stats::rnorm(n)
  if (include_factor) {
    response_1 <- response_1 + ifelse(predictors$site == "B", 1.5, 0)
  }
  if (parameters == 1L) {
    responses <- matrix(
      response_1,
      ncol = 1L,
      dimnames = list(NULL, "theta_1")
    )
  } else {
    response_2 <- -1.5 * predictors$signal_1 +
      2 * predictors$signal_2 + 0.15 * stats::rnorm(n)
    if (include_factor) {
      response_2 <- response_2 - ifelse(predictors$site == "C", 1, 0)
    }
    responses <- cbind(theta_1 = response_1, theta_2 = response_2)
  }

  model_data <- predictors
  model_data$.ipcr_response <- I(responses)
  mlm <- stats::lm(
    stats::reformulate(".", response = ".ipcr_response"),
    data = model_data,
    na.action = stats::na.fail
  )
  result <- list(
    info = list(
      ipcr_type = "standard",
      name = "synthetic regularized-IPCR example",
      class = "synthetic",
      parameters = colnames(responses),
      predictors = names(predictors),
      analytic = FALSE
    ),
    IPCs = responses,
    scores = matrix(0, nrow = n, ncol = parameters),
    predictors = predictors,
    complete_rows = rep(TRUE, n),
    mlm = mlm,
    mancova = NULL,
    output = list()
  )
  class(result) <- "ipcr"
  result
}

regularized_foldid <- function(n, folds = 5L) {
  rep(seq_len(folds), length.out = n)
}

test_that("ipcr_reg uses a multivariate group elastic net without changing IPCs", {
  source <- make_regularized_ipcr_example()
  result <- ipcr_reg(
    source,
    foldid = regularized_foldid(nrow(source$IPCs))
  )

  expect_s3_class(result, "ipcr_reg")
  expect_s3_class(result, "ipcr")
  expect_identical(class(source), "ipcr")
  expect_identical(result$IPCs, source$IPCs)
  expect_identical(result$scores, source$scores)
  expect_identical(result$regularization$family, "mgaussian")
  expect_identical(result$regularization$alpha, 0.5)
  expect_identical(result$regularization$s, "lambda.1se")
  expect_true(result$regularization$standardize)
  expect_true(result$regularization$standardize_response)
  expect_s3_class(result$regularization$fit, "cv.glmnet")
  expect_identical(
    dim(result$regularization$coefficients),
    c(4L, 2L)
  )
  expect_identical(
    rownames(result$regularization$coefficients),
    c("(Intercept)", "signal_1", "signal_2", "noise")
  )
  expect_true(all(c("signal_1", "signal_2") %in%
                    result$regularization$selected_predictors))
  expect_identical(
    result$regularization$importance$Selected,
    result$regularization$importance$Design_Column %in%
      result$regularization$selected_columns
  )
})

test_that("ipcr_reg uses the exact factor contrasts and records their terms", {
  source <- make_regularized_ipcr_example(include_factor = TRUE)
  result <- ipcr_reg(
    source,
    foldid = regularized_foldid(nrow(source$IPCs))
  )

  expect_identical(
    result$regularization$design_columns,
    c("signal_1", "signal_2", "noise", "siteB", "siteC")
  )
  expect_identical(
    unname(result$regularization$design_terms[c("siteB", "siteC")]),
    c("site", "site")
  )
  site_row <- result$regularization$predictor_terms[
    result$regularization$predictor_terms$Predictor == "site",
    ,
    drop = FALSE
  ]
  expect_identical(site_row$Design_Columns, "siteB, siteC")
  expect_true(site_row$Selected)
  expect_true(nzchar(site_row$Selected_Columns))

  newdata <- data.frame(
    signal_1 = c(0, 1),
    signal_2 = c(0, 0),
    noise = c(0, 0),
    site = factor(c("A", "C"), levels = c("A", "B", "C"))
  )
  expect_identical(dim(stats::predict(result, newdata = newdata)), c(2L, 2L))
})

test_that("single-parameter IPCR uses a Gaussian elastic net", {
  source <- make_regularized_ipcr_example(parameters = 1L)
  result <- ipcr_reg(
    source,
    foldid = regularized_foldid(nrow(source$IPCs))
  )

  expect_identical(result$regularization$family, "gaussian")
  expect_null(result$regularization$standardize_response)
  expect_identical(dim(stats::coef(result)), c(4L, 1L))
  expect_identical(dim(stats::fitted(result)), c(nrow(source$IPCs), 1L))
  expect_identical(dim(stats::residuals(result)), c(nrow(source$IPCs), 1L))
})

test_that("regularized accessors use glmnet rather than the ordinary mlm", {
  source <- make_regularized_ipcr_example()
  result <- ipcr_reg(
    source,
    foldid = regularized_foldid(nrow(source$IPCs))
  )

  expect_equal(stats::coef(result), result$regularization$coefficients)
  expect_equal(stats::coefficients(result), stats::coef(result))
  expect_identical(stats::nobs(result), nrow(source$IPCs))
  expect_identical(dim(stats::predict(result)), c(nrow(source$IPCs), 2L))
  expect_equal(
    stats::residuals(result),
    source$IPCs - stats::fitted(result)
  )
  expect_error(stats::confint(result), "not available for regularized IPCR")
  expect_error(stats::vcov(result), "not available for regularized IPCR")
  expect_error(stats::effects(result), "not available for regularized IPCR")
  expect_error(stats::sigma(result), "not available for regularized IPCR")
  expect_error(
    plot_differences(result),
    "not available for regularized IPCR"
  )
})

test_that("ipcr_reg preserves converged provenance and rejects legacy nonconvergence", {
  source <- make_regularized_ipcr_example()
  source$convergence <- list(
    converged = TRUE,
    status = "Converged after 2 iterations."
  )
  class(source) <- c("ipcr_it", "ipcr")
  result <- ipcr_reg(
    source,
    foldid = regularized_foldid(nrow(source$IPCs))
  )

  expect_identical(class(result), c("ipcr_reg", "ipcr_it", "ipcr"))
  expect_identical(result$convergence, source$convergence)

  source$convergence$converged <- FALSE
  expect_error(
    ipcr_reg(
      source,
      foldid = regularized_foldid(nrow(source$IPCs))
    ),
    "requires a converged ipcr_it object"
  )
})

test_that("print and summary identify regularization without inference", {
  source <- make_regularized_ipcr_example(include_factor = TRUE)
  result <- ipcr_reg(
    source,
    foldid = regularized_foldid(nrow(source$IPCs))
  )

  expect_output(print(result), "Multivariate group elastic net", fixed = TRUE)
  expect_output(print(result), "Selected predictors", fixed = TRUE)

  everything <- paste(
    utils::capture.output(summary(result)),
    collapse = "\n"
  )
  predictors <- paste(
    utils::capture.output(summary(result, what = "predictors")),
    collapse = "\n"
  )
  coefficients <- paste(
    utils::capture.output(summary(result, what = "coefficients")),
    collapse = "\n"
  )
  expect_match(everything, "standardized importance", fixed = TRUE)
  expect_match(everything, "Regularized IPC regression coefficients", fixed = TRUE)
  expect_match(everything, "p-values are not available", fixed = TRUE)
  expect_false(grepl("Regularized IPC regression coefficients", predictors, fixed = TRUE))
  expect_false(grepl("standardized importance", coefficients, fixed = TRUE))
  expect_error(summary(result, what = "unknown"), "'what' must be one of")
  expect_error(summary(result, digits = 1.5), "'digits'")
})

test_that("regularized summaries preserve iterated likelihood warnings", {
  source <- make_regularized_ipcr_example()
  source$convergence <- list(
    converged = TRUE,
    status = "Converged after 2 iterations.",
    C_l_rate_adjusted = 0.02,
    likelihood_warning_threshold = 0.001,
    likelihood_stable = FALSE
  )
  class(source) <- c("ipcr_it", "ipcr")
  result <- ipcr_reg(
    source,
    foldid = regularized_foldid(nrow(source$IPCs))
  )

  expect_output(
    summary(result),
    "Likelihood stability warning",
    fixed = TRUE
  )
  expect_output(summary(result), "C_l / learning_rate", fixed = TRUE)
})

test_that("ipcr_reg validates objects, tuning arguments, and folds", {
  source <- make_regularized_ipcr_example()
  n <- nrow(source$IPCs)
  folds <- regularized_foldid(n)

  expect_error(ipcr_reg(stats::lm(mpg ~ wt, data = mtcars)), "returned by ipcr")
  expect_error(ipcr_reg(source, alpha = -0.1), "'alpha'")
  expect_error(ipcr_reg(source, alpha = 1.1), "'alpha'")
  expect_error(ipcr_reg(source, s = "unknown"), "'s'")
  expect_error(ipcr_reg(source, nfolds = 2), "'nfolds'")
  expect_error(ipcr_reg(source, foldid = folds[-1L]), "'foldid'")
  expect_error(ipcr_reg(source, foldid = rep(1:2, length.out = n)), "three folds")
  expect_error(ipcr_reg(source, family = "gaussian"), "controlled by ipcr_reg")

  result <- ipcr_reg(source, foldid = folds)
  expect_error(ipcr_reg(result, foldid = folds), "already regularized")
  expect_warning(ipcr_reg(source, alpha = 0, foldid = folds), "ridge regression")

  one_predictor <- make_regularized_ipcr_example()
  one_predictor$predictors <- one_predictor$predictors["signal_1"]
  one_predictor$mlm <- stats::lm(
    one_predictor$IPCs ~ signal_1,
    data = one_predictor$predictors
  )
  expect_error(ipcr_reg(one_predictor, foldid = folds), "at least two predictor")
})

test_that("ipcr_reg works with a result returned by ipcr", {
  set.seed(1202)
  n <- 90L
  model_data <- data.frame(x = stats::rnorm(n))
  model_data$y <- 1 + model_data$x + stats::rnorm(n)
  predictors <- data.frame(
    group = rep(0:1, length.out = n),
    age = stats::rnorm(n),
    noise = stats::rnorm(n)
  )
  standard <- ipcr(stats::lm(y ~ x, data = model_data), predictors)
  result <- ipcr_reg(
    standard,
    foldid = regularized_foldid(n)
  )

  expect_s3_class(result, "ipcr_reg")
  expect_identical(result$IPCs, standard$IPCs)
  expect_identical(result$complete_rows, standard$complete_rows)
  expect_identical(result$predictors, standard$predictors)
  expect_identical(result$regularization$design_columns, names(predictors))
})
