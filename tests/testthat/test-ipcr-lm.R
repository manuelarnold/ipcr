make_lm_slope_heterogeneity_example <- function() {
  set.seed(4)

  n1 <- 100L
  n2 <- 100L
  n <- n1 + n2

  group_predictor <- rep(c(0, 1), times = c(n1, n2))
  x <- stats::rnorm(n)
  y <- 1 + 0.4 * x * group_predictor + stats::rnorm(n)

  list(
    fit = stats::lm(y ~ x),
    predictors = data.frame(group_predictor = group_predictor),
    n = n
  )
}

test_that("ipcr fits the lm slope-heterogeneity example", {
  example <- make_lm_slope_heterogeneity_example()

  result <- ipcr(
    x = example$fit,
    predictors = example$predictors
  )

  expect_s3_class(result, "ipcr")
  expect_named(
    result,
    c(
      "info", "IPCs", "scores", "predictors", "complete_rows",
      "mlm", "mancova", "output"
    ),
    ignore.order = FALSE
  )
  expect_s3_class(result$mlm, "mlm")
  expect_s3_class(result$mancova, "Anova.mlm")

  expect_identical(dim(result$IPCs), c(example$n, 3L))
  expect_identical(dim(result$scores), c(example$n, 3L))
  expect_identical(result$predictors, example$predictors)
  expect_identical(result$complete_rows, rep(TRUE, example$n))
  expect_identical(
    colnames(result$IPCs),
    c("(Intercept)", "x", "residual")
  )

  ipcr_coefficients <- stats::coef(result$mlm)
  expect_identical(
    rownames(ipcr_coefficients),
    c("(Intercept)", "group_predictor")
  )
})

test_that("lm IPC means reproduce the fitted model parameters", {
  example <- make_lm_slope_heterogeneity_example()

  result <- ipcr(
    x = example$fit,
    predictors = example$predictors
  )

  expected_parameters <- c(
    stats::coef(example$fit),
    residual = mean(stats::residuals(example$fit)^2)
  )

  expect_equal(
    colMeans(result$IPCs),
    expected_parameters,
    tolerance = 1e-10
  )
  expect_equal(
    result$IPCs,
    get_ipcs(example$fit),
    tolerance = 1e-10
  )
  expect_true(is.matrix(get_ipcs(example$fit)))
  expect_true(is.matrix(get_scores(example$fit)))
  expect_equal(result$scores, get_scores(example$fit), tolerance = 1e-12)
  expect_equal(
    result$IPCs[, "residual"],
    stats::residuals(example$fit)^2,
    tolerance = 1e-10
  )
})

test_that("ipcr requires one predictor row per lm observation", {
  example <- make_lm_slope_heterogeneity_example()
  incomplete_predictors <- example$predictors[-example$n, , drop = FALSE]

  expect_error(
    ipcr(
      x = example$fit,
      predictors = incomplete_predictors
    ),
    "number of rows in 'predictors'"
  )
})
