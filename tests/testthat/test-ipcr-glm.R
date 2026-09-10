test_that("glm IPCs contain regression coefficients but no residual parameter", {
  set.seed(815)
  n <- 120L
  data <- data.frame(
    x = stats::rnorm(n),
    predictor = rep(c(0, 1), each = n / 2L)
  )
  probability <- stats::plogis(-0.4 + 0.7 * data$x)
  data$y <- stats::rbinom(n, size = 1, prob = probability)
  fit <- stats::glm(y ~ x, data = data, family = stats::binomial())

  result <- ipcr(fit, predictors = data["predictor"])

  expect_s3_class(result, "ipcr")
  expect_identical(colnames(result$IPCs), names(stats::coef(fit)))
  expect_false("residual" %in% colnames(result$IPCs))
  expect_equal(
    colMeans(result$IPCs),
    stats::coef(fit),
    tolerance = 1e-6
  )
  expect_equal(
    result$IPCs,
    get_ipcs(fit),
    tolerance = 1e-10
  )
})
