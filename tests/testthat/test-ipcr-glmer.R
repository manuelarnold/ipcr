test_that("glmer IPCs use one row per cluster and omit residual dispersion", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("merDeriv")

  set.seed(816)
  n_clusters <- 20L
  observations_per_cluster <- 8L
  data <- data.frame(
    cluster = factor(rep(seq_len(n_clusters), each = observations_per_cluster)),
    x = stats::rnorm(n_clusters * observations_per_cluster)
  )
  random_intercept <- stats::rnorm(n_clusters, sd = 0.6)
  probability <- stats::plogis(
    -0.4 + 0.5 * data$x + random_intercept[data$cluster]
  )
  data$y <- stats::rbinom(nrow(data), size = 1, prob = probability)

  fit <- lme4::glmer(
    y ~ x + (1 | cluster),
    data = data,
    family = stats::binomial(),
    nAGQ = 5
  )
  predictors <- data.frame(
    group = rep(c(0, 1), length.out = n_clusters),
    row.names = levels(data$cluster)
  )

  result <- ipcr(fit, predictors)
  variance_parameters <- as.data.frame(lme4::VarCorr(fit))$vcov
  names(variance_parameters) <- paste(
    "cov",
    names(lme4::getME(fit, "theta")),
    sep = "_"
  )
  fitted_parameters <- c(lme4::fixef(fit), variance_parameters)

  expect_s3_class(result, "ipcr")
  expect_identical(dim(result$IPCs), c(n_clusters, length(fitted_parameters)))
  expect_identical(colnames(result$IPCs), names(fitted_parameters))
  expect_false("residual" %in% colnames(result$IPCs))
  expect_equal(
    colMeans(result$IPCs),
    fitted_parameters,
    tolerance = 1e-5
  )
  expect_equal(
    result$IPCs,
    get_ipcs(fit),
    tolerance = 1e-10
  )
})
