test_that("ipcr fits a simple lmer random-intercept model", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("merDeriv")

  data <- lme4::sleepstudy
  fit <- lme4::lmer(
    Reaction ~ Days + (1 | Subject),
    data = data,
    REML = FALSE
  )
  cluster_names <- levels(data$Subject)
  predictors <- data.frame(
    group = rep(c(0, 1), length.out = length(cluster_names)),
    row.names = cluster_names
  )

  result <- ipcr(x = fit, predictors = predictors)
  scores <- get_scores(fit)

  theta <- lme4::getME(fit, "theta")
  variance_parameters <- as.data.frame(lme4::VarCorr(fit))$vcov
  names(variance_parameters) <- c(
    paste("cov", names(theta), sep = "_"),
    "residual"
  )
  fitted_parameters <- c(
    lme4::fixef(fit),
    variance_parameters
  )

  expect_s3_class(result, "ipcr")
  expect_s4_class(fit, "lmerMod")
  expect_true(is.matrix(scores))
  expect_identical(
    dim(scores),
    c(length(cluster_names), length(fitted_parameters))
  )
  expect_identical(
    dim(result$IPCs),
    c(length(cluster_names), length(fitted_parameters))
  )
  expect_identical(colnames(result$IPCs), names(fitted_parameters))
  expect_identical(
    rownames(stats::coef(result$mlm)),
    c("(Intercept)", "group")
  )
  expect_lt(
    max(abs(colMeans(result$IPCs) - fitted_parameters)),
    1e-5
  )
  expect_equal(
    result$IPCs,
    get_ipcs(fit),
    tolerance = 1e-10
  )
  expect_identical(nrow(result$output$ipcr), 8L)
  expect_identical(nrow(result$output$F_tests), 4L)
})

test_that("lmer support currently requires ML and one grouping factor", {
  reml_fit <- lme4::lmer(
    Reaction ~ Days + (1 | Subject),
    data = lme4::sleepstudy,
    REML = TRUE
  )
  predictors <- data.frame(group = rep(c(0, 1), length.out = 18L))

  expect_error(
    ipcr(reml_fit, predictors),
    "REML = FALSE"
  )
})

test_that("lmer random-slope variances follow the merDeriv parameter order", {
  fit <- lme4::lmer(
    Reaction ~ Days + (Days | Subject),
    data = lme4::sleepstudy,
    REML = FALSE
  )
  cluster_names <- levels(lme4::sleepstudy$Subject)
  predictors <- data.frame(
    group = rep(c(0, 1), length.out = length(cluster_names)),
    row.names = cluster_names
  )
  variance_parameters <- as.data.frame(
    lme4::VarCorr(fit),
    order = "lower.tri"
  )$vcov
  names(variance_parameters) <- c(
    paste("cov", names(lme4::getME(fit, "theta")), sep = "_"),
    "residual"
  )
  fitted_parameters <- c(lme4::fixef(fit), variance_parameters)

  result <- ipcr(fit, predictors)

  expect_identical(colnames(result$IPCs), names(fitted_parameters))
  expect_equal(
    colMeans(result$IPCs),
    fitted_parameters,
    tolerance = 1e-4
  )
})
