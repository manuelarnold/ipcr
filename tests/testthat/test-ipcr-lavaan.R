make_lavaan_visual_cfa_example <- function() {
  data <- lavaan::HolzingerSwineford1939
  data$age_months <- 12 * data$ageyr + data$agemo
  data <- data[stats::complete.cases(data), ]

  predictors <- data.frame(
    sex = data$sex - 1,
    months = as.numeric(scale(data$age_months, center = TRUE, scale = FALSE)),
    school = as.numeric(data$school) - 1,
    grade = data$grade - 7
  )

  fit <- lavaan::cfa(
    model = "visual =~ x1 + x2 + x3",
    data = data
  )

  list(
    fit = fit,
    predictors = predictors,
    n = nrow(data)
  )
}

test_that("ipcr fits the lavaan visual CFA from Arnold et al. (2021)", {
  skip_if_not_installed("lavaan")
  example <- make_lavaan_visual_cfa_example()

  result <- ipcr(
    x = example$fit,
    predictors = example$predictors
  )
  fitted_parameters <- lavaan::coef(example$fit)

  expect_s3_class(result, "ipcr")
  expect_s4_class(example$fit, "lavaan")
  expect_length(fitted_parameters, 6L)
  expect_identical(
    dim(result$IPCs),
    c(example$n, length(fitted_parameters))
  )
  expect_identical(colnames(result$IPCs), names(fitted_parameters))
  expect_identical(
    rownames(stats::coef(result$mlm)),
    c("(Intercept)", "sex", "months", "school", "grade")
  )
  expect_lt(
    max(abs(colMeans(result$IPCs) - fitted_parameters)),
    1e-6
  )
  expect_equal(
    result$IPCs,
    get_ipcs(example$fit),
    tolerance = 1e-8
  )
  expect_identical(nrow(result$output$ipcr), 30L)
  expect_identical(nrow(result$output$F_tests), 6L)
})

test_that("lavaan simple equality constraints use reduced coordinates", {
  skip_if_not_installed("lavaan")
  data <- lavaan::HolzingerSwineford1939
  predictors <- data.frame(sex = data$sex - 1)
  fits <- list(
    shared_label = lavaan::cfa(
      "visual =~ x1 + a*x2 + a*x3",
      data = data
    ),
    explicit_equality = lavaan::cfa(
      paste(
        "visual =~ x1 + a*x2 + b*x3",
        "a == b",
        sep = "\n"
      ),
      data = data
    ),
    simple_constraint_storage = lavaan::cfa(
      "visual =~ x1 + a*x2 + a*x3",
      data = data,
      ceq.simple = TRUE
    )
  )

  for (fit in fits) {
    result <- ipcr(fit, predictors)
    fitted_parameters <- ipcr:::coef_ipcr.lavaan(fit)
    bread <- ipcr:::bread_ipcr.lavaan(fit)

    expect_true(
      isTRUE(fit@Model@eq.constraints) ||
        isTRUE(fit@Model@ceq.simple.only)
    )
    expect_length(fitted_parameters, 5L)
    expect_identical(dim(result$scores), c(nrow(data), 5L))
    expect_identical(dim(bread), c(5L, 5L))
    expect_identical(colnames(result$IPCs), names(fitted_parameters))
    expect_lt(
      max(abs(colMeans(result$IPCs) - fitted_parameters)),
      1e-6
    )
  }
})

test_that("non-simple lavaan constraints are rejected", {
  skip_if_not_installed("lavaan")
  data <- lavaan::HolzingerSwineford1939
  fit <- lavaan::cfa(
    paste(
      "visual =~ x1 + a*x2 + b*x3",
      "a == b^2",
      sep = "\n"
    ),
    data = data
  )

  expect_error(
    ipcr(fit, data.frame(sex = data$sex - 1)),
    "Only simple linear equality constraints"
  )
})

test_that("multigroup lavaan models retain all observations", {
  skip_if_not_installed("lavaan")
  data <- lavaan::HolzingerSwineford1939
  fit <- lavaan::cfa(
    "visual =~ x1 + x2 + x3",
    data = data,
    group = "school"
  )
  predictors <- data.frame(sex = data$sex - 1)
  result <- ipcr(fit, predictors)
  fitted_parameters <- lavaan::coef(fit)

  expect_identical(lavaan::lavInspect(fit, "ngroups"), 2L)
  expect_identical(dim(result$IPCs), c(nrow(data), 18L))
  expect_identical(dim(result$scores), c(nrow(data), 18L))
  expect_identical(colnames(result$IPCs), names(fitted_parameters))
  expect_lt(
    max(abs(colMeans(result$IPCs) - fitted_parameters)),
    1e-6
  )
  expect_equal(
    result$IPCs,
    get_ipcs(fit),
    tolerance = 1e-8
  )

  constrained_fit <- lavaan::cfa(
    "visual =~ x1 + x2 + x3",
    data = data,
    group = "school",
    group.equal = "loadings"
  )
  constrained_result <- ipcr(constrained_fit, predictors)
  constrained_parameters <- ipcr:::coef_ipcr.lavaan(constrained_fit)
  expect_identical(
    dim(constrained_result$IPCs),
    c(nrow(data), length(constrained_parameters))
  )
  expect_lt(
    max(abs(colMeans(constrained_result$IPCs) - constrained_parameters)),
    1e-6
  )
})

test_that("lavaan FIML models provide casewise IPCs", {
  skip_if_not_installed("lavaan")
  data <- lavaan::HolzingerSwineford1939
  data[seq(2L, nrow(data), by = 7L), "x1"] <- NA_real_
  data[seq(3L, nrow(data), by = 11L), "x3"] <- NA_real_
  fit <- lavaan::cfa(
    "visual =~ x1 + x2 + x3",
    data = data,
    missing = "fiml"
  )
  predictors <- data.frame(sex = data$sex - 1)
  result <- ipcr(fit, predictors)
  fitted_parameters <- lavaan::coef(fit)

  expect_identical(fit@Options$missing, "ml")
  expect_identical(dim(result$IPCs), c(nrow(data), 9L))
  expect_identical(dim(result$scores), c(nrow(data), 9L))
  expect_true(all(is.finite(result$scores)))
  expect_lt(max(abs(colSums(result$scores))), 1e-4)
  expect_lt(
    max(abs(colMeans(result$IPCs) - fitted_parameters)),
    1e-6
  )
})

test_that("continuous non-ML lavaan estimators use lavaan components", {
  skip_if_not_installed("lavaan")
  data <- lavaan::HolzingerSwineford1939
  predictors <- data.frame(sex = data$sex - 1)

  for (estimator in c("GLS", "ULS", "WLS")) {
    fit <- lavaan::cfa(
      "visual =~ x1 + x2 + x3",
      data = data,
      estimator = estimator
    )
    result <- ipcr(fit, predictors)
    fitted_parameters <- lavaan::coef(fit)

    expect_identical(fit@Options$estimator, estimator)
    expect_identical(dim(result$scores), c(nrow(data), 6L))
    expect_true(all(is.finite(result$scores)))
    expect_lt(
      max(abs(colMeans(result$IPCs) - fitted_parameters)),
      1e-6
    )
  }
})

make_ordered_lavaan_example <- function(n = 200L) {
  set.seed(927)
  factor_score <- stats::rnorm(n)
  make_indicator <- function(loading) {
    as.ordered(cut(
      loading * factor_score + stats::rnorm(n),
      breaks = c(-Inf, -0.5, 0.5, Inf),
      labels = c("low", "middle", "high")
    ))
  }

  data.frame(
    i1 = make_indicator(0.8),
    i2 = make_indicator(0.7),
    i3 = make_indicator(0.9),
    predictor = stats::rnorm(n)
  )
}

test_that("ordered lavaan indicators are supported with full WLS", {
  skip_if_not_installed("lavaan")
  data <- make_ordered_lavaan_example()
  fit <- lavaan::cfa(
    "factor =~ i1 + i2 + i3",
    data = data,
    ordered = c("i1", "i2", "i3"),
    estimator = "WLS"
  )
  result <- ipcr(fit, data["predictor"])
  fitted_parameters <- lavaan::coef(fit)

  expect_true(fit@Model@categorical)
  expect_identical(fit@Options$estimator, "WLS")
  expect_identical(dim(result$scores), c(nrow(data), 9L))
  expect_true(all(is.finite(result$scores)))
  expect_lt(
    max(abs(colMeans(result$IPCs) - fitted_parameters)),
    1e-6
  )
})

test_that("ordered lavaan indicators reject estimators without scores", {
  skip_if_not_installed("lavaan")
  data <- make_ordered_lavaan_example()
  fit <- lavaan::cfa(
    "factor =~ i1 + i2 + i3",
    data = data,
    ordered = c("i1", "i2", "i3")
  )

  expect_identical(fit@Options$estimator.orig, "WLSMV")
  expect_error(
    ipcr(fit, data["predictor"]),
    "require estimator = 'WLS'"
  )
})

test_that("ipcr_it recovers group differences in a lavaan location-scale model", {
  skip_if_not_installed("lavaan")
  set.seed(928)
  n <- 100L
  group <- rep(c(0, 1), each = n / 2L)
  y <- 0.7 * group + stats::rnorm(
    n,
    sd = ifelse(group == 0, 0.8, 1)
  )
  data <- data.frame(y = y)
  fit <- lavaan::sem(
    paste("y ~~ y", "y ~ 1", sep = "\n"),
    data = data,
    meanstructure = TRUE
  )
  named_lavaan_fit <- fit

  result <- ipcr_it(
    named_lavaan_fit,
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
  expect_identical(result$info$name, "named_lavaan_fit")
  expect_true(result$convergence$converged)
  expect_identical(result$convergence$attempts$learning_rate, 1)
  expect_equal(
    stats::coef(result)["group", "y~1"],
    unname(diff(group_means)),
    tolerance = 1e-6
  )
  expect_equal(
    stats::coef(result)["group", "y~~y"],
    unname(diff(group_variances)),
    tolerance = 1e-6
  )
  expect_identical(dim(result$IPCs), c(n, 2L))
  expect_identical(dim(result$scores), c(n, 2L))
  expect_identical(tail(result$iteration_history$status, 1L), "converged")
})

test_that("iterated lavaan RAM scores agree at the pooled estimates", {
  skip_if_not_installed("lavaan")
  data <- lavaan::HolzingerSwineford1939
  fit <- lavaan::cfa(
    "visual =~ x1 + a*x2 + a*x3",
    data = data,
    meanstructure = TRUE
  )
  predictors <- data.frame(sex = data$sex - 1)
  standard_result <- ipcr(fit, predictors)
  spec <- lavaan_ram_spec_ipcr(fit)
  parameter_estimates <- coef_ipcr(fit)
  parameter_values <- matrix(
    parameter_estimates,
    nrow = nrow(spec$data),
    ncol = length(parameter_estimates),
    byrow = TRUE
  )
  components <- analytical_iterated_lavaan_components_ipcr(
    spec = spec,
    rows = seq_len(nrow(spec$data)),
    parameter_values = parameter_values,
    learning_rate = 1
  )
  likelihood <- analytical_iterated_lavaan_components_ipcr(
    spec = spec,
    rows = seq_len(nrow(spec$data)),
    parameter_values = parameter_values,
    learning_rate = 1,
    likelihood_only = TRUE
  )

  expect_equal(components$scores, standard_result$scores, tolerance = 1e-10)
  expect_equal(components$IPCs, standard_result$IPCs, tolerance = 1e-10)
  expect_identical(names(likelihood), "log_likelihood")
  expect_equal(
    likelihood$log_likelihood,
    components$log_likelihood,
    tolerance = 1e-12
  )

  expect_warning(
    result <- ipcr_it(
      fit,
      predictors = predictors,
      conv = 1e6,
      learning_rate = 1,
      max_it = 1
    ),
    "learning-rate-adjusted likelihood change",
    fixed = TRUE
  )
  expect_true(result$convergence$converged)
  expect_false(result$convergence$likelihood_stable)
  expect_identical(NCOL(result$IPCs), 8L)
})

test_that("ipcr_it supports lavaan FIML models", {
  skip_if_not_installed("lavaan")
  data <- lavaan::HolzingerSwineford1939
  data[seq(2L, nrow(data), by = 7L), "x1"] <- NA_real_
  data[seq(3L, nrow(data), by = 11L), "x3"] <- NA_real_
  fit <- lavaan::cfa(
    "visual =~ x1 + x2 + x3",
    data = data,
    missing = "fiml"
  )
  predictors <- data.frame(sex = data$sex - 1)
  standard_result <- ipcr(fit, predictors)
  spec <- lavaan_ram_spec_ipcr(fit)
  parameter_estimates <- coef_ipcr(fit)
  parameter_values <- matrix(
    parameter_estimates,
    nrow = nrow(spec$data),
    ncol = length(parameter_estimates),
    byrow = TRUE
  )
  components <- analytical_iterated_lavaan_components_ipcr(
    spec = spec,
    rows = seq_len(nrow(spec$data)),
    parameter_values = parameter_values,
    learning_rate = 1
  )
  likelihood <- analytical_iterated_lavaan_components_ipcr(
    spec = spec,
    rows = seq_len(nrow(spec$data)),
    parameter_values = parameter_values,
    learning_rate = 1,
    likelihood_only = TRUE
  )

  expect_equal(components$scores, standard_result$scores, tolerance = 1e-10)
  expect_equal(components$IPCs, standard_result$IPCs, tolerance = 1e-10)
  expect_identical(names(likelihood), "log_likelihood")
  expect_equal(
    likelihood$log_likelihood,
    components$log_likelihood,
    tolerance = 1e-12
  )

  expect_warning(
    result <- ipcr_it(
      fit,
      predictors = predictors,
      conv = 1e6,
      learning_rate = 1,
      max_it = 1
    ),
    "learning-rate-adjusted likelihood change",
    fixed = TRUE
  )
  expect_true(result$convergence$converged)
  expect_false(result$convergence$likelihood_stable)
  expect_true(all(is.finite(result$IPCs)))
})

test_that("ipcr_it rejects unsupported lavaan iteration scopes", {
  skip_if_not_installed("lavaan")
  data <- lavaan::HolzingerSwineford1939
  predictors <- data.frame(sex = data$sex - 1)
  multigroup_fit <- lavaan::cfa(
    "visual =~ x1 + x2 + x3",
    data = data,
    group = "school"
  )
  gls_fit <- lavaan::cfa(
    "visual =~ x1 + x2 + x3",
    data = data,
    estimator = "GLS"
  )

  expect_error(
    ipcr_it(multigroup_fit, predictors),
    "only single-group lavaan models"
  )
  expect_error(
    ipcr_it(gls_fit, predictors),
    "only continuous lavaan models fitted by ML"
  )
  expect_error(
    ipcr_it(
      gls_fit,
      predictors,
      variance_parameterization = "log"
    ),
    "available only for analytical OpenMx MxRAMModel"
  )
})
