make_plotting_data <- function(n = 90L) {
  set.seed(820)
  data <- data.frame(x = stats::rnorm(n))
  data$y <- 1 + 0.7 * data$x + stats::rnorm(n)
  data
}

fit_plotting_ipcr <- function(predictors) {
  data <- make_plotting_data(nrow(predictors))
  ipcr(stats::lm(y ~ x, data = data), predictors)
}

test_that("heatmaps build for numeric, binary, and multi-level factors", {
  n <- 90L
  cases <- list(
    numeric = list(
      predictors = data.frame(age = seq_len(n)),
      design_names = "age"
    ),
    binary = list(
      predictors = data.frame(group = rep(c(0, 1), each = n / 2L)),
      design_names = "group"
    ),
    factor = list(
      predictors = data.frame(
        site = factor(rep(c("A", "B", "C"), length.out = n))
      ),
      design_names = c("siteB", "siteC")
    )
  )

  for (case in cases) {
    result <- fit_plotting_ipcr(case$predictors)
    correlation_plot <- plot(result)

    expect_s3_class(correlation_plot, "ggplot")
    expect_no_error(ggplot2::ggplot_build(correlation_plot))
    expect_setequal(
      unique(correlation_plot$data$predictor),
      case$design_names
    )
  }
})

test_that("heatmap annotations use named aesthetics", {
  result <- fit_plotting_ipcr(data.frame(age = seq_len(60L)))
  correlation_plot <- plot(result, print_corr = TRUE)

  expect_length(correlation_plot$layers, 2L)
  expect_no_error(ggplot2::ggplot_build(correlation_plot))
  expect_error(plot(result, print_corr = NA), "'print_corr'")
  expect_error(plot(result, print_corr = c(TRUE, FALSE)), "'print_corr'")
})

test_that("heatmaps retain and identify undefined correlations", {
  result <- suppressWarnings(
    fit_plotting_ipcr(data.frame(constant_predictor = rep(1, 60L)))
  )

  expect_warning(
    correlation_plot <- plot(result),
    "predictor columns: constant_predictor"
  )
  expect_true(all(is.na(correlation_plot$data$value)))
  built_plot <- ggplot2::ggplot_build(correlation_plot)
  expect_true(all(built_plot$data[[1L]]$fill == "grey80"))
})

test_that("plot_differences returns a continuous plot with t intervals", {
  n <- 90L
  result <- fit_plotting_ipcr(data.frame(age = seq(-2, 2, length.out = n)))

  difference_plot <- plot_differences(
    result,
    parameter = "x",
    predictor = "age",
    confidence_level = 0.95
  )
  expect_s3_class(difference_plot, "ggplot")
  expect_no_error(ggplot2::ggplot_build(difference_plot))
  expect_true(is.numeric(difference_plot$data$predictor_value))

  prediction_data <- difference_plot$layers[[2L]]$data
  positive_se <- prediction_data$std_error > 0
  observed_critical_values <- (
    prediction_data$upper[positive_se] -
      prediction_data$estimate[positive_se]
  ) / prediction_data$std_error[positive_se]
  expect_equal(
    observed_critical_values,
    rep(
      stats::qt(0.975, df = stats::df.residual(result$mlm)),
      sum(positive_se)
    ),
    tolerance = 1e-12
  )
})

test_that("plot_differences handles binary and multi-level predictors", {
  n <- 90L
  predictors <- data.frame(
    group = rep(c(0, 1), each = n / 2L),
    site = factor(rep(c("A", "B", "C"), length.out = n))
  )
  result <- fit_plotting_ipcr(predictors)

  group_plot <- plot_differences(
    result,
    parameter = "x",
    predictor = "group"
  )
  site_plot <- plot_differences(
    result,
    parameter = "x",
    predictor = "site"
  )

  expect_s3_class(group_plot, "ggplot")
  expect_true(is.factor(group_plot$data$predictor_value))
  expect_identical(
    as.character(group_plot$layers[[2L]]$data$predictor_value),
    c("0", "1")
  )
  expect_no_error(ggplot2::ggplot_build(group_plot))

  expect_s3_class(site_plot, "ggplot")
  expect_identical(
    as.character(site_plot$layers[[2L]]$data$predictor_value),
    c("A", "B", "C")
  )
  expect_no_error(ggplot2::ggplot_build(site_plot))
})

test_that("arbitrary two-valued numeric predictors remain continuous", {
  result <- fit_plotting_ipcr(data.frame(dose = rep(c(2, 5), each = 45L)))

  difference_plot <- plot_differences(
    result,
    parameter = "x",
    predictor = "dose"
  )

  expect_true(is.numeric(difference_plot$data$predictor_value))
  expect_identical(nrow(difference_plot$layers[[2L]]$data), 100L)
  expect_equal(
    range(difference_plot$layers[[2L]]$data$predictor_value),
    c(2, 5)
  )
  expect_no_error(ggplot2::ggplot_build(difference_plot))
})

test_that("plot_differences returns named plots for multiple requests", {
  n <- 90L
  result <- fit_plotting_ipcr(data.frame(
    age = seq(-2, 2, length.out = n),
    group = rep(c(0, 1), each = n / 2L)
  ))

  expect_silent(
    plots <- plot_differences(
      result,
      parameter = c("x", "residual"),
      predictor = c("age", "group")
    )
  )
  expect_type(plots, "list")
  expect_length(plots, 4L)
  expect_named(plots, c(
    "parameter: x | predictor: age",
    "parameter: x | predictor: group",
    "parameter: residual | predictor: age",
    "parameter: residual | predictor: group"
  ))
  expect_true(all(vapply(plots, inherits, logical(1), what = "ggplot")))
})

test_that("plot_differences uses the complete second-stage sample", {
  n <- 90L
  predictors <- data.frame(age = seq_len(n))
  predictors$age[c(4L, 21L)] <- NA_real_

  expect_warning(
    result <- fit_plotting_ipcr(predictors),
    "2 rows with missing predictor values"
  )
  difference_plot <- plot_differences(
    result,
    parameter = "x",
    predictor = "age"
  )

  expect_identical(nrow(difference_plot$data), sum(result$complete_rows))
  expect_false(anyNA(difference_plot$data))
})

test_that("plot_differences validates names and confidence levels", {
  result <- fit_plotting_ipcr(data.frame(age = seq_len(60L)))

  expect_error(
    plot_differences(result, parameter = "missing", predictor = "age"),
    "Unknown parameter name"
  )
  expect_error(
    plot_differences(result, parameter = "x", predictor = "missing"),
    "Unknown predictor name"
  )
  expect_error(
    plot_differences(result, parameter = c("x", "x"), predictor = "age"),
    "must not contain duplicate"
  )
  expect_error(
    plot_differences(result, parameter = "x", predictor = character()),
    "nonempty character vector"
  )
  expect_error(
    plot_differences(result, confidence_level = 0),
    "strictly between 0 and 1"
  )
  expect_error(
    plot_differences(result, confidence_level = 1),
    "strictly between 0 and 1"
  )
  expect_error(
    plot_differences(result, confidence_level = c(0.9, 0.95)),
    "confidence_level"
  )
  expect_error(plot_differences(list()), "class \"ipcr\"")
})

test_that("plot_differences rejects an unestimable predictor effect", {
  result <- suppressWarnings(
    fit_plotting_ipcr(data.frame(constant_predictor = rep(1, 60L)))
  )

  expect_error(
    plot_differences(
      result,
      parameter = "x",
      predictor = "constant_predictor"
    ),
    "not estimable"
  )
})

test_that("both plotting functions support a single IPC parameter", {
  set.seed(821)
  data <- data.frame(
    x = seq_len(30L),
    predictor = rep(c(0, 1), each = 15L)
  )
  data$y <- 2 * data$x + stats::rnorm(nrow(data))
  fit <- stats::nls(y ~ slope * x, data = data, start = list(slope = 1))

  expect_warning(
    result <- ipcr(fit, data["predictor"]),
    "not officially supported"
  )
  correlation_plot <- plot(result)
  difference_plot <- plot_differences(
    result,
    parameter = "slope",
    predictor = "predictor"
  )

  expect_s3_class(correlation_plot, "ggplot")
  expect_s3_class(difference_plot, "ggplot")
  expect_no_error(ggplot2::ggplot_build(correlation_plot))
  expect_no_error(ggplot2::ggplot_build(difference_plot))
})
