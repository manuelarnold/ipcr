make_predictor_handling_example <- function(n = 48L) {
  set.seed(914)
  model_data <- data.frame(
    model_predictor = stats::rnorm(n)
  )
  model_data$outcome <- 1 + 0.5 * model_data$model_predictor +
    stats::rnorm(n)

  list(
    fit = stats::lm(outcome ~ model_predictor, data = model_data),
    n = n
  )
}

test_that("ipcr uses x and predictors as its argument names", {
  expect_identical(
    names(formals(ipcr)),
    c("x", "predictors", "analytic")
  )
})

test_that("unnamed predictor columns receive deterministic names", {
  example <- make_predictor_handling_example()
  predictor_vector <- rep(c(0, 1), length.out = example$n)

  expect_warning(
    vector_result <- ipcr(example$fit, predictor_vector),
    "predictor1"
  )
  expect_identical(vector_result$info$predictors, "predictor1")
  expect_identical(
    rownames(stats::coef(vector_result)),
    c("(Intercept)", "predictor1")
  )

  predictor_matrix <- cbind(
    seq_len(example$n),
    rep(c(FALSE, TRUE), length.out = example$n)
  )
  expect_warning(
    matrix_result <- ipcr(example$fit, predictor_matrix),
    "predictor1, predictor2"
  )
  expect_identical(
    matrix_result$info$predictors,
    c("predictor1", "predictor2")
  )
})

test_that("partially unnamed predictor data retain their existing names", {
  example <- make_predictor_handling_example()
  predictors <- data.frame(
    age = seq_len(example$n),
    group = rep(c(0, 1), length.out = example$n),
    check.names = FALSE
  )
  names(predictors)[2] <- ""

  expect_warning(
    result <- ipcr(example$fit, predictors),
    "predictor1"
  )
  expect_identical(result$info$predictors, c("age", "predictor1"))
})

test_that("predictor names must be unique", {
  example <- make_predictor_handling_example()
  predictors <- data.frame(
    first = seq_len(example$n),
    second = rev(seq_len(example$n)),
    check.names = FALSE
  )
  names(predictors) <- c("duplicate", "duplicate")

  expect_error(
    ipcr(example$fit, predictors),
    "Predictor names must be unique"
  )
})

test_that("numeric, logical, and factor predictors are supported", {
  example <- make_predictor_handling_example()
  predictors <- data.frame(
    age = seq_len(example$n),
    treatment = rep(c(FALSE, TRUE), length.out = example$n),
    site = factor(rep(c("A", "B", "C"), length.out = example$n))
  )

  result <- ipcr(example$fit, predictors)

  expect_s3_class(result, "ipcr")
  expect_identical(
    rownames(stats::coef(result)),
    c("(Intercept)", "age", "treatmentTRUE", "siteB", "siteC")
  )

  correlation_plot <- plot(result)
  expect_s3_class(correlation_plot, "ggplot")
  expect_setequal(
    unique(correlation_plot$data$predictor),
    c("age", "treatmentTRUE", "siteB", "siteC")
  )
})

test_that("unsupported predictor columns are rejected clearly", {
  example <- make_predictor_handling_example()

  expect_error(
    ipcr(
      example$fit,
      data.frame(label = rep(c("A", "B"), length.out = example$n))
    ),
    "convert character columns to factors"
  )
  expect_error(
    ipcr(
      example$fit,
      data.frame(value = c(Inf, seq_len(example$n - 1L)))
    ),
    "infinite values"
  )
  expect_error(
    ipcr(example$fit, list(value = seq_len(example$n))),
    "vector, matrix, or data.frame"
  )
  expect_error(
    ipcr(example$fit, data.frame(row.names = seq_len(example$n))),
    "at least one column"
  )
})

test_that("missing predictors use an explicit complete-case regression", {
  example <- make_predictor_handling_example()
  predictors <- data.frame(
    age = seq_len(example$n),
    group = factor(rep(c("control", "treatment"), length.out = example$n))
  )
  predictors$age[c(3, 17)] <- NA_real_
  complete_rows <- stats::complete.cases(predictors)

  expect_warning(
    result <- ipcr(example$fit, predictors),
    "2 rows with missing predictor values"
  )

  expect_identical(nrow(result$IPCs), example$n)
  expect_identical(nrow(result$scores), example$n)
  expect_identical(stats::nobs(result), sum(complete_rows))
  expect_identical(result$predictors, predictors)
  expect_identical(result$complete_rows, complete_rows)
  expect_false("complete_predictor_rows" %in% names(result$info))

  expected <- stats::lm(
    as.matrix(result$IPCs[complete_rows, , drop = FALSE]) ~ age + group,
    data = predictors[complete_rows, , drop = FALSE]
  )
  expect_equal(stats::coef(result), stats::coef(expected))

  correlation_plot <- plot(result)
  expect_s3_class(correlation_plot, "ggplot")
  expect_false(anyNA(correlation_plot$data$value))
})

test_that("a usable complete-case predictor sample is required", {
  example <- make_predictor_handling_example()

  expect_error(
    ipcr(example$fit, data.frame(predictor = rep(NA_real_, example$n))),
    "no complete rows"
  )

  factor_predictor <- factor(rep(NA_character_, example$n),
                             levels = c("A", "B"))
  factor_predictor[1:2] <- "A"
  expect_warning(
    expect_error(
      ipcr(example$fit, data.frame(site = factor_predictor)),
      "at least two observed levels"
    ),
    "missing predictor values"
  )
})

test_that("Type-II MANCOVA is invariant to predictor order", {
  example <- make_predictor_handling_example(n = 120L)
  set.seed(915)
  first <- stats::rnorm(example$n)
  second <- 0.6 * first + stats::rnorm(example$n)
  predictors <- data.frame(first = first, second = second)

  result_first_second <- ipcr(example$fit, predictors[c("first", "second")])
  result_second_first <- ipcr(example$fit, predictors[c("second", "first")])

  expect_s3_class(result_first_second$mancova, "Anova.mlm")
  expect_s3_class(result_second_first$mancova, "Anova.mlm")
  expect_identical(result_first_second$mancova$type, "II")
  expect_identical(result_second_first$mancova$type, "II")

  for (predictor in names(predictors)) {
    expect_equal(
      result_first_second$mancova$SSP[[predictor]],
      result_second_first$mancova$SSP[[predictor]],
      tolerance = 1e-10
    )
    expect_equal(
      result_first_second$mancova$df[[predictor]],
      result_second_first$mancova$df[[predictor]]
    )
  }
  expect_equal(
    result_first_second$mancova$SSPE,
    result_second_first$mancova$SSPE,
    tolerance = 1e-10
  )
})
