test_that("shared second-stage helpers reproduce a standard IPCR result", {
  set.seed(1401)
  n <- 50L
  model_data <- data.frame(
    x = stats::rnorm(n),
    group = rep(0:1, length.out = n)
  )
  model_data$y <- 1 + model_data$x + stats::rnorm(n)
  result <- ipcr(
    stats::lm(y ~ x, data = model_data),
    predictors = model_data["group"]
  )

  rebuilt_mlm <- fit_second_stage_ipcr(
    IPCs = result$IPCs[result$complete_rows, , drop = FALSE],
    predictors = result$predictors[result$complete_rows, , drop = FALSE]
  )
  rebuilt_output <- second_stage_output_ipcr(
    mlm = rebuilt_mlm,
    parameter_names = result$info$parameters
  )

  expect_equal(stats::coef(rebuilt_mlm), stats::coef(result$mlm))
  expect_equal(rebuilt_output$output, result$output)
  expect_equal(rebuilt_output$mancova$SSPE, result$mancova$SSPE)
})

test_that("full-rank enforcement remains specific to iterated IPCR", {
  IPCs <- cbind(
    theta_1 = seq_len(10L),
    theta_2 = seq_len(10L)^2
  )
  predictors <- data.frame(constant = rep(1, NROW(IPCs)))

  standard_fit <- fit_second_stage_ipcr(IPCs, predictors)
  expect_true(anyNA(stats::coef(standard_fit)))
  expect_error(
    fit_second_stage_ipcr(
      IPCs,
      predictors,
      require_full_rank = TRUE
    ),
    "The iterated IPC regression design matrix is rank deficient.",
    fixed = TRUE
  )
})
