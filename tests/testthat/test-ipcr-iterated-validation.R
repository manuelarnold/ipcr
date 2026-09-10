test_that("common iterated arguments share one validator", {
  valid_arguments <- list(
    analytic = TRUE,
    conv = 0.01,
    learning_rate = c(1, 0.5, 0.1),
    max_it = 350L,
    iteration_info = FALSE
  )

  expect_invisible(
    do.call(validate_ipcr_it_common_arguments, valid_arguments)
  )

  invalid_arguments <- list(
    analytic = list(
      value = NA,
      message = "'analytic' must be either TRUE or FALSE."
    ),
    conv = list(
      value = 0,
      message = "'conv' must be one positive finite number."
    ),
    learning_rate = list(
      value = c(1, 0),
      message = paste0(
        "'learning_rate' must contain finite values greater than 0 ",
        "and no greater than 1."
      )
    ),
    max_it = list(
      value = 1.5,
      message = "'max_it' must be one positive whole number."
    ),
    iteration_info = list(
      value = NA,
      message = "'iteration_info' must be either TRUE or FALSE."
    )
  )

  for (argument in names(invalid_arguments)) {
    arguments <- valid_arguments
    arguments[[argument]] <- invalid_arguments[[argument]]$value
    expect_error(
      do.call(validate_ipcr_it_common_arguments, arguments),
      invalid_arguments[[argument]]$message,
      fixed = TRUE
    )
  }
})
