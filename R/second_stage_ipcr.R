fit_second_stage_ipcr <- function(
    IPCs, predictors, require_full_rank = FALSE) {
  response_name <- ".ipcr_response"
  while (response_name %in% names(predictors)) {
    response_name <- paste0(response_name, "_")
  }

  second_stage_data <- predictors
  second_stage_data[[response_name]] <- I(IPCs)
  second_stage_formula <- stats::reformulate(
    termlabels = ".",
    response = response_name
  )
  mlm <- stats::lm(
    formula = second_stage_formula,
    data = second_stage_data,
    na.action = stats::na.fail
  )

  if (isTRUE(require_full_rank) && any(!is.finite(stats::coef(mlm)))) {
    stop(
      "The iterated IPC regression design matrix is rank deficient.",
      call. = FALSE
    )
  }
  mlm
}

second_stage_output_ipcr <- function(mlm, parameter_names) {
  q <- length(parameter_names)
  mancova <- NULL
  mancova_summary <- NULL
  if (q > 1L) {
    mancova_components <- fit_mancova_ipcr(mlm, q = q)
    mancova <- mancova_components$model
    mancova_summary <- mancova_components$summary
  }

  ipcr_table <- data.frame()
  F_tests_table <- data.frame()
  missing_F_parameters <- character()
  mlm_summaries <- if (inherits(mlm, "mlm")) {
    summary(mlm)
  } else {
    list(summary(mlm))
  }

  for (parameter_index in seq_len(q)) {
    mlm_summary <- mlm_summaries[[parameter_index]]
    coefficient_table <- mlm_summary$coefficients
    ipcr_table <- rbind(
      ipcr_table,
      data.frame(
        Parameter = parameter_names[parameter_index],
        Predictor = rownames(coefficient_table),
        Estimate = coefficient_table[, "Estimate"],
        Std_Error = coefficient_table[, "Std. Error"],
        t_Value = coefficient_table[, "t value"],
        p_Value = coefficient_table[, "Pr(>|t|)"],
        row.names = NULL
      )
    )

    fstatistic <- mlm_summary$fstatistic
    has_fstatistic <- length(fstatistic) == 3L &&
      all(is.finite(fstatistic))
    if (!has_fstatistic) {
      missing_F_parameters <- c(
        missing_F_parameters,
        parameter_names[parameter_index]
      )
    }
    F_value <- if (has_fstatistic) unname(fstatistic[1L]) else NA_real_
    df1 <- if (has_fstatistic) unname(fstatistic[2L]) else NA_real_
    df2 <- if (has_fstatistic) unname(fstatistic[3L]) else NA_real_
    F_p_value <- if (has_fstatistic) {
      stats::pf(F_value, df1 = df1, df2 = df2, lower.tail = FALSE)
    } else {
      NA_real_
    }
    F_tests_table <- rbind(
      F_tests_table,
      data.frame(
        Parameter = parameter_names[parameter_index],
        Rsquared = mlm_summary$r.squared,
        F_Value = F_value,
        DF1 = df1,
        DF2 = df2,
        p_Value = F_p_value,
        row.names = NULL
      )
    )
  }

  if (length(missing_F_parameters) > 0L) {
    warning(
      "The overall F-statistic is unavailable for IPC parameters: ",
      paste(missing_F_parameters, collapse = ", "),
      ". This usually indicates that no predictor effect is estimable or ",
      "that the second-stage regression has insufficient residual variation ",
      "or degrees of freedom.",
      call. = FALSE
    )
  }

  list(
    mancova = mancova,
    output = list(
      ipcr = ipcr_table,
      mancova = mancova_summary,
      F_tests = F_tests_table
    )
  )
}
