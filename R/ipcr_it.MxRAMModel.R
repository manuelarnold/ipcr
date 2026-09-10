validate_ipcr_it_MxRAMModel_arguments <- function(
    x, analytic, conv, learning_rate, max_it, iteration_info,
    variance_parameterization) {
  if (!inherits(x, "MxRAMModel")) {
    stop("'x' must inherit from OpenMx class 'MxRAMModel'.", call. = FALSE)
  }
  validate_ipcr_it_common_arguments(
    analytic = analytic,
    conv = conv,
    learning_rate = learning_rate,
    max_it = max_it,
    iteration_info = iteration_info
  )
  if (!is.character(variance_parameterization) ||
      length(variance_parameterization) != 1L ||
      is.na(variance_parameterization) ||
      !variance_parameterization %in% c("original", "log")) {
    stop(
      "'variance_parameterization' must be either \"original\" or \"log\".",
      call. = FALSE
    )
  }
  if (identical(variance_parameterization, "log") && !isTRUE(analytic)) {
    stop(
      "Log-variance parameterization is available only with analytic = TRUE.",
      call. = FALSE
    )
  }
  if (!isTRUE(analytic)) {
    warning(
      "Setting analytic = FALSE uses repeated numerical OpenMx score and ",
      "Jacobian calculations and may lead to a very long run time.",
      call. = FALSE
    )
  }

  check_model_scope_ipcr(x)
  check_model_convergence_ipcr(x)
  if (isTRUE(analytic)) {
    analytical_issue <- openmx_analytical_score_issue_ipcr(x)
    if (!is.null(analytical_issue)) {
      stop(
        "Analytical iterated IPCR does not support ",
        analytical_issue,
        ".",
        call. = FALSE
      )
    }
  }
  if (identical(variance_parameterization, "log")) {
    ram_variance_parameterization_ipcr(x, variance_parameterization)
  }

  invisible(TRUE)
}

ram_variance_parameterization_ipcr <- function(
    x, variance_parameterization = c("original", "log")) {
  variance_parameterization <- match.arg(variance_parameterization)
  natural_estimates <- coef_ipcr(x)
  natural_names <- names(natural_estimates)
  transformed <- rep(FALSE, length(natural_names))

  if (identical(variance_parameterization, "log")) {
    parameter_labels <- ram_parameter_labels_ipcr(x)
    diagonal <- row(parameter_labels$S) == col(parameter_labels$S)
    diagonal_labels <- unique(
      parameter_labels$S[diagonal & !is.na(parameter_labels$S)]
    )
    transformed <- natural_names %in% diagonal_labels
    if (!any(transformed)) {
      stop(
        "Log-variance parameterization requires at least one free parameter on the diagonal of the RAM S matrix.",
        call. = FALSE
      )
    }

    for (parameter_name in natural_names[transformed]) {
      occurs_in_A <- any(
        parameter_labels$A == parameter_name & !is.na(parameter_labels$A)
      )
      occurs_in_M <- any(
        parameter_labels$M == parameter_name & !is.na(parameter_labels$M)
      )
      occurs_off_diagonal_S <- any(
        parameter_labels$S == parameter_name &
          !is.na(parameter_labels$S) & !diagonal
      )
      if (occurs_in_A || occurs_in_M || occurs_off_diagonal_S) {
        stop(
          "Free parameter '", parameter_name, "' occurs both on the RAM S ",
          "diagonal and at a non-variance RAM location. A parameter cannot ",
          "use both log-variance and original coordinates.",
          call. = FALSE
        )
      }
    }
    if (any(!is.finite(natural_estimates[transformed]) |
            natural_estimates[transformed] <= 0)) {
      offending <- natural_names[transformed][
        !is.finite(natural_estimates[transformed]) |
          natural_estimates[transformed] <= 0
      ]
      stop(
        "Log-variance parameterization requires positive finite fitted ",
        "variance estimates: ", paste(offending, collapse = ", "), ".",
        call. = FALSE
      )
    }
  }

  parameter_names <- natural_names
  parameter_names[transformed] <- paste0(
    "log(", natural_names[transformed], ")"
  )
  if (anyDuplicated(parameter_names)) {
    stop(
      "Log-transformed and original OpenMx parameter names are not unique.",
      call. = FALSE
    )
  }
  parameter_estimates <- natural_estimates
  parameter_estimates[transformed] <- log(
    natural_estimates[transformed]
  )
  names(parameter_estimates) <- parameter_names

  list(
    type = variance_parameterization,
    natural_names = natural_names,
    parameter_names = parameter_names,
    transformed = transformed,
    natural_estimates = natural_estimates,
    parameter_estimates = parameter_estimates
  )
}

iterated_coordinates_to_natural_ipcr <- function(parameter_values, spec) {
  parameter_values <- as.numeric(parameter_values)
  if (length(parameter_values) != length(spec$natural_names)) {
    stop(
      "Unable to align iterated parameter coordinates with the RAM parameters.",
      call. = FALSE
    )
  }
  natural_values <- parameter_values
  natural_values[spec$transformed] <- exp(
    parameter_values[spec$transformed]
  )
  names(natural_values) <- spec$natural_names
  if (any(!is.finite(natural_values))) {
    stop(
      "Back-transformed RAM parameter predictions are not finite.",
      call. = FALSE
    )
  }
  natural_values
}

transform_standard_result_ipcr <- function(
    standard_result, spec, analysis_predictors, complete_rows) {
  if (identical(spec$type, "original")) {
    return(standard_result)
  }

  n <- NROW(standard_result$IPCs)
  q <- length(spec$natural_names)
  pooled_natural <- matrix(
    spec$natural_estimates,
    nrow = n,
    ncol = q,
    byrow = TRUE
  )
  influence <- standard_result$IPCs - pooled_natural
  influence[, spec$transformed] <- sweep(
    influence[, spec$transformed, drop = FALSE],
    2L,
    spec$natural_estimates[spec$transformed],
    FUN = "/"
  )
  parameter_IPCs <- matrix(
    spec$parameter_estimates,
    nrow = n,
    ncol = q,
    byrow = TRUE
  ) + influence
  colnames(parameter_IPCs) <- spec$parameter_names

  parameter_scores <- standard_result$scores
  parameter_scores[, spec$transformed] <- sweep(
    parameter_scores[, spec$transformed, drop = FALSE],
    2L,
    spec$natural_estimates[spec$transformed],
    FUN = "*"
  )
  colnames(parameter_scores) <- spec$parameter_names

  standard_result$IPCs <- parameter_IPCs
  standard_result$scores <- parameter_scores
  standard_result$mlm <- fit_second_stage_ipcr(
    parameter_IPCs[complete_rows, , drop = FALSE],
    analysis_predictors,
    require_full_rank = TRUE
  )
  standard_result$info$parameters <- spec$parameter_names
  standard_result$info$natural_parameters <- spec$natural_names
  standard_result$info$log_variance_parameters <-
    spec$natural_names[spec$transformed]
  standard_result$info$variance_parameterization <- spec$type
  standard_result
}

backtransform_iterated_IPCs_ipcr <- function(IPCs, spec) {
  natural_IPCs <- as.matrix(IPCs)
  natural_IPCs[, spec$transformed] <- exp(
    natural_IPCs[, spec$transformed, drop = FALSE]
  )
  if (any(!is.finite(natural_IPCs))) {
    warning(
      "Some back-transformed individual parameter contributions are not finite.",
      call. = FALSE
    )
  }
  colnames(natural_IPCs) <- spec$natural_names
  natural_IPCs
}

iterated_coefficient_matrix_ipcr <- function(mlm, parameter_names) {
  coefficients <- stats::coef(mlm)
  if (is.null(dim(coefficients))) {
    coefficients <- matrix(
      coefficients,
      ncol = 1L,
      dimnames = list(names(coefficients), parameter_names)
    )
  } else {
    coefficients <- as.matrix(coefficients)
    colnames(coefficients) <- parameter_names
  }
  coefficients
}

standard_ipcr_coefficient_standard_errors <- function(
    mlm, parameter_names) {
  coefficients <- iterated_coefficient_matrix_ipcr(mlm, parameter_names)
  model_summaries <- if (inherits(mlm, "mlm")) {
    summary(mlm)
  } else {
    list(summary(mlm))
  }
  if (length(model_summaries) != length(parameter_names)) {
    stop(
      "Unable to align standard-IPCR standard errors with model parameters.",
      call. = FALSE
    )
  }

  standard_errors <- matrix(
    NA_real_,
    nrow = NROW(coefficients),
    ncol = NCOL(coefficients),
    dimnames = dimnames(coefficients)
  )
  for (parameter_index in seq_along(parameter_names)) {
    coefficient_table <- model_summaries[[parameter_index]]$coefficients
    if (!all(rownames(coefficients) %in% rownames(coefficient_table))) {
      stop(
        "Unable to align standard-IPCR standard errors with second-stage coefficients.",
        call. = FALSE
      )
    }
    standard_errors[, parameter_index] <- coefficient_table[
      rownames(coefficients),
      "Std. Error"
    ]
  }

  invalid <- !is.finite(standard_errors) | standard_errors <= 0
  if (any(invalid)) {
    indices <- which(invalid, arr.ind = TRUE)
    labels <- paste0(
      rownames(standard_errors)[indices[, "row"]],
      " / ",
      colnames(standard_errors)[indices[, "col"]]
    )
    stop(
      "The standard-error-scaled convergence rule requires positive finite ",
      "standard-IPCR standard errors. Invalid coefficients: ",
      paste(labels, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  standard_errors
}

iterated_coefficient_change_ipcr <- function(
    current, updated, standard_errors, learning_rate) {
  if (!identical(dim(current), dim(updated)) ||
      !identical(dim(current), dim(standard_errors))) {
    stop(
      "Unable to align coefficient changes with their convergence scales.",
      call. = FALSE
    )
  }

  max(abs(updated - current) / (learning_rate * standard_errors))
}

log_likelihood_change_per_observation_ipcr <- function(
    current, previous, n) {
  if (!is.finite(previous)) {
    return(NA_real_)
  }
  abs(current - previous) / n
}

likelihood_stability_ipcr <- function(
    C_l, learning_rate, threshold = 1e-3) {
  if (!is.numeric(C_l) || length(C_l) != 1L || !is.finite(C_l) || C_l < 0) {
    stop("'C_l' must be a single nonnegative finite number.", call. = FALSE)
  }
  if (!is.numeric(learning_rate) || length(learning_rate) != 1L ||
      !is.finite(learning_rate) || learning_rate <= 0) {
    stop(
      "'learning_rate' must be a single positive finite number.",
      call. = FALSE
    )
  }
  if (!is.numeric(threshold) || length(threshold) != 1L ||
      !is.finite(threshold) || threshold <= 0) {
    stop(
      "'threshold' must be a single positive finite number.",
      call. = FALSE
    )
  }

  C_l_rate_adjusted <- C_l / learning_rate
  list(
    C_l_rate_adjusted = C_l_rate_adjusted,
    threshold = threshold,
    stable = C_l_rate_adjusted <= threshold
  )
}

likelihood_stability_warning_ipcr <- function(
    C_l_rate_adjusted, threshold) {
  paste0(
    "Iterated IPCR met the C_SE convergence criterion, but its ",
    "learning-rate-adjusted likelihood change (C_l / learning_rate = ",
    format(C_l_rate_adjusted, digits = 4L),
    ") exceeds the diagnostic threshold of ",
    format(threshold, digits = 4L),
    ". The returned solution may be numerically unstable. Inspect the ",
    "convergence diagnostics and consider rerunning with ",
    "iteration_info = TRUE."
  )
}

iterated_parameter_predictions_ipcr <- function(mlm, parameter_names) {
  predictions <- stats::fitted(mlm)
  if (is.null(dim(predictions))) {
    predictions <- matrix(predictions, ncol = 1L)
  } else {
    predictions <- as.matrix(predictions)
  }
  if (NCOL(predictions) != length(parameter_names)) {
    stop(
      "Unable to align iterated IPC predictions with the OpenMx parameters.",
      call. = FALSE
    )
  }
  colnames(predictions) <- parameter_names
  predictions
}

natural_scale_effects_ipcr <- function(mlm, spec) {
  if (!any(spec$transformed)) {
    return(NULL)
  }
  coefficients <- iterated_coefficient_matrix_ipcr(
    mlm,
    spec$parameter_names
  )
  if (!"(Intercept)" %in% rownames(coefficients)) {
    stop(
      "Natural-scale variance effects require a second-stage intercept.",
      call. = FALSE
    )
  }
  predictor_names <- setdiff(rownames(coefficients), "(Intercept)")
  rows <- lapply(which(spec$transformed), function(parameter_index) {
    log_intercept <- coefficients["(Intercept)", parameter_index]
    reference_variance <- exp(log_intercept)
    do.call(rbind, lapply(predictor_names, function(predictor_name) {
      log_effect <- coefficients[predictor_name, parameter_index]
      variance_ratio <- exp(log_effect)
      data.frame(
        Parameter = spec$natural_names[parameter_index],
        Log_Parameter = spec$parameter_names[parameter_index],
        Predictor = predictor_name,
        Log_Effect = unname(log_effect),
        Variance_Ratio = unname(variance_ratio),
        Reference_Variance = unname(reference_variance),
        Variance_Change = unname(
          reference_variance * (variance_ratio - 1)
        ),
        row.names = NULL
      )
    }))
  })
  effects <- do.call(rbind, rows)
  rownames(effects) <- NULL
  effects
}

initial_iterated_attempt_ipcr <- function(
    standard_result, pooled_parameters, analysis_predictors, complete_rows,
    learning_rate) {
  analysis_pooled <- pooled_parameters[complete_rows, , drop = FALSE]
  analysis_standard_IPCs <- standard_result$IPCs[
    complete_rows,
    ,
    drop = FALSE
  ]
  analysis_IPCs <- analysis_pooled + learning_rate *
    (analysis_standard_IPCs - analysis_pooled)
  mlm <- if (identical(learning_rate, 1)) {
    standard_result$mlm
  } else {
    fit_second_stage_ipcr(
      analysis_IPCs,
      analysis_predictors,
      require_full_rank = TRUE
    )
  }
  list(mlm = mlm)
}

abort_ipcr_it_nonconvergence <- function(
    attempts, iteration_history, iteration_info,
    suggest_log_variance = FALSE) {
  message <- paste0(
    "Iterated IPCR did not converge for any supplied learning rate. ",
    "No result was returned because inference from a nonconverged ",
    "iteration is unreliable. Inspect the 'attempts' field of this ",
    "ipcr_it_nonconvergence condition for diagnostics."
  )
  if (isTRUE(iteration_info)) {
    message <- paste0(
      message,
      " The 'iteration_history' field contains the complete iteration history."
    )
  }
  if (isTRUE(suggest_log_variance)) {
    message <- paste0(
      message,
      " If predicted variance parameters cross their bounds, consider ",
      "variance_parameterization = \"log\"."
    )
  }

  condition_fields <- list(
    message = message,
    call = NULL,
    attempts = attempts
  )
  if (isTRUE(iteration_info)) {
    condition_fields$iteration_history <- iteration_history
  }
  stop(structure(
    condition_fields,
    class = c("ipcr_it_nonconvergence", "error", "condition")
  ))
}

validate_iterated_components_ipcr <- function(
    components, n_rows, n_parameters, likelihood_only = FALSE) {
  if (!is.list(components) ||
      length(components$log_likelihood) != 1L ||
      !is.numeric(components$log_likelihood) ||
      !is.finite(components$log_likelihood)) {
    stop(
      "The iterated-IPCR backend returned an invalid log-likelihood.",
      call. = FALSE
    )
  }
  if (isTRUE(likelihood_only)) {
    return(components)
  }

  expected_dimension <- c(n_rows, n_parameters)
  for (component_name in c("IPCs", "scores")) {
    component <- components[[component_name]]
    if (!is.matrix(component) || !is.numeric(component) ||
        !identical(dim(component), expected_dimension) ||
        any(!is.finite(component))) {
      stop(
        "The iterated-IPCR backend returned an invalid ",
        component_name,
        " matrix.",
        call. = FALSE
      )
    }
  }
  components
}

run_iterated_ipcr_backend <- function(
    fit, predictors, analytic, conv, learning_rate, max_it, iteration_info,
    component_evaluator, parameterization = NULL,
    model_name = deparse(substitute(fit))) {
  standard_result <- ipcr(
    x = fit,
    predictors = predictors,
    analytic = analytic
  )
  standard_result$info$name <- model_name
  if (is.null(parameterization)) {
    parameterization <- ram_variance_parameterization_ipcr(
      fit,
      variance_parameterization = "original"
    )
  }
  parameter_estimates <- parameterization$parameter_estimates
  parameter_names <- parameterization$parameter_names
  n <- nobs_ipcr(fit)
  q <- length(parameter_estimates)
  complete_rows <- standard_result$complete_rows
  analysis_rows <- which(complete_rows)
  analysis_predictors <- droplevels(
    standard_result$predictors[complete_rows, , drop = FALSE]
  )
  standard_result <- transform_standard_result_ipcr(
    standard_result = standard_result,
    spec = parameterization,
    analysis_predictors = analysis_predictors,
    complete_rows = complete_rows
  )
  standard_coefficient_standard_errors <-
    standard_ipcr_coefficient_standard_errors(
      standard_result$mlm,
      parameter_names
    )
  pooled_parameters <- matrix(
    parameter_estimates,
    nrow = n,
    ncol = q,
    byrow = TRUE,
    dimnames = list(NULL, parameter_names)
  )

  attempt_summaries <- data.frame(
    learning_rate = numeric(),
    iterations = integer(),
    C_SE = numeric(),
    C_l = numeric(),
    status = character(),
    reason = character(),
    stringsAsFactors = FALSE
  )
  iteration_history <- data.frame(
    learning_rate = numeric(),
    iteration = integer(),
    C_SE = numeric(),
    C_l = numeric(),
    log_likelihood = numeric(),
    status = character(),
    stringsAsFactors = FALSE
  )
  selected_result <- NULL

  for (alpha in learning_rate) {
    attempt <- tryCatch(
      initial_iterated_attempt_ipcr(
        standard_result = standard_result,
        pooled_parameters = pooled_parameters,
        analysis_predictors = analysis_predictors,
        complete_rows = complete_rows,
        learning_rate = alpha
      ),
      error = function(e) e
    )
    if (inherits(attempt, "error")) {
      attempt_summaries <- rbind(
        attempt_summaries,
        data.frame(
          learning_rate = alpha,
          iterations = 0L,
          C_SE = NA_real_,
          C_l = NA_real_,
          status = "initialization failure",
          reason = conditionMessage(attempt),
          stringsAsFactors = FALSE
        )
      )
      if (isTRUE(iteration_info)) {
        iteration_history <- rbind(
          iteration_history,
          data.frame(
            learning_rate = alpha,
            iteration = 0L,
            C_SE = NA_real_,
            C_l = NA_real_,
            log_likelihood = NA_real_,
            status = "initialization failure",
            stringsAsFactors = FALSE
          )
        )
      }
      next
    }
    current_mlm <- attempt$mlm
    attempt_status <- "maximum iterations reached"
    attempt_reason <- ""
    completed_iterations <- 0L
    previous_log_likelihood <- NA_real_
    previous_history_row <- NA_integer_
    attempt_C_SE <- NA_real_
    attempt_C_l <- NA_real_

    for (iteration in seq_len(as.integer(max_it))) {
      local_components <- tryCatch(
        {
          current_parameters <- iterated_parameter_predictions_ipcr(
            current_mlm,
            parameter_names
          )
          components <- component_evaluator(
            rows = analysis_rows,
            parameter_values = current_parameters,
            learning_rate = alpha,
            likelihood_only = FALSE
          )
          validate_iterated_components_ipcr(
            components,
            n_rows = length(analysis_rows),
            n_parameters = q
          )
        },
        error = function(e) e
      )
      if (inherits(local_components, "error")) {
        attempt_status <- "numerical failure"
        attempt_reason <- conditionMessage(local_components)
        if (isTRUE(iteration_info)) {
          iteration_history <- rbind(
            iteration_history,
            data.frame(
              learning_rate = alpha,
              iteration = iteration,
              C_SE = NA_real_,
              C_l = NA_real_,
              log_likelihood = NA_real_,
              status = attempt_status,
              stringsAsFactors = FALSE
            )
          )
        }
        break
      }

      current_log_likelihood <- local_components$log_likelihood
      current_C_l <- log_likelihood_change_per_observation_ipcr(
        current = current_log_likelihood,
        previous = previous_log_likelihood,
        n = length(analysis_rows)
      )
      if (is.finite(current_C_l)) {
        attempt_C_l <- current_C_l
      }
      if (isTRUE(iteration_info) && !is.na(previous_history_row)) {
        iteration_history$C_l[previous_history_row] <- current_C_l
        iteration_history$log_likelihood[previous_history_row] <-
          current_log_likelihood
      }
      updated_IPCs <- standard_result$IPCs
      updated_IPCs[complete_rows, ] <- local_components$IPCs
      updated_scores <- standard_result$scores
      updated_scores[complete_rows, ] <- local_components$scores
      next_mlm <- tryCatch(
        fit_second_stage_ipcr(
          local_components$IPCs,
          analysis_predictors,
          require_full_rank = TRUE
        ),
        error = function(e) e
      )
      if (inherits(next_mlm, "error")) {
        attempt_status <- "regression failure"
        attempt_reason <- conditionMessage(next_mlm)
        if (isTRUE(iteration_info)) {
          iteration_history <- rbind(
            iteration_history,
            data.frame(
              learning_rate = alpha,
              iteration = iteration,
              C_SE = NA_real_,
              C_l = NA_real_,
              log_likelihood = current_log_likelihood,
              status = attempt_status,
              stringsAsFactors = FALSE
            )
          )
        }
        break
      }

      coefficient_change <- tryCatch(
        {
          current_coefficients <- iterated_coefficient_matrix_ipcr(
            current_mlm,
            parameter_names
          )
          next_coefficients <- iterated_coefficient_matrix_ipcr(
            next_mlm,
            parameter_names
          )
          iterated_coefficient_change_ipcr(
            current = current_coefficients,
            updated = next_coefficients,
            standard_errors = standard_coefficient_standard_errors,
            learning_rate = alpha
          )
        },
        error = function(e) e
      )
      if (inherits(coefficient_change, "error")) {
        attempt_status <- "regression failure"
        attempt_reason <- conditionMessage(coefficient_change)
        if (isTRUE(iteration_info)) {
          iteration_history <- rbind(
            iteration_history,
            data.frame(
              learning_rate = alpha,
              iteration = iteration,
              C_SE = NA_real_,
              C_l = NA_real_,
              log_likelihood = current_log_likelihood,
              status = attempt_status,
              stringsAsFactors = FALSE
            )
          )
        }
        break
      }
      completed_iterations <- iteration
      attempt_C_SE <- coefficient_change
      iteration_status <- if (
        coefficient_change <= conv
      ) "converged" else "running"
      if (isTRUE(iteration_info)) {
        iteration_history <- rbind(
          iteration_history,
          data.frame(
            learning_rate = alpha,
            iteration = iteration,
            C_SE = coefficient_change,
            C_l = NA_real_,
            log_likelihood = current_log_likelihood,
            status = iteration_status,
            stringsAsFactors = FALSE
          )
        )
      }

      if (coefficient_change <= conv) {
        final_validation <- tryCatch(
          {
            final_parameters <- iterated_parameter_predictions_ipcr(
              next_mlm,
              parameter_names
            )
            components <- component_evaluator(
              rows = analysis_rows,
              parameter_values = final_parameters,
              learning_rate = alpha,
              likelihood_only = TRUE
            )
            validate_iterated_components_ipcr(
              components,
              n_rows = length(analysis_rows),
              n_parameters = q,
              likelihood_only = TRUE
            )
          },
          error = function(e) e
        )
        if (inherits(final_validation, "error")) {
          attempt_status <- "numerical failure"
          attempt_reason <- conditionMessage(final_validation)
          if (isTRUE(iteration_info)) {
            iteration_history$status[NROW(iteration_history)] <- attempt_status
          }
          break
        }

        final_C_l <- log_likelihood_change_per_observation_ipcr(
          current = final_validation$log_likelihood,
          previous = current_log_likelihood,
          n = length(analysis_rows)
        )
        if (isTRUE(iteration_info)) {
          iteration_history$C_l[NROW(iteration_history)] <- final_C_l
          iteration_history$log_likelihood[NROW(iteration_history)] <-
            final_validation$log_likelihood
        }
        attempt_C_l <- final_C_l

        attempt_status <- "converged"
        selected_result <- list(
          mlm = next_mlm,
          IPCs = updated_IPCs,
          scores = updated_scores,
          learning_rate = alpha,
          iterations = iteration,
          log_likelihood = final_validation$log_likelihood,
          C_SE = coefficient_change,
          C_l = final_C_l
        )
        break
      }

      if (iteration == as.integer(max_it)) {
        final_validation <- tryCatch(
          {
            final_parameters <- iterated_parameter_predictions_ipcr(
              next_mlm,
              parameter_names
            )
            components <- component_evaluator(
              rows = analysis_rows,
              parameter_values = final_parameters,
              learning_rate = alpha,
              likelihood_only = TRUE
            )
            validate_iterated_components_ipcr(
              components,
              n_rows = length(analysis_rows),
              n_parameters = q,
              likelihood_only = TRUE
            )
          },
          error = function(e) e
        )
        if (inherits(final_validation, "error")) {
          attempt_status <- "numerical failure"
          attempt_reason <- conditionMessage(final_validation)
          if (isTRUE(iteration_info)) {
            iteration_history$status[NROW(iteration_history)] <- attempt_status
          }
          break
        }

        final_C_l <- log_likelihood_change_per_observation_ipcr(
          current = final_validation$log_likelihood,
          previous = current_log_likelihood,
          n = length(analysis_rows)
        )
        if (isTRUE(iteration_info)) {
          iteration_history$C_l[NROW(iteration_history)] <- final_C_l
          iteration_history$log_likelihood[NROW(iteration_history)] <-
            final_validation$log_likelihood
          iteration_history$status[NROW(iteration_history)] <- attempt_status
        }
        attempt_C_l <- final_C_l

      }

      previous_log_likelihood <- current_log_likelihood
      if (isTRUE(iteration_info)) {
        previous_history_row <- NROW(iteration_history)
      }
      current_mlm <- next_mlm
    }

    attempt_summaries <- rbind(
      attempt_summaries,
      data.frame(
        learning_rate = alpha,
        iterations = completed_iterations,
        C_SE = attempt_C_SE,
        C_l = attempt_C_l,
        status = attempt_status,
        reason = attempt_reason,
        stringsAsFactors = FALSE
      )
    )
    if (identical(attempt_status, "converged")) {
      break
    }
  }

  if (is.null(selected_result)) {
    suggest_log_variance <-
      inherits(fit, "MxRAMModel") && isTRUE(analytic) &&
      identical(parameterization$type, "original")
    abort_ipcr_it_nonconvergence(
      attempts = attempt_summaries,
      iteration_history = iteration_history,
      iteration_info = iteration_info,
      suggest_log_variance = suggest_log_variance
    )
  }

  result <- standard_result
  result$info$ipcr_type <- "iterated"
  result$info$analytic <- analytic
  result$info$variance_parameterization <- parameterization$type
  result$info$natural_parameters <- parameterization$natural_names
  result$info$log_variance_parameters <-
    parameterization$natural_names[parameterization$transformed]
  result$info$converged <- TRUE
  result$info$learning_rate <- selected_result$learning_rate
  result$info$iterations <- selected_result$iterations
  result$IPCs <- selected_result$IPCs
  result$scores <- selected_result$scores
  result$mlm <- selected_result$mlm
  output_components <- second_stage_output_ipcr(
    mlm = result$mlm,
    parameter_names = parameter_names
  )
  result$mancova <- output_components$mancova
  result$output <- output_components$output
  result$parameterization <- list(
    type = parameterization$type,
    natural_parameters = parameterization$natural_names,
    reported_parameters = parameterization$parameter_names,
    log_variance_parameters =
      parameterization$natural_names[parameterization$transformed]
  )
  if (identical(parameterization$type, "log")) {
    result$natural_scale_IPCs <- backtransform_iterated_IPCs_ipcr(
      result$IPCs,
      parameterization
    )
    result$output$natural_scale_effects <- natural_scale_effects_ipcr(
      result$mlm,
      parameterization
    )
  } else {
    result$output$natural_scale_effects <- NULL
  }
  likelihood_stability <- likelihood_stability_ipcr(
    C_l = selected_result$C_l,
    learning_rate = selected_result$learning_rate
  )
  result$convergence <- list(
    converged = TRUE,
    status = paste0(
      "Converged after ", selected_result$iterations,
      " iterations with learning rate ", selected_result$learning_rate, "."
    ),
    criterion = conv,
    C_SE = selected_result$C_SE,
    C_l = selected_result$C_l,
    C_l_rate_adjusted = likelihood_stability$C_l_rate_adjusted,
    likelihood_warning_threshold = likelihood_stability$threshold,
    likelihood_stable = likelihood_stability$stable,
    standard_ipcr_standard_errors = standard_coefficient_standard_errors,
    max_iterations_per_rate = as.integer(max_it),
    learning_rate = selected_result$learning_rate,
    log_likelihood = selected_result$log_likelihood,
    attempts = attempt_summaries
  )
  if (isTRUE(iteration_info)) {
    result$iteration_history <- iteration_history
  }
  class(result) <- c("ipcr_it", "ipcr")
  if (!isTRUE(likelihood_stability$stable)) {
    warning(
      likelihood_stability_warning_ipcr(
        C_l_rate_adjusted = likelihood_stability$C_l_rate_adjusted,
        threshold = likelihood_stability$threshold
      ),
      call. = FALSE
    )
  }
  result
}

#' @export
ipcr_it.MxRAMModel <- function(
    x, predictors, analytic = TRUE, conv = 0.01,
    learning_rate = c(1, 0.5, 0.1),
    max_it = 350L, iteration_info = FALSE,
    variance_parameterization = c("original", "log")) {
  model_name <- deparse(substitute(x))
  variance_parameterization <- match.arg(variance_parameterization)
  validate_ipcr_it_MxRAMModel_arguments(
    x = x,
    analytic = analytic,
    conv = conv,
    learning_rate = learning_rate,
    max_it = max_it,
    iteration_info = iteration_info,
    variance_parameterization = variance_parameterization
  )

  parameterization <- ram_variance_parameterization_ipcr(
    x,
    variance_parameterization = variance_parameterization
  )
  analytical_context <- if (isTRUE(analytic)) {
    analytical_iterated_ram_context_ipcr(
      x,
      parameterization = parameterization
    )
  } else {
    NULL
  }
  component_evaluator <- if (isTRUE(analytic)) {
    function(
        rows, parameter_values, learning_rate, likelihood_only = FALSE) {
      analytical_iterated_ram_components_ipcr(
        x = x,
        rows = rows,
        parameter_values = parameter_values,
        learning_rate = learning_rate,
        context = analytical_context,
        likelihood_only = likelihood_only
      )
    }
  } else {
    function(
        rows, parameter_values, learning_rate, likelihood_only = FALSE) {
      numerical_iterated_ram_components_ipcr(
        x = x,
        rows = rows,
        parameter_values = parameter_values,
        learning_rate = learning_rate,
        likelihood_only = likelihood_only
      )
    }
  }

  run_iterated_ipcr_backend(
    fit = x,
    predictors = predictors,
    analytic = analytic,
    conv = conv,
    learning_rate = learning_rate,
    max_it = max_it,
    iteration_info = iteration_info,
    component_evaluator = component_evaluator,
    parameterization = parameterization,
    model_name = model_name
  )
}
