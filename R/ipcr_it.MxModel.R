validate_ipcr_it_MxModel_arguments <- function(
    x, analytic, conv, learning_rate, max_it, iteration_info) {
  if (!inherits(x, "MxModel")) {
    stop("'x' must inherit from OpenMx class 'MxModel'.", call. = FALSE)
  }
  if (inherits(x, "MxRAMModel")) {
    stop(
      "MxRAMModel objects must use ipcr_it.MxRAMModel().",
      call. = FALSE
    )
  }
  validate_ipcr_it_common_arguments(
    analytic = analytic,
    conv = conv,
    learning_rate = learning_rate,
    max_it = max_it,
    iteration_info = iteration_info
  )

  check_model_scope_ipcr(x)
  check_model_convergence_ipcr(x)

  if (length(x@constraints) > 0L) {
    stop(
      "Numerical iterated IPCR for general MxModels does not support explicit mxConstraint() objects.",
      call. = FALSE
    )
  }
  if (length(x@penalties) > 0L) {
    stop(
      "Numerical iterated IPCR for general MxModels does not support penalties.",
      call. = FALSE
    )
  }
  if (any(!is.na(x$data@weight)) || any(!is.na(x$data@frequency))) {
    stop(
      "Numerical iterated IPCR for general MxModels requires unweighted raw data.",
      call. = FALSE
    )
  }
  if (!inherits(
    x$fitfunction,
    c("MxFitFunctionML", "MxFitFunctionFIML")
  )) {
    stop(
      "Numerical iterated IPCR for general MxModels currently requires an ML fit function with a row-additive likelihood.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

iterated_openmx_parameter_groups_ipcr <- function(parameter_values) {
  key_columns <- lapply(
    seq_len(NCOL(parameter_values)),
    function(column) sprintf("%.17g", parameter_values[, column])
  )
  group_key <- do.call(paste, c(key_columns, list(sep = "\034")))
  split(seq_len(NROW(parameter_values)), match(group_key, unique(group_key)))
}

run_openmx_without_optimizer_ipcr <- function(x) {
  result <- NULL
  invisible(utils::capture.output(
    result <- suppressMessages(OpenMx::mxRun(
      x,
      useOptimizer = FALSE,
      silent = TRUE
    )),
    type = "output"
  ))
  result
}

subset_openmx_raw_model_ipcr <- function(x, rows) {
  observed <- openmx_raw_data_ipcr(x)
  subset_model <- x
  subset_model$data$observed <- observed[rows, , drop = FALSE]
  subset_model$data$numObs <- length(rows)
  subset_model
}

# Evaluate a damped IPC update for a general raw-data MxModel. At each distinct
# predicted parameter vector, imxRowGradients() supplies all casewise scores.
# Their outer product estimates the per-row information. Re-running the model
# on the corresponding data subset supplies the row-additive likelihood used
# to validate iterations and record likelihood diagnostics.
numerical_iterated_mxmodel_components_ipcr <- function(
    x, rows, parameter_values, learning_rate, likelihood_only = FALSE) {
  rows <- as.integer(rows)
  parameter_names <- names(OpenMx::omxGetParameters(x))
  q <- length(parameter_names)
  parameter_values <- as.matrix(parameter_values)

  if (!is.numeric(parameter_values) ||
      !identical(dim(parameter_values), c(length(rows), q))) {
    stop(
      "The predicted OpenMx parameter matrix has incompatible dimensions.",
      call. = FALSE
    )
  }
  if (any(!is.finite(parameter_values))) {
    stop("Predicted OpenMx parameters must be finite.", call. = FALSE)
  }
  colnames(parameter_values) <- parameter_names

  n <- nobs_ipcr(x)
  scores <- corrections <- NULL
  if (!isTRUE(likelihood_only)) {
    scores <- matrix(
      0,
      nrow = length(rows),
      ncol = q,
      dimnames = list(NULL, parameter_names)
    )
    corrections <- scores
  }
  log_likelihood <- 0
  groups <- iterated_openmx_parameter_groups_ipcr(parameter_values)

  for (group_positions in groups) {
    parameters <- parameter_values[group_positions[1L], ]
    names(parameters) <- parameter_names
    local_model <- OpenMx::omxSetParameters(
      model = x,
      labels = parameter_names,
      values = parameters
    )
    group_rows <- rows[group_positions]
    if (!isTRUE(likelihood_only)) {
      local_model <- run_openmx_without_optimizer_ipcr(local_model)
      local_scores <- as.matrix(numerical_scores_MxModel_ipcr(local_model))

      if (!identical(dim(local_scores), c(n, q))) {
        stop(
          "OpenMx numerical row gradients have incompatible dimensions at a predicted parameter vector.",
          call. = FALSE
        )
      }
      if (!is.null(colnames(local_scores)) &&
          all(parameter_names %in% colnames(local_scores))) {
        local_scores <- local_scores[, parameter_names, drop = FALSE]
      } else {
        colnames(local_scores) <- parameter_names
      }
      if (any(!is.finite(local_scores))) {
        stop(
          "OpenMx returned non-finite numerical row gradients at a predicted parameter vector.",
          call. = FALSE
        )
      }

      information <- crossprod(local_scores) / n
      if (qr(information)$rank < q || !is.finite(rcond(information))) {
        stop(
          "The outer-product information matrix is singular at a predicted OpenMx parameter vector.",
          call. = FALSE
        )
      }
      group_scores <- local_scores[group_rows, , drop = FALSE]
      group_corrections <- tryCatch(
        t(solve(information, t(group_scores))),
        error = function(e) {
          stop(
            "Unable to solve the outer-product information matrix: ",
            conditionMessage(e),
            call. = FALSE
          )
        }
      )
      scores[group_positions, ] <- group_scores
      corrections[group_positions, ] <- group_corrections
    }

    subset_model <- subset_openmx_raw_model_ipcr(local_model, group_rows)
    subset_model <- run_openmx_without_optimizer_ipcr(subset_model)
    subset_fit <- subset_model$output$fit
    if (length(subset_fit) != 1L || !is.finite(subset_fit)) {
      stop(
        "OpenMx returned a non-finite likelihood for a raw-data subset at a predicted parameter vector.",
        call. = FALSE
      )
    }
    log_likelihood <- log_likelihood - 0.5 * subset_fit
  }

  if (isTRUE(likelihood_only)) {
    return(list(log_likelihood = as.numeric(log_likelihood)))
  }

  IPCs <- parameter_values + learning_rate * corrections
  colnames(IPCs) <- parameter_names
  list(
    IPCs = IPCs,
    scores = scores,
    log_likelihood = as.numeric(log_likelihood)
  )
}

#' @export
ipcr_it.MxModel <- function(
    x, predictors, analytic = FALSE, conv = 0.01,
    learning_rate = c(1, 0.5, 0.1),
    max_it = 350L, iteration_info = FALSE,
    variance_parameterization = "original") {
  model_name <- deparse(substitute(x))
  validate_original_variance_parameterization_ipcr(
    variance_parameterization,
    "general MxModel iteration"
  )
  validate_ipcr_it_MxModel_arguments(
    x = x,
    analytic = analytic,
    conv = conv,
    learning_rate = learning_rate,
    max_it = max_it,
    iteration_info = iteration_info
  )

  warning(
    if (isTRUE(analytic)) {
      paste0(
        "Analytical scores are unavailable for general MxModels. ",
        "Iterated IPCR will use repeated numerical OpenMx row gradients, "
      )
    } else {
      paste0(
        "Iterated IPCR for a general MxModel uses repeated numerical ",
        "OpenMx row gradients, "
      )
    },
    "which can lead to a very long run time. Use an MxRAMModel with ",
    "analytic = TRUE when possible.",
    call. = FALSE
  )
  analytic <- FALSE

  run_iterated_ipcr_backend(
    fit = x,
    predictors = predictors,
    analytic = analytic,
    conv = conv,
    learning_rate = learning_rate,
    max_it = max_it,
    iteration_info = iteration_info,
    component_evaluator = function(
        rows, parameter_values, learning_rate, likelihood_only = FALSE) {
      numerical_iterated_mxmodel_components_ipcr(
        x = x,
        rows = rows,
        parameter_values = parameter_values,
        learning_rate = learning_rate,
        likelihood_only = likelihood_only
      )
    },
    model_name = model_name
  )
}
