validate_ipcr_it_lavaan_arguments <- function(
    x, analytic, conv, learning_rate, max_it, iteration_info) {
  if (!inherits(x, "lavaan")) {
    stop("'x' must inherit from class 'lavaan'.", call. = FALSE)
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

  if (lavaan::lavInspect(x, what = "ngroups") != 1L) {
    stop(
      "Iterated IPCR currently supports only single-group lavaan models.",
      call. = FALSE
    )
  }
  if (!identical(toupper(x@Options$estimator), "ML") ||
      isTRUE(x@Model@categorical)) {
    stop(
      "Iterated IPCR currently supports only continuous lavaan models fitted by ML.",
      call. = FALSE
    )
  }
  if (!identical(tolower(x@Options$likelihood), "normal")) {
    stop(
      "Iterated IPCR currently requires lavaan's normal-theory likelihood.",
      call. = FALSE
    )
  }
  if (isTRUE(x@Options$conditional.x)) {
    stop(
      "Iterated IPCR does not currently support conditional.x lavaan models.",
      call. = FALSE
    )
  }
  if (isTRUE(x@Options$.sampling.weights)) {
    stop(
      "Iterated IPCR does not currently support sampling-weighted lavaan models.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

lavaan_ram_spec_ipcr <- function(x) {
  parameter_table <- lavaan::parameterTable(x)
  observed_variables <- lavaan::lavNames(x, type = "ov")
  latent_variables <- lavaan::lavNames(x, type = "lv")
  all_variables <- c(observed_variables, latent_variables)
  variable_count <- length(all_variables)
  observed_count <- length(observed_variables)
  parameter_estimates <- coef_ipcr(x)
  parameter_names <- names(parameter_estimates)
  q <- length(parameter_estimates)

  free_rows <- parameter_table$free > 0L
  supported_operations <- c("=~", "~", "~~", "~1")
  unsupported_rows <- free_rows &
    !parameter_table$op %in% supported_operations
  if (any(unsupported_rows)) {
    stop(
      "Iterated IPCR cannot translate free lavaan parameters with operators: ",
      paste(unique(parameter_table$op[unsupported_rows]), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  full_parameter_count <- x@Model@nx.free
  reduction <- lavaan_constraint_reduction_ipcr(x@Model)
  reduction_matrix <- if (is.null(reduction)) {
    diag(full_parameter_count)
  } else {
    reduction$K
  }
  if (!identical(dim(reduction_matrix), c(full_parameter_count, q))) {
    stop(
      "Unable to align the lavaan free parameters with their reduced coordinates.",
      call. = FALSE
    )
  }

  A <- matrix(
    0,
    nrow = variable_count,
    ncol = variable_count,
    dimnames = list(all_variables, all_variables)
  )
  S <- A
  m <- matrix(
    0,
    nrow = variable_count,
    ncol = 1L,
    dimnames = list(all_variables, NULL)
  )
  F <- matrix(
    0,
    nrow = observed_count,
    ncol = variable_count,
    dimnames = list(observed_variables, all_variables)
  )
  F[cbind(seq_len(observed_count), match(observed_variables, all_variables))] <- 1

  empty_derivative <- function() {
    list(
      A = matrix(0, variable_count, variable_count),
      S = matrix(0, variable_count, variable_count),
      m = matrix(0, variable_count, 1L)
    )
  }
  derivatives <- lapply(seq_len(q), function(index) empty_derivative())

  add_derivative <- function(matrix_name, row, column, free_index) {
    coefficients <- reduction_matrix[free_index, ]
    nonzero <- which(coefficients != 0)
    for (parameter_index in nonzero) {
      derivatives[[parameter_index]][[matrix_name]][row, column] <<-
        derivatives[[parameter_index]][[matrix_name]][row, column] +
        coefficients[parameter_index]
    }
  }

  structural_rows <- parameter_table$op %in% supported_operations
  for (table_row in which(structural_rows)) {
    row <- parameter_table[table_row, ]
    value <- row$est
    if (!is.finite(value)) {
      stop(
        "The fitted lavaan model contains a non-finite RAM parameter.",
        call. = FALSE
      )
    }

    matrix_name <- NULL
    matrix_row <- matrix_column <- NA_integer_
    if (identical(row$op, "=~") &&
        row$lhs %in% all_variables && row$rhs %in% all_variables) {
      matrix_name <- "A"
      matrix_row <- match(row$rhs, all_variables)
      matrix_column <- match(row$lhs, all_variables)
    } else if (identical(row$op, "~") &&
               row$lhs %in% all_variables && row$rhs %in% all_variables) {
      matrix_name <- "A"
      matrix_row <- match(row$lhs, all_variables)
      matrix_column <- match(row$rhs, all_variables)
    } else if (identical(row$op, "~~") &&
               row$lhs %in% all_variables && row$rhs %in% all_variables) {
      matrix_name <- "S"
      matrix_row <- match(row$lhs, all_variables)
      matrix_column <- match(row$rhs, all_variables)
    } else if (identical(row$op, "~1") && row$lhs %in% all_variables) {
      matrix_name <- "m"
      matrix_row <- match(row$lhs, all_variables)
      matrix_column <- 1L
    }

    if (is.null(matrix_name)) {
      if (row$free > 0L) {
        stop(
          "Unable to place the free lavaan parameter '",
          row$lhs, row$op, row$rhs,
          "' in RAM matrices.",
          call. = FALSE
        )
      }
      next
    }

    if (identical(matrix_name, "A")) {
      A[matrix_row, matrix_column] <- value
    } else if (identical(matrix_name, "S")) {
      S[matrix_row, matrix_column] <- value
      S[matrix_column, matrix_row] <- value
    } else {
      m[matrix_row, matrix_column] <- value
    }

    if (row$free > 0L) {
      add_derivative(
        matrix_name = matrix_name,
        row = matrix_row,
        column = matrix_column,
        free_index = row$free
      )
      if (identical(matrix_name, "S") && matrix_row != matrix_column) {
        add_derivative(
          matrix_name = matrix_name,
          row = matrix_column,
          column = matrix_row,
          free_index = row$free
        )
      }
    }
  }

  data <- as.matrix(lavaan::lavInspect(x, what = "data"))
  storage.mode(data) <- "double"
  data <- data[, observed_variables, drop = FALSE]
  missingness <- apply(is.na(data), 1L, paste0, collapse = "")
  missingness_groups <- split(
    seq_len(NROW(data)),
    missingness,
    drop = TRUE
  )
  information_patterns <- lapply(missingness_groups, function(positions) {
    list(
      observed = which(!is.na(data[positions[1L], ])),
      proportion = length(positions) / NROW(data)
    )
  })

  list(
    matrices = list(A = A, S = S, m = m, F = F),
    derivatives = derivatives,
    data = data,
    parameter_estimates = parameter_estimates,
    parameter_names = parameter_names,
    meanstructure = isTRUE(x@Options$meanstructure),
    fixed_means = colMeans(data, na.rm = TRUE),
    missingness = missingness,
    information_patterns = information_patterns
  )
}

lavaan_ram_matrices_at_ipcr <- function(spec, parameters) {
  parameter_change <- parameters - spec$parameter_estimates
  matrices <- spec$matrices

  for (parameter_index in seq_along(parameter_change)) {
    change <- parameter_change[parameter_index]
    if (change == 0) {
      next
    }
    matrices$A <- matrices$A +
      change * spec$derivatives[[parameter_index]]$A
    matrices$S <- matrices$S +
      change * spec$derivatives[[parameter_index]]$S
    matrices$m <- matrices$m +
      change * spec$derivatives[[parameter_index]]$m
  }

  matrices
}

iterated_lavaan_groups_ipcr <- function(
    data, rows, parameter_values, missingness = NULL) {
  if (is.null(missingness)) {
    missingness <- apply(is.na(data), 1L, paste0, collapse = "")
  }
  grouping_data <- cbind(
    data.frame(
      missingness = missingness[rows]
    ),
    as.data.frame(parameter_values, optional = TRUE)
  )
  key_columns <- lapply(grouping_data, function(column) {
    if (is.numeric(column)) {
      sprintf("%.17g", column)
    } else {
      as.character(column)
    }
  })
  group_key <- do.call(paste, c(key_columns, list(sep = "\034")))
  group <- match(group_key, unique(group_key))
  split(seq_along(rows), group)
}

lavaan_ram_pattern_components_ipcr <- function(
    spec, matrices, moments, observed, compute_information = TRUE) {
  covariance <- moments$covariance[observed, observed, drop = FALSE]
  covariance <- 0.5 * (covariance + t(covariance))
  covariance_chol <- chol(covariance)
  covariance_inverse <- chol2inv(covariance_chol)
  observed_count <- length(observed)

  if (!isTRUE(compute_information)) {
    return(list(
      covariance = covariance,
      covariance_chol = covariance_chol,
      covariance_inverse = covariance_inverse,
      observed_count = observed_count
    ))
  }

  covariance_count <- observed_count * (observed_count + 1L) / 2L
  duplication <- lavaan::lav_matrix_duplication(observed_count)
  covariance_weight <- 0.5 * t(duplication) %*%
    kronecker(covariance_inverse, covariance_inverse) %*%
    duplication
  full_jacobian <- ram_jacobian_ipcr(
    matrices = matrices,
    moments = moments,
    derivatives = spec$derivatives,
    observed = observed
  )

  if (isTRUE(spec$meanstructure)) {
    moment_count <- covariance_count + observed_count
    weight <- matrix(0, nrow = moment_count, ncol = moment_count)
    covariance_rows <- seq_len(covariance_count)
    mean_rows <- covariance_count + seq_len(observed_count)
    weight[covariance_rows, covariance_rows] <- covariance_weight
    weight[mean_rows, mean_rows] <- covariance_inverse
    jacobian <- full_jacobian
  } else {
    moment_count <- covariance_count
    weight <- covariance_weight
    jacobian <- full_jacobian[seq_len(covariance_count), , drop = FALSE]
  }

  list(
    covariance = covariance,
    covariance_chol = covariance_chol,
    covariance_inverse = covariance_inverse,
    observed_count = observed_count,
    covariance_count = covariance_count,
    moment_count = moment_count,
    weight = weight,
    jacobian = jacobian,
    information = t(jacobian) %*% weight %*% jacobian
  )
}

analytical_iterated_lavaan_components_ipcr <- function(
    spec, rows, parameter_values, learning_rate, likelihood_only = FALSE) {
  rows <- as.integer(rows)
  parameter_values <- as.matrix(parameter_values)
  parameter_names <- spec$parameter_names
  q <- length(parameter_names)
  if (!is.numeric(parameter_values) ||
      !identical(dim(parameter_values), c(length(rows), q))) {
    stop(
      "The predicted lavaan parameter matrix has incompatible dimensions.",
      call. = FALSE
    )
  }
  if (any(!is.finite(parameter_values))) {
    stop("Predicted lavaan parameters must be finite.", call. = FALSE)
  }
  colnames(parameter_values) <- parameter_names

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
  groups <- iterated_lavaan_groups_ipcr(
    data = spec$data,
    rows = rows,
    parameter_values = parameter_values,
    missingness = spec$missingness
  )

  for (group_positions in groups) {
    representative_position <- group_positions[1L]
    representative_row <- rows[representative_position]
    parameters <- parameter_values[representative_position, ]
    names(parameters) <- parameter_names
    observed <- which(!is.na(spec$data[representative_row, ]))
    if (length(observed) == 0L) {
      stop(
        "Iterated IPCR cannot update a row with no observed variables.",
        call. = FALSE
      )
    }

    matrices <- lavaan_ram_matrices_at_ipcr(spec, parameters)
    moments <- ram_moments_ipcr(matrices)
    pattern_components <- lavaan_ram_pattern_components_ipcr(
      spec = spec,
      matrices = matrices,
      moments = moments,
      observed = observed,
      compute_information = !isTRUE(likelihood_only)
    )
    covariance <- pattern_components$covariance
    covariance_chol <- pattern_components$covariance_chol
    covariance_inverse <- pattern_components$covariance_inverse
    means <- if (isTRUE(spec$meanstructure)) {
      moments$means[observed]
    } else {
      spec$fixed_means[observed]
    }
    observed_count <- pattern_components$observed_count

    original_rows <- rows[group_positions]
    centered <- sweep(
      spec$data[original_rows, observed, drop = FALSE],
      2L,
      means,
      FUN = "-"
    )
    if (!isTRUE(likelihood_only)) {
      moment_count <- pattern_components$moment_count
      weight <- pattern_components$weight
      jacobian <- pattern_components$jacobian
      information <- matrix(0, nrow = q, ncol = q)
      for (information_pattern in spec$information_patterns) {
        expected_components <- lavaan_ram_pattern_components_ipcr(
          spec = spec,
          matrices = matrices,
          moments = moments,
          observed = information_pattern$observed
        )
        information <- information + information_pattern$proportion *
          expected_components$information
      }
      information <- 0.5 * (information + t(information))
      information_chol <- tryCatch(chol(information), error = function(e) e)
      if (inherits(information_chol, "error")) {
        stop(
          "The lavaan information matrix is singular at the predicted parameters.",
          call. = FALSE
        )
      }
      deviations <- t(vapply(
        seq_len(NROW(centered)),
        function(position) {
          centered_row <- centered[position, ]
          covariance_deviation <- lavaan::lav_matrix_vech(
            tcrossprod(centered_row) - covariance
          )
          if (isTRUE(spec$meanstructure)) {
            c(covariance_deviation, centered_row)
          } else {
            covariance_deviation
          }
        },
        numeric(moment_count)
      ))
      group_scores <- deviations %*% weight %*% jacobian
      scores[group_positions, ] <- group_scores
      corrections[group_positions, ] <- t(solve(information, t(group_scores)))
    }

    quadratic <- rowSums((centered %*% covariance_inverse) * centered)
    log_determinant <- 2 * sum(log(diag(covariance_chol)))
    log_likelihood <- log_likelihood - 0.5 * sum(
      observed_count * log(2 * pi) + log_determinant + quadratic
    )
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
ipcr_it.lavaan <- function(
    x, predictors, analytic = TRUE, conv = 0.01,
    learning_rate = c(1, 0.5, 0.1),
    max_it = 350L, iteration_info = FALSE,
    variance_parameterization = "original") {
  model_name <- deparse(substitute(x))
  validate_original_variance_parameterization_ipcr(
    variance_parameterization,
    "lavaan iteration"
  )
  validate_ipcr_it_lavaan_arguments(
    x = x,
    analytic = analytic,
    conv = conv,
    learning_rate = learning_rate,
    max_it = max_it,
    iteration_info = iteration_info
  )
  spec <- lavaan_ram_spec_ipcr(x)

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
      analytical_iterated_lavaan_components_ipcr(
        spec = spec,
        rows = rows,
        parameter_values = parameter_values,
        learning_rate = learning_rate,
        likelihood_only = likelihood_only
      )
    },
    model_name = model_name
  )
}
