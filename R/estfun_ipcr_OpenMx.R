# Compute OpenMx row gradients without printing an optimizer-style progress
# line for every score extraction.
numerical_scores_MxModel_ipcr <- function(x) {
  scores <- NULL
  invisible(utils::capture.output(
    scores <- suppressMessages(-0.5 * OpenMx::imxRowGradients(model = x)),
    type = "output"
  ))

  parameter_names <- names(OpenMx::omxGetParameters(x))
  if (NCOL(scores) == length(parameter_names)) {
    colnames(scores) <- parameter_names
  }
  scores
}

# Return a logical matrix identifying OpenMx definition-variable labels.
definition_label_mask_ipcr <- function(labels) {
  mask <- vapply(
    as.vector(labels),
    function(label) {
      !is.na(label) && nzchar(label) &&
        isTRUE(OpenMx::imxIsDefinitionVariable(label))
    },
    logical(1)
  )
  matrix(mask, nrow = NROW(labels), ncol = NCOL(labels))
}

# Add OpenMx's generated names for unlabelled free A, S, and M parameters to
# local copies of the RAM label matrices. The fitted model itself is unchanged.
ram_parameter_labels_ipcr <- function(x) {
  labels <- list(
    A = x$A$labels,
    S = x$S$labels,
    M = x$M$labels
  )
  parameter_names <- names(x$output$estimate)
  model_prefix <- paste0(x$name, ".")

  for (parameter_name in parameter_names) {
    is_already_labelled <- any(vapply(
      labels,
      function(matrix_labels) {
        any(matrix_labels == parameter_name, na.rm = TRUE)
      },
      logical(1)
    ))
    if (is_already_labelled || !startsWith(parameter_name, model_prefix)) {
      next
    }

    generated_name <- substring(parameter_name, nchar(model_prefix) + 1L)
    match <- regexec(
      "^([ASM])\\[([0-9]+),([0-9]+)\\]$",
      generated_name
    )
    pieces <- regmatches(generated_name, match)[[1]]
    if (length(pieces) != 4L) {
      next
    }

    matrix_name <- pieces[2]
    row <- as.integer(pieces[3])
    column <- as.integer(pieces[4])
    if (row > NROW(labels[[matrix_name]]) ||
        column > NCOL(labels[[matrix_name]])) {
      next
    }

    labels[[matrix_name]][row, column] <- parameter_name
    if (identical(matrix_name, "S")) {
      labels[[matrix_name]][column, row] <- parameter_name
    }
  }

  labels
}

openmx_definition_labels_ipcr <- function(x) {
  matrix_labels <- unlist(
    lapply(x@matrices, function(matrix) as.vector(matrix$labels)),
    use.names = FALSE
  )
  matrix_labels <- matrix_labels[!is.na(matrix_labels)]
  unique(matrix_labels[vapply(
    matrix_labels,
    OpenMx::imxIsDefinitionVariable,
    logical(1)
  )])
}

ram_definition_labels_ipcr <- function(x) {
  labels <- unlist(
    lapply(
      list(x$A$labels, x$S$labels, x$M$labels),
      as.vector
    ),
    use.names = FALSE
  )
  labels <- labels[!is.na(labels)]
  unique(labels[vapply(
    labels,
    OpenMx::imxIsDefinitionVariable,
    logical(1)
  )])
}

# Return a concise reason why the fast RAM score calculation is unavailable.
# A NULL result means that the model is within the documented analytical scope.
openmx_analytical_score_issue_ipcr <- function(x) {
  if (length(x@constraints) > 0L) {
    return("models containing explicit mxConstraint() objects")
  }
  if (length(x@algebras) > 0L) {
    return("models containing OpenMx algebras")
  }
  if (length(x@submodels) > 0L) {
    return("RAM models containing submodels")
  }
  if (length(x@penalties) > 0L) {
    return("penalized OpenMx models")
  }
  if (any(x$F$free)) {
    return("RAM models with free parameters in the F matrix")
  }

  manifest_data <- as.data.frame(
    x$data$observed[, x$manifestVars, drop = FALSE]
  )
  continuous_manifest <- vapply(
    manifest_data,
    function(variable) is.numeric(variable) && !is.ordered(variable),
    logical(1)
  )
  if (!all(continuous_manifest)) {
    return("ordinal or non-numeric manifest variables")
  }

  if (any(!is.na(x$data@weight)) || any(!is.na(x$data@frequency))) {
    return("weighted or frequency-weighted raw data")
  }

  all_definition_labels <- openmx_definition_labels_ipcr(x)
  ram_definition_labels <- ram_definition_labels_ipcr(x)
  unsupported_definition_labels <- setdiff(
    all_definition_labels,
    ram_definition_labels
  )
  if (length(unsupported_definition_labels) > 0L) {
    return("definition variables outside the A, S, and M matrices")
  }

  if (length(ram_definition_labels) > 0L) {
    definition_names <- sub("^data\\.", "", ram_definition_labels)
    missing_columns <- setdiff(definition_names, colnames(x$data$observed))
    if (length(missing_columns) > 0L) {
      return("definition variables that are absent from the raw data")
    }
    definition_data <- as.data.frame(
      x$data$observed[, definition_names, drop = FALSE]
    )
    numeric_definition <- vapply(definition_data, is.numeric, logical(1))
    if (!all(numeric_definition) || anyNA(definition_data)) {
      return("non-numeric or missing definition-variable values")
    }
  }

  labels <- ram_parameter_labels_ipcr(x)
  parameter_names <- names(x$output$estimate)
  matched_parameters <- vapply(
    parameter_names,
    function(parameter_name) {
      any(vapply(
        labels,
        function(matrix_labels) {
          any(matrix_labels == parameter_name, na.rm = TRUE)
        },
        logical(1)
      ))
    },
    logical(1)
  )
  if (!all(matched_parameters)) {
    return("free parameters outside the A, S, and M matrices")
  }

  NULL
}

# Construct dA/dtheta, dS/dtheta, and dm/dtheta once. Repeated labels mark all
# matrix locations governed by the same free parameter and therefore implement
# equality restrictions specified through shared labels.
ram_derivatives_ipcr <- function(x) {
  labels <- ram_parameter_labels_ipcr(x)
  parameter_names <- names(x$output$estimate)
  latent_dimension <- NROW(x$A$values)

  lapply(parameter_names, function(parameter_name) {
    dA <- matrix(0, latent_dimension, latent_dimension)
    dS <- matrix(0, latent_dimension, latent_dimension)
    dm <- matrix(0, latent_dimension, 1L)
    dA[labels$A == parameter_name & !is.na(labels$A)] <- 1
    dS[labels$S == parameter_name & !is.na(labels$S)] <- 1
    mean_labels <- t(labels$M)
    dm[mean_labels == parameter_name & !is.na(mean_labels)] <- 1
    list(A = dA, S = dS, m = dm)
  })
}

# Substitute one row of definition-variable values and, when supplied, a
# candidate free-parameter vector into local RAM matrices. Keeping this as a
# pure matrix operation is important for iterated IPCR: no OpenMx optimizer or
# expectation recalculation is invoked for the observation-specific models.
ram_matrices_for_row_ipcr <- function(x, row, parameters = NULL,
                                      parameter_labels = NULL) {
  matrices <- list(
    A = x$A$values,
    S = x$S$values,
    M = x$M$values
  )
  labels <- list(
    A = x$A$labels,
    S = x$S$labels,
    M = x$M$labels
  )

  for (matrix_name in names(matrices)) {
    definition_positions <- which(
      definition_label_mask_ipcr(labels[[matrix_name]]),
      arr.ind = TRUE
    )
    if (NROW(definition_positions) == 0L) {
      next
    }

    for (position in seq_len(NROW(definition_positions))) {
      matrix_row <- definition_positions[position, 1]
      matrix_column <- definition_positions[position, 2]
      definition_label <- labels[[matrix_name]][matrix_row, matrix_column]
      definition_name <- sub("^data\\.", "", definition_label)
      matrices[[matrix_name]][matrix_row, matrix_column] <-
        x$data$observed[row, definition_name]
    }
  }

  if (!is.null(parameters)) {
    if (is.null(parameter_labels)) {
      parameter_labels <- ram_parameter_labels_ipcr(x)
    }
    parameter_names <- names(x$output$estimate)
    if (is.null(names(parameters))) {
      names(parameters) <- parameter_names
    }
    if (!identical(names(parameters), parameter_names)) {
      stop(
        "Candidate OpenMx parameter names or order do not match the fitted model.",
        call. = FALSE
      )
    }

    for (parameter_name in parameter_names) {
      for (matrix_name in names(matrices)) {
        parameter_positions <-
          parameter_labels[[matrix_name]] == parameter_name &
          !is.na(parameter_labels[[matrix_name]])
        matrices[[matrix_name]][parameter_positions] <-
          parameters[[parameter_name]]
      }
    }
  }

  list(
    A = matrices$A,
    S = matrices$S,
    m = t(matrices$M),
    F = x$F$values
  )
}

# Obtain the most restrictive OpenMx bounds attached to every free RAM
# parameter. Shared labels may occur at several matrix locations; all of their
# bounds must be respected by an observation-specific parameter prediction.
ram_parameter_bounds_ipcr <- function(x, parameter_labels = NULL) {
  if (is.null(parameter_labels)) {
    parameter_labels <- ram_parameter_labels_ipcr(x)
  }
  parameter_names <- names(x$output$estimate)
  lower <- stats::setNames(rep(-Inf, length(parameter_names)), parameter_names)
  upper <- stats::setNames(rep(Inf, length(parameter_names)), parameter_names)

  for (parameter_name in parameter_names) {
    for (matrix_name in names(parameter_labels)) {
      positions <- parameter_labels[[matrix_name]] == parameter_name &
        !is.na(parameter_labels[[matrix_name]])
      if (!any(positions)) {
        next
      }

      matrix_lower <- x[[matrix_name]]$lbound[positions]
      matrix_upper <- x[[matrix_name]]$ubound[positions]
      matrix_lower <- matrix_lower[!is.na(matrix_lower)]
      matrix_upper <- matrix_upper[!is.na(matrix_upper)]
      if (length(matrix_lower) > 0L) {
        lower[[parameter_name]] <- max(lower[[parameter_name]], matrix_lower)
      }
      if (length(matrix_upper) > 0L) {
        upper[[parameter_name]] <- min(upper[[parameter_name]], matrix_upper)
      }
    }
  }

  list(lower = lower, upper = upper)
}

ram_moments_ipcr <- function(matrices) {
  identity <- diag(NROW(matrices$A))
  inverse_path <- solve(identity - matrices$A)
  filtered_inverse <- matrices$F %*% inverse_path
  latent_covariance <- inverse_path %*% matrices$S %*% t(inverse_path)

  list(
    inverse_path = inverse_path,
    filtered_inverse = filtered_inverse,
    latent_covariance = latent_covariance,
    covariance = filtered_inverse %*% matrices$S %*%
      t(filtered_inverse),
    means = as.vector(filtered_inverse %*% matrices$m)
  )
}

ram_jacobian_ipcr <- function(matrices, moments, derivatives, observed) {
  vapply(
    derivatives,
    function(derivative) {
      covariance_half <- moments$filtered_inverse %*% derivative$A %*%
        moments$latent_covariance %*% t(matrices$F)
      covariance_derivative <- covariance_half + t(covariance_half) +
        moments$filtered_inverse %*% derivative$S %*%
        t(moments$filtered_inverse)
      mean_derivative <- moments$filtered_inverse %*% derivative$A %*%
        moments$inverse_path %*% matrices$m +
        moments$filtered_inverse %*% derivative$m

      c(
        lavaan::lav_matrix_vech(
          covariance_derivative[observed, observed, drop = FALSE]
        ),
        as.vector(mean_derivative[observed, , drop = FALSE])
      )
    },
    numeric(length(observed) * (length(observed) + 3L) / 2L)
  )
}

openmx_score_groups_ipcr <- function(x, manifest_data) {
  missingness <- apply(
    is.na(manifest_data),
    1L,
    paste0,
    collapse = ""
  )
  grouping_data <- data.frame(missingness = missingness)

  definition_labels <- ram_definition_labels_ipcr(x)
  if (length(definition_labels) > 0L) {
    definition_names <- unique(sub("^data\\.", "", definition_labels))
    grouping_data <- cbind(
      grouping_data,
      x$data$observed[, definition_names, drop = FALSE]
    )
  }

  group <- interaction(grouping_data, drop = TRUE, lex.order = TRUE)
  split(seq_len(NROW(manifest_data)), group, drop = TRUE)
}

# Construct the row-level quantities that do not change during analytical RAM
# iteration. Keeping this separate lets ipcr_it.MxRAMModel() prepare them once.
iterated_ram_grouping_data_ipcr <- function(x, manifest_data) {
  missingness <- apply(
    is.na(manifest_data),
    1L,
    paste0,
    collapse = ""
  )
  grouping_data <- data.frame(missingness = missingness)

  definition_labels <- ram_definition_labels_ipcr(x)
  if (length(definition_labels) > 0L) {
    definition_names <- unique(sub("^data\\.", "", definition_labels))
    grouping_data <- cbind(
      grouping_data,
      x$data$observed[, definition_names, drop = FALSE]
    )
  }
  grouping_data
}

# Group rows that share every quantity needed to construct local RAM moments
# and expected information. Scores still vary with the observed data, but a
# group can reuse its covariance inverse, Jacobian, and Fisher information.
iterated_ram_groups_ipcr <- function(
    x, manifest_data, rows, parameter_values, grouping_data = NULL) {
  if (is.null(grouping_data)) {
    grouping_data <- iterated_ram_grouping_data_ipcr(x, manifest_data)
  }
  grouping_data <- grouping_data[rows, , drop = FALSE]
  grouping_data <- cbind(
    grouping_data,
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

analytical_iterated_ram_context_ipcr <- function(
    x,
    parameterization = ram_variance_parameterization_ipcr(
      x,
      variance_parameterization = "original"
    )) {
  manifest_data <- as.matrix(
    x$data$observed[, x$manifestVars, drop = FALSE]
  )
  storage.mode(manifest_data) <- "double"
  parameter_labels <- ram_parameter_labels_ipcr(x)

  list(
    manifest_data = manifest_data,
    parameter_names = parameterization$parameter_names,
    parameterization = parameterization,
    parameter_labels = parameter_labels,
    parameter_bounds = ram_parameter_bounds_ipcr(x, parameter_labels),
    derivatives = ram_derivatives_ipcr(x),
    grouping_data = iterated_ram_grouping_data_ipcr(x, manifest_data)
  )
}

# Evaluate the damped IPC update at observation-specific RAM parameter values.
# This is the computational core used by ipcr_it.MxRAMModel().
analytical_iterated_ram_components_ipcr <- function(
    x, rows, parameter_values, learning_rate, context = NULL,
    likelihood_only = FALSE) {
  if (is.null(context)) {
    context <- analytical_iterated_ram_context_ipcr(x)
  }
  manifest_data <- context$manifest_data
  rows <- as.integer(rows)
  parameter_names <- context$parameter_names
  parameterization <- context$parameterization
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

  parameter_labels <- context$parameter_labels
  parameter_bounds <- context$parameter_bounds
  derivatives <- context$derivatives
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

  groups <- iterated_ram_groups_ipcr(
    x = x,
    manifest_data = manifest_data,
    rows = rows,
    parameter_values = parameter_values,
    grouping_data = context$grouping_data
  )

  for (group_positions in groups) {
    representative_position <- group_positions[1L]
    representative_row <- rows[representative_position]
    coordinate_parameters <- parameter_values[representative_position, ]
    names(coordinate_parameters) <- parameter_names
    parameters <- iterated_coordinates_to_natural_ipcr(
      coordinate_parameters,
      parameterization
    )

    below_bound <- parameters < parameter_bounds$lower
    above_bound <- parameters > parameter_bounds$upper
    if (any(below_bound | above_bound)) {
      offending <- parameter_names[below_bound | above_bound]
      stop(
        "Predicted OpenMx parameters violate fitted-model bounds: ",
        paste(offending, collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    observed <- which(!is.na(manifest_data[representative_row, ]))
    if (length(observed) == 0L) {
      stop(
        "Iterated IPCR cannot update a row with no observed manifest variables.",
        call. = FALSE
      )
    }

    matrices <- ram_matrices_for_row_ipcr(
      x = x,
      row = representative_row,
      parameters = parameters,
      parameter_labels = parameter_labels
    )
    moments <- ram_moments_ipcr(matrices)
    covariance <- moments$covariance[
      observed,
      observed,
      drop = FALSE
    ]
    covariance <- 0.5 * (covariance + t(covariance))
    covariance_chol <- chol(covariance)
    covariance_inverse <- chol2inv(covariance_chol)
    means <- moments$means[observed]
    observed_count <- length(observed)

    original_rows <- rows[group_positions]
    centered <- sweep(
      manifest_data[original_rows, observed, drop = FALSE],
      2L,
      means,
      FUN = "-"
    )
    if (!isTRUE(likelihood_only)) {
      covariance_count <- observed_count * (observed_count + 1L) / 2L
      moment_count <- covariance_count + observed_count
      duplication <- lavaan::lav_matrix_duplication(observed_count)
      weight <- matrix(0, nrow = moment_count, ncol = moment_count)
      covariance_rows <- seq_len(covariance_count)
      mean_rows <- covariance_count + seq_len(observed_count)
      weight[covariance_rows, covariance_rows] <-
        0.5 * t(duplication) %*%
        kronecker(covariance_inverse, covariance_inverse) %*%
        duplication
      weight[mean_rows, mean_rows] <- covariance_inverse

      jacobian <- ram_jacobian_ipcr(
        matrices = matrices,
        moments = moments,
        derivatives = derivatives,
        observed = observed
      )
      coordinate_derivative <- rep(1, q)
      coordinate_derivative[parameterization$transformed] <-
        parameters[parameterization$transformed]
      jacobian <- sweep(
        jacobian,
        2L,
        coordinate_derivative,
        FUN = "*"
      )
      information <- t(jacobian) %*% weight %*% jacobian
      information <- 0.5 * (information + t(information))
      chol(information)
      deviations <- t(vapply(
        seq_len(NROW(centered)),
        function(position) {
          centered_row <- centered[position, ]
          c(
            lavaan::lav_matrix_vech(
              tcrossprod(centered_row) - covariance
            ),
            centered_row
          )
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

# Group rows for the numerical iterated calculation. In contrast to the fast
# analytical RAM path, definition variables can occur anywhere in an OpenMx
# model, so all definition-variable labels must contribute to the key.
iterated_openmx_numerical_groups_ipcr <- function(
    x, manifest_data, rows, parameter_values) {
  grouping_data <- data.frame(
    missingness = apply(
      is.na(manifest_data[rows, , drop = FALSE]),
      1L,
      paste0,
      collapse = ""
    )
  )

  definition_labels <- openmx_definition_labels_ipcr(x)
  if (length(definition_labels) > 0L) {
    definition_names <- unique(sub("^data\\.", "", definition_labels))
    grouping_data <- cbind(
      grouping_data,
      x$data$observed[rows, definition_names, drop = FALSE]
    )
  }
  grouping_data <- cbind(
    grouping_data,
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

openmx_manifest_moment_names_ipcr <- function(observed, manifest_count) {
  covariance_names <- matrix(
    NA_character_,
    nrow = manifest_count,
    ncol = manifest_count
  )
  for (column in seq_len(manifest_count)) {
    for (row in seq_len(manifest_count)) {
      covariance_names[row, column] <- paste0(
        "cov", max(row, column), "_", min(row, column)
      )
    }
  }

  c(
    lavaan::lav_matrix_vech(
      covariance_names[observed, observed, drop = FALSE]
    ),
    paste0("mean", observed)
  )
}

# Evaluate iterated IPC updates through OpenMx's numerical manifest-moment
# Jacobian. This is intentionally a slow fallback: each distinct predicted
# parameter vector requires an expectation recalculation and a numerical
# Jacobian, while the initial scores also require imxRowGradients().
numerical_iterated_ram_components_ipcr <- function(
    x, rows, parameter_values, learning_rate, likelihood_only = FALSE) {
  manifest_data <- as.matrix(
    x$data$observed[, x$manifestVars, drop = FALSE]
  )
  storage.mode(manifest_data) <- "double"
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

  parameter_bounds <- ram_parameter_bounds_ipcr(x)
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
  manifest_count <- NCOL(manifest_data)

  groups <- iterated_openmx_numerical_groups_ipcr(
    x = x,
    manifest_data = manifest_data,
    rows = rows,
    parameter_values = parameter_values
  )

  for (group_positions in groups) {
    representative_position <- group_positions[1L]
    representative_row <- rows[representative_position]
    parameters <- parameter_values[representative_position, ]
    names(parameters) <- parameter_names

    below_bound <- parameters < parameter_bounds$lower
    above_bound <- parameters > parameter_bounds$upper
    if (any(below_bound | above_bound)) {
      offending <- parameter_names[below_bound | above_bound]
      stop(
        "Predicted OpenMx parameters violate fitted-model bounds: ",
        paste(offending, collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    observed <- which(!is.na(manifest_data[representative_row, ]))
    if (length(observed) == 0L) {
      stop(
        "Iterated IPCR cannot update a row with no observed manifest variables.",
        call. = FALSE
      )
    }

    local_model <- OpenMx::omxSetParameters(
      model = x,
      labels = parameter_names,
      values = parameters
    )
    invisible(utils::capture.output(
      local_model <- suppressMessages(OpenMx::mxRun(
        local_model,
        useOptimizer = FALSE,
        silent = TRUE
      )),
      type = "output"
    ))

    covariance <- OpenMx::mxGetExpected(
      model = local_model,
      component = "covariance",
      defvar.row = representative_row
    )[observed, observed, drop = FALSE]
    covariance <- 0.5 * (covariance + t(covariance))
    covariance_chol <- chol(covariance)
    covariance_inverse <- chol2inv(covariance_chol)
    means <- as.vector(OpenMx::mxGetExpected(
      model = local_model,
      component = "means",
      defvar.row = representative_row
    ))[observed]

    original_rows <- rows[group_positions]
    centered <- sweep(
      manifest_data[original_rows, observed, drop = FALSE],
      2L,
      means,
      FUN = "-"
    )
    observed_count <- length(observed)
    if (!isTRUE(likelihood_only)) {
      full_jacobian <- OpenMx::omxManifestModelByParameterJacobian(
        model = local_model,
        defvar.row = representative_row
      )
      moment_names <- openmx_manifest_moment_names_ipcr(
        observed = observed,
        manifest_count = manifest_count
      )
      moment_rows <- match(moment_names, rownames(full_jacobian))
      if (anyNA(moment_rows)) {
        stop(
          "Unable to align the numerical OpenMx Jacobian with the observed moments.",
          call. = FALSE
        )
      }
      jacobian <- full_jacobian[moment_rows, , drop = FALSE]
      if (!is.null(colnames(jacobian)) &&
          all(parameter_names %in% colnames(jacobian))) {
        jacobian <- jacobian[, parameter_names, drop = FALSE]
      } else if (NCOL(jacobian) == q) {
        colnames(jacobian) <- parameter_names
      } else {
        stop(
          "Unable to align the numerical OpenMx Jacobian with the free parameters.",
          call. = FALSE
        )
      }

      covariance_count <- observed_count * (observed_count + 1L) / 2L
      moment_count <- covariance_count + observed_count
      duplication <- lavaan::lav_matrix_duplication(observed_count)
      weight <- matrix(0, nrow = moment_count, ncol = moment_count)
      covariance_rows <- seq_len(covariance_count)
      mean_rows <- covariance_count + seq_len(observed_count)
      weight[covariance_rows, covariance_rows] <-
        0.5 * t(duplication) %*%
        kronecker(covariance_inverse, covariance_inverse) %*%
        duplication
      weight[mean_rows, mean_rows] <- covariance_inverse
      information <- t(jacobian) %*% weight %*% jacobian
      information <- 0.5 * (information + t(information))
      chol(information)
      deviations <- t(vapply(
        seq_len(NROW(centered)),
        function(position) {
          centered_row <- centered[position, ]
          c(
            lavaan::lav_matrix_vech(
              tcrossprod(centered_row) - covariance
            ),
            centered_row
          )
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

analytical_scores_MxRAMModel_ipcr <- function(x) {
  manifest_data <- as.matrix(
    x$data$observed[, x$manifestVars, drop = FALSE]
  )
  storage.mode(manifest_data) <- "double"
  n <- NROW(manifest_data)
  parameter_names <- names(x$output$estimate)
  derivatives <- ram_derivatives_ipcr(x)
  scores <- matrix(
    0,
    nrow = n,
    ncol = length(parameter_names),
    dimnames = list(NULL, parameter_names)
  )

  groups <- openmx_score_groups_ipcr(x, manifest_data)
  for (group_rows in groups) {
    representative <- group_rows[1]
    observed <- which(!is.na(manifest_data[representative, ]))

    # A row with no observed manifest variable has no likelihood contribution
    # for the RAM parameters. OpenMx ordinarily removes or rejects such rows,
    # but retaining a zero score is the appropriate limiting result.
    if (length(observed) == 0L) {
      next
    }

    matrices <- ram_matrices_for_row_ipcr(x, representative)
    moments <- ram_moments_ipcr(matrices)
    covariance <- moments$covariance[
      observed,
      observed,
      drop = FALSE
    ]
    means <- moments$means[observed]
    covariance_inverse <- solve(covariance)
    observed_count <- length(observed)
    covariance_count <- observed_count * (observed_count + 1L) / 2L
    moment_count <- covariance_count + observed_count
    duplication <- lavaan::lav_matrix_duplication(observed_count)

    weight <- matrix(0, nrow = moment_count, ncol = moment_count)
    covariance_rows <- seq_len(covariance_count)
    mean_rows <- covariance_count + seq_len(observed_count)
    weight[covariance_rows, covariance_rows] <-
      0.5 * t(duplication) %*%
      kronecker(covariance_inverse, covariance_inverse) %*%
      duplication
    weight[mean_rows, mean_rows] <- covariance_inverse

    jacobian <- ram_jacobian_ipcr(
      matrices = matrices,
      moments = moments,
      derivatives = derivatives,
      observed = observed
    )
    deviations <- t(vapply(
      group_rows,
      function(row) {
        centered <- as.numeric(manifest_data[row, observed]) - means
        c(
          lavaan::lav_matrix_vech(
            tcrossprod(centered) - covariance
          ),
          centered
        )
      },
      numeric(moment_count)
    ))

    scores[group_rows, ] <- deviations %*% weight %*% jacobian
  }

  scores
}

#' @noRd
#' @export
estfun_ipcr.MxModel <- function(x, analytic = FALSE, ...) {
  numerical_scores_MxModel_ipcr(x)
}

#' @noRd
#' @export
estfun_ipcr.MxRAMModel <- function(x, analytic = FALSE, ...) {
  if (isTRUE(analytic)) {
    analytical_scores_MxRAMModel_ipcr(x)
  } else {
    numerical_scores_MxModel_ipcr(x)
  }
}
