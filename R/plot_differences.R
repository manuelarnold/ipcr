#' Plot Estimated Conditional Differences in Model Parameters
#'
#' Generates plots of estimated model parameters as a function of one or more
#' predictors.
#'
#' @param x A standard or iterated object inheriting from class \code{"ipcr"}.
#'   Regularized \code{"ipcr_reg"} objects are not supported because ordinary
#'   coefficient confidence intervals are unavailable.
#' @param parameter A character vector specifying model-parameter names as
#'   displayed by \code{print}. By default, all model parameters are plotted.
#' @param predictor A character vector specifying original predictor names. By
#'   default, the effects of all predictors are plotted.
#' @param confidence_level A single numeric value strictly between zero and one
#'   indicating the confidence level. The default is \code{0.95}.
#'
#' @details
#' Continuous predictors are shown with a fitted regression line and a
#' confidence band. Factors, logical predictors, and numeric predictors coded
#' with exactly the values zero and one are shown using estimated group values.
#' Multi-level factors are displayed in a single plot.
#'
#' Predictions and confidence intervals are computed from the design matrix of
#' the second-stage IPC regression. Other numeric predictors are set to zero,
#' logical predictors to \code{FALSE}, and factors to their reference level.
#' Confidence intervals use regression t critical values and the residual
#' degrees of freedom of the second-stage regression.
#'
#' @return If exactly one parameter-predictor combination is requested, a
#'   \code{ggplot} object. Otherwise, a named list of \code{ggplot} objects.
#'
#' @examples
#' # Generate data
#' set.seed(1)
#' ## predictors
#' group <- rep(0:1, each = 50)
#' age <- rnorm(n = 100)
#' predictors <- data.frame(group = group, age = age)
#' ## Model data
#' x <- rnorm(n = 100)
#' y <- 0.5 + (0.25 + 0.75 * group) * x +
#'   rnorm(n = 100, sd = sqrt(0.75))
#' d <- data.frame(x = x, y = y)
#'
#' # Fit a linear regression
#' m <- lm(y ~ x, data = d)
#'
#' # Investigate model with IPC regression
#' res <- ipcr(m, predictors = predictors)
#'
#' # Return all plots in a named list
#' plots <- plot_differences(res)
#'
#' # Plot the values of the regression slope x by group
#' plot_differences(res, parameter = "x", predictor = "group")
#'
#' @seealso \code{\link{plot.ipcr}}
#' @export
plot_differences <- function(x, parameter = NULL, predictor = NULL,
                             confidence_level = 0.95) {
  if (!inherits(x, "ipcr")) {
    stop("'x' must be an object of class \"ipcr\".", call. = FALSE)
  }
  if (inherits(x, "ipcr_reg")) {
    stop(
      "plot_differences() is not available for regularized IPCR because ",
      "ordinary coefficient confidence intervals are unavailable.",
      call. = FALSE
    )
  }
  if (!is.numeric(confidence_level) || length(confidence_level) != 1L ||
      is.na(confidence_level) || !is.finite(confidence_level) ||
      confidence_level <= 0 || confidence_level >= 1) {
    stop(
      "'confidence_level' must be a single finite number strictly between 0 and 1.",
      call. = FALSE
    )
  }

  parameter <- validate_plot_names_ipcr(
    requested = parameter,
    available = x$info$parameters,
    argument = "parameter"
  )
  predictor <- validate_plot_names_ipcr(
    requested = predictor,
    available = x$info$predictors,
    argument = "predictor"
  )

  model_frame <- stats::model.frame(x$mlm)
  analysis_IPCs <- normalize_ipcr_response_matrix(
    stats::model.response(model_frame),
    parameters = x$info$parameters
  )
  analysis_predictors <- model_frame[, x$info$predictors, drop = FALSE]
  design_matrix <- stats::model.matrix(x$mlm)
  coefficient_matrix <- normalize_ipcr_coefficient_matrix(
    stats::coef(x$mlm),
    parameters = x$info$parameters
  )

  if (!identical(rownames(coefficient_matrix), colnames(design_matrix))) {
    stop(
      "Unable to align the second-stage coefficients with its design matrix.",
      call. = FALSE
    )
  }
  if (!all(parameter %in% colnames(analysis_IPCs))) {
    stop(
      "Unable to align requested parameters with the stored IPC responses.",
      call. = FALSE
    )
  }

  coefficient_covariance <- stats::vcov(x$mlm, complete = TRUE)
  number_terms <- nrow(coefficient_matrix)
  number_parameters <- ncol(coefficient_matrix)
  expected_covariance_size <- number_terms * number_parameters
  if (!is.matrix(coefficient_covariance) ||
      !identical(
        dim(coefficient_covariance),
        c(expected_covariance_size, expected_covariance_size)
      )) {
    stop(
      "Unable to align the second-stage coefficient covariance matrix.",
      call. = FALSE
    )
  }

  residual_df <- stats::df.residual(x$mlm)
  if (length(residual_df) != 1L || !is.finite(residual_df) || residual_df <= 0) {
    stop(
      "A positive residual degrees of freedom is required for confidence intervals.",
      call. = FALSE
    )
  }
  critical_value <- stats::qt(
    p = 0.5 + confidence_level / 2,
    df = residual_df
  )

  model_terms <- stats::delete.response(stats::terms(x$mlm))
  reference_data <- predictor_reference_data_ipcr(analysis_predictors)
  plots <- vector("list", length(parameter) * length(predictor))
  plot_names <- character(length(plots))
  plot_index <- 1L

  for (parameter_name in parameter) {
    parameter_index <- match(parameter_name, colnames(coefficient_matrix))
    covariance_index <- (parameter_index - 1L) * number_terms +
      seq_len(number_terms)
    parameter_covariance <- coefficient_covariance[
      covariance_index,
      covariance_index,
      drop = FALSE
    ]
    parameter_coefficients <- coefficient_matrix[, parameter_index]

    for (predictor_name in predictor) {
      predictor_values <- analysis_predictors[[predictor_name]]
      categorical <- is.factor(predictor_values) ||
        is.logical(predictor_values) ||
        is_zero_one_ipcr(predictor_values)

      grid_values <- if (is.factor(predictor_values)) {
        levels(predictor_values)
      } else if (is.logical(predictor_values)) {
        c(FALSE, TRUE)
      } else if (categorical) {
        c(0, 1)
      } else {
        seq(
          from = min(predictor_values),
          to = max(predictor_values),
          length.out = 100L
        )
      }

      newdata <- reference_data[rep(1L, length(grid_values)), , drop = FALSE]
      if (is.factor(predictor_values)) {
        newdata[[predictor_name]] <- factor(
          grid_values,
          levels = levels(predictor_values),
          ordered = is.ordered(predictor_values)
        )
      } else {
        newdata[[predictor_name]] <- grid_values
      }

      prediction_design <- stats::model.matrix(
        model_terms,
        data = newdata,
        contrasts.arg = x$mlm$contrasts,
        xlev = x$mlm$xlevels
      )
      prediction_design <- prediction_design[
        , rownames(coefficient_matrix), drop = FALSE
      ]
      prediction <- ipcr_plot_prediction(
        design = prediction_design,
        coefficients = parameter_coefficients,
        covariance = parameter_covariance,
        critical_value = critical_value,
        parameter = parameter_name,
        predictor = predictor_name
      )

      if (categorical) {
        plots[[plot_index]] <- categorical_difference_plot_ipcr(
          predictor_values = predictor_values,
          IPC_values = analysis_IPCs[, parameter_name],
          grid_values = grid_values,
          prediction = prediction,
          predictor_name = predictor_name,
          parameter_name = parameter_name
        )
      } else {
        plots[[plot_index]] <- continuous_difference_plot_ipcr(
          predictor_values = predictor_values,
          IPC_values = analysis_IPCs[, parameter_name],
          grid_values = grid_values,
          prediction = prediction,
          predictor_name = predictor_name,
          parameter_name = parameter_name
        )
      }

      plot_names[plot_index] <- paste0(
        "parameter: ", parameter_name, " | predictor: ", predictor_name
      )
      plot_index <- plot_index + 1L
    }
  }

  if (length(plots) == 1L) {
    return(plots[[1L]])
  }
  names(plots) <- make.unique(plot_names)
  plots
}

validate_plot_names_ipcr <- function(requested, available, argument) {
  if (is.null(requested)) {
    return(available)
  }
  if (!is.character(requested) || length(requested) == 0L ||
      anyNA(requested) || any(!nzchar(requested))) {
    stop(
      "'", argument, "' must be NULL or a nonempty character vector.",
      call. = FALSE
    )
  }
  if (anyDuplicated(requested)) {
    stop("'", argument, "' must not contain duplicate names.", call. = FALSE)
  }
  unavailable <- requested[!requested %in% available]
  if (length(unavailable) > 0L) {
    stop(
      "Unknown ", argument, " name", if (length(unavailable) == 1L) "" else "s",
      ": ", paste(unavailable, collapse = ", "), ". Available names are: ",
      paste(available, collapse = ", "), ".",
      call. = FALSE
    )
  }
  requested
}

normalize_ipcr_coefficient_matrix <- function(coefficients, parameters) {
  if (is.null(dim(coefficients))) {
    coefficients <- matrix(
      coefficients,
      ncol = 1L,
      dimnames = list(names(coefficients), parameters)
    )
  }
  coefficients
}

predictor_reference_data_ipcr <- function(predictors) {
  reference_values <- lapply(predictors, function(column) {
    if (is.factor(column)) {
      return(factor(
        levels(column)[1L],
        levels = levels(column),
        ordered = is.ordered(column)
      ))
    }
    if (is.logical(column)) {
      return(FALSE)
    }
    0
  })
  as.data.frame(reference_values, optional = TRUE, check.names = FALSE)
}

is_zero_one_ipcr <- function(x) {
  if (!is.numeric(x)) {
    return(FALSE)
  }
  values <- sort(unique(x))
  length(values) == 2L && all(values == c(0, 1))
}

ipcr_plot_prediction <- function(design, coefficients, covariance,
                                 critical_value, parameter, predictor) {
  active_columns <- apply(abs(design) > sqrt(.Machine$double.eps), 2L, any)
  if (any(!is.finite(coefficients[active_columns])) ||
      any(!is.finite(
        covariance[active_columns, active_columns, drop = FALSE]
      ))) {
    stop(
      "The effect of predictor '", predictor, "' on parameter '", parameter,
      "' is not estimable from the second-stage regression.",
      call. = FALSE
    )
  }

  active_design <- design[, active_columns, drop = FALSE]
  active_coefficients <- coefficients[active_columns]
  active_covariance <- covariance[
    active_columns,
    active_columns,
    drop = FALSE
  ]
  estimate <- drop(active_design %*% active_coefficients)
  variance <- rowSums((active_design %*% active_covariance) * active_design)
  tolerance <- sqrt(.Machine$double.eps)
  if (any(variance < -tolerance)) {
    stop(
      "A negative prediction variance was obtained for predictor '", predictor,
      "' and parameter '", parameter, "'.",
      call. = FALSE
    )
  }
  standard_error <- sqrt(pmax(variance, 0))

  data.frame(
    estimate = estimate,
    std_error = standard_error,
    lower = estimate - critical_value * standard_error,
    upper = estimate + critical_value * standard_error
  )
}

categorical_difference_plot_ipcr <- function(predictor_values, IPC_values,
                                             grid_values, prediction,
                                             predictor_name, parameter_name) {
  display_levels <- as.character(grid_values)
  observed_data <- data.frame(
    predictor_value = factor(
      as.character(predictor_values),
      levels = display_levels
    ),
    IPCs = IPC_values
  )
  prediction_data <- cbind(
    data.frame(
      predictor_value = factor(display_levels, levels = display_levels)
    ),
    prediction
  )

  ggplot2::ggplot(
    observed_data,
    ggplot2::aes(x = predictor_value, y = IPCs, colour = predictor_value)
  ) +
    ggplot2::geom_jitter(width = 0.1, height = 0, alpha = 0.6) +
    ggplot2::geom_errorbar(
      data = prediction_data,
      ggplot2::aes(
        x = predictor_value,
        ymin = lower,
        ymax = upper,
        colour = predictor_value
      ),
      inherit.aes = FALSE,
      width = 0.2
    ) +
    ggplot2::geom_point(
      data = prediction_data,
      ggplot2::aes(
        x = predictor_value,
        y = estimate,
        colour = predictor_value
      ),
      inherit.aes = FALSE,
      size = 3
    ) +
    ggplot2::labs(x = predictor_name, y = parameter_name) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none")
}

continuous_difference_plot_ipcr <- function(predictor_values, IPC_values,
                                            grid_values, prediction,
                                            predictor_name, parameter_name) {
  observed_data <- data.frame(
    predictor_value = predictor_values,
    IPCs = IPC_values
  )
  prediction_data <- cbind(
    data.frame(predictor_value = grid_values),
    prediction
  )

  ggplot2::ggplot(
    observed_data,
    ggplot2::aes(x = predictor_value, y = IPCs)
  ) +
    ggplot2::geom_point() +
    ggplot2::geom_ribbon(
      data = prediction_data,
      ggplot2::aes(x = predictor_value, ymin = lower, ymax = upper),
      inherit.aes = FALSE,
      alpha = 0.2,
      fill = "#0072B2"
    ) +
    ggplot2::geom_line(
      data = prediction_data,
      ggplot2::aes(x = predictor_value, y = estimate),
      inherit.aes = FALSE,
      colour = "#0072B2",
      linewidth = 1.05
    ) +
    ggplot2::labs(x = predictor_name, y = parameter_name) +
    ggplot2::theme_bw()
}

utils::globalVariables(c(
  "estimate", "IPCs", "lower", "parameter", "predictor", "predictor_value",
  "upper", "value"
))
