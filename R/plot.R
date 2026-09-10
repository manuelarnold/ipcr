#' @title Plot Correlations Between Predictors and IPCs
#'
#' @description Generates a heatmap visualizing the correlations between
#' predictors and individual parameter contributions (IPCs).
#'
#' @param x An object of class \code{"ipcr"}.
#' @param print_corr Logical; if \code{TRUE}, correlation coefficients are
#' displayed on the heatmap. Default is \code{FALSE}.
#' @param ... Additional arguments (currently unused).
#'
#' @details This function is a wrapper for \code{\link[ggplot2]{ggplot}}.
#' Factor and logical predictors are represented by the columns of the design
#' matrix used for the IPC regression. Correlations involving a constant IPC or
#' predictor column are undefined and are displayed in grey.
#'
#' @return A \code{ggplot} object.
#'
#' @export

plot.ipcr <- function(x, print_corr = FALSE, ...) {
  if (!is.logical(print_corr) || length(print_corr) != 1L || is.na(print_corr)) {
    stop("'print_corr' must be either TRUE or FALSE.", call. = FALSE)
  }

  # Use the exact complete-case sample and factor coding from the IPC
  # regression. The intercept is not a predictor and is removed.
  analysis_IPCs <- normalize_ipcr_response_matrix(
    stats::model.response(stats::model.frame(x$mlm)),
    parameters = x$info$parameters
  )
  predictors <- stats::model.matrix(x$mlm)
  predictors <- predictors[, colnames(predictors) != "(Intercept)", drop = FALSE]

  ipc_standard_deviations <- apply(analysis_IPCs, 2L, stats::sd)
  predictor_standard_deviations <- apply(predictors, 2L, stats::sd)
  valid_ipcs <- is.finite(ipc_standard_deviations) &
    ipc_standard_deviations > 0
  valid_predictors <- is.finite(predictor_standard_deviations) &
    predictor_standard_deviations > 0

  if (any(!valid_ipcs) || any(!valid_predictors)) {
    undefined_columns <- c(
      if (any(!valid_ipcs)) {
        paste0(
          "IPC parameters: ",
          paste(colnames(analysis_IPCs)[!valid_ipcs], collapse = ", ")
        )
      },
      if (any(!valid_predictors)) {
        paste0(
          "predictor columns: ",
          paste(colnames(predictors)[!valid_predictors], collapse = ", ")
        )
      }
    )
    warning(
      "Correlations are undefined for constant columns (",
      paste(undefined_columns, collapse = "; "),
      "). Undefined heatmap cells are shown in grey.",
      call. = FALSE
    )
  }

  COR <- matrix(
    NA_real_,
    nrow = ncol(analysis_IPCs),
    ncol = ncol(predictors),
    dimnames = list(colnames(analysis_IPCs), colnames(predictors))
  )
  if (any(valid_ipcs) && any(valid_predictors)) {
    COR[valid_ipcs, valid_predictors] <- stats::cor(
      x = analysis_IPCs[, valid_ipcs, drop = FALSE],
      y = predictors[, valid_predictors, drop = FALSE]
    )
  }

  # Transform data into long format
  number_parameters <- nrow(COR)
  number_predictors <- ncol(COR)
  long_data <- data.frame(
    parameter = rep(rownames(COR), times = number_predictors),
    predictor = rep(colnames(COR), each = number_parameters),
    value = c(COR)
  )

  # Heatmap
  res <- ggplot2::ggplot(
    data = long_data,
    ggplot2::aes(x = predictor, y = parameter, fill = value)
  ) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                                  midpoint = 0, limit = c(-1, 1), space = "Lab",
                                  name = "Corr.", na.value = "grey80") +
    ggplot2::labs(x = "Predictors", y = "Parameters") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        vjust = 1,
        size = 12,
        hjust = 1
      )
    ) +
    ggplot2::coord_fixed()

  if (print_corr) {
    res <- res + ggplot2::geom_text(
      ggplot2::aes(label = round(value, digits = 2)),
      color = "black"
    )
  }

  res
}

normalize_ipcr_response_matrix <- function(response, parameters) {
  response <- as.matrix(response)
  if (ncol(response) != length(parameters)) {
    stop(
      "Unable to align the stored IPC responses with the model parameters.",
      call. = FALSE
    )
  }
  colnames(response) <- parameters
  response
}
