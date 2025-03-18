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
#' Currently, arguments passed to \code{ggplot} cannot be modified.
#'
#' @export

plot.ipcr <- function(x, print_corr = FALSE, ...) {

  # Re-code characters and factors into dummy variables
  predictors <- x$mlm$model[, -1]
  if (any(unlist(lapply(predictors, function(x) {is.character(x) | is.factor(x)})))) {
    predictors_strings <- Filter(function(x) {is.character(x) | is.factor(x)}, predictors)
    string_formula <- paste(colnames(predictors_strings), collapse = "+")
    dummies <- stats::model.matrix(stats::formula(paste("~", string_formula)), data = predictors_strings)[, -1]
    predictors <- predictors[, !unlist(lapply(predictors, function(x) {is.character(x) | is.factor(x)}))]
    predictors <- cbind(predictors, dummies)
  }

  # Correlation matrix
  COR <- stats::cor(x = x$ipc, y = predictors)

  # Transform data into long format
  p <- nrow(COR)
  q <- ncol(COR)


  long_data <- data.frame(parameter = rep(x$info$parameters, times = q),
                          covariate = rep(x$info$predictors, each = p),
                          value = c(COR))

  # Heatmap
  res <- ggplot2::ggplot(data = long_data,
                         ggplot2::aes(x = long_data[, 2],
                                      y = long_data[, 1],
                                      fill = long_data[, 3])) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                                  midpoint = 0, limit = c(-1, 1), space = "Lab",
                                  name = "Corr.") +
    ggplot2::labs(x = "Predictors", y = "Parameters") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, size = 12,
                                                       hjust = 1)) +
    ggplot2::coord_fixed()

  if (print_corr) {
    res +
      ggplot2::geom_text(ggplot2:: aes(x = long_data[, 2], y =  long_data[, 1],
                                       label = round(long_data[, 3], digits = 2)),
                         color = "black")
  } else {
    res
  }
}
