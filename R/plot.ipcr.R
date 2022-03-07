#' @title Plot Correlations between Covariates and IPCs
#' @description This functions plots a heat map that visualizes the correlation between
#' covariates and IPCs.
#' @param x an ipcr object.
#' @param print_corr a logical value; if TRUE correlation coefficients are added on the heatmap.
#' @param ... other arguments.
#' @details This function is a wrapper for \code{\link[ggplot2]{ggplot}}. Currently,
#' arguments passed to \code{ggplot} cannot be changed.
#' @export

plot.ipcr <- function(x, print_corr = FALSE, ...) {

  # check for covariates
  if (x$info$covariates[1] == "no covariates found") {
    stop("no covariates found")
  }

  # Re-code characters and factors into dummy variables
  covariates <- x$regression_list[[1]]$model[, -1]
  if (any(unlist(lapply(covariates, function(x) {is.character(x) | is.factor(x)})))) {
    covariates_strings <- Filter(function(x) {is.character(x) | is.factor(x)}, covariates)
    string_formula <- paste(colnames(covariates_strings), collapse = "+")
    dummies <- stats::model.matrix(stats::formula(paste("~", string_formula)), data = covariates_strings)[, -1]
    covariates <- covariates[, !unlist(lapply(covariates, function(x) {is.character(x) | is.factor(x)}))]
    covariates <- cbind(covariates, dummies)
  }

  # Correlation matrix
  COR <- stats::cor(x = x$IPCs, y = covariates)

  # Transform data into long format
  p <- nrow(COR)
  q <- ncol(COR)


  long_data <- data.frame(parameter = rep(x$info$parameters, times = q),
                          covariate = rep(x$info$covariates, each = p),
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
    ggplot2::labs(x = "Covariates", y = "Parameters") +
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
