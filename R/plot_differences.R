#' Plot Estimated Conditional Differences in Model Parameters
#'
#' Generates plots of estimated model parameters as a function of
#' one or more predictors.
#'
#' @param x An object of class \code{"ipcr"}.
#' @param parameter A character string specifying the name of a model parameter
#'   as displayed by \code{print}. By default, all model parameters are plotted.
#' @param predictor A character string specifying the name of a predictor.
#'   By default, the effects of all predictors are plotted.
#' @param confidence_level A numeric value indicating the confidence level
#'   for the confidence intervals. The default is \code{0.95}, corresponding
#'   to a 95 percent confidence interval.
#'
#' @details
#' The function plots regression lines for continuous predictors and group means for
#' dummy variables. Confidence intervals are included based on the model parameter
#' estimates and their corresponding variances and covariances. The plots visualize
#' the estimated parameter values as a function of a predictor, assuming all other
#' predictors are set to zero. Mean-centering predictors may improve interpretability.
#' This function is a wrapper for \code{\link[ggplot2]{ggplot}}.
#'
#' @examples
#' # Generate data
#' ## predictors
#' z1 <- rep(0:1, each = 50)
#' z2 <- rnorm(n = 100)
#' predictors <- data.frame(z1 = z1, z2 = z2)
#' ## Model data
#' x <- rnorm(n = 100)
#' y <- 0.5 + 0.75*z1 + rnorm(n = 100, sd = sqrt(0.75))
#' d <- data.frame(x = x, y = y)
#'
#' # Fit a linear regression
#' m <- lm(y ~ x, data = d)
#'
#' # Investigate model with IPC regression
#' res <- ipcr(m, predictors = predictors)
#'
#' # Plotting parameter values as a function of the predictors
#' ## All plots
#' plot_differences(res)
#' ## Plot the values of the regression slope x as a function of the predictor z1
#' plot_differences(res, parameter = "x", predictor = "z1")
#'
#' @seealso \code{\link{plot.ipcr}}
#' @export
#'
plot_differences <- function(x, parameter = NULL, predictor = NULL,
                             confidence_level = 0.95) {

  # z-value for confidence interval
  z <- stats::qnorm(p = 0.5 + confidence_level / 2)

  if (is.null(parameter)) {parameter <- x$info$parameters}
  if (is.null(predictor)) {
    predictor_id <- seq_along(x$info$predictors)
  } else {
    predictor_id <- which(x$info$predictor %in% predictor)
  }

  # loop through plots
  for (q in parameter) {
    for (k in predictor_id) {

      # model parameter estimates
      intercept <- coef(x)[1, q]
      slope <- coef(x)[k + 1, q] # position 1 is the intercept
      VCOV <- stats::vcov(x)[grep(pattern = q, x = rownames(vcov(x))),
                             grep(pattern = q, x = colnames(vcov(x)))]
      VCOV <- VCOV[c(1, k + 1), c(1, k + 1)]
      var_intercept <- VCOV[1, 1]
      cov_intercept_slope <- VCOV[2, 1]
      var_slope <- VCOV[2, 2]

      # plotting data
      df <- data.frame(k = x$mlm$model[, k + 1],
                       IPCs = x$IPCs[, q])

      # for dummy variables
      if (length(unique(df[, "k"])) == 2) {

        print(ggplot2::ggplot(data = df,
                              mapping = ggplot2::aes(x = factor(k),
                                                     y = IPCs,
                                                     col = factor(k))) +
                ggplot2::geom_jitter() +
                ggplot2::annotate(geom = "point", x = 1, y = intercept,
                                  color = "#F8766D", size = 5) +
                ggplot2::annotate(geom = "point", x = 2, y = intercept + slope,
                                  color = "#00BFC4", size = 5) +
                ggplot2::geom_errorbar(ggplot2::aes(x = 1,
                                                    ymin = intercept - z*sqrt(var_intercept),
                                                    ymax = intercept + z*sqrt(var_intercept)),
                                       width = 0.25, col = "#F8766D") +
                ggplot2::geom_errorbar(ggplot2::aes(x = 2,
                                                    ymin = intercept + slope - z*sqrt(var_intercept + 2*cov_intercept_slope + var_slope),
                                                    ymax = intercept + slope + z*sqrt(var_intercept + 2*cov_intercept_slope + var_slope)),
                                       width = 0.25, col = "#00BFC4") +
                ggplot2::xlab(x$info$predictors[k]) +
                ggplot2::ylab(q) +
                ggplot2::theme_bw() +
                ggplot2::theme(legend.position = "none"))

      } else { # for other variables (not dummys)

        # regression line and confidence interval functions
        regression_line <- function(x) {
          intercept + slope*x
        }

        upper_CI <- function(x) {
          intercept + slope*x + z*sqrt(var_intercept + 2*x*cov_intercept_slope + x^2*var_slope)
        }

        lower_CI <- function(x) {
          intercept + slope*x - z*sqrt(var_intercept + 2*x*cov_intercept_slope + x^2*var_slope)
        }

        # plot
        print(ggplot2::ggplot(data = df,
                              mapping = ggplot2::aes(x = k, y = IPCs)) +
                ggplot2::geom_point() +
                ggplot2::geom_function(fun = regression_line, col = "#0072B2", size = 1.05) +
                ggplot2::geom_function(fun = upper_CI, col = "#0072B2", linetype="dashed") +
                ggplot2::geom_function(fun = lower_CI, col = "#0072B2", linetype="dashed") +
                ggplot2::xlab(x$info$predictors[k]) +
                ggplot2::ylab(q) +
                ggplot2::theme_bw())
      }
    }
  }
}
