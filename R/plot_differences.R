#' @title Plot Estimated Conditional Differences in a Model Parameter
#' @description This functions plots the conditional values of the model parameters
#' as a function of the covariates.
#' @param x an ipcr object.
#' @param parameter a string. The name of a model parameter as in \code{coef}.
#' Per default all model parameters are plotted
#' @param covariate a string. The name of a covariate. Per default, all the effects of all
#' covariates are plotted.
#' @param confidence_level a numeric. The confidence interval plotted. 0.95 (resulting) in
#' 95\% confidence interval is the default.
#' @param ... other arguments.
#' @details
#' Note that regression lines (for continuous covariates), means (for dummy variables),
#' and confidence intervals are plotted, using the estimates of the model parameters and
#' the corresponding variances and covariances of the parameter estimates. The plots show
#' the estimated parameter value as a function of a covariate with all other covariates
#' set equal to zero. (Mean) centering of the covariates may increase the interpretability
#' of the plots.
#'
#' This function is a wrapper for \code{\link[ggplot2]{ggplot}}.
#' @examples
#' # Generate data
#' ## Covariates
#' z1 <- rep(0:1, each = 50)
#' z2 <- rnorm(n = 100)
#' covariates <- data.frame(z1 = z1, z2 = z2)
#' ## Model data
#' x <- rnorm(n = 100)
#' y <- 0.5 + 0.75*z1 + rnorm(n = 100, sd = sqrt(0.75))
#' d <- data.frame(x = x, y = y)
#'
#' # Fit a linear regression
#' m <- lm(y ~ x, data = d)
#'
#' # Investigate model with IPC regression
#' ipc <- ipcr(m, covariates = covariates)
#'
#' # Plotting parameter values as a function of the covariates
#' ## All plots
#' plot_differences(ipc)
#' ## Plot the values of the regression slope x as a function of the covariate z1
#' plot_differences(ipc, parameter = "x", covariate = "z1")
#' @seealso \code{\link{plot.ipcr}},
#' @export

plot_differences <- function(x, parameter = NULL, covariate = NULL,
                             confidence_level = 0.95, ...) {

  # check for covariates
  if (x$info$covariates[1] == "no covariates found") {
    stop("no covariates found")
  }

  # z-value for confidence interval
  z <- stats::qnorm(p = 0.5 + confidence_level / 2)

  if (is.null(parameter)) {parameter <- x$info$parameters}
  if (is.null(covariate)) {
    covariate_id <- seq_along(x$info$covariates)
  } else {
    covariate_id <- which(x$info$covariate %in% covariate)
  }

  # loop through plots
  for (q in parameter) {
    for (k in covariate_id) {

      # model parameter estimates
      m <- x$regression_list[[q]]
      intercept <- stats::coef(m)[1]
      slope <- stats::coef(m)[k + 1] # position 1 is the intercept
      VCOV <- stats::vcov(m)[c(1, k + 1), c(1, k + 1)]
      var_intercept <- VCOV[1, 1]
      cov_intercept_slope <- VCOV[2, 1]
      var_slope <- VCOV[2, 2]

      # plotting data
      df <- data.frame(k = x$regression_list[[1]]$model[, k + 1],
                       IPCs = x$IPCs[, q])

      # for dummy variables
      if (length(unique(df[, "k"])) == 2) {

        print(ggplot2::ggplot(data = df,
                              mapping = ggplot2::aes(x = factor(k),
                                                     y = IPCs,
                                                     col = factor(k))) +
                ggplot2::geom_jitter() +
                ggplot2::geom_point(ggplot2::aes(x = 1, y = intercept), col = "#F8766D", size = 5) +
                ggplot2::geom_point(ggplot2::aes(x = 2, y = intercept + slope), col = "#00BFC4", size = 5) +
                ggplot2::geom_errorbar(ggplot2::aes(x = 1,
                                                    ymin = intercept - z*sqrt(var_intercept),
                                                    ymax = intercept + z*sqrt(var_intercept)),
                                       width = 0.25, col = "#F8766D") +
                ggplot2::geom_errorbar(ggplot2::aes(x = 2,
                                                    ymin = intercept + slope - z*sqrt(var_intercept + 2*cov_intercept_slope + var_slope),
                                                    ymax = intercept + slope + z*sqrt(var_intercept + 2*cov_intercept_slope + var_slope)),
                                       width = 0.25, col = "#00BFC4") +
                ggplot2::xlab(x$info$covariates[k]) +
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
                ggplot2::xlab(x$info$covariates[k]) +
                ggplot2::ylab(q) +
                ggplot2::theme_bw())
      }
    }
  }
}
