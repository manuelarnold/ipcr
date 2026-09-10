#' @title Extract Scores
#' @description Extracts the casewise or clusterwise scores of a fitted model
#' using the score implementation for its backend. Supported lavaan, OpenMx,
#' and mixed-model objects use their dedicated score extractors; compatible
#' unsupported model classes are attempted through
#' \code{\link[sandwich]{estfun}}.
#' @param x a fitted model object.
#' @param analytic Logical. If \code{FALSE} (default), functions of
#' \pkg{lavaan}, \pkg{OpenMx}, or \pkg{sandwich} will be used to compute scores.
#' If \code{TRUE}, fast analytical scores are requested for a continuous
#' raw-data \code{MxRAMModel}. Unsupported OpenMx models use numerical scores
#' with a warning. See \code{\link{ipcr}} for the supported analytical scope.
#' @param ... Additional arguments passed to the backend-specific score method.
#' @return A numeric matrix containing the empirical estimating functions.
#' Typically, this is an \eqn{n * k} matrix corresponding to \eqn{n} IPC units
#' (clusters for mixed models) and \eqn{k} parameters. Its columns are named
#' according to the fitted model parameters.
#'
#' The estimating function (or score function) for a model is the derivative of the
#' objective function with respect to the parameter vector. The empirical estimating
#' functions are the evaluations of the estimating function at the observed
#' data (\eqn{n} observations) and the estimated parameters (of dimension
#' \eqn{k}).
#' @references
#' Zeileis, A.  (2006). Object-oriented computation of sandwich estimators. \emph{Journal of Statistical Software, 16}(9),
#' 1-16. doi: \href{https://doi.org/10.18637/jss.v016.i09}{10.18637/jss.v016.i09}
#'
#' Zeileis, A., Köll, S., Graham, N. (2020). Various versatile variances: An object-
#' oriented implementation of clustered covariances in R. \emph{Journal of Statistical Software, 95}(1),
#' 1-36. doi \href{https://doi.org/10.18637/jss.v095.i01}{10.18637/jss.v095.i01}
#' @seealso \code{\link[sandwich]{estfun}}
#' @export

get_scores <- function(x, analytic = FALSE, ...) {
  analytic <- prepare_model_ipcr(x, analytic)
  check_model_convergence_ipcr(x)
  validate_scores_ipcr(
    param_estimates = coef_ipcr(x),
    scores = estfun_ipcr(x, analytic = analytic, ...),
    n = nobs_ipcr(x)
  )
}
