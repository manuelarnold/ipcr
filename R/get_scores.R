#' @title Extract Scores
#' @description This functions extracts the scores of a fitted model. Arguments are passed
#' to the generic function \code{\link[sandwich]{estfun}}.
#' @param x a fitted model object.
#' @param ... arguments passed to methods.
#' @return A \code{data.frame} containing the empirical estimating functions. Typically,
#' this should be an \eqn{n * k} matrix corresponding to \eqn{n} observations and \eqn{k}
#' parameters. The columns should be named as in \code{\link[stats]{coef}} or
#' \code{\link[stats]{terms}}, respectively.
#'
#' The estimating function (or score function) for a model is the derivative of the
#' objective function with respect to the parameter vector. The empirical estimating
#' functions is the evaluation of the estimating function at the observed data (\eqn{n}
#' observations) and the estimated parameters (of dimension \eqn{k}).
#' @references
#' Zeileis, A.  (2006). Object-oriented computation of sandwich estimators. \emph{Journal of Statistical Software, 16}(9),
#' 1-16. doi: \href{https://doi.org/10.18637/jss.v016.i09}{10.18637/jss.v016.i09}
#'
#' Zeileis, A., Köll, S., Graham, N. (2020). Various versatile variances: An object-
#' oriented implementation of clustered covariances in R. \emph{Journal of Statistical Software, 95}(1),
#' 1-36. doi \href{https://doi.org/10.18637/jss.v095.i01}{10.18637/jss.v095.i01}
#' @seealso \code{\link[sandwich]{estfun}}
#' @export

get_scores <- function(x, ...) {
  data.frame(sandwich::estfun(x, ...))
}
