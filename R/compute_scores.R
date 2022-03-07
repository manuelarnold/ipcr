#' @title Compute Scores
#' @description This functions calculates the scores of a fitted model. The
#' scores are either computed analytically or numerically.
#' @param x a fitted model object.
#' @param analytic a logical value; if TRUE scores are computed analytically
#' (default), else scores are computed numerically.
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
#' @seealso \code{\link[sandwich]{estfun}}

compute_scores <- function(x, analytic = TRUE, ...) {

  if (analytic == TRUE) {

    scores <- data.frame(estfun_ipcr(x, ...))

  }

  if (analytic == FALSE) {

    scores <- compute_numeric_scores(x, ...)

  }

  scores
}
