#' @title Calculate Variance-Covariance Matrices for the IPC Regression Models
#' @description Returns the variance-covariance matrix of the main parameters of the IPC regression models. The “main” parameters of the model correspond to those returned by coef, and do not contain a nuisance scale parameter (sigma).
#' @param object an ipcr object.
#' @param parameter a character vector with target model parameters that are to be inspected.
#' @param method a character string indicating whether to show results iterated IPC regression (iterated) or standard IPC regression (standard).
#' @param ... further arguments passed to and from methods.
#' @return A single integer.
#' @export

vcov.ipcr <- function(object, parameter = NULL, method = "iterated", ...) {

  if (method == "iterated" & is.null(object$iterated_regression)) {
    method <- "standard"
  }

  if (method == "standard" & is.null(object$regression)) {
    stop("No IPC regression found")
  }

  if (is.null(parameter)) {
    parameter <- object$info$parameters
  }

  if (method == "standard") {
    res <- lapply(X = object$regression, FUN = vcov)
  } else {
    res <- lapply(X = object$iterated_regression, FUN = vcov)
  }

  res[parameter]

}
