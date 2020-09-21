#' @title Extract Residual Standard Deviation 'Sigma'
#' @description Extract the estimated standard deviation of the errors from the IPC regression models.
#' @param object an ipcr object.
#' @param parameter a character vector with target model parameters that are to be inspected.
#' @param method a character string indicating whether to show results iterated IPC regression (iterated) or standard IPC regression (standard).
#' @param ... further arguments passed to and from methods.
#' @return A single integer.
#' @export sigma.ipcr

sigma.ipcr <- function(object, parameter = NULL, method = "iterated", ...) {

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
    res <- sapply(X = object$regression, FUN = sigma)
  } else {
    res <- sapply(X = object$iterated_regression, FUN = sigma)
  }

  res[parameter]

}
