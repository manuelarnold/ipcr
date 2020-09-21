#' @title Extract Model Fitted Values
#' @description This function extracts the fitted values from IPC regression models.
#' @param object an ipcr object.
#' @param parameter a character vector with target model parameters that are to be inspected.
#' @param method a character string indicating whether to show results iterated IPC regression (iterated) or standard IPC regression (standard).
#' @... further arguments passed to and from methods.
#' @return Returns a list.
#' @export fitted.ipcr

fitted.ipcr <- function(object, parameter = NULL, method = "iterated", ...) {

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
    res <- sapply(X = object$regression[parameter], FUN = function(x) {x$fitted.values})
  } else {
    res <- sapply(X = object$iterated_regression[parameter], FUN = function(x) {x$fitted.values})
  }

  res

}
