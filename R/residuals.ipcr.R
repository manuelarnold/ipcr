#' @title Extract Model Resiudals
#' @description Extracts model residuals from IPC regression models.
#' @param object an ipcr object.
#' @param parameter a character vector with target model parameters that are to be inspected.
#' @param method a character string indicating whether to show results iterated IPC regression (iterated) or standard IPC regression (standard).
#' @param ... further arguments passed to and from methods.
#' @return Returns a list.
#' @export residuals.ipcr

residuals.ipcr <- function(object, parameter = NULL, method = "iterated", ...) {

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
    res <- sapply(X = object$regression[parameter], FUN = residuals)
  } else {
    res <- sapply(X = object$iterated_regression[parameter], FUN = residuals)
  }

  res

}
