#' @export
AIC.ipcr <- function(object, parameter = NULL, k = 2, ...) {
  if (is.null(parameter)) {
    parameter <- object$info$parameters
  }
  res <- sapply(X = object$regression_list, FUN = AIC, ...)
  res[parameter]
}

#' @export
BIC.ipcr <- function(object, parameter = NULL, ...) {
  if (is.null(parameter)) {
    parameter <- object$info$parameters
  }
  res <- sapply(X = object$regression_list, ..., FUN = BIC, ...)
  res[parameter]
}

#' @export
effects.ipcr <- function(object, parameter = NULL, ...) {
  if (is.null(parameter)) {
    parameter <- object$info$parameters
  }
  res <- sapply(X = object$regression_list, FUN = effects, ...)
  res[, parameter]
}

#' @export
vcov.ipcr <- function(object, parameter = NULL, ...) {

  if (is.null(parameter)) {
    parameter <- object$info$parameters
  }
  res <- lapply(X = object$regression_list, FUN = vcov, ...)
  res[parameter]
}
