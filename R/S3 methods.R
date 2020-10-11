#' @export

effects.ipcr <- function(object, parameter = NULL, ...) {
  if (is.null(parameter)) {
    parameter <- object$info$parameters
  }
  res <- sapply(X = object$regression_list, FUN = effects, ...)
  res[, parameter]
}

vcov.ipcr <- function(object, parameter = NULL, ...) {

  if (is.null(parameter)) {
    parameter <- object$info$parameters
  }
  res <- lapply(X = object$regression, FUN = vcov, ...)
  res[[, parameter]]
}
