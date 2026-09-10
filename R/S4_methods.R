# Accessor methods delegate to the second-stage regression stored in an ipcr
# object. ipcr objects are S3 lists, so these must be S3 methods as well.

#' @exportS3Method stats::coef
coef.ipcr <- function(object, ...) {
  stats::coef(object$mlm, ...)
}

#' @export
coefficients.ipcr <- function(object, ...) {
  stats::coefficients(object$mlm, ...)
}

#' @exportS3Method stats::fitted
fitted.ipcr <- function(object, ...) {
  stats::fitted(object$mlm, ...)
}

#' @export
nobs.ipcr <- function(object, ...) {
  stats::nobs(object$mlm, ...)
}

#' @export
predict.ipcr <- function(object, ...) {
  stats::predict(object$mlm, ...)
}

#' @export
residuals.ipcr <- function(object, ...) {
  stats::residuals(object$mlm, ...)
}

#' @export
sigma.ipcr <- function(object, ...) {
  stats::sigma(object$mlm, ...)
}
