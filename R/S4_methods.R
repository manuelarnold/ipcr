#' @export
setOldClass("ipcr")

#' @export
setMethod(f = "coef", signature = signature(object = "ipcr"),
          definition = function(object, ...) {
            coef(object$mlm, ...)
          })

#' @export
setMethod(f = "coefficients", signature = signature(object = "ipcr"),
          definition = function(object, ...) {
            coefficients(object$mlm, ...)
          })

#' @export
setMethod(f = "fitted", signature = signature(object = "ipcr"),
          definition = function(object, ...) {
            fitted(object$mlm, ...)
          })

#' @export
setMethod(f = "nobs", signature = signature(object = "ipcr"),
          definition = function(object, ...) {
            nobs(object$mlm, ...)
          })

#' @export
setMethod(f = "predict", signature = signature(object = "ipcr"),
          definition = function(object, ...) {
            predict(object$mlm, ...)
          })

#' @export
setMethod(f = "residuals", signature = signature(object = "ipcr"),
          definition = function(object, ...) {
            residuals(object$mlm, ...)
          })

#' @export
setMethod(f = "sigma", signature = signature(object = "ipcr"),
          definition = function(object, ...) {
            sigma(object$mlm, ...)
          })
