#' @export
setOldClass("ipcr")

#' @export
setMethod(f = "coef", signature = signature(object = "ipcr"),
          definition = function(object, parameter = NULL, ...) {
            if (is.null(parameter)) {
              parameter <- object$info$parameters
            }
            res <- sapply(X = object$regression_list, FUN = coef)
            res[, parameter]
          })

#' @export
setMethod(f = "fitted", signature = signature(object = "ipcr"),
          definition = function(object, parameter = NULL, ...) {
            if (is.null(parameter)) {
              parameter <- object$info$parameters
            }
            res <- sapply(X = object$regression_list, FUN = fitted, ...)
            res[, parameter]
          })

#' @export
setMethod(f = "logLik", signature = signature(object = "ipcr"),
          definition = function(object, parameter = NULL, ...) {
            if (is.null(parameter)) {
              parameter <- object$info$parameters
            }
            res <- sapply(X = object$regression_list, FUN = logLik, ...)
            res[parameter]
          })

#' @export
setMethod(f = "nobs", signature = signature(object = "ipcr"),
          definition = function(object, parameter = NULL, ...) {
            if (is.null(parameter)) {
              parameter <- object$info$parameters
            }
            res <- sapply(X = object$regression_list, FUN = nobs, ...)
            res[parameter]
          })

#' @export
setMethod(f = "predict", signature = signature(object = "ipcr"),
          definition = function(object, parameter = NULL, ...) {
            if (is.null(parameter)) {
              parameter <- object$info$parameters
            }
            res <- sapply(X = object$regression_list, FUN = predict, ...)
            res[, parameter]
          })

#' @export
setMethod(f = "residuals", signature = signature(object = "ipcr"),
          definition = function(object, parameter = NULL, ...) {
            if (is.null(parameter)) {
              parameter <- object$info$parameters
            }
            res <- sapply(X = object$regression_list, FUN = residuals, ...)
            res[, parameter]
          })

#' @export
setMethod(f = "sigma", signature = signature(object = "ipcr"),
          definition = function(object, parameter = NULL, ...) {
            if (is.null(parameter)) {
              parameter <- object$info$parameters
            }
            res <- sapply(X = object$regression_list, FUN = sigma, ...)
            res[parameter]
          })
