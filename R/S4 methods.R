#' @export

setOldClass("ipcr")

setMethod(f = "AIC", signature = signature(object = "ipcr"),
          definition = function(object, parameter = NULL, ..., k = 2) {
            if (is.null(parameter)) {
              parameter <- object$info$parameters
            }
            res <- sapply(X = object$regression_list, FUN = AIC, ...)
            res[parameter]
          })

setMethod(f = "BIC", signature = signature(object = "ipcr"),
          definition = function(object, parameter = NULL, ...) {
            if (is.null(parameter)) {
              parameter <- object$info$parameters
            }
            res <- sapply(X = object$regression_list, ..., FUN = BIC, ...)
            res[parameter]
          })

setMethod(f = "coef", signature = signature(object = "ipcr"),
          definition = function(object, parameter = NULL, ...) {
            if (is.null(parameter)) {
              parameter <- object$info$parameters
            }
            res <- sapply(X = object$regression_list, FUN = coef)
            res[, parameter]
          })

setMethod(f = "fitted", signature = signature(object = "ipcr"),
          definition = function(object, parameter = NULL, ...) {
            if (is.null(parameter)) {
              parameter <- object$info$parameters
            }
            res <- sapply(X = object$regression_list, FUN = fitted, ...)
            res[, parameter]
          })

setMethod(f = "logLik", signature = signature(object = "ipcr"),
          definition = function(object, parameter = NULL, ...) {
            if (is.null(parameter)) {
              parameter <- object$info$parameters
            }
            res <- sapply(X = object$regression_list, FUN = logLik, ...)
            res[parameter]
          })

setMethod(f = "nobs", signature = signature(object = "ipcr"),
          definition = function(object, parameter = NULL, ...) {
            if (is.null(parameter)) {
              parameter <- object$info$parameters
            }
            res <- sapply(X = object$regression_list, FUN = nobs, ...)
            res[parameter]
          })

setMethod(f = "predict", signature = signature(object = "ipcr"),
          definition = function(object, parameter = NULL, ...) {
            if (is.null(parameter)) {
              parameter <- object$info$parameters
            }
            res <- sapply(X = object$regression_list, FUN = predict, ...)
            res[, parameter]
          })

setMethod(f = "residuals", signature = signature(object = "ipcr"),
          definition = function(object, parameter = NULL, ...) {
            if (is.null(parameter)) {
              parameter <- object$info$parameters
            }
            res <- sapply(X = object$regression_list, FUN = residuals, ...)
            res[, parameter]
          })

setMethod(f = "sigma", signature = signature(object = "ipcr"),
          definition = function(object, parameter = NULL, ...) {
            if (is.null(parameter)) {
              parameter <- object$info$parameters
            }
            res <- sapply(X = object$regression_list, FUN = sigma, ...)
            res[parameter]
          })
