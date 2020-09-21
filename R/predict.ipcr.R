#' @title Predict Method for Linear Model Fits
#' @description Predicted values based on the IPC regression models.
#' @param object an ipcr object.
#' @param parameter a character vector with target model parameters that are to be inspected.
#' @param method a character string indicating whether to show results iterated IPC regression (iterated) or standard IPC regression (standard).
#' @param ... further arguments passed to and from methods.
#' @return A single integer.
#' @export predict.ipcr

predict.ipcr <- function(object, parameter = NULL, method = "iterated", newdata,
                         se.fit = FALSE, scale = NULL, df = Inf,
                         interval = c("none", "confidence", "prediction"),
                         level = 0.95, type = c("response", "terms"), terms = NULL,
                         na.action = na.pass, pred.var = res.var/weights,
                         weights = 1, ...) {

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
    res <- sapply(X = object$regression[parameter], FUN = predict, newdata, se.fit, scale, df,
                  interval, level, type, terms, na.action, pred.var, weights)
  } else {
    res <- sapply(X = object$iterated_regression[parameter], FUN = predict, newdata, se.fit,
                  scale, df, interval, level, type, terms, na.action, pred.var,
                  weights)
  }

  res

}
