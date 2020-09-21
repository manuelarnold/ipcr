#' @title Extract the Number of Observations from a MxRAMModel Fit
#' @description Extract the number of observations from a MxRAMModel fit.
#' @param object a fitted MXRAMModel-type model.
#' @... arguments passed to methods.
#' @return A single integer.
#' @export

nobs.MxRAMModel <- function(object, ...) {
  object$data$numObs
}
