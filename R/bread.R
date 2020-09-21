#' @title Bread Matrix of the Sandwich Matrix
#' @description These functions extract an estimator for the bread of sandwiches.
#' @param x a fitted model object.
#' @... arguments passed to methods.
#' @return Returns a matrix.
#' @export

bread.lavaan <- function(x, ...) {
  solve(lavaan::lavInspect(x, what = "information.expected"))
}

bread.MxRAMModel <- function(x, ...) {
  x$output$vcov * x$data$numObs
}
