#' @title Iterated Individual Parameter Contribution Regression
#' @description Generic function for computing iterated individual parameter contribution regression.
#' @param x a fitted model object.
#' @... arguments passed to methods.
#' @return Returns a list.
#' @export

iterated_ipcr <- function(x, ...) {
  UseMethod("iterated_ipcr")
}

iterated_ipcr.lm <- function(x, ...) {
  stop("Iterated IPC regression has not been implemented yet for lm-type models.")
}

iterated_ipcr.glm <- function(x, ...) {
  stop("Iterated IPC regression has not been implemented yet for glm-type models.")
}
