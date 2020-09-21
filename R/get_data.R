#' @title Get the data used to fit the model
#' @description This function gets the data used to fit the model and returns it as data frame.
#' @param x a fitted model object
#' @... arguments passed to methods.
#' @return Returns a data frame.
#' @export

get_data <- function(x, ...) {
  UseMethod("get_data")
}

get_data.glm <- function(x, ...) {
  x$data[, all.vars(formula(x)), drop = FALSE]
}

get_data.lavaan <- function(x, ...) {
  lavaan::lavInspect(object = x, what = "data")
}

get_data.lm <- function(x, ...) {
  x$model
}

get_data.MxRAMModel <- function(x, ...) {
  x$data$observed[, x$manifestVars, drop = FALSE]
}
