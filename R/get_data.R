get_data <- function(x, ...) {
  UseMethod("get_data")
}

get_data.glm <- function(x, ...) {
  x$data[, all.vars(stats::formula(x)), drop = FALSE]
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
