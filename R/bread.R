bread.lavaan <- function(x, ...) {
  solve(lavaan::lavInspect(x, what = "information.expected"))
}

bread.MxRAMModel <- function(x, ...) {
  x$output$vcov * x$data$numObs
}
