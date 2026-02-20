#' @noRd
get_RAM_LDL_matrices <- function(x, ...) {
  UseMethod("get_RAM_LDL_matrices")
}

#' @noRd
get_RAM_LDL_matrices.lavaan <- function(x, ...) {
  sandwich::estfun(x, ...)
}
