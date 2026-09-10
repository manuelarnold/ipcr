#' @export
confint.ipcr <- function(object, ...) {
  stats::confint(object$mlm, ...)
}

#' @export
effects.ipcr <- function(object, ...) {
  effects(object$mlm, ...)
}

#' @export
vcov.ipcr <- function(object, ...) {
  vcov(object$mlm, ...)
}
