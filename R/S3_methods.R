#' @export
coefci.ipcr <- function(x, ...) {
  lmtest::coefci(x$mlm, ...)
}

#' @export
coeftest.ipcr <- function(x, ...) {
  lmtest::coeftest(x$mlm, ...)
}

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
