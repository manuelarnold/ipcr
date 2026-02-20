#' @noRd
estfun_ipcr <- function(x, analytic = FALSE, ...) {
  UseMethod("estfun_ipcr")
}

#' @noRd
estfun_ipcr.default <- function(x, ...) {
  sandwich::estfun(x, ...)
}

#' @noRd
estfun_ipcr.lavaan <- function(x, ...) {
  lavaan::estfun.lavaan(x)
}

#' @noRd
estfun_ipcr.lm <- function(x, ...) {
  res <- sandwich::estfun(x)
  sigma2_ML <- mean(x$residuals^2)
  res <- cbind(res, "residual" = x$residuals^2 / (2 * sigma2_ML^2) -
                 1 / (2 * sigma2_ML))
}

#' @noRd
estfun_ipcr.lmerMod <- function(x, ...) {
  merDeriv::estfun.lmerMod(x, level = 1)
}

#' @noRd
estfun_ipcr.glmerMod <- function(x, ...) {
  merDeriv::estfun.glmerMod(x, level = 1)
}
