estfun_ipcr <- function(x, ...) {
  UseMethod("estfun_ipcr")
}

estfun_ipcr.default <- function(x, ...) {
  sandwich::estfun(x, ...)
}

estfun_ipcr.lm <- function(x, ...) {
  res <- sandwich::estfun(x)
  sigma2_ML <- mean(x$residuals^2)
  res <- cbind(res, "residual" = x$residuals^2 / (2 * sigma2_ML^2) - 1 / (2 * sigma2_ML))
}

estfun_ipcr.lmerMod <- function(x, ...) {
  merDeriv::estfun.lmerMod(x, level = 1)
}

estfun_ipcr.glmerMod <- function(x, ...) {
  merDeriv::estfun.glmerMod(x, level = 1)
}
