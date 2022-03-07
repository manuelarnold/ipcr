coef_ipcr <- function(x, ...) {
  UseMethod("coef_ipcr")
}

coef_ipcr.default <- function(x, ...) {
  coef(x, ...)
}

coef_ipcr.lm <- function(x, ...) {
  res <- stats::coef(x)
  res <- c(res, "residual" = mean(x$residuals^2))
  res
}

coef_ipcr.lmerMod <- function(x, ...) {
  var_corr <- as.data.frame(lme4::VarCorr(x))$vcov
  parts <- lme4::getME(x, "theta")
  names(var_corr) <- c(paste("cov", names(parts), sep = "_"),
                       "residual")
  res <- c(lme4::fixef(x), var_corr)
  res
}

coef_ipcr.glmerMod <- function(x, ...) {
  var_corr <- as.data.frame(lme4::VarCorr(x))$vcov
  parts <- lme4::getME(x, "theta")
  names(var_corr) <- c(paste("cov", names(parts), sep = "_"),
                       "residual")
  res <- c(lme4::fixef(x), var_corr)
  res
}

coef_ipcr.lavaan <- function(x, ...) {
  if (x@Model@eq.constraints) {
    res <- lavaan::coef(x)[!duplicated(names(lavaan::coef(x)))]
  } else {
    res <- lavaan::coef(x)
  }
  res
}
