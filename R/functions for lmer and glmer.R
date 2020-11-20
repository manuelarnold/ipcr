coef.lmerMod <- function(x, ...) {
  var_corr <- as.data.frame(lme4::VarCorr(x))$vcov
  parts <- lme4::getME(x, "theta")
  names(var_corr) <- c(paste("cov", names(parts), sep = "_"),
                       "residual")
  res <- c(lme4::fixef(x), var_corr)
  res
}

coef.glmerMod <- function(x, ...) {
  var_corr <- as.data.frame(lme4::VarCorr(x))$vcov
  parts <- lme4::getME(x, "theta")
  names(var_corr) <- c(paste("cov", names(parts), sep = "_"),
                       "residual")
  res <- c(lme4::fixef(x), var_corr)
  res
}

estfun.lmerMod <- function(x, ...) {
  merDeriv::estfun.lmerMod(x, level = 1)
}

estfun.glmerMod <- function(x, ...) {
  merDeriv::estfun.glmerMod(x, level = 1)
}

bread.lmerMod <- function(x, ...) {
  as.matrix(merDeriv::bread.lmerMod(x, full = TRUE))
}

bread.glmerMod <- function(x, ...) {
  as.matrix(merDeriv::bread.glmerMod(x, full = TRUE))
}
