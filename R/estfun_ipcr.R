#' @noRd
estfun_ipcr <- function(x, analytic = FALSE, ...) {
  UseMethod("estfun_ipcr")
}

#' @noRd
#' @export
estfun_ipcr.default <- function(x, ...) {
  sandwich::estfun(x, ...)
}

#' @noRd
#' @export
estfun_ipcr.lavaan <- function(x, ...) {
  scores <- lavaan::lavScores(x)
  parameter_names <- names(coef_ipcr(x))
  if (NCOL(scores) == length(parameter_names)) {
    colnames(scores) <- parameter_names
  }
  scores
}

#' @noRd
#' @export
estfun_ipcr.lm <- function(x, ...) {
  res <- sandwich::estfun(x)
  sigma2_ML <- mean(x$residuals^2)
  res <- cbind(res, "residual" = x$residuals^2 / (2 * sigma2_ML^2) -
                 1 / (2 * sigma2_ML))
}

#' @noRd
#' @export
estfun_ipcr.glm <- function(x, ...) {
  sandwich::estfun(x, ...)
}

#' @noRd
#' @export
estfun_ipcr.lmerMod <- function(x, ...) {
  merDeriv::estfun.lmerMod(x, level = 2)
}

#' @noRd
#' @export
estfun_ipcr.glmerMod <- function(x, ...) {
  merDeriv::estfun.glmerMod(x)
}
