# Internal coef() that behaves differently to the standard coef()


#' @noRd
coef_ipcr <- function(x, ...) {
  UseMethod("coef_ipcr")
}

#' @noRd
#' @export
coef_ipcr.default <- function(x, ...) {
  coef(x, ...)
}

#' @noRd
#' @export
coef_ipcr.lm <- function(x, ...) {
  res <- stats::coef(x)
  res <- c(res, "residual" = mean(x$residuals^2))
  res
}

mixed_variance_parameters_ipcr <- function(x, include_residual) {
  variance_parameters <- as.data.frame(
    lme4::VarCorr(x),
    order = "lower.tri"
  )$vcov
  random_effect_names <- paste(
    "cov",
    names(lme4::getME(x, "theta")),
    sep = "_"
  )
  parameter_names <- if (include_residual) {
    c(random_effect_names, "residual")
  } else {
    random_effect_names
  }

  if (length(variance_parameters) != length(parameter_names)) {
    stop(
      "Unable to align mixed-model variance estimates with merDeriv parameters.",
      call. = FALSE
    )
  }
  stats::setNames(variance_parameters, parameter_names)
}

#' @noRd
#' @export
coef_ipcr.glm <- function(x, ...) {
  stats::coef(x)
}

#' @noRd
#' @export
coef_ipcr.lmerMod <- function(x, ...) {
  variance_parameters <- mixed_variance_parameters_ipcr(
    x,
    include_residual = TRUE
  )
  res <- c(lme4::fixef(x), variance_parameters)
  res
}

#' @noRd
#' @export
coef_ipcr.glmerMod <- function(x, ...) {
  variance_parameters <- mixed_variance_parameters_ipcr(
    x,
    include_residual = FALSE
  )
  res <- c(lme4::fixef(x), variance_parameters)
  res
}

#' @noRd
#' @export
coef_ipcr.MxModel <- function(x, ...) {
  OpenMx::omxGetParameters(x)
}

#' @noRd
#' @export
coef_ipcr.lavaan <- function(x, ...) {
  estimates <- lavaan::coef(x)
  reduction <- lavaan_constraint_reduction_ipcr(x@Model)
  if (!is.null(reduction)) {
    estimates <- estimates[reduction$keep]
  }
  estimates
}
