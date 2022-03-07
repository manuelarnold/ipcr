#' @export

confint.ipcr <- function(object, parm, level = 0.95, parameter = NULL, ...) {

  if (missing(parm)) {
    cov_names <- c("(Intercept)", object$info$covariates)
  } else {
    cov_names <- parm
  }

  if (is.null(parameter)) {
    parameter <- object$info$parameters
  }

  q <- length(parameter)
  k <- length(cov_names)

  res <- do.call(rbind, lapply(object$regression_list[parameter],
                               FUN = confint, parm = parm, level = level, ...))
  res <- as.data.frame(res)

  res <- cbind(
    Parameter = rep(parameter, each = k),
    Covariate = rep(cov_names, times = q),
    res
  )

  rownames(res) <- NULL

  res
}
