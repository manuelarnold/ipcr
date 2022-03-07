#' @export

coefci.ipcr <- function(x, parm = NULL, level = 0.95, vcov. = NULL, df = NULL,
                        parameter = NULL, ...) {

  if (is.null(parm)) {
    cov_names <- c("(Intercept)", x$info$covariates)
  } else {
    cov_names <- parm
  }

  if (is.null(parameter)) {
    parameter <- x$info$parameters
  }

  q <- length(parameter)
  k <- length(cov_names)

  res <- do.call(rbind, lapply(x$regression_list[parameter],
                               FUN = coefci, parm = parm, level = level, vcov. = vcov.,
                               df = df, ...))
  res <- as.data.frame(res)

  res <- cbind(
    Parameter = rep(parameter, each = k),
    Covariate = rep(cov_names, times = q),
    res
  )

  rownames(res) <- NULL

  res
}
