#' @export

coeftest.ipcr <- function(x, vcov. = NULL, df = NULL, parameter = NULL, ...) {

  if (is.null(parameter)) {
    parameter <- x$info$parameters
  }

  q <- length(parameter)
  cov_names <- c("(Intercept)", x$info$covariates)
  k <- length(cov_names)

  res <- do.call(rbind, lapply(x$regression[parameter],
                               FUN = coeftest, vcov. = vcov., df = df, ...))
  res <- as.data.frame(res)

  res <- cbind(
    Parameter = rep(parameter, each = k),
    Covariate = rep(cov_names, times = q),
    res
  )

  rownames(res) <- NULL

  res
}
