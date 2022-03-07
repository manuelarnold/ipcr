coefficients_matrix <- function(x) {
  q <- length(x$info$parameters)
  cov_names <- c("(Intercept)", x$info$covariates)
  k <- length(cov_names)

  res <- do.call(rbind, lapply(x$regression, FUN = function(X) summary(X)$coefficients))
  res <- as.data.frame(res)

  res <- cbind(
    Parameter = rep(x$info$parameters, each = k),
    Covariate = rep(cov_names, times = q),
    res
  )

  rownames(res) <- NULL

  res
}
