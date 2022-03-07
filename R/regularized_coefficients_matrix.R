regularized_coefficents_matrix <- function(x) {
  q <- length(x$info$parameters)
  cov_names <- c("(Intercept)", x$info$covariates)
  k <- length(cov_names)

  Estimate <- do.call(rbind, lapply(x$regularized_regression_list,
                               FUN = function(X) coef(X, s = x$info$s)))
  Estimate <- Estimate[, 1]

  res <- data.frame(
    Parameter = rep(x$info$parameters, each = k),
    Covariate = rep(cov_names, times = q),
    Estimate = Estimate
  )

  res
}
