#' @title Calculate Individual Deviations from the Sample Moments
#' @description Internal function that calculates the individual deviations from the sample mean and sample covariance.
#' @param x list with information about the model.
#' @return Matrix with individual deviations from the sample mean and sample covariance.
#' @export

get_moment_deviations <- function(x) {
  cd <- scale(x = x$data_obs, center = TRUE, scale = FALSE)
  mc <- matrix(data = apply(X = cd, MARGIN = 1,
                            FUN = function (x) {matrixcalc::vech(x %*% t(x))}),
               nrow = x$N, ncol = x$p_star, byrow = TRUE)
  vech_cov <- matrix(data = rep(x = matrixcalc::vech(x$exp_cov), times = x$N),
                     byrow = TRUE, nrow = x$N, ncol = x$p_star)
  md <- mc - vech_cov
  if (x$mean_structure) {
    exp_means <- mxGetExpected(model = x$fit, component = "means")
    means <- matrix(data = rep(x = exp_means, times = x$N), byrow = TRUE,
                    nrow = x$N, ncol = x$p)
    mean_dev <- x$data_obs - means
    md <- cbind(md, mean_dev)
  }
  md <- as.matrix(x = md)
}
