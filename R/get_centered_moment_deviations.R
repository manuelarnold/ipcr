#' @title Calculate Centered Individual Deviations
#' @description Internal function that calculates the individual deviations from the predicted individual- or group-specific mean and covariance.
#' @param x list with information about the model.
#' @return Matrix with centered individual devations from the sample moments.
#' @export

get_centered_moment_deviations <- function(x) {
  reg_list <- apply(X = x$data_obs, MARGIN = 2, FUN = function(X) {
    lm(X ~ x$covariates)
  })
  data_centered <- x$data_obs -
    sapply(X = reg_list, FUN = function(x) {x$fitted.values})
  cent_mc <- matrix(data = apply(X = data_centered, MARGIN = 1,
                                 FUN = function(x) {matrixcalc::vech(x %*% t(x))}),
                    nrow = x$N, ncol = x$p_star, byrow = TRUE)
  if (x$mean_structure) {
    cent_mc <- cbind(cent_mc, x$data_obs)
  }
  return(cent_mc)
}
