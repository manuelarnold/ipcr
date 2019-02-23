#' @title Calculate the Weight Matrix
#' @description Internal function that calculates the estimated weight matrix.
#' @param x list with information about the model.
#' @return Estimated weight matrix.
#' @export

get_weight_matrix2 <- function(x, exp_cov_inv_up) {
  V <- 0.5 * t(x$D) %*% kronecker(X = exp_cov_inv_up, Y = exp_cov_inv_up) %*% x$D
  if (x$mean_structure) {
    V_m_cov <- matrix(data = 0, nrow = x$p_star_means, ncol = x$p_star_means)
    V_m_cov[1:x$p_star, 1:x$p_star] <- V
    V_m_cov[(x$p_star + 1):x$p_star_means, (x$p_star + 1):x$p_star_means] <- exp_cov_inv
    V <- V_m_cov
  }
  return(V)
}
