#' @title Calculate the Weight Matrix
#' @description Internal function that calculates the estimated weight matrix.
#' @param x list with information about the model.
#' @return Estimated weight matrix.
#' @export

get_weight_matrix <- function(x) {
  V <- 0.5 * t(x$D) %*% kronecker(X = x$exp_cov_inv, Y = x$exp_cov_inv) %*% x$D
  if (x$mean_structure) {
    V_m_cov <- matrix(data = 0, nrow = x$p_star_means, ncol = x$p_star_means)
    V_m_cov[1:x$p_star, 1:x$p_star] <- V
    V_m_cov[(x$p_star + 1):x$p_star_means, (x$p_star + 1):x$p_star_means] <- x$exp_cov_inv
    V <- V_m_cov
  }
  return(V)
}
