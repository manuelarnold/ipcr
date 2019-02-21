#' @title Calculate Initial Individual Parameter Contributions
#' @description Internal function that calculates the initial parameter contributions.
#' @param x list with information about the model.
#' @return Matrix of initial individual parameter contributions.
#' @export

get_initial_IPCs2 <- function(x) {
  jacobian <- get_analytic_jac(x = x)
  if (x$mean_structure == FALSE) {
    jacobian <- jacobian[1:x$p_star, ]
  }
  V <- get_weight_matrix(x = x)
  W <- solve(t(jacobian) %*% V %*% jacobian) %*% t(jacobian) %*% V
  IPCs <- data.frame(matrix(data = rep(x = x$param_estimates, times = x$N),
                            byrow = TRUE, nrow = x$N, ncol = x$q) +
                       x$moment_deviations %*% t(W))
  colnames(IPCs) <- x$param_names
  return(IPCs)
}
