#' @title Prepare Input for the Analytic Gradient
#' @description Internal function that calculates the partial derivative of the RAM matrices with respect to the parameters.
#' @param x list with information about the model.
#' @return List with partial derivatives.
#' @export

get_analytic_jac_input <- function(x) {

  Zero <- matrix(0, nrow = x$p_unf, ncol = x$p_unf)
  A <- lapply(seq_len(x$q), function(x) {Zero})
  S <- A
  zero <- matrix(0, nrow = x$p_unf, ncol = 1)
  m <- lapply(seq_len(x$q), function(x) {zero})


  for (i in 1:x$q) {
    A[[i]][which(x$fit$A$labels == x$param_names[i], arr.ind = TRUE)] <- 1
  }

  for (i in 1:x$q) {
    S[[i]][which(x$fit$S$labels == x$param_names[i], arr.ind = TRUE)] <- 1
  }

  for (i in 1:x$q) {
    m[[i]][which(x$fit$M$labels == x$param_names[i])] <- 1
  }

  RAM_deriv <- list(A = A, S = S, m = m)

}
