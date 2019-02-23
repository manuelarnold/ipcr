#' @title Calculate the Analyitc Jacobian Matrix
#' @description Internal function that anayltically calculates the Jacobian matrix.
#' @param x list with information about the model.
#' @param deriv list with partial derivatives of the RAM matrices.
#' @return Jacobian matrix.
#' @export

get_analytic_jac_up <- function(x, A_up, S_up, M_up, F_up, B_up = B_up) {

  E <- B_up %*% S_up %*% t(B_up)
  FB <- F_up %*% B_up

  ###
  jac <- matrix(0, nrow = x$p_star_means, ncol = x$q)

  for (i in 1:x$q) {
    symm <- FB %*% x$RAM_deriv$A[[i]] %*% E %*% t(F_up)
    jac[1:x$p_star, i] <- vech(symm + t(symm) + FB %*% x$RAM_deriv$S[[i]] %*% t(FB))
  }

  for (i in 1:x$q) {
    jac[(x$p_star+1):x$p_star_means, i] <- FB %*% x$RAM_deriv$A[[i]] %*% B_up %*% M_up +
      FB %*% x$RAM_deriv$m[[i]]
  }

  return(jac)

}
