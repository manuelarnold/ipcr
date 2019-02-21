#' @title Calculate the Analyitc Jacobian Matrix
#' @description Internal function that anayltically calculates the Jacobian matrix.
#' @param x list with information about the model.
#' @param deriv list with partial derivatives of the RAM matrices.
#' @return Jacobian matrix.
#' @export

get_analytic_jac <- function(x) {

  F_RAM <- x$fit$F$values
  m <- t(x$fit$M$values)
  B <- solve(diag(x = 1, nrow = x$p_unf) - x$fit$A$values)
  E <- B %*% x$fit$S$values %*% t(B)
  FB <- F_RAM %*% B

  ###
  jac <- matrix(0, nrow = x$p_star_means, ncol = x$q)

  for (i in 1:x$q) {
    symm <- FB %*% x$RAM_deriv$A[[i]] %*% E %*% t(F_RAM)
    jac[1:x$p_star, i] <- vech(symm + t(symm) + FB %*% x$RAM_deriv$S[[i]] %*% t(FB))
  }


  for (i in 1:x$q) {
    jac[(x$p_star+1):x$p_star_means, i] <- FB %*% x$RAM_deriv$A[[i]] %*% B %*% m +
      FB %*% x$RAM_deriv$m[[i]]
  }

  return(jac)

}
