#' @title Update the Individual Contributions to the Parameter Estimates
#' @description Internal function that updates the individual contributions to the parameter estimates.
#' @param x list with information about the model.
#' @param A_up object.
#' @param S_up object.
#' @param M_up object.
#' @param F_up object.
#' @param Ident object.
#' @param RAM_params object.
#' @param RAM_coord object.
#' @return Matrix of initial individual parameter contributions.
#' @export

get_updated_IPCs2 <- function(x, A_up, S_up, M_up, F_up, Ident, RAM_params,
                              RAM_coord) {

  # Storing objects
  updated_IPCs <- matrix(NA, nrow = x$N, ncol = x$q)
  md <- x$cent_md
  param_estimates <- rep(NA, x$q)

  # Update the IPCs of individual or group i
  for (i in unique(x$group[, 2])) {
    ID_group <- x$group[which(x$group[, 2] == i), 1]
    n_group <- length(ID_group)
    IPC_pred <- c(1, x$covariates[which(x$group[, 2] == i)[1],])



    # Update RAM matrices
    for (j in 1:x$q){

      param_estimates[j] <- sum(coef(x$IPC_reg[[j]]) * IPC_pred)

      if (RAM_params[j] == "A") {
        A_up[RAM_coord[[j]]] <- sum(coef(x$IPC_reg[[j]]) * IPC_pred)
      }

      if (RAM_params[j] == "S") {
        S_up[RAM_coord[[j]]] <- sum(coef(x$IPC_reg[[j]]) * IPC_pred)
      }

      if (RAM_params[j] == "M") {
        M_up[RAM_coord[[j]]] <- sum(coef(x$IPC_reg[[j]]) * IPC_pred)
      }
    }

    # Update sample covariance
    B_up <- solve(Ident - A_up)
    exp_cov_up <- F_up %*% B_up %*% S_up %*% t(B_up) %*% t(F_up)
    exp_cov_inv_up <- solve(exp_cov_up)

    # Update weight matrix and jacobian matrix
    jacobian <- get_analytic_jac_up(x = x, A_up = A_up, S_up = S_up,
                                    M_up = M_up, F_up = F_up, B_up = B_up)
    if (x$mean_structure == FALSE) {
      jacobian <- jacobian[1:x$p_star, ]
    }
    V <- get_weight_matrix2(x = x, exp_cov_inv_up = exp_cov_inv_up)
    W <- solve(t(jacobian) %*% V %*% jacobian) %*% t(jacobian) %*% V

    # Update the centered contributions to the sample moments
    md <- as.data.frame(md)
    md[ID_group, 1:x$p_star] <- x$cent_md[ID_group, 1:x$p_star] -
      matrix(rep(x = matrixcalc::vech(exp_cov_up), times = n_group), byrow = TRUE,
             nrow = n_group, ncol = x$p_star)
    if (x$mean_structure) {
      ### !!!! Orientation could be wrong
      exp_means <- F_up %*% B_up %*% M_up
      means_matrix <- matrix(rep(x = exp_means, times = x$N), byrow = TRUE,
                             nrow = x$N, ncol = x$p)
      means_dev <- x$data_obs - means_matrix
      md[ID_group, (x$p_star + 1):x$p_star_means] <- means_dev[ID_group, ]
    }
    md <- as.matrix(md)

    updated_IPCs[ID_group, ] <- md[ID_group, ] %*% t(W) +
      matrix(rep(x = param_estimates, times = n_group), byrow = TRUE,
             nrow = n_group, ncol = x$q)
  }

  colnames(updated_IPCs) <- x$param_names
  return(updated_IPCs)
}
