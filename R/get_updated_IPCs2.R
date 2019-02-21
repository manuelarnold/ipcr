#' @title Update the Individual Contributions to the Parameter Estimates
#' @description Internal function that updates the individual contributions to the parameter estimates.
#' @param x list with information about the model.
#' @return Matrix of initial individual parameter contributions.
#' @export

get_updated_IPCs2 <- function(x) {

  # Storing object for the updated deviations from the sample moments
  md <- x$cent_md

  # Storing object for the updated IPCs
  updated_IPCs <- matrix(NA, nrow = x$N, ncol = x$q)

  # Update the IPCs of individual or group i
  for (i in unique(x$group[, 2])) {
    ID_group <- x$group[which(x$group[, 2] == i), 1]
    n_group <- length(ID_group)
    IPC_pred <- c(1, x$covariates[which(x$group[, 2] == i)[1],])

    for (j in 1:x$q) {
      x$param_estimates[j] <- sum(coef(x$IPC_reg[[j]]) * IPC_pred)
    }
    x$fit <- omxSetParameters(model = x$fit, labels = x$param_names,
                              values = x$param_estimates)
    x$fit <- suppressMessages(mxRun(model = x$fit, useOptimizer = FALSE))

    # Update expected covariance and inverse
    x$exp_cov <- mxGetExpected(model = x$fit, component = "covariance")
    x$exp_cov_inv <- solve(x$exp_cov)

    jacobian <- get_analytic_jac(x = x)
    if (x$mean_structure == FALSE) {
      jacobian <- jacobian[1:x$p_star, ]
    }
    V <- get_weight_matrix(x = x)
    W <- solve(t(jacobian) %*% V %*% jacobian) %*% t(jacobian) %*% V

    md <- as.data.frame(md)
    md[ID_group, 1:x$p_star] <- x$cent_md[ID_group, 1:x$p_star] -
      matrix(rep(x = matrixcalc::vech(x$exp_cov), times = n_group), byrow = TRUE,
             nrow = n_group, ncol = x$p_star)
    if (x$mean_structure) {
      exp_means <- mxGetExpected(model = x$fit, component = "means")
      means_matrix <- matrix(rep(x = exp_means, times = x$N), byrow = TRUE,
                             nrow = x$N, ncol = x$p)
      means_dev <- x$data_obs - means_matrix
      md[ID_group, (x$p_star + 1):x$p_star_means] <- means_dev[ID_group, ]
    }
    md <- as.matrix(md)

    updated_IPCs[ID_group, ] <- md[ID_group, ] %*% t(W) +
      matrix(rep(x = x$param_estimates, times = n_group), byrow = TRUE,
             nrow = n_group, ncol = x$q)
  }

  colnames(updated_IPCs) <- x$param_names
  return(updated_IPCs)
}
