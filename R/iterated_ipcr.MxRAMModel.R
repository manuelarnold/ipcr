iterated_ipcr.MxRAMModel <- function(x, IPC, iteration_info, covariates, conv, max_it,
                                     linear_MxModel, ...) {

  # Preparations --------
  ## Model properties
  data_obs <- as.matrix(x$data$observed[, x$manifestVars, drop = FALSE])
  covariates_matrix <- as.matrix(covariates)
  covariates_design_matrix <- cbind(1, covariates_matrix)
  n <- x$data$numObs
  p <- length(x$manifestVars)
  p_unf = nrow(x$A$values)
  p_star = (p * (p + 1)) / 2
  p_star_means = p * (p + 3) / 2
  ms <- any(x$M$free)
  param_estimates<- x$output$estimate
  param_names <- names(param_estimates)
  q <- length(param_estimates)
  exp_cov <- OpenMx::mxGetExpected(model = x, component = "covariance")
  exp_cov_inv <- solve(exp_cov)
  if(ms) {exp_mean <- OpenMx::mxGetExpected(model = x, component = "means")}
  k <- NCOL(covariates)
  Ident <- diag(x = 1, nrow = p_unf)
  Dup <- lavaan::lav_matrix_duplication(n = p)
  indices_n <- seq_len(n)
  indices_p_star <- seq_len(p_star)
  indices_param <- seq_len(q)
  indices_p_star_p_means <- (p_star + 1):p_star_means


  ## RAM matrices
  F_RAM <- x$F$values
  A <- x$A$values
  S <- x$S$values
  m <- t(x$M$values)
  B <- solve(Ident - A)
  FB <- F_RAM %*% B
  E <- B %*% S %*% t(B)


  ## Jacobian matrix
  if (linear_MxModel) { # Analytic Jacobian matrix

    ### Derivative Matrices
    Zero <- matrix(0, nrow = p_unf, ncol = p_unf)
    A_deriv <- lapply(indices_param, function(x) {Zero})
    S_deriv <- A_deriv
    zero <- matrix(0, nrow = p_unf, ncol = 1)
    m_deriv <- lapply(indices_param, function(x) {zero})

    for (i in indices_param) {
      A_deriv[[i]][which(x$A$labels == param_names[i], arr.ind = TRUE)] <- 1
    }

    for (i in indices_param) {
      S_deriv[[i]][which(x$S$labels == param_names[i], arr.ind = TRUE)] <- 1
    }

    for (i in indices_param) {
      m_deriv[[i]][which(x$M$labels == param_names[i])] <- 1
    }

    ## Analytic Jacobian
    jac <- matrix(0, nrow = p_star_means, ncol = q)

    for (i in indices_param) {
      symm <- FB %*% A_deriv[[i]] %*% E %*% t(F_RAM)
      jac[seq_len(p_star), i] <- lavaan::lav_matrix_vech(symm + t(symm) + FB %*% S_deriv[[i]] %*% t(FB))
    }

    if (ms) {
      for (i in indices_param) {
        jac[indices_p_star_p_means, i] <- FB %*% A_deriv[[i]] %*% B %*% m +
          FB %*% m_deriv[[i]]
      }
    }

  } else { # Numeric Jacobian matrix
    jac <- OpenMx::omxManifestModelByParameterJacobian(model = x)
  }

  if (!ms) {jac <- jac[indices_p_star, , drop = FALSE]}



  # Initial IPC regression --------
  ## Individual deviations from the sample moments
  data_obs_c <- scale(x = data_obs, center = TRUE, scale = FALSE)
  mc <- matrix(data = apply(X = data_obs_c, MARGIN = 1,
                            FUN = function (x) {lavaan::lav_matrix_vech(x %*% t(x))}),
               nrow = n, ncol = p_star, byrow = TRUE)
  vech_cov <- matrix(data = rep(x = lavaan::lav_matrix_vech(exp_cov), times = n),
                     byrow = TRUE, nrow = n, ncol = p_star)
  md <- mc - vech_cov
  if (ms) {
    means <- matrix(data = rep(x = exp_mean, times = n), byrow = TRUE,
                    nrow = n, ncol = p)
    mean_dev <- data_obs - means
    md <- cbind(md, mean_dev)
  }

  ## Weight matrix V
  V <- 0.5 * t(Dup) %*% kronecker(X = exp_cov_inv, Y = exp_cov_inv) %*% Dup
  if (ms) {
    V_m_cov <- matrix(data = 0, nrow = p_star_means, ncol = p_star_means)
    V_m_cov[indices_p_star, indices_p_star] <- V
    V_m_cov[indices_p_star_p_means, indices_p_star_p_means] <- exp_cov_inv
    V <- V_m_cov
  }

  ## Initial IPCs
  W <- solve(t(jac) %*% V %*% jac) %*% t(jac) %*% V
  IPCs <- matrix(data = rep(x = param_estimates, times = n),
                 byrow = TRUE, nrow = n, ncol = q) +
    md %*% t(W)
  IPCs <- as.data.frame(IPCs)
  colnames(IPCs) <- param_names

  ## Initial IPC regression
  ipcr_data <- cbind(IPCs, covariates)
  param_names_ipcr <- paste0("IPCs_", gsub("\\(|\\)", "", param_names))
  IV <- paste(colnames(covariates), collapse = " + ")
  colnames(ipcr_data)[seq_len(q)] <- param_names_ipcr
  ipcr_list <- lapply(param_names_ipcr, FUN = function(x) {
    do.call(what = "lm",
            args = list(formula = paste(x, "~", IV), data = as.name("ipcr_data")))
  })
  names(ipcr_list) <- param_names



  # Start iterated IPC regression --------
  ## Storing objects for the updating procedure
  it_est <- matrix(sapply(X = ipcr_list, FUN = function(x) {coef(x)}),
                   nrow = 1, ncol = q * (k + 1))
  it_se <- matrix(sapply(X = ipcr_list, FUN = function(x) {sqrt(diag(vcov(x)))}),
                  nrow = 1, ncol = q * (k + 1))

  ## Center moment deviations at the covariate
  center_reg_list <- apply(X = data_obs, MARGIN = 2, FUN = function(y) {
    stats::lm(y ~ covariates_matrix)})
  data_centered <- data_obs -
    sapply(X = center_reg_list, FUN = function(x) {x$fitted.values})
  cent_md <- matrix(data = apply(X = data_centered, MARGIN = 1,
                                 FUN = function(x) {lavaan::lav_matrix_vech(x %*% t(x))}),
                    nrow = n, ncol = p_star, byrow = TRUE)
  if (ms) {
    cent_md <- cbind(cent_md, data_obs)
  }

  ## Calculate model fit
  if (iteration_info) {
    log_lik_individual <- rep(NA, times = n)
    for (i in indices_n) {
      param_estimates <- sapply(X = ipcr_list, FUN = function(x) {sum(coef(x) * covariates_design_matrix[i, ])})
      x <- OpenMx::omxSetParameters(model = x, labels = param_names,
                                    values = param_estimates)
      x <- suppressMessages(OpenMx::mxRun(model = x, useOptimizer = FALSE))
      data_individual <- t(data_obs[i, , drop = FALSE])
      sigma_individual <- OpenMx::mxGetExpected(model = x, component = "covariance")
      sigma_inv_individual <- solve(sigma_individual)
      if (ms) {
        mu_individual <- t(OpenMx::mxGetExpected(model = x, component = "means"))
        log_lik_individual[i] <- t(data_individual - mu_individual) %*% sigma_inv_individual %*%
          (data_individual - mu_individual) + log(det(sigma_individual))
      }  else {
        log_lik_individual[i] <- t(data_individual) %*% sigma_inv_individual %*%
          data_individual + log(det(sigma_individual))
      }
    }
    log_lik <- -0.5 * sum(log_lik_individual) + n * p * log(2 * pi)
  }

  ## Groups of observations with same covariate values in the covariate
  group <- transform(covariates, group_ID = as.numeric(interaction(covariates, drop = TRUE)))$group_ID

  ## Assign SEM parameters to the corresponding RAM matrices
  RAM_params <- rep(NA, times = q)
  RAM_params[which(param_names %in% x$A$labels)] <- "A"
  RAM_params[which(param_names %in% x$S$labels)] <- "S"
  RAM_params[which(param_names %in% x$M$labels)] <- "M"

  ## Get coordinates of the parameters in the corresponding RAM matrices
  RAM_coord <- list()
  for (i in indices_param) {
    if (RAM_params[i] == "A") {
      RAM_coord[[i]] <- which(x$A$labels == param_names[i], arr.ind = TRUE)
    }
    if (RAM_params[i] == "S") {
      RAM_coord[[i]] <- which(x$S$labels == param_names[i], arr.ind = TRUE)
    }
    if (RAM_params[i] == "M") {
      RAM_coord[[i]] <- which(t(x$M$labels) == param_names[i], arr.ind = TRUE)
    }
  }



  # Start the iteration process --------
  nr_iterations <- 0
  difference <- rep(x = conv + 1, times = NCOL(it_est))
  breakdown <- FALSE
  updated_IPCs <- matrix(NA, nrow = n, ncol = q)
  colnames(updated_IPCs) <- param_names
  unique_groups <- unique(group)
  cent_md_up <- cent_md

  ## while loop
  while(nr_iterations < max_it & isFALSE(all(abs(difference) < conv)) &
        breakdown == FALSE) {

    ### Update IPCs ###
    # Update the IPCs of individual or group i
    for (i in unique_groups) { # Start loop with index i
      ID_group <- which(group == i) # HERE, THIS IS WRONG
      n_group <- length(ID_group)
      IPC_pred <- covariates_design_matrix[group == i, , drop = FALSE][1, ]

      # Update RAM matrices
      for (j in indices_param){

        param_estimates[j] <- sum(coef(ipcr_list[[j]]) * IPC_pred)

        if (RAM_params[j] == "A") {
          A[RAM_coord[[j]]] <- param_estimates[j]
        }

        if (RAM_params[j] == "S") {
          S[RAM_coord[[j]]] <- param_estimates[j]
        }

        if (RAM_params[j] == "M") {
          m[RAM_coord[[j]]] <- param_estimates[j]
        }
      } # end loop with index j: updating paramaters

      # Update sample covariance
      B <- solve(Ident - A)
      FB <- F_RAM %*% B
      E <- B %*% S %*% t(B)
      exp_cov <- F_RAM %*% E %*% t(F_RAM)
      exp_cov_inv <- solve(exp_cov)


      # Update Jacobian matrix
      if (linear_MxModel) { # Analytic Jacobian matrix

        for (j in indices_param) {
          symm <- FB %*% A_deriv[[j]] %*% E %*% t(F_RAM)
          jac[indices_p_star, j] <- lavaan::lav_matrix_vech(symm + t(symm) + FB %*% S_deriv[[j]] %*% t(FB))
        }

        if (ms) {
          for (j in indices_param) {
            jac[indices_p_star_p_means, j] <- FB %*% A_deriv[[j]] %*% B %*% m +
              FB %*% m_deriv[[j]]
          }
        }


      } else { # Numeric Jacobian matrix
        x <- OpenMx::omxSetParameters(model = x, labels = param_names,
                              values = param_estimates)
        x <- suppressMessages(OpenMx::mxRun(model = x, useOptimizer = FALSE))
        jac <- OpenMx::omxManifestModelByParameterJacobian(model = x)
      }

      if (!ms) {
        jac <- jac[indices_p_star, , drop = FALSE]
      }

      # Update weight matrix and W matrix
      V <- 0.5 * t(Dup) %*% kronecker(X = exp_cov_inv, Y = exp_cov_inv) %*% Dup
      if (ms) {
        V_m_cov <- matrix(data = 0, nrow = p_star_means, ncol = p_star_means)
        V_m_cov[indices_p_star, indices_p_star] <- V
        V_m_cov[indices_p_star_p_means, indices_p_star_p_means] <- exp_cov_inv
        V <- V_m_cov
      }
      W <- solve(t(jac) %*% V %*% jac) %*% t(jac) %*% V

      # Update the centered contributions to the sample moments
      cent_md_up[ID_group, indices_p_star] <- cent_md[ID_group, indices_p_star] -
        matrix(rep(x = lavaan::lav_matrix_vech(exp_cov), times = n_group), byrow = TRUE,
               nrow = n_group, ncol = p_star)
      if (ms) {
        ### !!!! Orientation could be wrong
        exp_means <- FB %*% m
        means_matrix <- matrix(rep(x = exp_means, times = n), byrow = TRUE,
                               nrow = n, ncol = p)
        means_dev <- data_obs - means_matrix
        cent_md_up[ID_group, indices_p_star_p_means] <- means_dev[ID_group, ]
      }
      cent_md_up <- as.matrix(cent_md_up)

      updated_IPCs[ID_group, ] <- cent_md_up[ID_group, ] %*% t(W) +
        matrix(rep(x = param_estimates, times = n_group), byrow = TRUE,
               nrow = n_group, ncol = q)
    } # end loop with index i: Find observations with identical covariates


    if (class(updated_IPCs)[1] == "try-error") {
      breakdown <- TRUE
    }

    if (breakdown == FALSE) {

      # Estimate updated IPC regression parameter
      ipcr_data <- cbind(updated_IPCs, covariates)
      colnames(ipcr_data)[indices_param] <- param_names_ipcr
      for (j in indices_param) {
        ipcr_list[[j]] <- do.call(what = "lm",
                                  args = list(formula = paste(param_names_ipcr[j], "~", IV),
                                              data = as.name("ipcr_data")))
      }


      ## Store results
      it_est <- rbind(it_est, c(sapply(X = ipcr_list,
                                       FUN = function(x) {coef(x)})))
      it_se <- rbind(it_se, c(sapply(X = ipcr_list,
                                     FUN = function(x) {sqrt(diag(vcov(x)))})))

      ## Calculate model fit
      if (iteration_info) {
        for (i in indices_n) {
          param_estimates <- sapply(X = ipcr_list, FUN = function(x) {sum(coef(x) * covariates_design_matrix[i, ])})
          x <- OpenMx::omxSetParameters(model = x, labels = param_names,
                                        values = param_estimates)
          x <- suppressMessages(OpenMx::mxRun(model = x, useOptimizer = FALSE))
          data_individual <- t(data_obs[i, , drop = FALSE])
          sigma_individual <- OpenMx::mxGetExpected(model = x, component = "covariance")
          sigma_inv_individual <- solve(sigma_individual)
          if (ms) {
            mu_individual <- t(OpenMx::mxGetExpected(model = x, component = "means"))
            log_lik_individual[i] <- t(data_individual - mu_individual) %*% sigma_inv_individual %*%
              (data_individual - mu_individual) + log(det(sigma_individual))
          }  else {
            log_lik_individual[i] <- t(data_individual) %*% sigma_inv_individual %*%
              data_individual + log(det(sigma_individual))
          }
        }
        log_lik <- c(log_lik, -0.5 * sum(log_lik_individual) + n * p * log(2 * pi))
      }

      # Covergence criteria
      difference <- it_est[nrow(it_est), ] -
        it_est[nrow(it_est) - 1, ]
      nr_iterations <- nr_iterations + 1
      cat("Iteration:", nr_iterations, "\n")

    }
  }



  # Catch errors --------
  if (breakdown == TRUE) {
    cat("The iterated IPC regression aborted prematurely.")
    IPC$info$iterated_status <- "Iterated IPC regression aborted prematurely."
  }

  if(breakdown == FALSE & nr_iterations == max_it &
     isFALSE(all(abs(difference) <= conv))) {
    cat(paste("The iterated IPC regression algorithm did not converge after", max_it, "iterations. Consider to increase the maximum number of iterations."))
    IPC$info$iterated_status <- paste("Iterated IPC regression reached the maximum number of", max_it, "iterations without converging.")
  }

  if(breakdown == FALSE & nr_iterations < max_it &
     isFALSE(all(abs(difference) > conv))) {
    cat("Iterated IPC regression converged.")
    IPC$info$iterated_status <- paste("Iterated IPC regression converged succesfully after", nr_iterations, "iterations.")
    IPC$IPCs <- as.data.frame(updated_IPCs)
    IPC$regression_list <- ipcr_list
  }



  # Prepare output --------
  if(iteration_info) {IPC$iteration_matrix <- cbind(log_lik, it_est)}

  IPC

}
