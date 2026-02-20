#' @noRd
ipcr_it.MxRAMModel <- function(
    fit, predictors, analytic, conv = NULL,
    learning_rate = c(seq(from = 1, to = 0.1, by = -0.1), 0.05,
                      0.01, 0.005, 0.001),
    max_it = 350, iteration_info = FALSE) {

#  ipcr_it.MxRAMModel <- function(x, IPC, iteration_info, covariates, conv,
#                                 max_it, linear_MxModel, ...) {

  # Preprocess predictors ----
  if (is.null(predictors)) {
    stop("No predictors were provided.")
  }

  predictors <- as.data.frame(predictors)

  if (is.null(names(predictors)) | any(is.na(names(predictors)))) {
    warning("Some predictor are not named. Renaming all predictors using the
            order of the data.frame.")
    pred_names <- paste0("predictor", seq_len(NCOL(predictors)))
    colnames(predictors) <- pred_names
  } else {
    pred_names <- names(predictors)
  }



  # Preparations --------
  ## Model properties
  data_obs <- as.matrix(fit$data$observed[, fit$manifestVars, drop = FALSE])
  predictors_matrix <- as.matrix(predictors)
  predictors_design_matrix <- cbind(1, predictors_matrix)
  n <- fit$data$numObs
  p <- length(fit$manifestVars)
  p_unf = nrow(fit$A$values)
  p_star = (p * (p + 1)) / 2
  p_star_means = p * (p + 3) / 2
  ms <- any(fit$M$free)
  param_estimates<- fit$output$estimate
  param_names <- names(param_estimates)
  q <- length(param_estimates)
  exp_cov <- OpenMx::mxGetExpected(model = fit, component = "covariance")
  exp_cov_inv <- solve(exp_cov)
  if(ms) {exp_mean <- OpenMx::mxGetExpected(model = fit, component = "means")}
  k <- NCOL(predictors)
  Ident <- diag(1, nrow = p_unf)
  Dup <- lavaan::lav_matrix_duplication(n = p)
  indices_n <- seq_len(n)
  indices_p_star <- seq_len(p_star)
  indices_param <- seq_len(q)
  indices_p_star_p_means <- (p_star + 1):p_star_means


  ## RAM matrices
  F_RAM <- fit$F$values
  A <- fit$A$values
  S <- fit$S$values
  m <- t(fit$M$values)
  B <- solve(Ident - A)
  FB <- F_RAM %*% B
  E <- B %*% S %*% t(B)


  ## Jacobian matrix
  if (analytic) { # Analytic Jacobian matrix

    ### Derivative Matrices
    Zero <- matrix(0, nrow = p_unf, ncol = p_unf)
    A_deriv <- lapply(indices_param, function(x) {Zero})
    S_deriv <- A_deriv
    zero <- matrix(0, nrow = p_unf, ncol = 1)
    m_deriv <- lapply(indices_param, function(x) {zero})

    for (i in indices_param) {
      A_deriv[[i]][which(fit$A$labels == param_names[i], arr.ind = TRUE)] <- 1
    }

    for (i in indices_param) {
      S_deriv[[i]][which(fit$S$labels == param_names[i], arr.ind = TRUE)] <- 1
    }

    for (i in indices_param) {
      m_deriv[[i]][which(fit$M$labels == param_names[i])] <- 1
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
    jac <- OpenMx::omxManifestModelByParameterJacobian(model = fit)
  }

  if (!ms) {jac <- jac[indices_p_star, , drop = FALSE]}



  # Initial IPC regression --------
  ## Individual deviations from the sample moments
  data_obs_c <- scale(data_obs, center = TRUE, scale = FALSE)
  mc <- matrix(data = apply(X = data_obs_c, MARGIN = 1,
                            FUN = function (x) {lavaan::lav_matrix_vech(x %*% t(x))}),
               nrow = n, ncol = p_star, byrow = TRUE)
  vech_cov <- matrix(data = rep(lavaan::lav_matrix_vech(exp_cov), times = n),
                     byrow = TRUE, nrow = n, ncol = p_star)
  md <- mc - vech_cov
  if (ms) {
    means <- matrix(data = rep(exp_mean, times = n), byrow = TRUE,
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
  IPCs <- matrix(data = rep(param_estimates, times = n), byrow = TRUE, nrow = n,
                 ncol = q) + md %*% t(W)
  colnames(IPCs) <- param_names

  ## Initial IPC regression
  mlm <- lm(IPCs ~ ., data = predictors)



  # Start iterated IPC regression --------
  ## Storing objects for the updating procedure
  it_est <- matrix(coef(mlm), nrow = 1, ncol = q * (k + 1))
  it_se <- matrix(sqrt(diag(vcov(mlm))), nrow = 1, ncol = q * (k + 1))

  ## Center moment deviations at the covariate
  data_centered <- as.matrix(resid(lm(data_obs ~ ., data = predictors)))
  cent_md <- matrix(data = apply(X = data_centered, MARGIN = 1,
                                 FUN = function(x) {lavaan::lav_matrix_vech(x %*% t(x))}),
                    nrow = n, ncol = p_star, byrow = TRUE)
  if (ms) {
    cent_md <- cbind(cent_md, data_obs)
  }

  ## Calculate model fit
  if (iteration_info) {
    log_lik_individual <- rep(NA, times = n)
    param_estimates_ind <- predict(mlm)
    for (i in indices_n) {
      fit <- OpenMx::omxSetParameters(model = fit, labels = param_names,
                                    values = param_estimates_ind[i])
      fit <- suppressMessages(OpenMx::mxRun(model = fit, useOptimizer = FALSE))
      data_individual <- t(data_obs[i, , drop = FALSE])
      sigma_individual <- OpenMx::mxGetExpected(model = fit, component = "covariance")
      sigma_inv_individual <- solve(sigma_individual)
      if (ms) {
        mu_individual <- t(OpenMx::mxGetExpected(model = fit, component = "means"))
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
  group <- transform(predictors, group_ID = as.numeric(interaction(predictors, drop = TRUE)))$group_ID

  ###
  ### HERE HERE HERE
  ###


  ## Assign SEM parameters to the corresponding RAM matrices
  RAM_params <- rep(NA, times = q)
  RAM_params[which(param_names %in% fit$A$labels)] <- "A"
  RAM_params[which(param_names %in% fit$S$labels)] <- "S"
  RAM_params[which(param_names %in% fit$M$labels)] <- "M"

  ## Get coordinates of the parameters in the corresponding RAM matrices
  RAM_coord <- list()
  for (i in indices_param) {
    if (RAM_params[i] == "A") {
      RAM_coord[[i]] <- which(fit$A$labels == param_names[i], arr.ind = TRUE)
    }
    if (RAM_params[i] == "S") {
      RAM_coord[[i]] <- which(fit$S$labels == param_names[i], arr.ind = TRUE)
    }
    if (RAM_params[i] == "M") {
      RAM_coord[[i]] <- which(t(fit$M$labels) == param_names[i], arr.ind = TRUE)
    }
  }



  # Start the iteration process --------
  nr_iterations <- 0
  difference <- rep(conv + 1, times = NCOL(it_est))
  updated_IPCs <- matrix(NA, nrow = n, ncol = q)
  colnames(updated_IPCs) <- param_names
  unique_groups <- unique(group)
  cent_md_up <- cent_md

  ## while loop
  while(nr_iterations < max_it & isFALSE(all(abs(difference) < conv))) {

    ### Try to update the IPCs of individuals and/or groups
    updated_IPCs <- try(expr = {

      for (i in unique_groups) { # Start loop with index i

        ID_group <- which(group == i)
        n_group <- length(ID_group)
        IPC_pred <- predictors_design_matrix[group == i, , drop = FALSE][1, ]

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
        if (analytic) { # Analytic Jacobian matrix

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
          fit <- OpenMx::omxSetParameters(model = fit, labels = param_names,
                                        values = param_estimates)
          fit <- suppressMessages(OpenMx::mxRun(model = fit, useOptimizer = FALSE))
          jac <- OpenMx::omxManifestModelByParameterJacobian(model = fit)
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
          matrix(rep(lavaan::lav_matrix_vech(exp_cov), times = n_group), byrow = TRUE,
                 nrow = n_group, ncol = p_star)
        if (ms) {
          ### !!!! Orientation could be wrong
          exp_means <- FB %*% m
          means_matrix <- matrix(rep(exp_means, times = n), byrow = TRUE,
                                 nrow = n, ncol = p)
          means_dev <- data_obs - means_matrix
          cent_md_up[ID_group, indices_p_star_p_means] <- means_dev[ID_group, ]
        }
        cent_md_up <- as.matrix(cent_md_up)

        updated_IPCs[ID_group, ] <- cent_md_up[ID_group, ] %*% t(W) +
          matrix(rep(param_estimates, times = n_group), byrow = TRUE,
                 nrow = n_group, ncol = q)
      } # end loop with index i: Find observations with identical predictors

      updated_IPCs
    }, # end expr of try()

    outFile = stop("Iterated IPC regression aborted prematurely.\n", call. = FALSE)

    )
    # end try


    # Estimate updated IPC regression parameter
    ipcr_data <- cbind(updated_IPCs, predictors)
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
      param_estimates_ind <- predict(mlm)
      for (i in indices_n) {
        fit <- OpenMx::omxSetParameters(model = fit, labels = param_names,
                                      values = param_estimates_ind[i])
        fit <- suppressMessages(OpenMx::mxRun(model = fit, useOptimizer = FALSE))
        data_individual <- t(data_obs[i, , drop = FALSE])
        sigma_individual <- OpenMx::mxGetExpected(model = fit, component = "covariance")
        sigma_inv_individual <- solve(sigma_individual)
        if (ms) {
          mu_individual <- t(OpenMx::mxGetExpected(model = fit, component = "means"))
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



  # Catch errors --------
  if(nr_iterations == max_it & isFALSE(all(abs(difference) <= conv))) {
    warning("The iterated IPC regression algorithm did not converge after ", max_it,
            " iterations. Consider to increase the maximum number of iterations.",
            call. = FALSE)
    IPC$info$iterated_status <- paste("Iterated IPC regression reached the maximum number of", max_it, "iterations without converging.")
  }

  if(nr_iterations < max_it & isFALSE(all(abs(difference) > conv))) {
    cat("Iterated IPC regression converged.")
    IPC$info$iterated_status <- paste("Iterated IPC regression converged succesfully after", nr_iterations, "iterations.")
    IPC$IPCs <- as.data.frame(updated_IPCs)
    IPC$regression_list <- ipcr_list
  }



  # Prepare output --------
  if(iteration_info) {IPC$iteration_matrix <- cbind(log_lik, it_est)}

  IPC

}
