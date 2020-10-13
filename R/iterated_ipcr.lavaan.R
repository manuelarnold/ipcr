iterated_ipcr.lavaan <- function(x, IPC, iteration_info, covariates, conv, max_it, linear, ...) {

  # Preparations --------
  ## Model properties
  data_obs <- lavaan::lavInspect(object = x, what = "data")
  covariates_matrix <- as.matrix(covariates)
  covariates_design_matrix <- cbind(1, covariates_matrix)
  n <- x@SampleStats@ntotal
  p <- x@pta$nvar[[1]]
  p_star <- (p * (p + 1)) / 2
  p_star_means = p * (p + 3) / 2
  ms <- x@Model@meanstructure
  param_estimates <- coef(x)
  param_names <- names(param_estimates)
  q <- length(param_estimates)
  exp_cov <- x@Fit@Sigma.hat[[1]]
  if(ms) {exp_mean <- x@Fit@Mu.hat[[1]]}
  k <- NCOL(covariates)
  indices_n <- seq_len(n)
  #  indices_p_star <- seq_len(p_star)
  #  indices_p_star_p_means <- (p_star + 1):p_star_means
  if (ms) {
    indices_p_star <- (p + 1):(p_star_means)
    indices_p_star_p_means <- 1:p
  } else {
    indices_p_star <- seq_len(p_star)
  }
  indices_param <- seq_len(q)

  ## Prepare stuff for lavaan model syntax
  ParTable <- as.data.frame(x@ParTable)
  ParTable[ParTable$op == "~1", "rhs"] <- 1
  ParTable[ParTable$op == "~1", "op"] <- "~"
  fixed_params <- ParTable$free == 0
  values <- rep(NA, times = NROW(ParTable))
  values[fixed_params] <- paste0(ParTable[fixed_params, "est"], "*")



  # Initial IPC regression --------
  ## Calculate initial IPCs
  scores <- lavaan::estfun.lavaan(x)
  bread_matrix <- bread.lavaan(x)
  IPCs <- data.frame(matrix(param_estimates, nrow = n, ncol = q, byrow = TRUE) +
                       scores %*% t(bread_matrix))
  colnames(IPCs) <- param_names


  ## Initial IPC regression
  ipcr_data <- cbind(IPCs, covariates_matrix)
  param_names_ipcr <- paste0("IPCs_", gsub("\\(|\\)", "", param_names))
  param_names_ipcr <- gsub(pattern = "~~",replacement = ".WITH.", x = param_names_ipcr)
  param_names_ipcr <- gsub(pattern = "~",replacement = ".ON.", x = param_names_ipcr)
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

  ## Calculate model fit
  if (iteration_info) {
    log_lik_individual <- rep(NA, times = n)
    for (i in indices_n) {
      param_estimates <- sapply(X = ipcr_list, FUN = function(x) {sum(coef(x) * covariates_design_matrix[i, ])})
      values[!fixed_params] <- paste0("start(", param_estimates, ")*")
      model_syntax <- paste0(ParTable$lhs, ParTable$op, values, ParTable$rhs, collapse = "\n")
      x <- lavaan::lavaan(model = model_syntax, data = data_obs, do.fit = FALSE)
      data_individual <- t(data_obs[i, , drop = FALSE])
      sigma_individual <- x@Fit@Sigma.hat[[1]]
      sigma_inv_individual <- solve(sigma_individual)
      if (ms) {
        mu_individual <- x@Fit@Mu.hat[[1]]
        log_lik_individual[i] <- t(data_individual - mu_individual) %*% sigma_inv_individual %*%
          (data_individual - mu_individual) + log(det(sigma_individual))
      }  else {
        log_lik_individual[i] <- t(data_individual) %*% sigma_inv_individual %*%
          data_individual + log(det(sigma_individual))
      }
    }
    log_lik <- -0.5 * sum(log_lik_individual) + n * p * log(2 * pi)
  }

  ## Center moment deviations at the covariates
  center_reg_list <- apply(X = data_obs, MARGIN = 2, FUN = function(y) {
    stats::lm(y ~ covariates_matrix)})
  data_centered <- data_obs -
    sapply(X = center_reg_list, FUN = function(x) {x$fitted.values})
  cent_md <- matrix(data = apply(X = data_centered, MARGIN = 1,
                                 FUN = function(x) {lavaan::lav_matrix_vech(x %*% t(x))}),
                    nrow = n, ncol = p_star, byrow = TRUE)
  if (ms) {
    cent_md <- cbind(data_obs, cent_md)
  }

  ## Groups of observations with same covariate values in the covariate
  group <- transform(covariates, group_ID = as.numeric(interaction(covariates, drop = TRUE)))$group_ID

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

    ### Update the IPCs of individual or group i
    for (i in unique_groups) { # Start loop with index i
      ID_group <- which(group == i) # HERE, THIS IS WRONG
      n_group <- length(ID_group)
      IPC_pred <- covariates_design_matrix[group == i, , drop = FALSE][1, ]

      ### Update Parameters
      param_estimates <- sapply(X = ipcr_list, FUN = function(x) {sum(coef(x) * IPC_pred)})

      ### Update model synatx
      values[!fixed_params] <- paste0("start(", param_estimates, ")*")
      model_syntax <- paste0(ParTable$lhs, ParTable$op, values, ParTable$rhs, collapse = "\n")

      ### Update lavaan model
      x <- lavaan::lavaan(model = model_syntax, data = data_obs, do.fit = FALSE)

      ### Updated Jacobian matrix
      jac <- lavaan:::computeDelta(x@Model)[[1]]

      ### Updated weight matrix
      V <- lavaan::lavInspect(object = x, what = "wls.v")

      ### Update transformation matrix W
      W <- solve(t(jac) %*% V %*% jac) %*% t(jac) %*% V

      ### Update the centered contributions to the sample moments
      cent_md_up[ID_group, indices_p_star] <- cent_md[ID_group, indices_p_star] -
        matrix(rep(x = lavaan::lav_matrix_vech(x@Fit@Sigma.hat[[1]]), times = n_group), byrow = TRUE,
               nrow = n_group, ncol = p_star)
      if (ms) {
        means_matrix <- matrix(rep(x@Fit@Mu.hat[[1]], times = n), byrow = TRUE,
                               nrow = n, ncol = p)
        means_dev <- data_obs - means_matrix
        cent_md_up[ID_group, indices_p_star_p_means] <- means_dev[ID_group, ]
      }
      cent_md_up <- as.matrix(cent_md_up)

      ### Calculate the updated the IPCs of the selected group
      updated_IPCs[ID_group, ] <- cent_md_up[ID_group, ] %*% t(W) +
        matrix(rep(x = param_estimates, times = n_group), byrow = TRUE,
               nrow = n_group, ncol = q)
    } # end loop with index i: Find observations with identical covariates

    ## Check for failed attempt
    if (class(updated_IPCs)[1] == "try-error") {
      breakdown <- TRUE
    }

    if (breakdown == FALSE) {

      ## Estimate updated IPC regression parameter
      ipcr_data <- cbind(updated_IPCs, covariates)
      colnames(ipcr_data)[indices_param] <- param_names_ipcr
      updated_IPCs <- as.data.frame(updated_IPCs)
      ipcr_list <- lapply(param_names_ipcr, FUN = function(x) {
        do.call(what = "lm",
                args = list(formula = paste(x, "~", IV), data = as.name("ipcr_data")))
      })
      names(ipcr_list) <- param_names

      ## Store results
      it_est <- rbind(it_est, c(sapply(X = ipcr_list,
                                       FUN = function(x) {coef(x)})))
      it_se <- rbind(it_se, c(sapply(X = ipcr_list,
                                     FUN = function(x) {sqrt(diag(vcov(x)))})))

      ## Calculate model fit
      if (iteration_info) {
        for (i in indices_n) {
          param_estimates <- sapply(X = ipcr_list, FUN = function(x) {sum(coef(x) * covariates_design_matrix[i, ])})
          values[!fixed_params] <- paste0("start(", param_estimates, ")*")
          model_syntax <- paste0(ParTable$lhs, ParTable$op, values, ParTable$rhs, collapse = "\n")
          x <- lavaan::lavaan(model = model_syntax, data = data_obs, do.fit = FALSE)
          data_individual <- t(data_obs[i, , drop = FALSE])
          sigma_individual <- x@Fit@Sigma.hat[[1]]
          sigma_inv_individual <- solve(sigma_individual)
          if (ms) {
            mu_individual <- x@Fit@Mu.hat[[1]]
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
