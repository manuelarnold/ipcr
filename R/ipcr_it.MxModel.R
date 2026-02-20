#' @noRd
ipcr_it.MxRAMModel <- function(
    fit, predictors, analytic, conv = NULL,
    learning_rate = c(seq(from = 1, to = 0.1, by = -0.1), 0.05,
                      0.01, 0.005, 0.001),
    max_it = 350, iteration_info = FALSE) {

  check_arguments_ipcr_it_MxModel(fit = fit, predictors = predictors,
                                  analytic = analytic, conv = conv,
                                  learning_rate = learning_rate,
                                  max_it = max_it,
                                  iteration_info = iteration_info)

  if (analytic) {
    warning("There is no analytical method for models with the MxModel class
            (excluding models with the MxRAMModel class). Iterated IPCR will be
            performed numerically.")
    analytic <- FALSE
  }

  # Preparations --------
  ## Model properties
  manifest_vars_names <- get_manifest_vars_from_MxModel(fit = fit)
  dat <- as.matrix(fit$data$observed[, manifest_vars_names, drop = FALSE])
  dat_c <- scale(dat, center = TRUE, scale = FALSE)
  predictors_design_matrix_df <- cbind(1, predictors)
  n <- fit$data$numObs
  p <- length(manifest_vars_names)
  p_star = (p * (p + 1)) / 2
  p_star_means = p * (p + 3) / 2
  ms <- !is.null(fit$expectation$means)
  param_estimates <- fit$output$estimate
  param_names <- names(param_estimates)
  q <- length(param_estimates)
  k <- NCOL(predictors)
  #Dup <- lavaan::lav_matrix_duplication(n = p)
  #indices_n <- seq_len(n)
  #indices_p_star <- seq_len(p_star)
  #indices_p_star_p_means <- (p_star + 1):p_star_means

  # Group observations with identical predictors, missings and definition variables
  grp <- get_mx_groups(fit = fit, predictors = predictors)
  unique_groups <- unique(grp$groups)
  unique_NA_def <- unique(grp$NA_def_comb_key)

  ## ipcr object
  IPCR <- list("info" = list(ipcr_type = "iterated",
                             name = deparse(substitute(fit)),
                             class = class(fit),
                             parameters = param_names,
                             predictors = colnames(predictors),
                             analytic = analytic))

  # Start computing the initial IPCs ###########################################
  IPCs <- matrix(NA, nrow = n, ncol = q, dimnames = list(NULL, param_names))

  # Loop over unique definition variables and missing data patterns
  for (i in unique_NA_def) {
    # grp: object with the partitioning information (predictors, def, NA)
    # idx: indices
    # NA: missings
    # def: definition variables
    # comb: combination
    # sel: select (changes with each iteration)
    idx_NA_def_comb <- which(grp$NA_def_comb_key == i)
    NA_def_comb_sel <- idx_NA_def_comb[1]
    NA_pattern_sel <- grp$NA_pattern_ID[idx_NA_def_comb[1]]

    # Select centered data and moments
    dat_sel <- dat[
      idx_NA_def_comb,
      grp$NA_patterns$vars_index[[grp$NA_pattern_ID[NA_def_comb_sel]]],
      drop = FALSE]
    dat_c_sel <- dat_c[
      idx_NA_def_comb,
      grp$NA_patterns$vars_index[[grp$NA_pattern_ID[NA_def_comb_sel]]],
      drop = FALSE]
    exp_cov_sel <- OpenMx::mxGetExpected(
      model = fit, component = "covariance", defvar.row = NA_def_comb_sel
    )[grp$NA_patterns$vars_index[[grp$NA_pattern_ID[NA_def_comb_sel]]],
      grp$NA_patterns$vars_index[[grp$NA_pattern_ID[NA_def_comb_sel]]],
      drop = FALSE]
    exp_cov_inv_sel <- solve(exp_cov_sel)
    exp_means_sel <- OpenMx::mxGetExpected(
      model = fit, component = "means", defvar.row = NA_def_comb_sel)[
        , grp$NA_patterns$vars_index[[grp$NA_pattern_ID[NA_def_comb_sel]]],
        drop = FALSE]
    n_sel <- length(idx_NA_def_comb)
    p_sel <- grp$NA_patterns$n_vars[[i]]
    p_star_sel <- (p_sel * (p_sel + 1)) / 2
    Dup_sel <- lavaan::lav_matrix_duplication(n = p_sel)

    # Jacobian: missing data and ordering (cov first, then means)
    jac <- OpenMx::omxManifestModelByParameterJacobian(
      model = fit, defvar.row = NA_def_comb_sel)
    jac <- jac[grp$NA_patterns$jacobian_idx[[NA_pattern_sel]], ,
               drop = FALSE]
    if (!ms) {jac <- jac[seq_len(p_star_sel), , drop = FALSE]}

    # Check for unidentified parameters in jac due to missing data. The
    # corresponding columns in jac are deleted, and the corresponding IPCs
    # are filled with zeros later.
    nid_par <- apply(X = jac, MARGIN = 2, FUN = function(col) all(col == 0))
    jac <- jac[, !nid_par]



    # Initial IPC regression --------
    ## Individual deviations from the sample moments
    mc <- matrix(data = apply(X = dat_c_sel, MARGIN = 1,
                              FUN = function (x) {lavaan::lav_matrix_vech(x %*% t(x))}),
                 nrow = n_sel, ncol = p_star_sel, byrow = TRUE)
    vech_cov <- matrix(data = rep(lavaan::lav_matrix_vech(exp_cov_sel),
                                  times = n_sel),
                       byrow = TRUE, nrow = n_sel, ncol = p_star_sel)
    md <- mc - vech_cov
    if (ms) {
      means <- matrix(data = rep(exp_means_sel, times = n_sel), byrow = TRUE,
                      nrow = n_sel, ncol = p_sel)
      mean_dev <- dat_sel - means
      md <- cbind(md, mean_dev)
    }

    ## Weight matrix V
    V <- 0.5 * t(Dup_sel) %*% kronecker(X = exp_cov_inv_sel, Y = exp_cov_inv_sel) %*% Dup_sel
    if (ms) {
      p_star_means_sel <- p_sel * (p_sel + 3) / 2
      V_m_cov <- matrix(data = 0, nrow = p_star_means_sel, ncol = p_star_means_sel)
      V_m_cov[seq_len(p_star_sel), seq_len(p_star_sel)] <- V
      V_m_cov[(p_star_sel + 1):p_star_means_sel, (p_star_sel + 1):p_star_means_sel] <- exp_cov_inv_sel
      V <- V_m_cov
    }

    ## Initial IPCs
    W <- solve(t(jac) %*% V %*% jac) %*% t(jac) %*% V
    IPCs[idx_NA_def_comb, ] <- matrix(
      data = rep(param_estimates, times = n_sel), byrow = TRUE,
      nrow = n_sel, ncol = q, dimnames = list(NULL, param_names)) + md %*% t(W)
  }

  ## Initial IPC regression
  mlm <- lm(IPCs ~ ., data = predictors)



  # Start iterated IPC regression --------
  ## Storing objects for the updating procedure
  it_est <- matrix(coef(mlm), nrow = 1, ncol = q * (k + 1))
  it_se <- matrix(sqrt(diag(vcov(mlm))), nrow = 1, ncol = q * (k + 1))

  ## Center moment deviations at the covariate
  data_centered <- as.matrix(resid(lm(dat ~ ., data = predictors)))
  cent_md <- matrix(data = apply(X = data_centered, MARGIN = 1,
                                 FUN = function(x) {lavaan::lav_matrix_vech(x %*% t(x))}),
                    nrow = n, ncol = p_star, byrow = TRUE)
  if (ms) {
    cent_md <- cbind(cent_md, dat)
  }

  ## Calculate model fit
  if (iteration_info) {
    log_lik_individual <- rep(NA, times = n)
    param_estimates_ind <- predict(mlm)
    for (i in indices_n) {
      fit <- OpenMx::omxSetParameters(model = fit, labels = param_names,
                                      values = param_estimates_ind[i])
      fit <- suppressMessages(OpenMx::mxRun(model = fit, useOptimizer = FALSE))
      data_individual <- t(dat[i, , drop = FALSE])
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

  ###
  ### HERE HERE HERE
  ##

  # Start the iteration process --------
  nr_iterations <- 0
  difference <- rep(conv + 1, times = NCOL(it_est))
  updated_IPCs <- matrix(NA, nrow = n, ncol = q)
  colnames(updated_IPCs) <- param_names

  cent_md_up <- cent_md

  ## while loop
  while(nr_iterations < max_it & isFALSE(all(abs(difference) < conv))) {

    ### Try to update the IPCs of individuals and/or groups
    updated_IPCs <- try(expr = {

      for (i in unique_groups) { # Start loop with index i

        ID_group <- which(group == i)
        n_group <- length(ID_group)
        IPC_pred <- predictors_design_matrix_df[group == i, , drop = FALSE][1, ]
        param_estimates <- predict(mlm, newdata = IPC_pred )

        fit <- OpenMx::omxSetParameters(model = fit, labels = param_names,
                                        values = param_estimates)
        fit <- suppressMessages(OpenMx::mxRun(model = fit, useOptimizer = FALSE))
        jac <- OpenMx::omxManifestModelByParameterJacobian(model = fit)
        exp_cov <- OpenMx::mxGetExpected(model = fit, component = "covariance")
        exp_cov_inv <- solve(exp_cov)

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
          exp_means <- OpenMx::mxGetExpected(model = fit, component = "means")
          means_matrix <- matrix(rep(exp_means, times = n), byrow = TRUE, # orientation does not matter here
                                 nrow = n, ncol = p)
          means_dev <- dat - means_matrix
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


    ## Estimate updated IPC regression parameter
    mlm <- lm(updated_IPCs ~ ., data = predictors)

    ## Store results
    it_est <- rbind(it_est, c(coef(mlm)))
    it_se <- rbind(it_se, sqrt(diag(vcov(mlm))))

    ## Calculate model fit
    if (iteration_info) {
      param_estimates_ind <- predict(mlm)
      for (i in indices_n) {
        fit <- OpenMx::omxSetParameters(model = fit, labels = param_names,
                                        values = param_estimates_ind[i])
        fit <- suppressMessages(OpenMx::mxRun(model = fit, useOptimizer = FALSE))
        data_individual <- t(dat[i, , drop = FALSE])
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
    IPCR$status <- paste("Iterated IPC regression reached the maximum number of", max_it, "iterations without converging.")
  }

  if(nr_iterations < max_it & isFALSE(all(abs(difference) > conv))) {
    cat("Iterated IPC regression converged.")
    # Prepare output (from ipcr())

    IPCR$status <- paste("Iterated IPC regression converged succesfully after", nr_iterations, "iterations.")
    IPCR$IPCs <- as.data.frame(updated_IPCs)
    IPCR$mlm <- mlm

    ## data.frame with coefficients of the IPCR equations
    ### Initialize an empty data.frame for the regression coefficients
    ipcr_table <- data.frame()

    ### Initialize empty data.frame for the F tests
    F_tests_table <- data.frame()

    # Loop through models and extract coefficients
    for (i in seq_len(NCOL(mlm$coefficients))) {
      mlm_summary <- summary(mlm)[[i]]

      # Create a data frame for this model's coefficients
      temp_df_coef <- data.frame(
        Parameter = colnames(mlm$coefficients)[i],
        Predictor = rownames(mlm_summary$coefficients),
        Estimate = mlm_summary$coefficients[, "Estimate"],
        Std_Error = mlm_summary$coefficients[, "Std. Error"],
        t_Value = mlm_summary$coefficients[, "t value"],
        p_Value = mlm_summary$coefficients[, "Pr(>|t|)"],
        row.names = NULL
      )

      # Create a data.frame for this model's F-test
      F_Value <- mlm_summary$fstatistic[1]
      DF1 <- mlm_summary$fstatistic[2]
      DF2 <- mlm_summary$fstatistic[3]
      temp_df_F_test <- data.frame(
        Parameter = colnames(mlm$coefficients)[i],
        Rsquared = mlm_summary$r.squared,
        F_Value = F_Value,
        DF1 = DF1,
        DF2 = DF2,
        p_Value = pf(q = F_Value, df1 = DF1, df2 = DF2, lower.tail = FALSE),
        row.names = NULL
      )

      # Combine with the main tables
      ipcr_table <- rbind(ipcr_table, temp_df_coef)
      F_tests_table <- rbind(F_tests_table, temp_df_F_test)
    }

    ## Table with IPCR coefficients
    IPCR$output$ipcr <- ipcr_table

    ## MANCOVA table
    IPCR$output$mancova <- summary(IPCR$mancova)$stats

    ## F-test table
    IPCR$output$F_tests <- F_tests_table

    ## MANCOVA
    if (q > 1) {
      IPCR$mancova <- manova(mlm)
    } else {
      IPCR$mancova <- NULL
    }
  }



  # Prepare output --------
  if(iteration_info) {IPCR$iteration_matrix <- cbind(log_lik, it_est)}

  ## Assign class
  class(IPCR) <- "ipcr_it"

  IPCR

}
