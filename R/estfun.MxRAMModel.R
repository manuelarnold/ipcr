estfun.MxRAMModel <- function(x, ...) {

  # Model information
  N <- x$data$numObs
  data_obs <- x$data$observed[, x$manifestVars, drop = FALSE]
  p <- length(x$manifestVars)
  mean_structure <- any(x$M$free)
  p_star <- (p * (p + 1)) / 2
  p_star_seq <- seq_len(p_star)
  p_star_means <- p * (p + 3) / 2

  Dup <- lavaan::lav_matrix_duplication(n = p)

  # Does the model have definition variables?
  if (OpenMx::imxHasDefinitionVariable(x)) {

    p_star_p_means_seq <- (p_star + 1):p_star_means
    param_names <- names(x$output$estimate)
    q <- length(param_names)
    q_seq <- seq_len(q)
    p_unf <- NROW(x$A$values)

    Ident <- diag(x = 1, nrow = p_unf)
    Zero <- matrix(0, nrow = p_unf, ncol = p_unf)
    A_deriv <- lapply(q_seq, function(x) {Zero})
    S_deriv <- A_deriv
    zero <- matrix(0, nrow = p_unf, ncol = 1)
    m_deriv <- lapply(q_seq, function(x) {zero})

    for (i in q_seq) {
      A_deriv[[i]][which(x$A$labels == param_names[i], arr.ind = TRUE)] <- 1
    }

    for (i in q_seq) {
      S_deriv[[i]][which(x$S$labels == param_names[i], arr.ind = TRUE)] <- 1
    }

    for (i in q_seq) {
      m_deriv[[i]][which(x$M$labels == param_names[i])] <- 1
    }

    # Scores
    scores <- matrix(NA, nrow = N, ncol = q)
    colnames(scores) <- names(x$output$estimate)

    # Individual sample moments
    cd <- scale(x = data_obs, center = TRUE, scale = FALSE)
    if (p == 1) {
      mc <- matrix(apply(X = cd, MARGIN = 1,
                         FUN = function(x) {lavaan::lav_matrix_vech(x %*% t(x))}))
    } else {
      mc <- t(apply(X = cd, MARGIN = 1,
                    FUN = function(x) {lavaan::lav_matrix_vech(x %*% t(x))}))
    }

    # Prepare empty object for the individual deviations from the sample moments
    if (mean_structure) {
      md <- matrix(NA, nrow = N, ncol = p_star_means)
    } else {
      md <- matrix(NA, nrow = N, ncol = p_star)
    }

    # Get definition variables
    df <- identify_definition_variables(x)
    df_labels <- df[[1]]
    df_nr <- length(df_labels)
    df_data <- x$data$observed[, df[[2]], drop = FALSE]
    df_indices <- seq_along(df_labels)

    ## RAM matrices
    F_RAM <- x$F$values
    A <- x$A$values
    S <- x$S$values
    m <- t(x$M$values)
    B <- solve(Ident - A)
    FB <- F_RAM %*% B
    E <- B %*% S %*% t(B)

    if (!imxHasConstraint(x)) {

      jac <- matrix(0, nrow = p_star_means, ncol = q)

    }

    ## Assign definition variables to the corresponding RAM matrices
    RAM_df <- rep(NA, times = df_nr)
    RAM_df[which(df_labels %in% x$A$labels)] <- "A"
    RAM_df[which(df_labels %in% x$S$labels)] <- "S"
    RAM_df[which(df_labels %in% x$M$labels)] <- "M"

    # Get coordinates of the definition variables in the corresponding RAM matrices
    RAM_coord <- list()
    for (j in df_indices) {
      if (RAM_df[j] == "A") {
        RAM_coord[[j]] <- which(x$A$labels == df_labels[j], arr.ind = TRUE)
      }
      if (RAM_df[j] == "S") {
        RAM_coord[[j]] <- which(x$S$labels == df_labels[j], arr.ind = TRUE)
      }
      if (RAM_df[j] == "M") {
        RAM_coord[[j]] <- which(t(x$M$labels) == df_labels[j], arr.ind = TRUE)
      }
    }

    # Get groups of individuals with identical definition variables
    group <- transform(df_data,
                       group_ID = as.numeric(interaction(df_data,
                                                         drop = TRUE)))
    unique_groups <- unique(group$group_ID)

    # Loop over the all the different definition variable values
    for (i in unique_groups) { # Start loop with index i

      group_rows <- which(group$group_ID == i) # CHECK IF NEEDED
      group_n <- NROW(group_rows)
      df_values <- as.numeric(group[group$group_ID == i, ][1, ])

      # Update the definition variable values in the RAM matrices
      for (j in df_indices){

        if (RAM_df[j] == "A") {
          A[RAM_coord[[j]]] <- df_values[j]
        }

        if (RAM_df[j] == "S") {
          S[RAM_coord[[j]]] <- df_values[j]
        }

        if (RAM_df[j] == "M") {
          m[RAM_coord[[j]]] <- df_values[j]
        }
      } # end loop with index j

      # Update sample covariance
      B <- solve(Ident - A)
      FB <- F_RAM %*% B
      E <- B %*% S %*% t(B)
      exp_cov <- F_RAM %*% E %*% t(F_RAM)
      exp_cov_inv <- solve(exp_cov)

      # Update Jacobian matrix
      if (!imxHasConstraint(x)) { # Analytic Jacobian matrix

        for (j in seq_len(q)) {
          symm <- FB %*% A_deriv[[j]] %*% E %*% t(F_RAM)
          jac[p_star_seq, j] <- lavaan::lav_matrix_vech(symm + t(symm) + FB %*% S_deriv[[j]] %*% t(FB))
        }

        if (mean_structure) {
          for (j in seq_len(q)) {
            jac[(p_star+1):p_star_means, j] <- FB %*% A_deriv[[j]] %*% B %*% m +
              FB %*% m_deriv[[j]]
          }
        }


      } else { # Numeric Jacobian matrix
        x <- OpenMx::omxSetParameters(model = x, labels = param_names,
                                      values = param_estimates)
        x <- suppressMessages(OpenMx::mxRun(model = x, useOptimizer = FALSE))
        jac <- OpenMx::omxManifestModelByParameterJacobian(model = x)
      }

      if (mean_structure == FALSE) {jac <- jac[p_star_seq, , drop = FALSE]}


      V <- 0.5 * t(Dup) %*% kronecker(X = exp_cov_inv, Y = exp_cov_inv) %*% Dup
      if (mean_structure) {
        V_m_cov <- matrix(data = 0, nrow = p_star_means, ncol = p_star_means)
        V_m_cov[p_star_seq, p_star_seq] <- V
        V_m_cov[p_star_p_means_seq, p_star_p_means_seq] <- exp_cov_inv
        V <- V_m_cov
      }

      # Individual deviations from the sample moments
      vech_cov <- matrix(data = rep(x = lavaan::lav_matrix_vech(exp_cov),
                                    times = group_n),
                         byrow = TRUE, nrow = group_n, ncol = p_star)
      md[group_rows, 1:p_star] <- mc[group_rows, ] - vech_cov
      if (mean_structure) {
        means <- matrix(data = rep(x = FB %*% m, times = group_n), byrow = TRUE,
                        nrow = group_n, ncol = p)
        means_dev <- as.matrix(data_obs[group_rows, ]) - means
        md[group_rows, (p_star+1):p_star_means] <- means_dev
      }

      scores[group_rows, ]  <- md[group_rows, ] %*% V %*% jac
    }

  } else {

    # no definition variables

    # Model implied covariance matrix
    exp_cov <- OpenMx::mxGetExpected(model = x, component = "covariance")
    exp_cov_inv <- solve(exp_cov)

    # Partial derivaitive of model-implied moments w.r.t. model parameter estimates
    jac <- OpenMx::omxManifestModelByParameterJacobian(model = x)

    # Weght matrix
    if (mean_structure == FALSE) {jac <- jac[p_star_seq, ]}
    V <- 0.5 * t(Dup) %*% kronecker(X = exp_cov_inv, Y = exp_cov_inv) %*% Dup
    if (mean_structure) {
      V_m_cov <- matrix(data = 0, nrow = p_star_means, ncol = p_star_means)
      V_m_cov[p_star_seq, p_star_seq] <- V
      V_m_cov[(p_star + 1):p_star_means, (p_star + 1):p_star_means] <- exp_cov_inv
      V <- V_m_cov
    }

    # Individual deviations from the sample moments
    cd <- scale(x = data_obs, center = TRUE, scale = FALSE)
    if (p == 1) {
      mc <- matrix(apply(X = cd, MARGIN = 1,
                         FUN = function(x) {lavaan::lav_matrix_vech(x %*% t(x))}))
    } else {
      mc <- t(apply(X = cd, MARGIN = 1,
                    FUN = function(x) {lavaan::lav_matrix_vech(x %*% t(x))}))
    }
    vech_cov <- matrix(data = rep(x = lavaan::lav_matrix_vech(exp_cov), times = N),
                       byrow = TRUE, nrow = N, ncol = p_star)
    md <- mc - vech_cov
    if (mean_structure) {
      exp_means <- OpenMx::mxGetExpected(model = x, component = "means")
      means <- matrix(data = rep(x = exp_means, times = N), byrow = TRUE,
                      nrow = N, ncol = p)
      mean_dev <- data_obs - means
      md <- as.matrix(cbind(md, mean_dev))
    }

    # Calculates scores
    scores <- md %*% V %*% jac
  }

  scores

}
