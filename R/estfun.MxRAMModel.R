estfun.MxRAMModel <- function(x, ...) {

  # Model information
  N <- x$data$numObs
  data_obs <- x$data$observed[, x$manifestVars, drop = FALSE]
  p <- length(x$manifestVars)
  mean_structure <- any(x$M$free)
  p_star <- (p * (p + 1)) / 2
  p_star_seq <- seq_len(p_star)
  p_star_means <- p * (p + 3) / 2

  # Model implied covariance matrix
  exp_cov <- OpenMx::mxGetExpected(model = x, component = "covariance")
  exp_cov_inv <- solve(exp_cov)

  # Partial derivaitive of model-implied moments w.r.t. model parameter estimates
  jac <- OpenMx::omxManifestModelByParameterJacobian(model = x)

  # Weght matrix
  if (mean_structure == FALSE) {jac <- jac[p_star_seq, ]}
  Dup <- lavaan::lav_matrix_duplication(n = p)
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
  md %*% V %*% jac

}
