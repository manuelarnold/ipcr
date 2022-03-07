compute_numeric_scores <- function(x, ...) {
  UseMethod("compute_numeric_scores")
}


compute_numeric_scores.default <- function(x, ...) {
  stop(paste(c("Numeric scores are currently only implemented for lavaan",
               "andMxRAMModels-type models. Request analytic scores for other",
               "model types."), collapse = ""))
}


compute_numeric_scores.MxRAMModel <- function(x, ...) {

  parameters <- OpenMx::omxGetParameters(x)
  param_names <- names(parameters)
  n_param <- length(parameters)
  scores <- matrix(NA, nrow = x$data$numObs, ncol = n_param)

  for (j in seq_len(n_param)) {

    epsilon <- 1e-5 # !!!change to adjustable argument late on!!!

    model_left <- OpenMx::omxSetParameters(x,
                                           labels = param_names[j],
                                           values = parameters[j] - epsilon)
    model_left <- OpenMx::mxRun(model_left, useOptimizer = FALSE, silent = TRUE)
    ind_logLik_left <- individual_FIML(
      data = x$data$observed,
      exp_m = OpenMx::mxGetExpected(model_left, component = "mean"),
      exp_cov = OpenMx::mxGetExpected(model_left, component = "covariance")
    )

    model_right <- omxSetParameters(x,
                                    labels = param_names[j],
                                    values = parameters[j] + epsilon)
    model_right <- OpenMx::mxRun(model_right, useOptimizer = FALSE, silent = TRUE)
    ind_logLik_right <- individual_FIML(
      data = x$data$observed,
      exp_m = OpenMx::mxGetExpected(model_right, component = "mean"),
      exp_cov = OpenMx::mxGetExpected(model_right, component = "covariance")
    )

    scores[, j] <- (ind_logLik_right-ind_logLik_left)/(2*epsilon)

  }

  scores

}

