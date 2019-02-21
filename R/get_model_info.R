#' @title Information from the OpenMx model object
#' @description Internal function that retrieves information from an OpenMx model object.
#' @param fit fitted single-group RAM-type model.
#' @param covariates IPC regression predictor variables.
#' @return List with model information.
#' @export

get_model_info <- function(fit, covariates) {

  p <- length(fit$manifestVars)
  ms <- any(fit$M$free)
  exp_cov <- mxGetExpected(model = fit, component = "covariance")

  info <- list(fit = fit,
               N = fit$output$data[[1]][[1]],
               data_obs = fit$data$observed[, fit$manifestVars, drop = FALSE],
               p = p,
               p_unf = nrow(fit$A$values),
               p_star = (p * (p + 1)) / 2,
               p_star_means = p * (p + 3) / 2,
               mean_structure = ms,
               param_estimates = fit$output$estimate,
               param_names = names(fit$output$estimate),
               q = length(fit$output$estimate))
  if(ms) {
    info <- c(info, list(exp_mean = mxGetExpected(model = fit,
                                                  component = "means")))
  }
  info <- c(info, list(exp_cov = exp_cov,
                       exp_cov_inv = solve(exp_cov),
                       covariates = covariates,
                       z = NCOL(covariates)))
  if (p == 1) {
    info <- c(info, list(D = matrix(data = 1)))
  } else {
    info <- c(info, list(D = duplication.matrix(n = p)))
  }
  return(info)
}
