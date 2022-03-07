print_info <- function(x) {
  output <- paste("IPC regression")
  output <- c(output, paste("Model name:", x$name))
  output <- c(output, paste("Model name:", x$name))
  output <- c(output, paste("Model class:", x$class))
  output <- c(output, paste("Model parameters:", paste(x$parameters, collapse = ", ")))
  output <- c(output, paste("Covariates:", paste(x$covariates, collapse = ", ")))
  output <- c(output, paste("Iterated IPC regression:", x$iterate))
  if(x$iterate) {
    output <- c(output, paste("Iterated IPC regression status:", x$iterated_status))
    output <- c(output, paste("Iterated IPC regression information matrix:", x$iteration_info))
    output <- c(output, paste("Convergence criterion:", x$conv))
    output <- c(output, paste("Maximum number of iterations:", x$max_it))
  }
  output <- c(output, paste("Regularizations:", x$regularization))
  if (x$regularization) {
    output <- c(output, paste("Penality term (s):", x$s))
    output <- c(output, paste("Elastic net mixing parameter (alpha):", x$alpha))
    output <- c(output, paste("Observation weights:", x$weights))
    output <- c(output, paste("Number of penality terms:", x$nlambda))
    output <- c(output, paste("Standardize covariates:", x$standardize))
    output <- c(output, paste("Number of folds:", x$nfolds))
  }
  if (x$class == "MxRAMModel") {
    output <- c(output, paste("Assumed MxRAMModel:", x$linear_OpenMx))
  }
  output
}
