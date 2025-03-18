#' @export

print.ipcr <- function(x, ...) {
  #
  output <- paste("Individual Parameter Contribution Regression\n")
  output <- c(output, paste("Model name:", x$info$name))
  output <- c(output, paste("Model class:", x$info$class))
  output <- c(output, paste("Model parameters:", paste(x$info$parameters,
                                                       collapse = ", ")))
  output <- c(output, paste("Predictors:", paste(x$info$predictors,
                                                 collapse = ", ")))
  if (x$info$class %in% c("MxModel")) {
    output <- c(output, paste("Linear OpenMx Model:", x$info$linear_MxModel))
  }
  writeLines(output)
}

# Stuff for print.ipcr_it

# output <- c(output, paste("Iterated IPC regression status:", x$iterated_status))
# output <- c(output, paste("Iterated IPC regression information matrix:", x$iteration_info))
# output <- c(output, paste("Convergence criterion:", x$conv))
# output <- c(output, paste("Maximum number of iterations:", x$max_it))


# Stuff for print.ipcr_reg and print.ipcr_it_reg

# output <- c(output, paste("Penality term (s):", x$s))
# output <- c(output, paste("Elastic net mixing parameter (alpha):", x$alpha))
# output <- c(output, paste("Observation weights:", x$weights))
# output <- c(output, paste("Number of penality terms:", x$nlambda))
# output <- c(output, paste("Standardize covariates:", x$standardize))
# output <- c(output, paste("Number of folds:", x$nfolds))
