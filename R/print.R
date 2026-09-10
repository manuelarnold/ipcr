#' @export

print.ipcr <- function(x, ...) {
  output <- paste("Individual Parameter Contribution Regression\n")
  output <- c(output, paste("Model name:", x$info$name))
  output <- c(
    output,
    paste("Model class:", paste(x$info$class, collapse = ", "))
  )
  output <- c(output, paste("Model parameters:", paste(x$info$parameters,
                                                       collapse = ", ")))
  output <- c(output, paste("Predictors:", paste(x$info$predictors,
                                                 collapse = ", ")))
  writeLines(output)
  invisible(x)
}

#' @export
print.ipcr_it <- function(x, ...) {
  print.ipcr(x, ...)
  if (identical(x$info$variance_parameterization, "log")) {
    cat(
      "Variance parameterization: Log scale for free RAM S-diagonal parameters\n"
    )
  }
  cat("Iterated IPCR status:", x$convergence$status, "\n")
  if (identical(x$convergence$likelihood_stable, FALSE)) {
    cat(
      "Likelihood stability warning:",
      likelihood_stability_warning_ipcr(
        C_l_rate_adjusted = x$convergence$C_l_rate_adjusted,
        threshold = x$convergence$likelihood_warning_threshold
      ),
      "\n"
    )
  }
  invisible(x)
}
