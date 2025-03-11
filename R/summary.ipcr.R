#' @title Inidividual Parameter Contribution Regression Summary
#' @description This functions returns the coefficients of the individual parameter
#' contribution (IPC) regression equations.
#' @param object an ipcr object.
#' @param regularization a logical value; whether to show results for regularized IPC
#' regression (if computed) or non-regularized results.
#' @param digits integer indicating the number of decimal places to be used.
#' @param verbose a logical value; if TRUE ipcr settings, non-regularized, and regularized
#' results are shown.
#' @param ... further arguments passed to and from methods.
#' @return NULL
#' @export

summary.ipcr <- function(object, ...) {
  cat("Individual Parameter Contribution Regression\n\n")

  cat("Marginal Linear Effects\n\n")

  cat("MANCOVA:\n")
  print(object$output$mancova)
  cat("\n")

  cat("F-Tests:\n")
  print(object$output$F_test)
  cat("\n")

  cat("Individual Parameter Contribution Regression Coefficients:\n")
  print(object$output$ipcr)
}
