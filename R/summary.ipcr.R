#' @title Inidividual Parameter Contribution Regression Summary
#' @description This functions returns the coefficients of the individual parameter
#' contribution (IPC) regression equations.
#' @param object an ipcr object.
#' @param regularization a logical value. Wheter to show results for regularized IPC
#' regression (if computed) or non-regularized results.
#' @param digits integer indicating the number of decimal places to be used.
#' @param vebose a logical value; if TRUE ipcr settings, non-regularized, and regularized
#' results are shown.
#' @param ... ... further arguments passed to and from methods.
#' @return NULL
#' @export

summary.ipcr <- function(object, regularization = TRUE, digits = 3, verbose = FALSE, ...) {

  # Verbose summary
  if (verbose) {
    output <- "ipcr settings:"
    output <- c(output, "")
    output <- c(output, print_info(object$info)[-1])
    if (is.list(object$regression_list)) {
      output <- c(output, "")
      if(object$info$iterate) {
        output <- c(output, "Iterated IPC regression coefficients:")
      } else {
        output <- c(output, "Standard IPC regression coefficients:")
      }
      rounded_coefficients <- data.frame(
        lapply(object$output$coefficients_matrix,
               function(y) {if(is.numeric(y)) round(y, digits = digits) else y}))
      output <- c(output, capture.output(rounded_coefficients))
    }
    if (is.list(object$regularized_regression_list)) {
      output <- c(output, "")
      output <- c(output, "Regularized IPC regression coefficients")
      zero_coefficients <- which(object$output$regularized_coefficients_matrix$Estimate == 0)
      rounded_coefficients <- data.frame(
        lapply(object$output$regularized_coefficients_matrix,
               function(y) {if(is.numeric(y)) round(y, digits = digits) else y}))
      rounded_coefficients[zero_coefficients, "Estimate"] <- "."
      output <- c(output, capture.output(rounded_coefficients))
    }
  } else {

    ## Standard summary
    # Regularized IPC regression was performed and requested
    if(is.list(object$regularized_regression_list) & regularization) {
      if(object$info$iterate) {
        output <- "Regularized iterated IPC regression coefficients:"
      } else {
        output <- "Regularized standard IPC regression coefficients:"
      }
      zero_coefficients <- which(object$output$regularized_coefficients_matrix$Estimate == 0)
      rounded_coefficients <- data.frame(
        lapply(object$output$regularized_coefficients_matrix,
               function(y) {if(is.numeric(y)) round(y, digits = digits) else y}))
      rounded_coefficients[zero_coefficients, "Estimate"] <- "."
      output <- c(output, capture.output(rounded_coefficients))
    } else if (is.list(object$regression_list)) {
      if(object$info$iterate) {
        output <- "Iterated IPC regression coefficients:"
      } else {
        output <- "Standard IPC regression coefficients:"
      }
      rounded_coefficients <- data.frame(
        lapply(object$output$coefficients_matrix,
               function(y) {if(is.numeric(y)) round(y, digits = digits) else y}))
      output <- c(output, capture.output(rounded_coefficients))
    } else {
      output <- "IPCs:"
      output <- c(output, object$output$info)
    }
  }

  writeLines(output)

}
