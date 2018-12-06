#' @title Inidividual Parameter Contribution Regression Summary
#' @description This functions returns the estimated coefficients of the individual parameter contribution (IPC) regression equations as well as the corresponding standard errors, t-values, and p-values.
#' @param object an ipcr object.
#' @param parameter a character vector with the target parameters of the structural equation model (SEM) that are to be inspected.
#' @param method a character string indicating the method to be used; for method = "both" static and iterated IPC regression is used (default), for method = "static" only static IPC coefficients are displayed, and for method = "iterate" only iterated coefficients are shown.
#' @return NULL
#' @examples
#' # Display the iterated individual contribution regression parameters of all SEM parameters.
#' summary(object = IPC_reg, method = "iterate")
#' @export summary.ipcr

summary.ipcr <- function(object, parameter = NULL, method = "both") {

  # Check if IPC regression parameters have been estimated
  if (any(c("StandardIPCRegression", "IteratedIPCRegression") %in% names(object))) {

    # Select target parameters
    if (is.null(parameter)) {
      parameter <- names(object$StandardIPCRegression)
    }

    # Number of parameters
    q <- length(parameter)
    param_position <- which(names(object$StandardIPCRegression) %in% parameter)
    significant <- FALSE
    output <- c()

    # Create output for static IPC regression estimates
    if (method %in% c("static", "both")) {

    output[1] <- "Static IPC regression parameter estimates:"

      for (i in 1:q) {
        output[length(output) + 1] <- paste("Parameter:", parameter[i])
        lm_summary <- capture.output(summary(object$StandardIPCRegression[[param_position[i]]]))
        start_row <- grep(pattern = "Coefficients:", x = lm_summary, fixed = TRUE) + 1
        if (identical(grep(pattern = "---", x = lm_summary, fixed = TRUE), integer(0))) {
          # no significant coefficents
          end_row <- grep(pattern = "Residual standard error", x = lm_summary, fixed = TRUE) - 2
        } else {
          significant <- TRUE
          end_row <- grep(pattern = "---", x = lm_summary, fixed = TRUE) - 1
        }
        number_row <- end_row - start_row
        output[(length(output) + 1):(length(output) + 1 + number_row)] <- lm_summary[start_row:end_row]
      }
      if (isTRUE(significant)) {
        output[length(output) + 1] <- "---"
        output[length(output) + 1] <- "Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1"
      }
    }

    # Create output for iterated IPC regression estimates
    if (method %in% c("iterated", "both")) {

      significant <- FALSE

      if(method == "both") {
        output[length(output) + 1] <- ""
      }

      output[length(output) + 1] <- "Iterated IPC regression parameter estimates:"

      for (i in 1:q) {
        output[length(output) + 1] <- paste("Parameter:", parameter[i])
        lm_summary <- capture.output(summary(object$IteratedIPCRegression[[param_position[i]]]))
        start_row <- grep(pattern = "Coefficients:", x = lm_summary, fixed = TRUE) + 1
        if (identical(grep(pattern = "---", x = lm_summary, fixed = TRUE), integer(0))) {
          # no significant coefficents
          end_row <- grep(pattern = "Residual standard error", x = lm_summary, fixed = TRUE) - 2
        } else {
          significant <- TRUE
          end_row <- grep(pattern = "---", x = lm_summary, fixed = TRUE) - 1
        }
        number_row <- end_row - start_row
        output[(length(output) + 1):(length(output) + 1 + number_row)] <- lm_summary[start_row:end_row]
      }
      if (isTRUE(significant)) {
        output[length(output) + 1] <- "---"
        output[length(output) + 1] <- "Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ' 1"
      }
    }

    writeLines(output)
    } else {
    cat("No IPC regression parameters estimated")
  }
}
