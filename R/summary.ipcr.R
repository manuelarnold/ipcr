#' @title Individual Parameter Contribution Regression Summary
#'
#' @description Summary method for objects of class \code{"ipcr"}.
#'
#' @param object An object of class \code{"ipcr"}.
#' @param what A character string specifying what to display. Options are:
#'   \itemize{
#'     \item \code{"everything"} (default): Displays all summary components.
#'     \item \code{"marginal"}: Shows only marginal linear effects.
#'     \item \code{"ipcr"}: Shows only individual parameter contribution regression coefficients.
#'   }
#' @param digits An integer indicating the number of decimal places to use in the output.
#'   Default is 3.
#' @param ... Additional arguments passed to or from other methods.
#'
#' @return Invisibly returns \code{NULL}.
#'
#' @export

summary.ipcr <- function(object, what = "everything", digits = 3, ...) {
  cat("Individual Parameter Contribution Regression\n\n")

  if (what %in% c("everything", "marginal")) {

  cat("Marginal Linear Effects:\n\n")

  cat("  MANCOVA:\n")
  mancova <- as.data.frame(object$output$mancova)
  mancova$stars <- sapply(mancova[, "Pr(>F)"], add_stars)
  mancova[, 1:6] <- round(mancova[, 1:6], digits = digits)
  mancova["Residuals", 2:7] <- ""
  mancova <- cbind(Predictor = rownames(mancova), mancova)
  mancova <- rbind(c(colnames(mancova)[1:7], ""), mancova)
  colnames(mancova) <- NULL
  mancova[, 1] <- format(mancova[, 1], justify = "left")
  mancova_output <- capture.output(print(mancova, row.names = FALSE))
  cat(paste0("  ", mancova_output), sep = "\n") # Add 2 spaces for indentation
  cat("\n")

  cat("  F-Tests:\n")

  F_tests_output <- object$output$F_tests
  F_tests_output$stars <- sapply(F_tests_output[, "p_Value"], add_stars)
  F_tests_output <- F_tests_output[, -c(4, 5)]
  F_tests_output[, 2:4] <- round(F_tests_output[, 2:4], digits = digits)
  # Align first column with heading to the left
  F_tests_output <- rbind(c(colnames(F_tests_output)[1:4], ""), F_tests_output)
  colnames(F_tests_output) <- NULL
  F_tests_output[, 1] <- format(F_tests_output[, 1], justify = "left")
  F_tests_output <- capture.output(print(F_tests_output, row.names = FALSE, ))
  cat(paste0("  ", F_tests_output), sep = "\n")# Add 2 spaces for indentation
  cat(paste("   Numerator degrees of freedom (df1):",  object$output$F_tests[1, "DF1"]))
  cat("\n")
  cat(paste("   Denominator degrees of freedom (df2):",  object$output$F_tests[1, "DF2"]))
  cat("\n\n")
  }


  cat("Individual Parameter Contribution Regression Coefficients:\n")
  ipcr_output <- object$output$ipcr
  ipcr_output$stars <- sapply(ipcr_output[, "p_Value"], add_stars)
  ipcr_output[, 3:6] <- round(ipcr_output[, 3:6], digits = digits)
  ipcr_output <- rbind(c(colnames(ipcr_output)[1:6], ""), ipcr_output)
  colnames(ipcr_output) <- NULL
  ipcr_output[, 1:2] <- format(ipcr_output[, 1:2], justify = "left")
  ipcr_output <- capture.output(print(ipcr_output, row.names = FALSE))
  cat(paste0("  ", ipcr_output), sep = "\n") # Add 2 spaces for indentation
}
