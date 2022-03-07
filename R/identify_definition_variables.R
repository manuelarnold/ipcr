#' Identifies definition variables in OpenMx models
#'
#' This is an internal function that identifies which variables in the
#' data set stored in an \pkg{OpenMx} model are definition variables.
#'
#' @param x An OpenMx model.

identify_definition_variables <- function(x) {
  definition_variables <- c()
  for (i in 1:length(x@matrices)) {
    definition_variables <- c(definition_variables,
                              sapply(x@matrices[[i]]$labels,
                                     OpenMx::imxIsDefinitionVariable))
  }
  definition_variables <- names(which(definition_variables))
  list(labels = definition_variables, # OpenMx labels
       data = sub(".*\\.", "", definition_variables)) # column names
}
