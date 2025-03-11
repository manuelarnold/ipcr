#' Add together two numbers
#'
#' @param x A number.
#' @param y A number.
#' @return A number.
#' @examples
#' add(1, 1)
#' add(10, 1)

ipcr_reg <- function(fit, covariates = NULL, iterate = FALSE, iteration_info = FALSE,
                     conv = 0.0001, max_it = 50, regularization = FALSE, s = "lambda.min",
                     alpha = 1, weights = NULL, nlambda = 100, standardize = TRUE,
                     nfolds = 10, linear_MxModel = TRUE) {

#   # Checks --------
#
# ## Get model data
# #model_data <- get_data(fit)
#
# # This needs to be checked and enabled
# #check_ipcr_arguments(fit = fit, iterate = iterate, conv = conv,
# #                    max_it = max_it, linear = linear, model_data = model_data)
#
# ## Check covariates and transform into data.frame
# ### Covariate status:
# ### 0: everything is fine (required for iterated IPCR)
# ### 1: missing data in covariates (run standard IPCR with complete observations)
# ### 2: no covariates found (only compute IPCs)
# if (is.null(covariates)) {
#   covariate_status <- 2
#   info_covariates <- "no covariates found"
# } else {
#   covariates <- as.data.frame(covariates)
#   if (all(stats::complete.cases(covariates))) {
#     covariate_status <- 0
#   } else {
#     covariate_status <- 1
#     warning("Missing data in covariates detected. Standard IPC regression is performed using complete observations (rows).")
#   }
#   if (is.null(names(covariates)) | any(is.na(names(covariates)))) {
#     warning("Some covariate are not named. Renaming all covariates using the order of the data.frame.")
#     colnames(covariates) <- paste0("covariate", seq_len(NCOL(covariates)))
#   }
#   info_covariates <- names(covariates)
# }
#
# ## Check for definition variables in a OpenMx model
# ### iterated IPCR is not yet implemented for OpenMx models with definition variables
# ### Gives a warning and switches iterated from TRUE to FALSE
# if (class(fit) == "MxRAMModel") {
#   if (OpenMx::imxHasDefinitionVariable(fit) & iterate) {
#     warning("Iterated IPC regression is not available for OpenMx models with definition variables. Standard IPC regression is carried out instead of iterated IPC regression")
#     iterated <- FALSE
#   }
# }
#
#
#
# # Storing object for output --------
#
# ## Model parameters
# param_estimates <- coef_ipcr(fit)
# q <- length(param_estimates)
# param_names <- names(param_estimates)
#
# ## ipcr object
# IPC <- list("info" = list(name = deparse(substitute(fit)),
#                           class = class(fit)[[1]],
#                           parameters = param_names,
#                           covariates = info_covariates,
#                           iterate = iterate,
#                           iterated_status = NULL,
#                           iteration_info = iteration_info,
#                           conv = conv,
#                           max_it = max_it,
#                           regularization = regularization,
#                           s = s,
#                           alpha = alpha,
#                           weights = weights,
#                           nlambda = nlambda,
#                           standardize = standardize,
#                           nfolds = nfolds,
#                           linear_MxModel = linear_MxModel))
#
#
#
# # Iterated IPC regression --------
#
# if (iterate) {
#   if(covariate_status != 0) {
#     stop("Error: Iterated IPC regression requires complete covariates.")
#   }
#   IPC <- iterated_ipcr(fit, IPC = IPC, covariates = covariates,
#                        iteration_info = iteration_info, conv = conv, max_it = max_it,
#                        linear_MxModel = linear_MxModel)
# } else { ### Begin: perform standard IPC regression ###
#
#
#
#   # Standard IPC regression --------
#
#   ## Information from fitted object
#   n <- nobs(fit)
#   q <- length(param_estimates)
#
#   ## Compute score
#   scores <- estfun_ipcr(fit)
#   bread_matrix <- bread_ipcr(fit)
#   IPCs <- data.frame(matrix(param_estimates, nrow = n, ncol = q, byrow = TRUE) +
#                        scores %*% t(bread_matrix))
#   colnames(IPCs) <- param_names
#   IPC$IPCs <- IPCs
#
#   ## Perform IPC regression
#   if (covariate_status %in% c(0, 1)) {
#
#     ipcr_data <- cbind(IPCs, covariates)
#     param_names_ipcr <- paste0("IPCs_", gsub("\\(|\\)|\\]|\\[", "", param_names))
#     param_names_ipcr <- gsub(pattern = "~~", replacement = ".WITH.", x = param_names_ipcr)
#     param_names_ipcr <- gsub(pattern = "=~", replacement = ".BY.", x = param_names_ipcr)
#     param_names_ipcr <- gsub(pattern = "~", replacement = ".ON.", x = param_names_ipcr)
#     param_names_ipcr <- gsub(pattern = ",", replacement = ".", x = param_names_ipcr)
#     colnames(ipcr_data)[seq_len(q)] <- param_names_ipcr
#     IV <- paste(colnames(covariates), collapse = " + ")
#     ipcr_list <- lapply(param_names_ipcr, FUN = function(x) {
#       do.call(what = "lm",
#               args = list(formula = paste(x, "~", IV), data = as.name("ipcr_data")))
#     })
#     names(ipcr_list) <- param_names
#     IPC$regression_list <- ipcr_list
#   }
# } # End: perform standard IPC regression
#
#
#
# # Regularization --------
# if (regularization) {
#
#   # Re-code characters and factors into dummy variables
#   if (any(unlist(lapply(covariates, function(x) {is.character(x) | is.factor(x)})))) {
#     covariates_strings <- Filter(function(x) {is.character(x) | is.factor(x)}, covariates)
#     string_formula <-  paste(colnames(covariates_strings), collapse = "+")
#     dummies <- stats::model.matrix(stats::formula(paste("~", string_formula)), data = covariates_strings)[, -1]
#     covariates <- covariates[, !unlist(lapply(covariates, function(x) {is.character(x) | is.factor(x)}))]
#     covariates <- cbind(covariates, dummies)
#   }
#
#   IPC$regularized_regression_list <- lapply(IPC$IPCs, FUN = function(y) {
#     glmnet::cv.glmnet(x = as.matrix(covariates), y = y, alpha = alpha,
#                       weights = weights, nlambda = nlambda, standardize = standardize,
#                       nfolds = nfolds)})
# }
#
#
#
# # Prepare output --------
# IPC$output$info <- print_info(IPC$info)
# if (is.list(IPC$regression_list)) {
#   IPC$output$coefficients_matrix <- coefficients_matrix(IPC)
# }
# if (is.list(IPC$regularized_regression_list)) {
#   IPC$output$regularized_coefficients_matrix <- regularized_coefficents_matrix(IPC)
# }
#
#
#
#   # class(IPCR) <- "ipcr_reg"
#   # class(IPCR) <- "ipcr_it_reg"

}
