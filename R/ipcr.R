#' @title Individual Parameter Contribution Regression
#' @description Explain parameter heterogeneity in a fitted model with individual parameter contribution (IPC) regression. IPC Regression allows to test if parameter values differ across values of covariates.
#' @param fit a fitted model object.
#' @param covariates a data frame or matrix with covariates.
#' @param iterated a logical value; if TRUE iterated IPC regression is performed after initial IPC regression.
#' @param convergence_criterion an integer used as stopping criterion for iterated IPC regression. Criterion is the largest difference in all parameter estimate between iterations.
#' @param max_it the maximum number of iterations for iterated IPC regressions.
#' @param linear for MxRAMModel-type models only. A logical value indicating if the structural equation model is linear (default) or non linear.
#' @return Returns an ipcr object.
#' @examples
#' # Specify an OpenMx model
#' m <- mxModel(model = "CFA",
#'             manifestVars = c("x1", "x2", "x3"),
#'             latentVars = "f",
#'             type = "RAM",
#'             mxData(observed = ipcr_data, type = "raw"),
#'             mxPath(from = "f", to = c("x1", "x2", "x3"), arrows = 1,
#'                    free = c(FALSE, TRUE, TRUE), values = 1,
#'                    labels = c("l1", "l2", "l3")),
#'             mxPath(from = "f", arrows = 2, free = TRUE, values = 0.75,
#'                    labels = "var_f"),
#'             mxPath(from = c("x1", "x2", "x3"), arrows = 2, free = TRUE,
#'                    values = 0.25, labels = c("e1", "e2", "e3")),
#'             mxPath(from = "one", to = c("x1", "x2", "x3"), arrows = 1,
#'                    free = FALSE, values = 0))
#'
#' # Fit the model
#' fit <- mxTryHard(model = m)
#'
#' # Investigate the parameter estiamtes with IPC regression
#' IPC_reg <- ipcr(fit = fit, covariates = ipcr_covariates, iterated = TRUE)
#'
#' # Get an overview about parameter differences
#' plot(IPC_reg)
#'
#' # Show IPC regression output
#' summary(object = IPC_reg)
#' @export ipcr3

ipcr <- function(fit, covariates = NULL, iterated = FALSE, conv = 0.0001,
                  max_it = 50, linear = TRUE) {

  #### Checks #####
  # Get model data
  model_data <- get_data(fit)

  # Transform covariates into data frame
  # covariates = NULL is preserved
  covariates <- as.data.frame(covariates)

  # Check if the arguments are correct
  check_ipcr_arguments(fit = fit, iterated = iterated, conv = conv,
                       max_it = max_it, linear = linear, model_data = model_data)

  # Storing object for output
  IPC <- list("info" = list(name = deparse(substitute(fit)),
                            class = class(fit)))

  # Status code
  # 0: everything is fine (required for iterated IPCR)
  # 1: missing data in covariates (run standard IPCR with complete observations)
  # 2: no covariates found (only compute IPCs)

  # Check covariates
  if (is.null(covariates)) {
    status_code <- 2
    IPC$info$covariates <- "no covariates found"
  } else {
    if (all(complete.cases(covariates))) {
      status_code <- 0
    } else {
      status_code <- 1
      warning("Missing data in covariates detected. Standard IPC regression is performed using complete observations (rows).")
    }
    if (is.null(names(covariates)) | any(is.na(names(covariates)))) {
      warning("Some covariate are not named. Renaming all covariates using the order of the data.frame.")
      colnames(covariates) <- paste0("covariate", seq_len(NCOL(covariates)))
    }
    IPC$info$covariates <- names(covariates)
  }

  # Iterated IPC regression
  if (all(iterated, status_code == 0)) {
    IPC <- iterated_ipcr(fit, IPC = IPC, covariates = covariates, conv = conv,
                         max_it = max_it, linear = linear)
  } else { ### Begin: perform standard IPC regression ###

    # Check covariates for missing and add information to the output
    if (is.null(covariates)) {
      IPC$info$covariates <- "no covariates specified"
    } else {
      if (all(complete.cases(covariates))) {
        covariates <- as.matrix(covariates)
        if (is.null(colnames(covariates)) | any(is.na(colnames(covariates)))) {
          colnames(covariates) <- paste0("z", 1:NCOL(covariates))
        }
        IPC$info$covariates <- colnames(covariates)
        IPC$covariates <- covariates

      } else {
        IPC$info$covariates <- "missing data in covariates"
      }
    }

    # Name the covariates if needed
    if (!is.null(covariates)) {
      if (is.null(colnames(covariates)) | any(is.na(colnames(covariates)))) {
        colnames(covariates) <- paste0("z", 1:NCOL(covariates))
        IPC$info$covariates <- colnames(covariates)
      }
    }


    #### Compute IPCs ####
    # Information from fitted object
    n <- nobs(fit)
    param_estimates <- coef(fit)[-c(2,3,4, 7:9)]
    q <- length(param_estimates)
    param_names <- IPC$info$parameters <- names(param_estimates)

    # Compute score
    scores <- sandwich::estfun(fit)
    bread_matrix <- sandwich::bread(fit)[-c(2,3,4, 7:9), -c(2,3,4, 7:9)]
    IPCs <- data.frame(matrix(param_estimates, nrow = n, ncol = q, byrow = TRUE) +
                         scores %*% t(bread_matrix))
    colnames(IPCs) <- param_names
    IPC$IPCs <- IPCs


    ##### IPC Regression ####
    if (IPC$info$covariates[1] != "no covariates specified" &
        IPC$info$covariates[1] != "missing data in covariates") {

      ipcr_data <- cbind(IPCs, covariates)
      param_names_ipcr <- paste0("IPCs_", gsub("\\(|\\)", "", param_names))
      param_names_ipcr <- gsub(pattern = "~~",replacement = ".WITH.", x = param_names_ipcr)
      param_names_ipcr <- gsub(pattern = "~",replacement = ".ON.", x = param_names_ipcr)
      colnames(ipcr_data)[seq_len(q)] <- param_names_ipcr
      ipcr_list <- list()
      IV <- paste(colnames(covariates), collapse = " + ")
      for (i in seq_len(q)) {
        ipcr_list[[i]] <- do.call(what = "lm",
                                  args = list(formula = paste(param_names_ipcr[i], "~", IV),
                                              data = as.name("ipcr_data")))
      }
      names(ipcr_list) <- param_names
      IPC$regression <- ipcr_list
    }
  } # End: perform standard IPC regression

  class(IPC) <- "ipcr"
  return(IPC)

}
