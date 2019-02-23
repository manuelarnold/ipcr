#' @title Individual Parameter Contribution Regression (With Analytic Jacobian)
#' @description Explain parameter heterogeneity in a structural equation model with individual parameter contribution (IPC) regression. IPC Regression allows to test if parameter values differ across values of covariates.
#' @param fit a structural equation model estimated with OpenMx. Only single-group RAM-type OpenMx model are allowed. The model must contain an estimated covariance structure. The mean structure is optional.
#' @param covariates a data frame or matrix with covariates.
#' @param iterated a logical value; if TRUE iterated IPC regression is performed after initial IPC regression.
#' @param convergence_criterion an integer used as stopping criterion for iterated IPC regression. Criterion is the largest difference in all parameter estimate between iterations.
#' @param max_it the maximum number of iterations for iterated IPC regressions.
#' @return Returns a new ipcr object.
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
#' @export

ipcr2 <- function(fit, covariates = NULL, iterated = FALSE, conv = 0.0001,
                 max_it = 50) {

  ### Preparations

  # Check if the arguments are correct
  get_arguments_checked(fit = fit, covariates = covariates, iterated = iterated,
                        conv = conv, max_it = max_it)

  # Storing object for output
  IPC <- list("Info" = list(name = fit$name))

  # Check covariates for missing and add information to the output
  if (is.null(covariates)) {
    IPC$Info$covariates <- "no covariates specified"
  } else {
    if (all(complete.cases(covariates))) {
      covariates <- as.matrix(covariates)
      if (is.null(colnames(covariates)) | any(is.na(colnames(covariates)))) {
        colnames(covariates) <- paste0("z", 1:NCOL(covariates))
      }
      IPC$Info$covariates <- colnames(covariates)
      IPC$Covariates <- covariates

    } else {
      IPC$Info$covariates <- "missing data in covariates"
    }
  }

  # Name the covariates if needed
  if (!is.null(covariates)) {
    if (is.null(colnames(covariates)) | any(is.na(colnames(covariates)))) {
      colnames(covariates) <- paste0("z", 1:NCOL(covariates))
      IPC$Info$covariates <- colnames(covariates)
    }
  }

  # Basic information about the data set and model
  X <- get_model_info(fit = fit, covariates = covariates)

  # If linear model is true
  X <- c(X, RAM_deriv = list(get_analytic_jac_input(x = X)))

  # Save model parameter names
  IPC$Info$parameters <- X$param_names


  ### Calculate initial IPCs

  # Get individual deviations from the sample moments
  X$moment_deviations <- get_moment_deviations(x = X)

  # Get initial IPCs
  X$initial_IPCs <- get_initial_IPCs(x = X)

  # Store initial IPCs
  IPC$InitialIPCs <- X$initial_IPCs

  # Status report
  cat("Initial IPCs computed\n")


  ### Initial IPC regression

  if (IPC$Info$covariates[1] != "no covariates specified" &
      IPC$Info$covariates[1] != "missing data in covariates") {

    IPC_reg_data <- cbind(X$initial_IPCs, X$covariates)
    IPC_reg_list <- list()
    IV <- paste(colnames(covariates), collapse = " + ")

    for (i in 1:X$q) {
      IPC_reg_list[[i]] <- do.call(what = "lm",
                                   args = list(formula = paste(X$param_names[i], "~", IV),
                                               data = as.name("IPC_reg_data")))
    }

    names(IPC_reg_list) <- X$param_names
    IPC$InitialIPCRegression <- IPC_reg_list


    # Status report
    cat("Initial IPC regression parameters estimated\n")

  } else {
    if (IPC$Info$covariates[1] == "no covariates specified") {
      warning("No covariates specified.", call. = FALSE)
    }
    if(IPC$Info$covariates[1] == "missing data in covariates") {
      warning("Missing data in covariates. IPC regression parameters cannot be estimated."
              , call. = FALSE)
    }
  }


  ### Iterated IPC Regression

  if (iterated == TRUE &
      IPC$Info$covariates[1] != "no covariates specified" &
      IPC$Info$covariates[1] != "missing data in covariates") {

    # Add convergence criterion to info
    IPC$Info$convergence_criterion <- conv

    # Storing objects for the updating procedure
    it_est <- matrix(sapply(X = IPC_reg_list,
                            FUN = function(x) {coef(x)}),
                     nrow = 1, ncol = X$q * (X$z + 1))
    it_se <- matrix(sapply(X = IPC_reg_list,
                           FUN = function(x) {sqrt(diag(vcov(x)))}),
                    nrow = 1, ncol = X$q * (X$z + 1))

    # Add new information to list
    X <- c(X, list(IPC_reg = IPC_reg_list,
                   cent_md = get_centered_moment_deviations(x = X),
                   group = get_groups(x = X)))

    ### Prepare input for get_updated_IPCs2.

    # These matrices change across groups
    A_up <- fit$A$values
    S_up <- fit$S$values
    M_up <- t(fit$M$values)
    F_up <- fit$F$values

    # Identity matrix for (I - A)^(-1)
    Ident <- diag(x = 1, nrow = X$p_unf)

    # Assign SEM parameters to the corresponding RAM matrices
    RAM_params <- rep(NA, X$q)
    RAM_params[which(X$param_names %in% X$fit$A$labels)] <- "A"
    RAM_params[which(X$param_names %in% X$fit$S$labels)] <- "S"
    RAM_params[which(X$param_names %in% X$fit$M$labels)] <- "M"

    # Get coordinates of the parameters in the corresponding RAM matrices
    RAM_coord <- list()
    for (i in 1:X$q) {
      if (RAM_params[i] == "A") {
        RAM_coord[[i]] <- which(X$fit$A$labels == X$param_names[i], arr.ind = TRUE)
      }
      if (RAM_params[i] == "S") {
        RAM_coord[[i]] <- which(X$fit$S$labels == X$param_names[i], arr.ind = TRUE)
      }
      if (RAM_params[i] == "M") {
        RAM_coord[[i]] <- which(X$fit$M$labels == X$param_names[i], arr.ind = TRUE)
      }
    }

    # Start the iteration process
    nr_iterations <- 0
    difference <- rep(x = conv + 1, times = NCOL(it_est))
    breakdown <- FALSE

    while(nr_iterations < max_it & isFALSE(all(abs(difference) < conv)) &
          breakdown == FALSE) {

      # Update the IPCs for every individual in each group
      updated_IPCs <- try(get_updated_IPCs2(x = X, updated_IPCs, A_up, S_up,
                                            M_up, F_up, Ident, RAM_params,
                                            RAM_coord),
                          silent = TRUE)
      if (class(updated_IPCs) == "try-error") {
        breakdown <- TRUE
      }

      if (breakdown == FALSE) {

        # Estimate updated IPC regression parameter
        IPC_reg_data_it <- data.frame(cbind(updated_IPCs, X$covariates))
        for (i in 1:X$q) {
          X$IPC_reg[[i]] <- do.call(what = "lm",
                                    args = list(formula = paste(X$param_names[i], "~", IV),
                                                data = as.name("IPC_reg_data_it")))
        }

        # Store results
        it_est <- rbind(it_est, c(sapply(X = X$IPC_reg,
                                         FUN = function(x) {coef(x)})))
        it_se <- rbind(it_se, c(sapply(X = X$IPC_reg,
                                       FUN = function(x) {sqrt(diag(vcov(x)))})))

        # Covergence criteria
        difference <- it_est[nrow(it_est), ] -
          it_est[nrow(it_est) - 1, ]
        nr_iterations <- nr_iterations + 1
        cat("Iteration:", nr_iterations, "\n")

      }
    }

    # Error reports
    if (breakdown == TRUE) {
      cat("The iterated IPC regression aborted prematurely.")
      IPC$Info$iterated <- "Iterated IPC regression aborted prematurely."
    }

    if(breakdown == FALSE & nr_iterations == max_it &
       isFALSE(all(abs(difference) <= conv))) {
      cat(paste("The iterated IPC regression algorithm did not converge after", max_it, "iterations. Consider to increase the maximum number of iterations."))
      IPC$Info$iterated <- paste("Iterated IPC regression reached the maximum number of", max_it, "iterations without converging.")
    }

    if(breakdown == FALSE & nr_iterations < max_it &
       isFALSE(all(abs(difference) > conv))) {
      cat("Iterated IPC regression converged.")
      IPC$Info$iterate <- paste("Iterated IPC regression converged succesfully after", nr_iterations, "iterations.")
      IPC$IteratedIPCs <- as.data.frame(updated_IPCs)
      IPC$IteratedIPCRegression <- X$IPC_reg
    }

  }

  class(IPC) <- "ipcr"
  return(IPC)

}
