#' @title Individual Parameter Contribution Regression
#' @description Explain heterogeneity in a parametric model with individual parameter
#' contribution (IPC) regression. IPCs regression allows studying heterogeneity
#' in parameter estimates as a linear function of covariates. \code{\link[ipcr]{ipcr}} was
#' mainly written for investigating parameter heterogeneity in structural equation models
#' fitted with \pkg{lavaan} or \pkg{OpenMx} but works also for models estimated with
#' \code{\link[stats4]{lm}} and \code{\link[stats4]{glm}}.
#' @param fit a fitted model object.
#' @param covariates a data frame or matrix with covariates. Interaction and polynomial
#' terms can be included as new variables which may require centering.
#' @param iterate a logical value; if TRUE iterated IPC regression is performed.
#' Currently, iterated IPC regression is only available for models fitted with
#' \pkg{lavaan} or \pkg{OpenMx}.
#' @param iteration_info a logical value; if TRUE the parameter values for each iteration
#' with corresponding log-likelihood value are stored in a matrix. Requesting this matrix
#' increases the runtime.
#' @param conv an integer used as stopping criterion for iterated IPC regression. Criterion
#'  is the largest difference in any parameter estimate between iterations.
#' @param max_it the maximum number of iterations for iterated IPC regressions.
#' @param regularization a logical value; if TRUE regularized linear regression models are
#' fitted via penalized maximum likelihood using k-fold cross-validation.
#' @param s a character."lambda.min" (default) gives the minimum mean cross-validated
#' error. The other option is lambda.1se", which gives the most regularized model such
#' that the error is within one standard error of the minimum. For regularized IPC
#' regression only.
#' @param alpha The elastic net mixing parameter with \eqn{0 \le \alpha \le 1}. alpha = 1
#' is the lasso penalty (default) and alpha = 0 the rigde penalty. For regularized IPC
#'  regression only.
#' @param weights observation weights for regularization. Can be total counts if responses are
#' proportion matrices. Default is 1 for each observation. For regularized IPC regression
#'  only.
#' @param nlambda the number of penalty terms. The default is 100. For regularized IPC
#'  regression only.
#' @param standardize a logical value; if TRUE variable are standarized prior to
#' regularization. This only affects regularization; standard/iterated IPC regression
#' coeffiecients are not standardized.
#' @param nfolds number of folds - default is 10. Although nfolds can be as large as the
#' sample size (leave-one-out cross validation), it is not recommended for large datasets.
#' Smallest value allowable is 3. For regularized IPC regression only.
#' @param linear_MxModel a logical value indicating if the structural equation model
#' contains non-linear functions of model parameters (FALSE) or not (TRUE, default). TRUE
#' speeds up the runtime of linear model. Only relevant for models fitted with
#' \pkg{OpenMx}.
#' @details \code{regularization=TRUE} adds a another list with regularized
#' linear models to the \code{ipcr} object. All arguments related to regularization are
#' passed to \code{\link[glmnet]{cv.glmnet}}. If requested, iterated IPCs will be used
#' for the regularized IPC regression models.
#' @return \code{ipcr} returns an object of \code{\link[base]{class}} "\code{ipcr}".
#'
#' The function \code{summary} prints a summary of the IPC regression equations.
#' \code{print} shows the arguments specified in the \code{ipcr} function call.
#' \code{plot} visualizes the correlation between IPCs and covariates in the form of a
#' heatmap. The generic functions \code{AIC}, \code{BIC}, \code{coef}, \code{confint},
#' \code{effects}, \code{fitted}, \code{logLik}, \code{nobs}, \code{predict},
#' \code{residuals}, \code{sigma}, and \code{vcvo} extract various information from the
#' value returned by \code{ipcr}. Heteroskedastic robust IPC regression can be performed
#' with the functions \code{\link[lmtest]{coeftest}} and \code{\link[lmtest]{coefci}}.
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
#' IPC_reg <- ipcr(fit = fit, covariates = ipcr_covariates, iterate = TRUE)
#'
#' # Get an overview about parameter differences
#' plot(IPC_reg)
#'
#' # Show IPC regression output
#' summary(IPC_reg)
#' @references
#' Arnold, M., Oberski, D. L., Brandmaier, A. M., & Voelkle, M. C. (2019). Identifying
#' heterogeneity in dynamic panel models with individual parameter contribution
#' regression. *Structural Equation Modeling, 27*, 613-628.
#' @seealso \code{\link{get_ipcs}} and \code{\link{get_scores}}
#' @export

ipcr <- function(fit, covariates = NULL, iterate = FALSE, iteration_info = FALSE,
                 conv = 0.0001, max_it = 50, regularization = FALSE, s = "lambda.min",
                 alpha = 1, weights = NULL, nlambda = 100, standardize = TRUE,
                 nfolds = 10, linear_MxModel = TRUE) {

  # Checks --------

  ## Get model data
  model_data <- get_data(fit)

  # This needs to be checked and enabled
  #check_ipcr_arguments(fit = fit, iterate = iterate, conv = conv,
  #                    max_it = max_it, linear = linear, model_data = model_data)

  ## Check covariates and transform into data.frame
  ### Covariate status:
  ### 0: everything is fine (required for iterated IPCR)
  ### 1: missing data in covariates (run standard IPCR with complete observations)
  ### 2: no covariates found (only compute IPCs)
  if (is.null(covariates)) {
    covariate_status <- 2
    info_covariates <- "no covariates found"
  } else {
    covariates <- as.data.frame(covariates)
    if (all(complete.cases(covariates))) {
      covariate_status <- 0
    } else {
      covariate_status <- 1
      warning("Missing data in covariates detected. Standard IPC regression is performed using complete observations (rows).")
    }
    if (is.null(names(covariates)) | any(is.na(names(covariates)))) {
      warning("Some covariate are not named. Renaming all covariates using the order of the data.frame.")
      colnames(covariates) <- paste0("covariate", seq_len(NCOL(covariates)))
    }
    info_covariates <- names(covariates)
  }


  # Storing object for output --------

  ## Model parameters
  param_estimates <- coef(fit)
  q <- length(param_estimates)
  param_names <- names(param_estimates)

  ## ipcr object
  IPC <- list("info" = list(name = deparse(substitute(fit)),
                            class = class(fit)[[1]],
                            parameters = param_names,
                            covariates = info_covariates,
                            iterate = iterate,
                            iterated_status = NULL,
                            iteration_info = iteration_info,
                            conv = conv,
                            max_it = max_it,
                            regularization = regularization,
                            s = s,
                            alpha = alpha,
                            weights = weights,
                            nlambda = nlambda,
                            standardize = standardize,
                            nfolds = nfolds,
                            linear_MxModel = linear_MxModel))



  # Iterated IPC regression --------

  if (iterate) {
    if(covariate_status != 0) {
      stop("Error: Iterated IPC regression requires complete covariates.")
    }
    IPC <- iterated_ipcr(fit, IPC = IPC, covariates = covariates,
                         iteration_info = iteration_info, conv = conv, max_it = max_it,
                         linear = linear)
  } else { ### Begin: perform standard IPC regression ###



    # Standard IPC regression --------

    ## Information from fitted object
    n <- nobs(fit)
    q <- length(param_estimates)

    ## Compute score
    scores <- sandwich::estfun(fit)
    bread_matrix <- sandwich::bread(fit)
    IPCs <- data.frame(matrix(param_estimates, nrow = n, ncol = q, byrow = TRUE) +
                         scores %*% t(bread_matrix))
    colnames(IPCs) <- param_names
    IPC$IPCs <- IPCs

    ## Perform IPC regression
    if (covariate_status %in% c(0, 1)) {

      ipcr_data <- cbind(IPCs, covariates)
      param_names_ipcr <- paste0("IPCs_", gsub("\\(|\\)", "", param_names))
      param_names_ipcr <- gsub(pattern = "~~",replacement = ".WITH.", x = param_names_ipcr)
      param_names_ipcr <- gsub(pattern = "~",replacement = ".ON.", x = param_names_ipcr)
      colnames(ipcr_data)[seq_len(q)] <- param_names_ipcr
      IV <- paste(colnames(covariates), collapse = " + ")
      ipcr_list <- lapply(param_names_ipcr, FUN = function(x) {
        do.call(what = "lm",
                args = list(formula = paste(x, "~", IV), data = as.name("ipcr_data")))
      })
      names(ipcr_list) <- param_names
      IPC$regression_list <- ipcr_list
    }
  } # End: perform standard IPC regression



  # Regularization --------
  if (regularization) {
    IPC$regularized_regression_list <- lapply(IPC$IPCs, FUN = function(y) {
      glmnet::cv.glmnet(x = as.matrix(covariates), y = y, alpha = alpha,
                        weights = weights, nlambda = nlambda, standardize = standardize,
                        nfolds = nfolds)})
  }



  # Prepare output --------
  IPC$output$info <- print_info(IPC$info)
  if (is.list(IPC$regression_list)) {
    IPC$output$coefficients_matrix <- coefficients_matrix(IPC)
    }
  if (is.list(IPC$regularized_regression_list)) {
    IPC$output$regularized_coefficients_matrix <- regularized_coefficents_matrix(IPC)
    }

  class(IPC) <- "ipcr"

  return(IPC)

}
