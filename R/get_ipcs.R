#' @title Calculate Individual Parameter Contributions
#' @description This function calculates the individual parameter contributions (IPCs) of
#'  a fitted model.
#' @param x a fitted model object.
#' @param analytic Logical. If \code{FALSE} (default), functions of
#' \pkg{lavaan}, \pkg{OpenMx}, or \pkg{sandwich} will be used to compute scores.
#' If \code{TRUE}, fast analytical scores are requested for a continuous
#' raw-data \code{MxRAMModel}. Unsupported OpenMx models use numerical scores
#' with a warning. See \code{\link{ipcr}} for the supported analytical scope.
#' @param ... Arguments passed to score-extraction methods.
#' @details \code{get_ipcs} is a convenience function for obtaining IPCs. The more
#' powerful \code{\link[ipcr]{ipcr}} function also provides IPCs and can be used to predict
#' differences in model parameters by regressing the IPCs on predictors.
#' @return A numeric \eqn{n * k} matrix containing the IPCs of the fitted
#' model, where \eqn{n} is the number of IPC units (clusters for mixed models)
#' and \eqn{k} is the number of fitted parameters. Its columns are named
#' according to the fitted model parameters.
#' @references
#' Arnold, M., Oberski, D. L., Brandmaier, A. M., & Voelkle, M. C. (2019). Identifying
#' heterogeneity in dynamic panel models with individual parameter contribution
#' regression. \emph{Structural Equation Modeling, 27}, 613-628. doi:
#' \href{https://doi.org/10.1080/10705511.2019.1667240}{10.1080/10705511.2019.1667240}
#' @seealso \code{\link[ipcr]{ipcr}}
#' @export

get_ipcs <- function(x, analytic = FALSE, ...) {
  analytic <- prepare_model_ipcr(x, analytic)
  check_model_convergence_ipcr(x)
  param_estimates <- coef_ipcr(x)
  scores <- estfun_ipcr(x, analytic = analytic, ...)
  bread_matrix <- bread_ipcr(x)
  ipc_components <- validate_ipcr_components(
    param_estimates = param_estimates,
    scores = scores,
    bread_matrix = bread_matrix,
    n = nobs_ipcr(x)
  )
  IPCs <- matrix(
    param_estimates,
    nrow = ipc_components$n,
    ncol = ipc_components$q,
    byrow = TRUE
  ) + ipc_components$scores %*% t(ipc_components$bread_matrix)
  colnames(IPCs) <- names(param_estimates)
  IPCs
}
