#' @title Calcuate Individual Parameter Contributions
#' @description This functions calculates the individual parameter contributions (IPCs) of
#'  a fitted model.
#' @param x a fitted model object.
#' @details \code{get_ipcs} is a convenience function for obtaining IPCs. The more
#' powerful \code{\link[ipcr]{ipcr}} function also provides IPCs and can be used to predict
#' differences in model parameters by regressing the IPCs on covariates.
#' @return A \code{data.frame} containing the of a fitted model. The dimensions of this
#' \code{data.frame} are \eqn{n * k}, where \eqn{n} denotes the number of observations and
#' \eqn{k} the number of parameters. The columns should be named as in
#' \code{\link[stats]{coef}}.
#' @references
#' Arnold, M., Oberski, D. L., Brandmaier, A. M., & Voelkle, M. C. (2019). Identifying
#' heterogeneity in dynamic panel models with individual parameter contribution
#' regression. \emph{Structural Equation Modeling, 27}, 613-628. doi:
#' \href{https://doi.org/10.1080/10705511.2019.1667240}{10.1080/10705511.2019.1667240}
#' @seealso \code{\link[ipcr]{ipcr}}
#' @export

get_ipcs <- function(x) {
  param_estimates <- coef_ipcr(x)
  scores <- estfun_ipcr(x)
  bread_matrix <- bread_ipcr(x)
  IPCs <- data.frame(matrix(param_estimates, nrow = nobs(x),
                            ncol = length(param_estimates), byrow = TRUE) +
                       scores %*% t(bread_matrix))
  colnames(IPCs) <- names(param_estimates)
  IPCs
}
