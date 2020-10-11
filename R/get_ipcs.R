#' @title Calcuate Individual Parameter Contributions
#' @description This functions calculates the individual parameter contributions of a
#' fitted model.
#' @param x a fitted model object.
#' @return Returns a data.frame containing the individual parameter contributions of a
#' fitted model.
#' @examples
#' @export

get_ipcs <- function(x) {
  param_estimates <- coef(x)
  scores <- sandwich::estfun(x)
  bread_matrix <- sandwich::bread(x)
  IPCs <- data.frame(matrix(param_estimates, nrow = nobs(x),
                            ncol = length(param_estimates), byrow = TRUE) +
                       scores %*% t(bread_matrix))
  colnames(IPCs) <- names(param_estimates)
  IPCs
}
