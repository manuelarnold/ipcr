#' @title Extracts Scores
#' @description This functions extracts the scores of a fitted model. Arguments are passed
#' to \link[sandwich]{estfun}.
#' @param x a fitted model object.
#' @param ... arguments passed to methods.
#' @return Returns a matrix containing the scores of a fitted model.
#' @examples
#' @export get_scores

get_scores <- function(x, ...) {
  data.frame(sandwich::estfun(x, ...))
}
