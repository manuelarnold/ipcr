#' @title Partition the Samples Into Similar Groups
#' @description Internal function that partitions the samples in groups with the same values of covariates.
#' @param x list with information about the model.
#' @param md matrix with individual deviations from the sample moments.
#' @return Matrix with groupings.
#' @export

get_groups <- function(x) {
  group <- matrix(data = NA, nrow = x$N, ncol = 2)
  group[, 1] <- 1:x$N # ID
  counter_group <- 1
  for (i in 1:x$N) {
    if (is.na(group[i, 2])) {
      for (j in i:x$N) {
        if (all(x$covariates[j, , drop = FALSE] == x$covariates[i, , drop = FALSE])) {
          group[j, 2] <- counter_group
        }
      }
      counter_group <- counter_group + 1
    }
  }
  return(group)
}
