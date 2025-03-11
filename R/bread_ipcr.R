#' @noRd
bread_ipcr <- function(x, ...) {
  UseMethod("bread_ipcr")
}

#' @noRd
bread_ipcr.default <- function(x, ...) {
  sandwich::bread(x, ...)
}

#' @noRd
bread_ipcr.lm <- function(x, ...) {
  bread_matrix <- sandwich::bread(x)
  dims <- NROW(bread_matrix)
  res <- matrix(0, nrow = dims + 1, ncol = dims + 1)
  res[1:dims, 1:dims] <- bread_matrix
  res[dims + 1, dims + 1] <- 2 * mean(x$residuals^2)
  rownames(res) <- colnames(res) <- c(names(stats::coef(x)), "residual")
  res
}

#' @noRd
bread_ipcr.lmerMod <- function(x, ...) {
  as.matrix(merDeriv::bread.lmerMod(x, full = TRUE))
}

#' @noRd
bread_ipcr.glmerMod <- function(x, ...) {
  as.matrix(merDeriv::bread.glmerMod(x, full = TRUE))
}

#' @noRd
bread_ipcr.lavaan <- function(x, ...) {
  if (x@Model@eq.constraints) {
    K <- eval(parse(text = "lavaan:::lav_constraints_R2K(x@Model)"))
    res <- solve(t(K) %*% lavaan::lavInspect(x, what = "information.expected") %*% K)
  } else {
    res <- solve(lavaan::lavInspect(x, what = "information.expected"))
  }
  res
}

#' @noRd
bread_ipcr.MxRAMModel <- function(x, ...) {
  x$output$vcov * x$data$numObs
}
