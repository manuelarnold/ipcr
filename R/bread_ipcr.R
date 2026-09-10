#' @noRd
bread_ipcr <- function(x, ...) {
  UseMethod("bread_ipcr")
}

#' @noRd
#' @export
bread_ipcr.default <- function(x, ...) {
  sandwich::bread(x, ...)
}

#' @noRd
#' @export
bread_ipcr.lm <- function(x, ...) {
  bread_matrix <- sandwich::bread(x)
  dims <- NROW(bread_matrix)
  res <- matrix(0, nrow = dims + 1, ncol = dims + 1)
  res[1:dims, 1:dims] <- bread_matrix
  residual_variance <- mean(x$residuals^2)
  res[dims + 1, dims + 1] <- 2 * residual_variance^2
  rownames(res) <- colnames(res) <- c(names(stats::coef(x)), "residual")
  res
}

#' @noRd
#' @export
bread_ipcr.glm <- function(x, ...) {
  sandwich::bread(x, ...)
}

#' @noRd
#' @export
bread_ipcr.lmerMod <- function(x, ...) {
  as.matrix(merDeriv::bread.lmerMod(x, full = TRUE))
}

#' @noRd
#' @export
bread_ipcr.glmerMod <- function(x, ...) {
  as.matrix(merDeriv::bread.glmerMod(x, full = TRUE))
}

#' @noRd
#' @export
bread_ipcr.lavaan <- function(x, ...) {
  information <- as.matrix(
    lavaan::lavInspect(x, what = "information.expected")
  )
  reduction <- lavaan_constraint_reduction_ipcr(x@Model)

  if (is.null(reduction)) {
    res <- solve(information)
  } else {
    K <- reduction$K
    res <- solve(t(K) %*% information %*% K)
  }

  parameter_names <- names(coef_ipcr(x))
  dimnames(res) <- list(parameter_names, parameter_names)
  res
}

# Reduce lavaan's full free-parameter coordinates to the coordinates used for
# simple equality constraints. The first branch is adapted from lavaan's
# lav_con_r2k() (lavaan is GPL-2 or later; ipcr is GPL-3).
lavaan_constraint_reduction_ipcr <- function(lavmodel) {
  if (isTRUE(lavmodel@eq.constraints)) {
    m_r <- lavmodel@ceq.JAC
    if (is.null(m_r) || NROW(m_r) == 0L) {
      stop(
        "Unable to obtain the lavaan equality-constraint Jacobian.",
        call. = FALSE
      )
    }

    npar_full <- NCOL(m_r)
    m_k <- diag(npar_full)
    for (i in seq_len(NROW(m_r))) {
      idx1 <- which(m_r[i, ] == 1)
      idx2 <- which(m_r[i, ] == -1)
      if (length(idx1) != 1L || length(idx2) != 1L ||
          any(!m_r[i, ] %in% c(-1, 0, 1))) {
        stop(
          "Only simple linear equality constraints are supported for lavaan models.",
          call. = FALSE
        )
      }
      m_k[idx2, idx1] <- 1
    }
    neg_idx <- which(colSums(m_r) < 0)
    m_k <- m_k[, -neg_idx, drop = FALSE]
  } else if (isTRUE(lavmodel@ceq.simple.only)) {
    m_k <- lavmodel@ceq.simple.K
  } else {
    return(NULL)
  }

  keep <- vapply(
    seq_len(NCOL(m_k)),
    function(column) which(m_k[, column] != 0)[1],
    integer(1)
  )
  list(K = m_k, keep = keep)
}

#' @noRd
#' @export
bread_ipcr.MxModel <- function(x, ...) {
  x$output$vcov * nobs_ipcr(x)
}
