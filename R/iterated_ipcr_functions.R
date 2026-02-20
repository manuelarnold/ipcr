#' @noRd
# get_RAM_matrices ----
# manifest first, then latent
get_RAM_matrices <- function(x, ...) {
  UseMethod("get_RAM_matrices")
}

#' @noRd
get_RAM_matrices.lavaan <- function(x, ...) {
  pt <- lavaan::parameterTable(x)
  ov <- lavaan::lavNames(x, type = "ov")   # manifest variables
  lv <- lavaan::lavNames(x, type = "lv")   # latent variables
  all_vars <- c(ov, lv)               # manifest first, then latent
  n_all <- length(all_vars)
  n_obs <- length(ov)

  # Initialize matrices
  A_values <- A_labels <- matrix(0, nrow = n_all, ncol = n_all,
                                 dimnames = list(all_vars, all_vars))
  A_labels[, ] <- NA
  S_values <- S_labels <- matrix(0, nrow = n_all, ncol = n_all,
                                 dimnames = list(all_vars, all_vars))
  S_labels[, ] <- NA
  m_values <- m_labels <- matrix(0, nrow = 1, ncol = n_all,
                                 dimnames = list(NULL, all_vars))
  m_labels[, ] <- NA

  # Fill matrices
  for (i in seq_len(nrow(pt))) {
    row <- pt[i, ]
    val <- row$est
    label <- if (!is.na(row$label) && nzchar(row$label)) {
      row$label
    } else if (row$free > 0) {
      paste0(row$lhs, row$op, row$rhs)
    } else {
      NA
    }

    if (row$lhs %in% all_vars && row$rhs %in% all_vars) {
      if (row$op == "~") {
        A_values[row$lhs, row$rhs] <- val
        A_labels[row$lhs, row$rhs] <- label
      } else if (row$op == "=~") {
        A_values[row$rhs, row$lhs] <- val
        A_labels[row$rhs, row$lhs] <- label
      } else if (row$op == "~~") {
        S_values[row$lhs, row$rhs] <- val
        S_values[row$rhs, row$lhs] <- val
        S_labels[row$lhs, row$rhs] <- label
        S_labels[row$rhs, row$lhs] <- label
      }
    }

    if (row$op == "~1" && row$lhs %in% all_vars) {
      m_values[1, row$lhs] <- val
      m_labels[1, row$lhs] <- label
    }
  }

  # F matrix: selects manifest variables from full variable list
  F_values <- cbind(diag(1, nrow = n_obs),
                    matrix(0, nrow = n_obs,ncol = length(lv)))
  rownames(F_values) <- ov
  colnames(F_values) <- all_vars
  F_labels <- matrix(NA, nrow = n_obs, ncol = n_all,
                     dimnames = list(ov, all_vars))

  list(
    par = names(coef_ipcr.lavaan(x)),
    ov = ov,
    lv = lv,
    A = list(values = A_values, labels = A_labels),
    S = list(values = S_values, labels = S_labels),
    F = list(values = F_values, labels = F_labels),
    m = list(values = m_values, labels = m_labels)
  )
}

#' @noRd
get_RAM_matrices.MxModel <- function(x, ...) {

  # Give pseudo-labels to matrices if unlabeled parameters are given
  # These labels should follow the standard OpenMx naming convention
  # If not, derivatives will not be correct
  candidate_param_id <- which(startsWith(x = names(x$output$estimate),
                                         prefix = x$name))
  if (length(candidate_param_id) > 0) {
    for (k in candidate_param_id) {
      candidate_param_name <- names(x$output$estimate)[k]
      cplen <- nchar(x$name)
      candidate_matrix <- substr(candidate_param_name, start = cplen + 2,
                                 stop = cplen + 2)
      candidate_pos <- as.integer(strsplit(substr(
        candidate_param_name,
        start = cplen + 4,
        stop = nchar(candidate_param_name) - 1), ",")[[1]])
      if (candidate_matrix == "A") {
        x$A$labels[candidate_pos[1], candidate_pos[2]] <- candidate_param_name
      } else if (candidate_matrix == "S") {
        x$S$labels[candidate_pos[1], candidate_pos[2]] <- candidate_param_name
      } else if (candidate_matrix == "M") {
        x$M$labels[candidate_pos] <- candidate_param_name
      }
    }
  }

  list(
    par = names(x$output$estimate),
    ov = x$manifestVars,
    lv = x$latentVars,
    A = list(values = x$A$values, labels = x$A$labels),
    S = list(values = x$S$values, labels = x$S$labels),
    F = list(values = x$F$labels, labels = x$F$labels),
    m = list(values = x$m$labels, labels = x$m$labels)
  )
}


#' @noRd
# reorder_RAM_matrices ----
reorder_RAM_matrices <- function(x) {
  # Reorder manifest and latent variables with respect of the number of
  # associated covariance parameters.
  # Re-ordering increases the chance that that 1 to 1 mapping is not lost due to
  # the Cholesky decomposition
  # Count non-NA covariance entries in S$labels for each variable
  S_labels <- x$S$labels
  is_cov <- !is.na(S_labels) & upper.tri(S_labels, diag = TRUE)

  cov_counts <- rowSums(is_cov, na.rm = TRUE)

  # Separate and sort variables
  ov_counts <- cov_counts[x$ov]
  lv_counts <- cov_counts[x$lv]

  ov_sorted <- names(sort(ov_counts, decreasing = TRUE))
  lv_sorted <- names(sort(lv_counts, decreasing = TRUE))

  new_order <- c(ov_sorted, lv_sorted)

  # Reorder A
  A_values <- x$A$values[new_order, new_order, drop = FALSE]
  A_labels <- x$A$labels[new_order, new_order, drop = FALSE]

  # Reorder S
  S_values <- x$S$values[new_order, new_order, drop = FALSE]
  S_labels <- x$S$labels[new_order, new_order, drop = FALSE]

  # Reorder F
  F_values <- x$F$values[, new_order, drop = FALSE]
  F_labels <- x$F$labels[, new_order, drop = FALSE]

  # Reorder m
  m_values <- x$m$values[new_order]
  m_labels <- x$m$labels[new_order]

  list(
    A = list(values = A_values, labels = A_labels),
    S = list(values = S_values, labels = S_labels),
    F = list(values = F_values, labels = F_labels),
    m = list(values = m_values, labels = m_labels),
    par = x$par,
    ov = ov_sorted,
    lv = lv_sorted
  )
}


#' @noRd
# compute_LDL_matrices ----
# LDL decomposition
# S = LDL^T
compute_LDL_matrices <- function(S) {
  n <- nrow(S)
  L <- diag(1, nrow = n)
  D <- numeric(length = n)

  for (j in 1:n) {
    if (j == 1) {
      D[j] <- S[j, j]
    } else {
      sum_ld <- sum(L[j, 1:(j-1)]^2 * D[1:(j-1)])
      D[j] <- S[j, j] - sum_ld
    }
    for (i in (j+1):n) {
      if (i > n) break
      sum_l <- sum(L[i, 1:(j-1)] * L[j, 1:(j-1)] * D[1:(j-1)])
      L[i, j] <- (S[i, j] - sum_l) / D[j]
    }
  }
  list(L = L, D = diag(D))
}


# LDL_transform_RAM_matrices ----
LDL_transform_RAM_matrices <- function(x) {
  LDL_matrices <- compute_LDL_matrices(x$S$values)

  # Check if there is a 1 to 1 mapping
  if (!identical(LDL_matrices$L[lower.tri(LDL_matrices$L)] == 0,
                 is.na(x$S$labels[lower.tri(x$S$labels)]))) {
    warning("Covariance structure is too complex for LDL transformation.
            Switching to log-constrained variance parameters.")
    return("log")
  }

  L_labels <- x$S$labels
  L_labels[!lower.tri(x$S$labels)] <- NA

  D_labels <- x$S$labels
  D_labels[lower.tri(D_labels)] <- NA
  D_labels[upper.tri(D_labels)] <- NA

  list(
    A = list(values = A_values, labels = A_labels),
    L = list(values = LDL_matrices$L, labels = L_labels),
    D = list(values = LDL_matrices$D, labels = D_labels),
    F = list(values = F_values, labels = F_labels),
    m = list(values = m_values, labels = m_labels),
    par = x$par,
    ov = x$ov,
    lv = x$lv
  )
}


# add_RAM_unc_deriv ----
add_RAM_unc_deriv <- function(x) {
  n_par <- length(x$par)
  par_seq <- seq_len(n_par)
  n_ov <- length(x$ov)

  Zero <- matrix(0, nrow = n_ov, ncol = n_ov)
  A_deriv <- lapply(par_seq, function(x) {Zero})
  S_deriv <- A_deriv
  zero <- matrix(0, nrow = n_ov, ncol = 1)
  m_deriv <- lapply(par_seq, function(x) {zero})

  for (i in par_seq) {
    A_deriv[[i]][which(x$A$labels == x$par[i], arr.ind = TRUE)] <- 1
  }

  for (i in par_seq) {
    S_deriv[[i]][which(x$S$labels == x$par[i], arr.ind = TRUE)] <- 1
  }

  for (i in par_seq) {
    m_deriv[[i]][which(x$m$labels == x$par[i])] <- 1
  }

  x$A$deriv <- A_deriv
  x$S$deriv <- S_deriv
  x$m$deriv <- m_deriv
  x
}


# initialize_RAM_log_deriv ----
add_RAM_log_deriv <- function(x) {
  n_par <- length(x$par)
  par_seq <- seq_len(n_par)
  n_ov <- length(x$ov)

  Zero <- matrix(0, nrow = n_ov, ncol = n_ov)
  A_deriv <- lapply(par_seq, function(x) {Zero})
  S_deriv <- A_deriv
  zero <- matrix(0, nrow = n_ov, ncol = 1)
  m_deriv <- lapply(par_seq, function(x) {zero})

  for (i in par_seq) {
    A_deriv[[i]][which(x$A$labels == x$par[i], arr.ind = TRUE)] <- 1
  }

  for (i in par_seq) {
    row_col <- which(x$S$labels == x$par[i], arr.ind = TRUE)
    if (row_col[1, 1] == row_col[1, 2]) {
      S_deriv[[i]][row_col] <- x$S$values[row_col]
    } else {
      S_deriv[[i]][row_col] <- 1
    }
  }

  for (i in par_seq) {
    m_deriv[[i]][which(x$m$labels == x$par[i])] <- 1
  }

  x$A$deriv <- A_deriv
  x$S$deriv <- S_deriv
  x$m$deriv <- m_deriv
  x
}


# add_RAM_LDL_deriv ----
add_RAM_LDL_deriv <- function(x) {
  n_par <- length(x$par)
  par_seq <- seq_len(n_par)
  n_ov <- length(x$ov)

  Zero <- matrix(0, nrow = n_ov, ncol = n_ov)
  A_deriv <- lapply(par_seq, function(x) {Zero})
  L_deriv <- D_deriv <- A_deriv
  zero <- matrix(0, nrow = n_ov, ncol = 1)
  m_deriv <- lapply(par_seq, function(x) {zero})

  for (i in par_seq) {
    A_deriv[[i]][which(x$A$labels == x$par[i], arr.ind = TRUE)] <- 1
  }

  for (i in par_seq) {
    L_deriv[[i]][which(x$L$labels == x$par[i], arr.ind = TRUE)] <- 1
  }

  for (i in par_seq) {
    row_col <- which(x$D$labels == x$par[i], arr.ind = TRUE)
    D_deriv[[i]][row_col] <- x$D$values[row_col]
  }

  for (i in par_seq) {
    m_deriv[[i]][which(x$m$labels == x$par[i])] <- 1
  }

  x$A$deriv <- A_deriv
  x$L$deriv <- L_deriv
  x$D$deriv <- D_deriv
  x$m$deriv <- m_deriv
  x
}
