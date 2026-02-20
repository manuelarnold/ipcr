get_mx_labels <- function(fit) {
  all_labels <- unlist(lapply(fit@matrices, FUN = function(m) c(m@labels)))
  if (length(fit@algebras)) {
    all_labels <- c(all_labels, unlist(lapply(fit@algebras, function(a){
      if (!is.null(a$labels)) {c(a@labels)}})))
  }
  if (length(fit@constraints))
    all_labels <- c(all_labels, unlist(lapply(fit@constraints, function(cn)
      c(cn@result@labels))))
  if (length(fit@submodels))
    all_labels <- c(all_labels, unlist(lapply(fit@submodels, get_mx_labels)))
}

get_mx_groups <- function(fit, predictors) {

  dat <- as.data.frame(fit$data$observed)
  N <- fit$data$numObs
  vars_names <- names(dat)
  n_vars <- length(vars_names)

  # Predictor key
  predictors_key <- apply(predictors, MARGIN = 1,
                          function(r) paste(r, collapse = "|"))

  # Missing pattern key
  NA_key <- apply(is.na(dat), MARGIN = 1,
                        function(r) paste(as.integer(r), collapse = ""))

  # Definition variable key
  all_labels <- get_mx_labels(fit = fit)
  def_names <- grep(pattern = "^data\\.", x = all_labels, value = TRUE)
  def_names <- unique(sub(pattern = "^data\\.", replacement = "",
                          x = def_names))
  if (length(def_names) == 0) {
    def_key <- rep("", times = N)
  } else {
    def_dat <- dat[, def_names, drop = FALSE]
    def_key  <- apply(def_dat , MARGIN = 1,
                      function(r) paste(r, collapse = "|"))
  }

  # combine missing data pattern and definition variables
  NA_def_comb <- paste(NA_key, def_key, sep = "||")
  NA_def_comb_key <- as.numeric(factor(NA_def_comb,
                                      levels = unique(NA_def_comb)))

  # combine and factor into group IDs
  comb_key <- paste(predictors_key, NA_key, def_key, sep = "||")
  groups <- as.numeric(factor(comb_key, levels = unique(comb_key)))

  # Store UNIQUE missing patterns
  unique_NA_key <- unique(NA_key)
  obs_vars_idx  <- lapply(unique_NA_key, FUN = function(k) {
    bits <- as.integer(strsplit(k, "")[[1]])
    which(bits == 0) # indexes of non-missing variables
  })
  obs_vars_names  <- lapply(obs_vars_idx, FUN = function(idx) vars_names[idx])
  n_obs_vars <- lapply(obs_vars_names, FUN = length)

  NA_patterns <- list(
    key        = unique_NA_key, # "010" pattern strings
    vars_index  = obs_vars_idx,      # integer indices
    vars_names  = obs_vars_names,    # corresponding variable names
    n_vars = n_obs_vars,
    jacobian_idx = vector("list", length(unique_NA_key))  # preallocate
  )

  cov_matrix <- matrix(paste0("cov", rep(1:n_vars, each = n_vars), "_",
                              rep(1:n_vars, times = n_vars)),
                       nrow = n_vars, ncol = n_vars, byrow = TRUE)
  mean_vector <- matrix(paste0("mean", 1:n_vars), nrow = n_vars, ncol = 1)

  for (i in seq_len(length(NA_patterns$key))) {
    NA_patterns$jacobian_idx[[i]] <- c(
      lavaan::lav_matrix_vech(
        cov_matrix[NA_patterns$vars_index[[i]],
                   NA_patterns$vars_index[[i]],
                   drop = FALSE]),
      mean_vector[NA_patterns$vars_index[[i]]])
  }

  # Output
  list(
    groups = groups, # length N vector
    NA_pattern_ID = match(NA_key, unique_NA_key), # pattern id per row
    NA_def_comb = NA_def_comb,
    NA_def_comb_key = NA_def_comb_key,
    NA_patterns = NA_patterns  # list described above
  )
}

get_manifest_vars_from_MxModel <- function(fit) {
  manifest_vars <- fit$manifestVars
  if (length(manifest_vars) > 0) {
    return(manifest_vars)
  }
  manifest_vars <- fit$expectation$dims
  if (length(manifest_vars) > 0) {
    return(manifest_vars)
  } else {
    stop("No manifest variables found.")
  }
}
