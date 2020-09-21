#' @title Iterated Individual Parameter Contribution Regression
#' @description This function performs iterated individual parameter contribution regression for fitted MxRAM-type models.
#' @param x a fitted model object
#' @... arguments passed to methods.
#' @return Returns a list.
#' @export

iterated_ipcr.lavaan <- function(x, IPC, covariates, conv, max_it, linear, ...) {

  obs_data <- lavaan::lavInspect(object = x, what = "data")
  covariates_matrix <- covariates
  n <- x@SampleStats@ntotal
  p <- x@pta$nvar
  ms <- x@Model@meanstructure
  param_estimates <- coef(x)
  param_names <- IPC$info$parameters <- names(param_estimates)
  q <- length(param_estimates)
  k <- NCOL(covariates)
  indices_param <- seq_len(q)


  # Prepare stuff for model syntax
  ParTable <- as.data.frame(x@ParTable)
  ParTable[ParTable$op == "~1", "rhs"] <- 1
  ParTable[ParTable$op == "~1", "op"] <- "~"
  fixed_params <- ParTable$free == 0
  values <- rep(NA, times = NROW(ParTable))
  values[fixed_params] <- paste0(ParTable[fixed_params, "est"], "*")


  # Initial IPCs
  scores <- lavaan::estfun.lavaan(x)
  bread_matrix <- bread.lavaan(x)
  IPCs <- data.frame(matrix(param_estimates, nrow = n, ncol = q, byrow = TRUE) +
                       scores %*% t(bread_matrix))
  colnames(IPCs) <- param_names
  IPC$IPCs <- IPCs


  # Initial IPC regression
  ipcr_data <- cbind(IPCs, covariates_matrix)
  param_names_ipcr <- paste0("IPCs_", gsub("\\(|\\)", "", param_names))
  param_names_ipcr <- gsub(pattern = "~~",replacement = ".WITH.", x = param_names_ipcr)
  param_names_ipcr <- gsub(pattern = "~",replacement = ".ON.", x = param_names_ipcr)
  IV <- paste(colnames(covariates), collapse = " + ")
  colnames(ipcr_data)[seq_len(q)] <- param_names_ipcr
  ipcr_list <- list()
  for (j in indices_param) {
    ipcr_list[[j]] <- do.call(what = "lm",
                              args = list(formula = paste(param_names_ipcr[j], "~", IV),
                                          data = as.name("ipcr_data")))
  }
  names(ipcr_list) <- param_names
  IPC$regression <- ipcr_list


  #### Iterated IPCR ####
  IPC$info$conv <- conv


  # Storing objects for the updating procedure
  it_est <- matrix(sapply(X = ipcr_list, FUN = function(x) {coef(x)}),
                   nrow = 1, ncol = q * (k + 1))
  it_se <- matrix(sapply(X = ipcr_list, FUN = function(x) {sqrt(diag(vcov(x)))}),
                  nrow = 1, ncol = q * (k + 1))

  # Center moment deviations at the covariate
  ### THIS DIFFERENT FROM OPENMX ###
  #center_reg_list <- apply(X = data_obs, MARGIN = 2, FUN = function(y) {
  #  lm(y ~ covariates_matrix)})
  #data_centered <- data_obs -
  #  sapply(X = center_reg_list, FUN = function(x) {x$fitted.values})
  #cent_md <- matrix(data = apply(X = data_centered, MARGIN = 1,
  #                               FUN = function(x) {lav_matrix_vech(x %*% t(x))}),
  #                  nrow = n, ncol = p_star, byrow = TRUE)
  #if (ms) {
  #  cent_md <- cbind(cent_md, data_obs)
  #}


  # Groups of observations with same covariate values in the covariate
  group <- transform(covariates, group_ID = as.numeric(interaction(covariates, drop = TRUE)))$group_ID


  # Start the iteration process
  nr_iterations <- 0
  difference <- rep(x = conv + 1, times = NCOL(it_est))
  breakdown <- FALSE
  updated_IPCs <- matrix(NA, nrow = n, ncol = q)
  colnames(updated_IPCs) <- param_names
  unique_groups <- unique(group)


  while(nr_iterations < max_it & isFALSE(all(abs(difference) < conv)) &
        breakdown == FALSE) {

    ### Update IPCs ###
    # Update the IPCs of individual or group i
    for (i in unique_groups) { # Start loop with index i
      ID_group <- which(group == i) # HERE, THIS IS WRONG
      n_group <- length(ID_group)
      IPC_pred <- as.numeric(c(1, covariates[group == i, , drop = FALSE][1, ])) #### Check here


      # Update Parameters
      for (j in indices_param){
        param_estimates[j] <- sum(coef(ipcr_list[[j]]) * IPC_pred)
      }


      # Update model synatx
      values[!fixed_params] <- paste0("start(", param_estimates, ")*")
      model_syntax <- paste0(ParTable$lhs, ParTable$op, values, ParTable$rhs, collapse = "\n")


      # Update lavaan model
      x <- lavaan::lavaan(model = model_syntax, data = obs_data, do.fit = FALSE)


      # Update IPCs
      scores_ID_group <- lavaan::estfun.lavaan(x)[ID_group, , drop = FALSE]
      bread_matrix <- bread.lavaan(x)
      IPCs_ID_group <- matrix(param_estimates, nrow = n_group, ncol = q, byrow = TRUE) +
        scores_ID_group %*% t(bread_matrix)
      updated_IPCs[ID_group, ] <- IPCs_ID_group
    } # end loop with index i: Find observations with identical covariates


    if (class(updated_IPCs)[1] == "try-error") {
      breakdown <- TRUE
    }

    if (breakdown == FALSE) {

      # Estimate updated IPC regression parameter
      ipcr_data <- cbind(updated_IPCs, covariates)
      colnames(ipcr_data)[indices_param] <- param_names_ipcr
      for (j in indices_param) {
        ipcr_list[[j]] <- do.call(what = "lm",
                                  args = list(formula = paste(param_names_ipcr[j], "~", IV),
                                              data = as.name("ipcr_data")))
      }


      # Store results
      it_est <- rbind(it_est, c(sapply(X = ipcr_list,
                                       FUN = function(x) {coef(x)})))
      it_se <- rbind(it_se, c(sapply(X = ipcr_list,
                                     FUN = function(x) {sqrt(diag(vcov(x)))})))

      # Covergence criteria
      difference <- it_est[nrow(it_est), ] -
        it_est[nrow(it_est) - 1, ]
      nr_iterations <- nr_iterations + 1
      cat("Iteration:", nr_iterations, "\n")

    }
  }

  # Error reports
  if (breakdown == TRUE) {
    cat("The iterated IPC regression aborted prematurely.")
    IPC$info$iterated <- "Iterated IPC regression aborted prematurely."
  }

  if(breakdown == FALSE & nr_iterations == max_it &
     isFALSE(all(abs(difference) <= conv))) {
    cat(paste("The iterated IPC regression algorithm did not converge after", max_it, "iterations. Consider to increase the maximum number of iterations."))
    IPC$info$iterated <- paste("Iterated IPC regression reached the maximum number of", max_it, "iterations without converging.")
  }

  if(breakdown == FALSE & nr_iterations < max_it &
     isFALSE(all(abs(difference) > conv))) {
    cat("Iterated IPC regression converged.")
    IPC$info$iterate <- paste("Iterated IPC regression converged succesfully after", nr_iterations, "iterations.")
    IPC$iterated_IPCs <- as.data.frame(updated_IPCs)
    IPC$iterated_regression <- ipcr_list
  }

  IPC

}
