#' @title Individual parameter contribution regression
#' @description Explain parameter heterogeneity in a structural equation model (SEM) with individual parameter contribution (IPC) regression.
#' @param fit a SEM estimated with OpenMx. The model must contain an estimated covariance structure. The mean structure is optional.
#' @param formula a regression equation formula. Interaction between regressors and polynominals are not implemented yet. You can add these terms as new regressors to the covariates data frame.
#' @param covariates a data frame or matrix with the IPC predictor variables.
#' @param iterate a logical value; if TRUE (default) iterated IPC regression is performed after static IPC regression.
#' @param convergence_criterion an integer used as stopping criterion for iterated IPC regression. Criterion is the largest difference in all parameter estimate between iterations.
#' @param max_it the maximum number of iterations.
#' @return Returns a new ipcr object.
#' @examples
#' # Specify an OpenMx model
#' Model <- mxModel(manifestVars = c("x1", "x2", "x3", "y1", "y2", "y3"),
#'                  latentVars = c("xi", "eta"),
#'                  type = "RAM",
#'                  mxData(observed = ipcr_data, type = "raw"),
#'                  mxPath(from = "xi", to = c("x1", "x2", "x3"),
#'                         connect = "single", arrows = 1,
#'                         free = c(FALSE, TRUE, TRUE), values = 1,
#'                         labels = c("l11", "l22", "l33")),
#'                  mxPath(from = "eta", to = c("y1", "y2", "y3"),
#'                  connect = "single", arrows = 1, free = c(FALSE, TRUE, TRUE),
#'                  values = 1, labels = c("l44", "l55", "l66")),
#'                  mxPath(from = "xi", to = "eta", arrows = 1, free = TRUE,
#'                  values = 0.5, labels = "beta"),
#'                  mxPath(from = c("xi", "eta", "x1", "x2", "x3", "y1", "y2", "y3"),
#'                         connect = "single", arrows = 2, free = TRUE,
#'                         values = c(0.75, 0.555, rep(0.25, 6)),
#'                         labels = c("phi11", "phi22", "e11", "e22", "e33", "e44", "e55", "e66")),
#'                  mxPath(from = "one", to = c("x1", "x2", "x3", "y1", "y2", "y3"),
#'                         connect = "single", arrows = 1, free = TRUE, values = 0,
#'                          labels = c("m_x1", "m_x2", "m_x3", "m_y1", "m_y2", "m_y3"))
#' )
#'
#' # Fit the model
#' Fit <- mxTryHard(model = Model)
#'
#' # Run the model
#' IPC_reg <- ipcr(fit = Fit, formula = ~ z1 + z2, covariates = ipcr_data, iterate = TRUE)
#'
#' # Display the iterated individual parameter contribution regression estimates of the
#' regression coefficient beta
#' summary(object = IPC_reg, parameter = "beta", method = "iterated")
#' @export ipcr

ipcr <- function(fit, formula = NULL, covariates = NULL, iterate = FALSE,
                 convergence_criterion = 0.0001, max_it = 50) {

  # Check for missing data SEM data
  # Compute the IPCs if the data is complete
  data_obs <- fit$data$observed[, fit$manifestVars, drop = FALSE]

  if (!any(is.na(data_obs))) {

    ### Calculate the initial IPCs

    # Basic information about the data set and model
    N <- fit$output$data[[1]][[1]]
    p <- length(fit$manifestVars)
    p_star <- (p * (p + 1)) / 2
    p_star_means <- p * (p + 3) / 2
    mean_structure <- any(fit$M$free)
    param_estimates <- omxGetParameters(model = fit)
    param_names <- names(param_estimates)
    q <- length(param_estimates)
    covariates <- as.matrix(covariates)

    # Calculate the individual deviations between observed and estimated moments
    centered_data <- scale(x = data_obs, center = TRUE, scale = FALSE)
    moment_contributions <- matrix(data = apply(X = centered_data, MARGIN = 1,
                                                FUN = function (x) {matrixcalc::vech(x %*% t(x))}),
                                   nrow = N, ncol = p_star, byrow = TRUE)
    expected_cov <- mxGetExpected(model = fit, component = "covariance")
    vech_cov_matrix <- matrix(data = rep(x = matrixcalc::vech(expected_cov), times = N),
                              byrow = TRUE, nrow = N, ncol = p_star)
    moment_deviations <- moment_contributions - vech_cov_matrix
    if (mean_structure) {
      expected_means <- mxGetExpected(model = fit, component = "means")
      means_matrix <- matrix(data = rep(x = expected_means, times = N),
                             byrow = TRUE, nrow = N, ncol = p)
      mean_dev <- data_obs - means_matrix
      moment_deviations <- cbind(moment_deviations, mean_dev)
    }
    moment_deviations <- as.matrix(x = moment_deviations)

    # Get Jacobian, weight matrix, and transformation matrix
    jacobian <- suppressMessages(expr = mxCheckIdentification(model = fit, details = FALSE)$jacobian)
    if (mean_structure == FALSE) {
      jacobian <- jacobian[1:p_star, ]
    }
    expected_cov_inv <- solve(expected_cov)

    # This can later be replaced with the lavaan duplication matrix (it can handle p = 1)
    if (p == 1) {
      D <- matrix(data = 1)
    } else {
      D <- duplication.matrix(n = p)
    }
    V <- 0.5 * t(D) %*% kronecker(X = expected_cov_inv, Y = expected_cov_inv) %*% D
    if (mean_structure) {
      V_m_cov <- matrix(data = 0, nrow = p_star_means, ncol = p_star_means)
      V_m_cov[1:p_star, 1:p_star] <- V
      V_m_cov[(p_star + 1):p_star_means, (p_star + 1):p_star_means] <- expected_cov_inv
      V <- V_m_cov
    }
    W <- solve(t(jacobian) %*% V %*% jacobian) %*% t(jacobian) %*% V

    # Calculate the static IPCs
    static_IPC <- data.frame(matrix(data = rep(x = param_estimates, times = N),
                                    byrow = TRUE, nrow = N, ncol = q) + moment_deviations %*% t(W))
    colnames(static_IPC) <- param_names

    # Store the results
    IPC <- list("StaticIPCs" = static_IPC)

    # Status report
    cat("Static IPCs computed\n")


    ### Estimate the initial individual IPC regression parameters

    #if(!is.null(formula) & !is.null(covariates)) { switch this on !!!!! bracket closes at the end

    # Get regression equations
    if (class(formula) == "formula") {
      formula <- Reduce(paste, deparse(formula))
    }
    formula <- gsub(pattern = " ", replacement = "", x = formula, fixed = TRUE)
    DV <- strsplit(x = formula, split = "~")[[1]][1]
    IV <- strsplit(x = formula, split = "~")[[1]][2]
    DV_separated <- strsplit(x = DV, split = "\\+")
    IV_separated <- strsplit(x = IV, split = "\\+")[[1]]
    if (identical(DV_separated[[1]], character(0))) {
      target_parameters <- param_names
    } else {
      target_parameters <- DV_separated
    }
    IPC_reg_eq <- c()
    IPC_reg_list <- list()

    # Select covariates from the provided data set
    covariates <- covariates[, IV_separated, drop = FALSE]

    # Check for missing data in the covariates
    # Compute the IPC regression parameter estimates if the data is complete
    if (!any(is.na(covariates))) {

      # Run initial IPC regression
      IPC_reg_data <- cbind(static_IPC, covariates)
      for (i in 1:length(target_parameters)) {
        IPC_reg_eq[i] <- paste(target_parameters[i], "~", IV)
        IPC_reg_list[[i]] <- do.call(what = "lm", args = list(formula = as.formula(IPC_reg_eq[i]),
                                                              data = as.name("IPC_reg_data")))
      }
      names(IPC_reg_list) <- target_parameters

      # Store the results
      IPC[[2]] <- IPC_reg_list
      names(IPC)[2] <- "StaticIPCRegression"

      # Status report
      cat("Static IPC regression parameters estimated\n")


      ### Estimate the iterated individual IPC regression parameters

      if(iterate) {

        # Calculate individual- or group-specific contributions to the observed moments
        igsmc_reg <- list()
        data_igsmc <- data.frame(cbind(data_obs, covariates))
        for (i in 1:p) {
          igsmc_reg[[i]] <- do.call(what = "lm",
                                    args = list(formula = as.formula(paste0(fit$manifestVars[i], "~", IV)),
                                                data = as.name("data_igsmc")))
        }
        data_igscentered <- data_obs - sapply(X = igsmc_reg, FUN = function(x){x$fitted.values})
        mom_cont_igsmc <- matrix(data = apply(X = data_igscentered, MARGIN = 1,
                                              FUN = function (x) {matrixcalc::vech(x %*% t(x))}),
                                 nrow = N, ncol = p_star, byrow = TRUE)
        if (mean_structure) {
          mom_cont_igsmc <- cbind(mom_cont_igsmc, data_obs)
        }

        # Storing matrices for the updating procedure
        estimates_matrix <- matrix(data = c(sapply(X = IPC_reg_list, FUN = function(x){coef(x)})),
                                   nrow = 1, ncol = q * (length(IV_separated) + 1))
        SE_matrix <- matrix(data = c(sapply(X = IPC_reg_list, FUN = function(x){sqrt(diag(vcov(x)))})),
                            nrow = 1, ncol = q * (length(IV_separated) + 1))

        # Partition the samples in groups with the same value of the covariates
        group <- matrix(data = NA, nrow = N, ncol = 2)
        group[, 1] <- 1:N # ID
        counter_group <- 1
        for (i in 1:N) {
          if (is.na(group[i, 2])) {
            for (j in i:N) {
              if (all(covariates[j, , drop = FALSE] == covariates[i, , drop = FALSE])) {
                group[j, 2] <- counter_group
              }
            }
            counter_group <- counter_group + 1
          }
        }

        # Updated objects for iterated IPC regression
        IPC_design_matrix <- cbind(1, covariates)
        fit_up <- list()
        IPC_reg_list_up <- IPC_reg_list
        param_up <- c()
        moment_deviations_up <- moment_deviations
        updated_IPC <- static_IPC
        nr_iterations <- 0
        difference <- rep(x = 1, times = ncol(estimates_matrix))

        # Start the iteration process
        while(nr_iterations < max_it & isFALSE(all(abs(difference) < convergence_criterion))){

          # Update the IPCs for every individual in each group
          for (i in unique(group[, 2])) {
            ID_group <- group[which(group[, 2] == i), 1]
            n_group <- length(ID_group)
            IPC_pred <- IPC_design_matrix[which(group[, 2] == i)[1],]
            for (j in 1:q) {
              param_up[j] <- sum(coef(IPC_reg_list_up[[j]]) * IPC_pred)
            }
            fit_up[[i]] <- omxSetParameters(model = fit, labels = param_names, values = param_up)
            fit_up[[i]] <- suppressMessages(mxRun(model = fit_up[[i]], useOptimizer = FALSE))
            suppressMessages(jacobian_up <- mxCheckIdentification(model = fit_up[[i]], details = FALSE)$jacobian)
            if (mean_structure == FALSE) {
              jacobian_up <- jacobian_up[1:p_star, ]
            }
            expected_cov_up <- mxGetExpected(model = fit_up[[i]], component = "covariance")
            expected_cov_inv_up <- solve(expected_cov_up)
            V_up <- 0.5 * t(D) %*% kronecker(X = expected_cov_inv_up, Y = expected_cov_inv_up) %*% D
            if (mean_structure) {
              V_m_cov_up <- matrix(data = 0, nrow = p_star_means, ncol = p_star_means)
              V_m_cov_up[1:p_star, 1:p_star] <- V_up
              V_m_cov_up[(p_star + 1):p_star_means, (p_star + 1):p_star_means] <- expected_cov_inv_up
              V_up <- V_m_cov_up
            }
            W_up <- solve(t(jacobian_up) %*% V_up %*% jacobian_up) %*% t(jacobian_up) %*% V_up

            # the individual deviations between observed and estimated moments
            moment_deviations_up <- data.frame(moment_deviations_up)
            moment_deviations_up[ID_group, 1:p_star] <- mom_cont_igsmc[ID_group, 1:p_star] -
              matrix(data = rep(x = matrixcalc::vech(expected_cov_up), times = n_group), byrow = TRUE,
                     nrow = n_group, ncol = p_star)
            if (mean_structure) {
              expected_means_up <- mxGetExpected(model = fit_up[[i]], component = "means")
              means_matrix_up <- matrix(data = rep(x = expected_means_up, times = N),
                                        byrow = TRUE, nrow = N, ncol = p)
              mean_dev_up <- data_obs - means_matrix_up
              moment_deviations_up[ID_group, (p_star + 1):p_star_means] <- mean_dev_up[ID_group, ]
            }
            moment_deviations_up <- as.matrix(moment_deviations_up)

            # Calculate the updated IPCs
            updated_IPC[ID_group, ] <- moment_deviations_up[ID_group, ] %*% t(W_up) +
              matrix(data = rep(x = param_up, times = n_group), byrow = TRUE, nrow = n_group, ncol = q)
          }

          # Estimate the iterated IPC regression parameters
          IPC_reg_data_up <- data.frame(cbind(updated_IPC, covariates))

          for (i in 1:length(target_parameters)) {
            IPC_reg_eq[i] <- paste(target_parameters[i], "~", IV)
            IPC_reg_list_up[[i]] <- do.call(what = "lm", args = list(formula = as.formula(IPC_reg_eq[i]),
                                                                     data = as.name("IPC_reg_data_up")))
          }
          # Store results
          estimates_up <- c(sapply(X = IPC_reg_list_up, FUN = function(x){coef(x)}))
          SE_up <- c(sapply(X = IPC_reg_list_up, FUN = function(x){sqrt(diag(vcov(x)))}))
          estimates_matrix <- rbind(estimates_matrix, estimates_up)
          SE_matrix <- rbind(SE_matrix, SE_up)

          # Covergence criteria
          difference <- estimates_matrix[nrow(estimates_matrix), ] -
            estimates_matrix[nrow(estimates_matrix) - 1, ]
          nr_iterations <- nr_iterations + 1
          cat("Iteration:", nr_iterations, "\n")
        }

        # Store the results
        IPC[[3]] <- updated_IPC
        names(IPC)[3] <- "IteratedIPCs"
        IPC[[4]] <- IPC_reg_list_up
        names(IPC)[4] <- "IteratedIPCRegression"

      }

    } else {
      warning("Incomplete covariates data. IPC regression parameter cannot be estimated.",
              call. = FALSE)
    }

    class(IPC) <- "ipcr"
    return(IPC)

  } else {
    warning("Incomplete SEM data. IPCs cannot be calculated",
            .Call = FALSE)
  }
}
