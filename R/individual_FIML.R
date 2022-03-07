individual_FIML <- function(data, exp_m, exp_cov){
  
  N <- NROW(data)
  ind_FIML <- rep(NA, length(N))
  is_missing <- is.na(data)
  
  for(i in seq_len(N)){
    
    n_missings <- sum(!is_missing[i, ])
    
    if(n_missings == 0) next
    
    # subset expected means
    filtered_m <- exp_m[!is_missing[i, ]]
    
    # subset expected covariance
    filtered_cov <- exp_cov[!is_missing[i, ], !is_missing[i, ]]
    
    ind_FIML[i] <- mvtnorm::dmvnorm(x = data[i, !is_missing[i, ]],
                                   mean = filtered_m,
                                   sigma = filtered_cov,
                                   log = TRUE)
    
  }
  
  ind_FIML
  
}
