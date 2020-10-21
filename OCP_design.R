# sample size calculation based on observed conditional power

OCP_design <- function(design, t1, delta){
  #---------------------------------------
  # Recalculate second stage sample size
  #---------------------------------------
  delta_obs <- t1 * sqrt(2 / design$n_1)   
  n_recalc <- design$n_1 * ((qnorm(0.2) - qnorm(1 - design$alpha_1) 
                             * sqrt(sqrt(design$n_1)^2 
                                    + sqrt(design$n_1)^2) / sqrt(design$n_2)) / t1 
                            + sqrt(design$n_1) / sqrt(design$n_2))^2 + design$n_1 
  n <- ifelse(t1 >= qnorm(1 - design$alpha_1)| t1 < qnorm(1 - design$alpha_0),
              design$n_1, ifelse(n_recalc < design$f * design$n_1, n_recalc, 
                                 design$f * design$n_1))
 
  #---------------------------------------
  # Calculate conditional power
  #---------------------------------------
  w_1 <- sqrt(design$n_1)
  w_2 <- sqrt(design$n_2)
  cp <- c()
  cp_true <- c()
  for(i in 1:length(n)){
    if(n[i]==design$n_1){
      cp[i] <- 0
      cp_true[i] <- 0
    }
    else{
      cp[i] <- 1-pnorm(qnorm(1 - design$alpha_1) * sqrt(w_1^2+w_2^2) / w_2 
                       - t1[i] * (w_1 / w_2 + sqrt((n[i] - design$n_1) / design$n_1)))
      cp_true[i] <- 1 - pnorm(qnorm(1 - design$alpha_1) * sqrt(w_1^2 + w_2^2)
                              /w_2 - t1[i] * (w_1 / w_2) 
                              - delta * sqrt((n[i] - design$n_1) / 2))
    }
  }
  
  res <- list(
    n = n,
    cp = cp,
    cp_true = cp_true
  )
  return(res) 
}