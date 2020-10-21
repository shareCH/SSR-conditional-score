# sample size calculation according to optimization function approach 
# by Jennison and Turnbull (2015)

#-----------------------------------------------#
#This defines the path to the program
#-----------------------------------------------#
setwd("PleaseDefinePathToProgram")
#-----------------------------------------------#

source("cond_power.r")

optfunc_design <- function(design, t1, delta){
  #---------------------------------------
  # Recalculate second stage sample size
  #---------------------------------------
  n_recalc <- c()
  tau <- c()
  for(j in 1:length(t1)){
    tau[j] <- 0.005
    f <- rep(0, design$n_1+design$n_2 - 1)
    cp <- rep(0, design$n_1+design$n_2 - 1)
    n <- rep(0, design$n_1+design$n_2 - 1)

    for(i in (design$n_1 + design$n_2):(design$f * design$n_1)){
      n[i] <- i
      cp[i] <- cond_power(t1[j], design$n_1, design$n_2, n[i], design$alpha_1)
      f[i] <- cp[i] - tau[j] / 4 * ((n[i]) - (design$n_1 + design$n_2))
    }
    n_recalc[j] <- which(f==max(f))
  }
  n <- ifelse(t1 >= qnorm(1 - design$alpha_1) | t1 < qnorm(1 - design$alpha_0),
              design$n_1, n_recalc)
  
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
      cp[i] <- 1 - pnorm(qnorm(1 - design$alpha_1) * sqrt(w_1^2 + w_2^2)
                         /w_2 - t1[i] * (w_1 / w_2 + sqrt((n[i] - design$n_1) / design$n_1)))
      cp_true[i] <- 1 - pnorm(qnorm(1 - design$alpha_1) * sqrt(w_1^2 + w_2^2) 
                              /w_2 - t1[i] * (w_1 / w_2) - delta * sqrt((n[i] - design$n_1) / 2))
    }
  }
  
  res <- list(
    n = n,
    cp = cp,
    cp_true = cp_true
  )
  return(res) 
}