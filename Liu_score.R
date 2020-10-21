# global sample size evaluation according to Liu et al. (2008)

Liu_score <- function(delta, n, cp, cp_true, design, n_sim){
  pow <- 0.8
  alpha <- design$alpha_glob
  
  f_s <- 2   #value chosen as proposed in Liu et al. (2008)
  f_p <- 0.2   #value chosen as proposed in Liu et al. (2008)
  
  set.seed(140)
  s1 <- rnorm(n_sim, delta * sqrt(design$n_1/2), 1)
  
  #-------------------------------------------------
  # Ideal sample size (paper: "N_{1-\beta}(\delta)")
  #-------------------------------------------------
  n_ideal_pre <- 2*(qnorm(pow) + qnorm(1 - alpha))^2 / delta^2
  n_ideal <- ifelse(is.infinite(n_ideal_pre), NA, n_ideal_pre)
  
  #----------------------------------------------------------------------------
  # create matrix of all n (include also the ones with efficacy/futility stops) 
  #----------------------------------------------------------------------------
  n_all <- ifelse(s1 < qnorm(1 - design$alpha_1) & 
                    s1 >= qnorm(1 - design$alpha_0), 0, design$n_1)
  n_all <- rbind(n_all, n_all, n_all, n_all, n_all)
  c <- 1
  for (i in 1:length(n_all[1,])){
    if (n_all[1,i] == 0){
      n_all[,i] <- n[,c]
      c <- c + 1
    }
    else 
      n_all[,i] = n_all[,i]
  }
  rownames(n_all) <- c("OCP", "restr OCP", "Promising", "OptFunc", "classicGS")
  
  
  #---------------------------------------------------------------------------
  # Ratio of study sample size to ideal sample size 
  # (paper: "SR(N \vert \delta, \beta)")
  #---------------------------------------------------------------------------
  sr <- matrix(data = NA, nrow = length(n_all[,1]), ncol = length(n_all[1,]))
  for (i in 1:length(n_all[,1])){
      sr[i,] <- n_all[i,] / n_ideal
  }
  
  #-----------------------------------------------------------------
  # ROS (relative oversize, paper: "ROS(\delta \vert f_s, \beta)")
  #-----------------------------------------------------------------
  ros <- ifelse(apply(sr - 1, 1, mean) > 0, apply(sr - 1, 1, mean) / (f_s - 1), 0)
  

  #------------------------------------
  # True power (paper: "pow")
  #------------------------------------
  pow_true <- c()

  for(i in 1:length(n_all[,1])){
    pow_true[i] <- (1 - pnorm(qnorm(1 - design$alpha_1), delta * sqrt(design$n_1/2), 1)) + (pnorm(qnorm(1-design$alpha_1), delta*sqrt(design$n_1/2), 1)-pnorm(qnorm(1-design$alpha_0), delta*sqrt(design$n_1/2), 1)) * mean(cp_true[i,])
  }
     
  #-------------------------------------------------
  # Fixed sample size at true power (paper: "N_pow")
  #-------------------------------------------------
  n_pow_true <- 2*(qnorm(pow_true) + qnorm(1 - alpha))^2 / delta^2
  
  #--------------------------------------------------------------------
  # Fixed sample size to achieve a scaled targeted power 
  # (paper: "N_{(1-f_p) \times (1-\beta)})
  #--------------------------------------------------------------------
  n_fpower <- 2*(qnorm(pow * (1 - f_p)) + qnorm(1 - alpha))^2 / delta^2
  
  #----------------------------------------------------------------
  # RUP (relative underpower, paper: "RUP(\delta \vert f_p, \beta))
  #----------------------------------------------------------------
  rup <- ifelse( (n_ideal - n_pow_true) > 0, (n_ideal - n_pow_true) / (n_ideal - n_fpower) , 0)
  
  #--------------------------------------------------
  # create matrix of all cp 
  # (include also the ones for the first stage 
  # which are defined as 0 and 1 respectively)
  #--------------------------------------------------
  cp_all <- ifelse(s1 < qnorm(1-design$alpha_1) & s1 >= qnorm(1-design$alpha_0), 100, ifelse(s1 >= qnorm(1-design$alpha_1), 1, 0))
  cp_all <- rbind(cp_all, cp_all, cp_all, cp_all, cp_all)
  cc <- 1
  for (i in 1:length(cp_all[1,])){
    if (cp_all[1,i] == 100){
      cp_all[,i] <- cp[,cc]
      cc <- cc + 1
    }
    else
      cp_all[,i] = cp_all[,i]
  }
  rownames(cp_all) <- c("OCP", "restr OCP", "Promising", "OptFunc", "classicGS")
   
  #-------------------------------------------------------------
  # total performance 
  # (paper: "R(\delta \vert f_p, f_s, \beta)")
  # (as average performance score)
  #-------------------------------------------------------------
  score <- rup + ros
  
  #------------------------------------
  # additional output
  #------------------------------------
  mean_n_all <- rowMeans(n_all)
  mean_cp_all <- rowMeans(cp_all)
  
  res <- list(
    ros = ros, # relative oversizing value (cf. Liu et al. (2008))
    rup = rup, # relative underpowering value (cf. Liu et al. (2008))
    score = score, # global score (cf. Liu et al. (2008))
    n_all = n_all, # global sample size values (for known interim test statistic)
    mean_n_all = mean_n_all, # expected global sample size
    cp_all = cp_all, # global power values (for known interim test statistic)
    pow_true = pow_true # global power for true effect size
  )
  
  return(res)  
}
