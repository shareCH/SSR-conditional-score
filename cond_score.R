# conditional performance measures and score according to Herrmann et al. (2020)

#--------------------------------------
#This defines the path to the program
#--------------------------------------
setwd("PleaseDefinePathToProgram")
#--------------------------------------

source("Design.r")

cond_score <- function(cp, #conditional power
                       ref_cp, #target conditional power
                       n, #total sample size (stage 1 and stage 2)
                       n_req, # target sample size
                       n_sim #number of simulations
                       ){
  #------------------------------------
  # Boxplots for n and cp
  #------------------------------------
  mean_n <- rowMeans(n)
  var_n <- apply(n, 1, var)
  rownames(n) <- c("OCP", "restrOCP", "Prom", "OptFunc", "GS")
  par(mar=c(6.1, 5.1, 4.1, 2.1), cex.axis = 1.2, cex.lab = 2, mgp = c(3.2, 1.2, 0))
  boxplot(n, use.cols = FALSE, col = "lightgray", ylab = "recalculated n per group")
  abline(h = n_req, col = "red", lty = 1)
  points(mean_n, pch = 4, cex = 2)
  
  mean_cp <- rowMeans(cp)
  var_cp <- apply(cp, 1, var)
  rownames(cp) <- c("OCP", "restrOCP", "Prom", "OptFunc", "GS")
  par(mar = c(6.1, 5.1, 4.1, 2.1), cex.axis = 1.2, cex.lab = 2, mgp = c(3.2, 1.2, 0))
  boxplot(cp, use.cols = FALSE, col = "lightgray", ylab = "second stage cp")
  abline(h = ref_cp, col = "red", lty = 1)
  points(mean_cp, pch = 4, cex = 2)
  
  #------------------------------------
  # Scores for n and cp
  #------------------------------------
  e_n <- c()
  v_n <- c()
  score_n <- c()
  e_cp <- c()
  v_cp <- c()
  score_cp <- c()
  var_n_max <- c()
  var_cp_max <- c()
  
  for(i in 1:nrow(n)){
    if(is.na(mean_n[i])){
      e_n[i] <- NA
      v_n[i] <- NA
      score_n[i] <- NA
      e_cp[i] <- NA
      v_cp[i] <- NA
      score_cp[i] <- NA
    }
    else if(mean_n[i] > design$n_1*design$f & var_n[i]==0){
      e_n[i] <- NA
      v_n[i] <- NA
      score_n[i] <- NA
      e_cp[i] <- NA
      v_cp[i] <- NA
      score_cp[i] <- NA
    }
    else{
    # location component n:
    e_n[i] <- 1 - abs(mean_n[i] - n_req) / ((design$f - 1) * design$n_1)
    # variation component n:
    var_n_max[i] <- ((design$n_1 * design$f - design$n_1) / 2)^2 * n_sim / (n_sim - 1)
    v_n[i] <- 1 - sqrt(var_n[i] / var_n_max[i])
    # total n subscore:
    score_n[i] <- apply(cbind(e_n[i], v_n[i]), 1, mean)
    # location component cp:
    e_cp[i] <- 1 - abs(mean_cp[i] - ref_cp) / (1 - design$alpha_glob)
    # variation component cp:
    var_cp_max[i] <- ((1 - 0) / 2)^2 * n_sim / (n_sim - 1)
    v_cp[i] <- 1 - sqrt(var_cp[i] / var_cp_max[i])
    # total cp subscore:
    score_cp[i] <- apply(cbind(e_cp[i], v_cp[i]), 1, mean)
    }
  }
  
  #------------------------------------
  # Common conditional score 
  #------------------------------------
  score_cond <- apply(cbind(score_cp, score_n), 1, mean)
  
  res <- list(
    mean_cp = mean_cp, # mean conditional pwer
    e_cp = e_cp, # location component conditional pwer
    var_cp = var_cp, # variance component conditional power
    v_cp = v_cp, # variation component conditional power
    score_cp = score_cp, # sub-score conditional power
    mean_n = mean_n, # mean conditional sample size
    e_n = e_n, # location component conditional sample size
    var_n = var_n, # variance conditional sample size
    v_n = v_n, # variation component conditional sample size
    score_n = score_n, # sub-score conditional sample size
    score_cond = score_cond # conditional score
  )
  return(res)  
}
