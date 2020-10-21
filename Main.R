####################################################################
# R Code for the article "A new conditional performance score      #
# for the evaluation of adaptive group sequential designs with     #
# sample size recalculation" by Carolin Herrmann, Maximilian Pilz, # 
# Meinhard Kieser and Geraldine Rauch (doi: 10.1002/sim.8534)      #
####################################################################

################################################################################
# To determine the performance please:
# - specify your specific path to the program for Main.R, cond_score.R, 
#   optfunc_design.R, promising_design.R
# - decide on the number of simulations n_sim in Main.R
# - determine the design (interim sample size per group, global one-sided 
#   significance level (one sided), initial incremental sample size for stage 2 
#   per group, futility stop bound, factor for maximal sample size) in Main.R
# - define the true underlying standardized effect size in Main.R
################################################################################


library(lattice)
library(RColorBrewer)
library(latticeExtra)

#-----------------------------------------------#
# Define path to the program:
#-----------------------------------------------#
setwd("PleaseDefinePathToProgram")
#-----------------------------------------------#

source("Design.r")
source("OCP_design.r")
source("restrOCP_design.r")
source("promising_design.r")
source("optfunc_design.r")
source("classicGS_design.r")
source("cond_power.r") 
source("cond_score.r")
source("Liu_score.r")


#---------------------------------------------#
# Specify the number of simulations n_sim
#---------------------------------------------#
n_sim <- 10000
#---------------------------------------------#

#-----------------------------------------------------------------------------#
# Determine the study design with
# interim sample size per group n_1, 
# global significance level (one sided) alpha_glob, 
# initial incremental sample size for stage 2 per group n_2, 
# futility stop bound alpha_0, 
# factor for maximal sample size f with n_max=f*n_1
#-----------------------------------------------------------------------------#
design <- design(n_1 = 50, alpha_glob = 0.025, n_2 = 50, alpha_0 = 0.5, f = 4)
#-----------------------------------------------------------------------------#

#----------------------------------------------------------------------------------
# Plot the sample size as a function of t1 and delta
# (here: observed conditional power approach "OCP", restricted observed conditional
# power approach "restrOCP", promising zone approach "Prom", optimization function 
# approach "OptFunc", classical group sequential approach "GS", fixed sample size
# approach "fixed")
#----------------------------------------------------------------------------------
set.seed(5)
t1 <- seq(-1,3,0.01)
delta <- t1*sqrt(2/design$n_1)

par(mar = c(5.1, 4.1, 4.1, 2.1), cex.axis = 1, cex.lab = 1, mgp = c(3,1,0))
plot(OCP_design(design, delta/sqrt(2/design$n_1), delta)$n ~ delta, 
     xlim = c(-1*sqrt(2/design$n_1), 3*sqrt(2/design$n_1)), ylim = c(0,300), 
     type = 'l', col = "red", lwd = 2, lty = 1, xlab = NA, ylab = NA, axes = F)
axis(side = 3)
mtext(side = 3, line = 3, 'Observed interim effect')
mtext(side = 1, line = 3, 'Observed value of the interim test statistic')
par(new = TRUE)
plot(restrOCP_design(design, t1, delta)$n ~ t1, xlim = c(-1,3), ylim = c(0,300),
     type = 'l', col = "orange", lwd = 2, lty = 2, xlab = NA, ylab = "Total n per group")
par(new = TRUE)
plot(promising_design(design, t1, delta)$n ~ t1, xlim = c(-1,3), ylim = c(0,300),
     type = 'l', col = "blue", lwd = 2, lty = 3, xlab = NA, ylab = "Total n per group")
par(new = TRUE)
plot(optfunc_design(design, t1, delta)$n ~ t1, xlim = c(-1,3), ylim = c(0,300),
     type = 'l', col = "green", lwd = 2, lty = 4, xlab = NA, ylab = "Total n per group")
par(new = TRUE)
plot(classicGS_design(design, t1, delta)$n ~ t1, xlim = c(-1,3), ylim = c(0,300),
     type = 'l', col = "purple", lwd = 2, lty = 5, xlab = NA, ylab = "Total n per group")
par(new = TRUE)
abline(v = qnorm(1-design$alpha_0), col="black", lty = 1)
abline(v = qnorm(1-design$alpha_1), col="black", lty = 1)
par(new = TRUE)
legend(-1, 300, c("OCP", "restrOCP", "Prom", "OptFunc", "GS"),
       lty = c(1, 2, 3, 4, 5), lwd = c(2, 2, 2, 2, 2), bty = "n", 
       col = c("red", "orange", "blue", "green", "purple"))



###################################################################
# Initiate Performance Simulation
###################################################################

#------------------------------------------------------------
# Decide on the underlying treatment effect Delta (e.g. 0.4)
#------------------------------------------------------------
delta <- 0.4
#------------------------------------------------------------

#------------------------------------------------------------
# Create values of test statistic
#------------------------------------------------------------
set.seed(140)
s1 <- rnorm(n_sim, delta*sqrt(design$n_1/2), 1)


#--------------------------------------------------------------------
# Determine required fixed sample size
# and define the target values for sample size and conditional power
#--------------------------------------------------------------------
if(delta != 0){
  n_req <- 1 * ((qnorm(0.975) + qnorm(0.8)) / delta)^2
  pow <- power.t.test(n = n_req, delta = delta, sd = 1, sig.level = 0.025,
                    power = NULL,
                    type = c("two.sample"),
                    alternative = "one.sided")$power
  while(pow < 0.8){
    n_req <- n_req + 1
    pow <- power.t.test(n = n_req, delta = delta, sd = 1, sig.level = 0.025,
                       power = NULL,
                       type = c("two.sample"),
                       alternative = "one.sided")$power
  }
  n_req <- n_req
  ref_cp <- 0.8
}
if(delta == 0){
  n_req <- design$n_1
  ref_cp <- 0.025
}

if(n_req > design$f * design$n_1){
  n_req <- design$n_1
  ref_cp <- 0.025
}

#---------------------------------------------------------------------------
# Restrict the simulated data to the recalculation area
#---------------------------------------------------------------------------
t1 <- s1[s1 < qnorm(1 - design$alpha_1) & s1 >= qnorm(1 - design$alpha_0)]

#------------------------------------------------------------------------
# Store the recalculated total sample sizes "n" per design in one matrix
#------------------------------------------------------------------------
n <- matrix(nrow = 5, ncol = length(t1))
rownames(n) <- c("OCP", "restrOCP", "Prom", "OptFunc", "classicGS")
n[1,] <- OCP_design(design, t1, delta)$n
n[2,] <- restrOCP_design(design, t1, delta)$n
n[3,] <- promising_design(design, t1, delta)$n
n[4,] <- optfunc_design(design, t1, delta)$n
n[5,] <- classicGS_design(design, t1, delta)$n

#--------------------------------------------------------------------------------
# Store the conditional power "cp" for the observed delta at the interim analysis
#--------------------------------------------------------------------------------
par(mfrow = c(1,1))
cp <- matrix(nrow = 5, ncol = length(t1))
rownames(cp) <- c("OCP", "restrOCP", "Prom", "OptFunc", "classicGS")
cp[1,] <- OCP_design(design, t1, delta)$cp
cp[2,] <- restrOCP_design(design, t1, delta)$cp
cp[3,] <- promising_design(design, t1, delta)$cp
cp[4,] <- optfunc_design(design, t1, delta)$cp
cp[5,] <- classicGS_design(design, t1, delta)$cp

#-----------------------------------------------------------
# Store the conditional power "cp_true" for the true delta 
#-----------------------------------------------------------
par(mfrow = c(1,1))
cp_true <- matrix(nrow = 5, ncol = length(t1))
rownames(cp_true) <- c("OCP", "restrOCP", "Prom", "OptFunc", "classicGS")
cp_true[1,] <- OCP_design(design, t1, delta)$cp_true
cp_true[2,] <- restrOCP_design(design, t1, delta)$cp_true
cp_true[3,] <- promising_design(design, t1, delta)$cp_true
cp_true[4,] <- optfunc_design(design, t1, delta)$cp_true
cp_true[5,] <- classicGS_design(design, t1, delta)$cp_true

#----------------------------------------------------------------------------
# Calculate the new conditional performance score "nscore" by Herrmann et al.
#----------------------------------------------------------------------------
nscore <- cond_score(cp, ref_cp, n, n_req, n_sim)

#--------------------------------------------------------------------
# Calculate the old global performance score "oscore" by Liu et al.
#--------------------------------------------------------------------
oscore <- Liu_score(delta, n, cp, cp_true, design, n_sim)

#----------------------------------------------------------------------------
# output for conditional and global performance:
# - out_cond with 
#   expected sample size per group within the recalculation area (RA) "Mean n", 
#   variance of sample size per group within the RA "Var n", conditional sample 
#   size sub-score "score n", expected conditional power within the RA "Mean cp", 
#   variance of conditional power within the RA "Var cp", conditional power 
#   sub-score "score cp", point-wise new conditional score "cond score"
# - out_cond_components with 
#   location component for conditional sample size "comp_e_n", variation component 
#   for conditional sample size "comp_v_n", conditional sample size sub-score 
#   "score n", location component for conditional power "comp_e_cp", variation 
#   component for conditional power "comp_v_cp", conditional power sub-score 
#   "score cp", point-wise new conditional score "cond score" 
# - out_Liu with 
#   expected sample size per group over both stages "mean_n_all", relative 
#   oversizing function "ROS", power of the design for true standardized treatment 
#   effect "pow_true", relative underpowering function "RUP", point-wise Liu score 
#   "LiuScore"
#----------------------------------------------------------------------------

# conditional performance evaluation
out_cond <- matrix(nrow = 5, ncol = 7)
rownames(out_cond) <- c("OCP", "restrOCP", "Promising", "Opt.Func.", "classicGS")
colnames(out_cond) <- c("Mean n", "Var n", "score n","Mean cp", "Var cp", 
                        "score cp", "cond score")  
for(i in 1:5){
  out_cond[i,] <- round(c(nscore$mean_n[i], nscore$var_n[i], nscore$score_n[i], 
                          nscore$mean_cp[i], nscore$var_cp[i], nscore$score_cp[i], 
                          nscore$score_cond[i]), digits = 3)
}
print(out_cond)

# components of conditional performance score
out_cond_components <- matrix(nrow = 5, ncol = 7)
rownames(out_cond_components) <- c("OCP", "restrOCP", "Promising", "Opt.Func.", 
                                   "classicGS")
colnames(out_cond_components) <- c("comp_e_n", "comp_v_n", "score n", "comp_e_cp", 
                                   "comp_v_cp", "score cp", "cond score")
for(i in 1:5){
  out_cond_components[i,] <- round(c(nscore$e_n[i], nscore$v_n[i], nscore$score_n[i], 
                                     nscore$e_cp[i], nscore$v_cp[i], nscore$score_cp[i], 
                                     nscore$score_cond[i]), digits = 3)
}
print(out_cond_components)

# global performance evaluation
out_Liu <- matrix(nrow=5, ncol=5)
rownames(out_Liu) <- c("OCP", "restrOCP", "Promising", "Opt.Func.", "classicGS")
colnames(out_Liu) <- c("mean_n_all", "ROS", "pow_true", "RUP", "LiuScore")
for(i in 1:5){
  out_Liu[i,] <- round(c(oscore$mean_n_all[i], oscore$ros[i], oscore$pow_true[i], 
                         oscore$rup[i], oscore$score[i]), digits = 3)
}
print(out_Liu)
