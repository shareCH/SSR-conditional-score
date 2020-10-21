# conditional power calculation

cond_power <- function(t1, # interim test statistic
                       n_1, # interim sample size per group
                       n_2, # initial incremental sample size for stage 2 per group
                       n_2new, # recalculated sample size for stage 2 per group
                       alpha_loc # local significance level
                       ){

w_1 <- sqrt(n_1)
w_2 <- sqrt(n_2)

cp <- 1 - pnorm(qnorm(1 - alpha_loc) * sqrt(w_1^2 + w_2^2) / w_2 
                - t1 * (w_1 / w_2 + sqrt(n_2new) / sqrt(n_1)))
}