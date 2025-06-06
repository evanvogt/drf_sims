###############
# title: simulating data for scenario 6 - two-way interaction with HTE + bigger ATE
# date started: 07/01/2025
# date finished:
# author: Ellie Van Vogt
###############
rm(list = ls(all=T))
set.seed(1998)
# libraries ----

# paths ----
path <- "/rds/general/user/evanvogt/projects/nihr_drf_simulations"
setwd(path)

# parameters -----
sims <- 1000

sizes <- c(250, 500, 1000, 5000) # sample sizes


X1_prob <- 0.4 # probability of being female
X3_prob <- 0.7 # mechvent prob

b0 <- 1 # creatinine increase
bW <- -1 # average treatment effect
b1 <- -0.05 # prognostic - female
b2 <- 2 # prognostic - APACHE ish
b3 <- 2
b4 <- 0.5
b34 <- -0.5 # two way interaction

# function for generating the data

generate_dataset <- function(n) {
  W <- rbinom(n, 1, 0.5)
  X1 <- rbinom(n, 1, X1_prob)
  X2 <- rnorm(n, 0, 1)
  X3 <- rbinom(n, 1, X3_prob)
  X4 <- rnorm(n, 0, 1)
  
  err <- rnorm(0, 0.5)
  
  Y <- b0 + b1*X1 + b2*X2 + W*(bW + b3*X3 + b4*X4 + b34*X3*X4) + err

  p0 <- b0 + b1*X1 + b2*X2
  p1 <- b0 + b1*X1 + b2*X2 + (bW + b3*X3 + b4*X4 + b34*X3*X4)
  tau <- p1 - p0
  
  dataset <- as.data.frame(cbind(Y, W, X1, X2, X3, X4))
  truth <- as.data.frame(cbind(p0, p1, tau))
  
  return(list(dataset = dataset, truth = truth))
}

# generating the data ----

for (size in sizes) {
  dataset <- lapply(1:sims, function(i) generate_dataset(size))
  saveRDS(dataset, file = paste0("live/data/competing_risk//scenario_6_", size, ".RDS"))
}

# save the true DGM function for the oracle DR learner
fmla <- "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b3*X$X3 + b4*X$X4 + b34*X$X3*X$X4)"
oracle_list <- list(fmla = fmla, b0 = b0, b1 = b1, b2 = b2, b3 = b3, b4 = b4, b34 = b34, bW = bW)
saveRDS(oracle_list, file = paste0("live/data/competing_risk//scenario_6_oracle.RDS"))

# true subgroup effects ----
# not sure about how to calculate the truth for this yet lol - might need to fiddle more with the DGM

# 
# # create a heat map of the variables to find the cut offs??
# thr <- -bW/b4 #1
# 
# 
# large <- generate_dataset(100000)
# 
# large <- cbind(large[[1]], large[[2]])
# 
# s1 <- mean(large$tau[large$X4 > thr])
# 
# s2 <- mean(large$tau[large$X4 < thr])
# 
# gates <- c(s1, s2)
# names(gates) <- c("X4>1", "X4<1")
# saveRDS(gates, paste0("live/data/competing_risk//scenario_5_true_GATEs", size, ".rds"))
