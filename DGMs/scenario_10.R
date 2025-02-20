###############
# title: simulating data for scenario 10 - exponential hte
# date started: 07/01/2025
# date finished:
# author: Ellie Van Vogt
###############
# doing this with the logistic link because otherwise the errors go weird
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

b0 <- -0.4 # baseline log odds risk
bW <- -0.2 # average treatment effect
b1 <- 0.5 # prognostic - female
b2 <- 0.5 # prognostic - APACHE ish
b3 <- 0.2
b4 <- -0.1


# function for generating the data


generate_dataset <- function(n) {
  W <- rbinom(n, 1, 0.5)
  X1 <- rbinom(n, 1, X1_prob)
  X2 <- rnorm(n, 0, 1)
  X4 <- rnorm(n, 0, 1)
  
  lp <- b0 + b1*X1 + b2*X2 + W*(bW + b4*exp(X4))
  prob <- plogis(lp)
  Y <- rbinom(n, 1, prob)
  
  p0 <- plogis(b0 + b1*X1 + b2*X2)
  p1 <- plogis(b0 + b1*X1 + b2*X2 + (bW + b4*exp(X4)))
  tau <- p1 - p0
  
  dataset <- as.data.frame(cbind(Y, W, X1, X2, X4))
  truth <- as.data.frame(cbind(p0, p1, tau))
  
  return(list(dataset = dataset, truth = truth))
}

# generating the data ----

for (size in sizes) {
  dataset <- lapply(1:sims, function(i) generate_dataset(size))
  saveRDS(dataset, file = paste0("live/data/scenario_10_", size, ".rds"))
}


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
# saveRDS(gates, paste0("live/data/scenario_5_true_GATEs", size, ".rds"))
