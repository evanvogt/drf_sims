###############
# title: simulating data for scenario 11 - larger ATE
# date started: 07/01/2025
# date finished:
# author: Ellie Van Vogt
###############
rm(list = ls(al = T))
set.seed(1998)
# libraries ----

# paths ----
path <- "/rds/general/user/evanvogt/projects/nihr_drf_simulations"
setwd(path)

# parameters -----
sims <- 1000

sizes <- c(250, 500, 1000, 5000) # sample sizes



X1_prob <- 0.4 # probability of being female

b0 <- -0.4 # baseline log odds risk
bW <- -1 # average treatment effect
b1 <- 0.5 # prognostic - female
b2 <- 0.5 # prognostic - APACHE ish
b4 <- 0.6 # predictive


# function for generating the data


generate_dataset <- function(n) {
  W <- rbinom(n, 1, 0.5)
  X1 <- rbinom(n, 1, X1_prob)
  X2 <- rnorm(n, 0, 1)
  X4 <- rnorm(n, 0, 1)
  
  lp <- b0 + b1*X1 + b2*X2 + W*(bW + b4*X4)
  prob <- plogis(lp)
  Y <- rbinom(n, 1, prob)
  
  p0 <- plogis(b0 + b1*X1 + b2*X2)
  p1 <- plogis(b0 + b1*X1 + b2*X2 + (bW + b4*X4))
  tau <- p1 - p0
  
  dataset <- as.data.frame(cbind(Y, W, X1, X2, X4))
  truth <- as.data.frame(cbind(p0, p1, tau))
  
  return(list(dataset = dataset, truth = truth))
}

# generating the data ----

for (size in sizes) {
  dataset <- lapply(1:sims, function(i) generate_dataset(size))
  saveRDS(dataset, file = paste0("live/data/scenario_11_", size, ".RDS"))
}

# save the true DGM function for the oracle DR learner
fmla <- "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b4*X$X4)"
oracle_list <- list(fmla = fmla, b0 = b0, b1 = b1, b2 = b2, b4 = b4, bW = bW)
saveRDS(oracle_list, file = paste0("live/data/scenario_11_oracle.RDS"))

# true subgroup effects ----
# generate the threshold for positive and negative treatment effect:

thr <- -bW/b4 #2.5

large <- generate_dataset(100000)

large <- cbind(large[[1]], large[[2]])

s1 <- mean(large$tau[large$X4 > thr])

s2 <- mean(large$tau[large$X4 < thr])

gates <- c(s1, s2)
names(gates) <- c("X4>1.67", "X4<1.67")
saveRDS(gates, paste0("live/data/scenario_11_true_GATEs", size, ".RDS"))
