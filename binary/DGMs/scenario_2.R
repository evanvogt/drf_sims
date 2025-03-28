###############
# title: simulating data for secnario 2 - simple HTE
# date started: 07/01/2025
# date finished:
# author: Ellie Van Vogt
###############
rm(list = ls(all = T))
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
b3 <- -0.4 # predictive



# function for generating the data


generate_dataset <- function(n) {
  W <- rbinom(n, 1, 0.5)
  X1 <- rbinom(n, 1, X1_prob)
  X2 <- rnorm(n, 0, 1)
  X3 <- rbinom(n, 1, X3_prob)
  
  lp <- b0 + b1*X1 + b2*X2 + W*(bW + b3*X3)
  prob <- plogis(lp)
  Y <- rbinom(n, 1, prob)
  
  p0 <- plogis(b0 + b1*X1 + b2*X2)
  p1 <- plogis(b0 + b1*X1 + b2*X2 + (bW + b3*X3))
  tau <- p1 - p0
  
  dataset <- as.data.frame(cbind(Y, W, X1, X2, X3))
  truth <- as.data.frame(cbind(p0, p1, tau))
  
  return(list(dataset = dataset, truth = truth))
}


# generating the data ----

for (size in sizes) {
  dataset <- lapply(1:sims, function(i) generate_dataset(size))
  saveRDS(dataset, file = paste0("live/data/scenario_2_", size, ".RDS"))
}

# save the true DGM function for the oracle DR learner
fmla <- "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b3*X$X3)"
oracle_list <- list(fmla = fmla, b0 = b0, b1 = b1, b2 = b2, b3 = b3, bW = bW)
saveRDS(oracle_list, file = paste0("live/data/scenario_2_oracle.RDS"))

# true subgroup effects ----
# generate the threshold for positive and negative treatment effect:


large <- generate_dataset(100000)

large <- cbind(large[[1]], large[[2]])

s1 <- mean(large$tau[large$X3 == 1])

s2 <- mean(large$tau[large$X3 == 0])

gates <- c(s1, s2)
names(gates) <- c("X3==1", "X3==0")
saveRDS(gates, paste0("live/data/scenario_2_true_GATEs", size, ".RDS"))

