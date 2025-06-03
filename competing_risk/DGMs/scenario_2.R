###############
# title: simulating data for secnario 2 - binary HTE
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

b0 <- 0.2 # creatinine increase
bW <- -0.1 # average treatment effect
b1 <- -0.05 # prognostic - female
b2 <- 2 # prognostic - APACHE ish
b3 <- 2 

# function for generating the data

generate_dataset <- function(n) {
  W <- rbinom(n, 1, 0.5)
  X1 <- rbinom(n, 1, X1_prob)
  X2 <- rnorm(n, 0, 1)
  X3 <- rbinom(n, 1, X3_prob)
  
  err <- rnorm(n, 0, 0.5)
  
  Y <- b0 + b1*X1 + b2*X2 + W*(bW + b3*X3) + err

  p0 <- b0 + b1*X1 + b2*X2
  p1 <- b0 + b1*X1 + b2*X2 + (bW + b3*X3)
  tau <- p1 - p0
  
  # add a bunch of variables with no relation to outcome or treatment
  X01 <- rnorm(n, 0, 1)
  X02 <- rnorm(n, 0, 1)
  X03 <- rnorm(n, 0, 1)
  cats <- sample(c("A", "B", "C"), size = n, replace = TRUE, prob = c(0.45, 0.3, 0.25))
  X04 <- as.integer(cats == "A")
  X05 <- as.integer(cats == "B")
  
  dataset <- as.data.frame(cbind(Y, W, X1, X2, X3, X01, X02, X03, X04, X05))
  truth <- as.data.frame(cbind(p0, p1, tau))
  
  return(list(dataset = dataset, truth = truth))
}


# generating the data ----

for (size in sizes) {
  dataset <- lapply(1:sims, function(i) generate_dataset(size))
  saveRDS(dataset, file = paste0("live/data/competing_risk//scenario_2_", size, ".RDS"))
}

# save the true DGM function for the oracle DR learner
fmla <- "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b3*X$X3)"
oracle_list <- list(fmla = fmla, b0 = b0, b1 = b1, b2 = b2, b3 = b3, bW = bW)
saveRDS(oracle_list, file = paste0("live/data/competing_risk//scenario_2_oracle.RDS"))

# true subgroup effects ----
# generate the threshold for positive and negative treatment effect:


large <- generate_dataset(100000)

large <- cbind(large[[1]], large[[2]])

s1 <- mean(large$tau[large$X3 == 1])

s2 <- mean(large$tau[large$X3 == 0])

gates <- c(s1, s2)
names(gates) <- c("X3==1", "X3==0")
saveRDS(gates, paste0("live/data/competing_risk//scenario_2_true_GATEs", size, ".RDS"))

