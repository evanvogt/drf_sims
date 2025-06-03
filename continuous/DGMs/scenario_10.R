###############
# title: simulating data for scenario 10 - exponential hte
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

sizes <- as.numeric(commandArgs(trailingOnly = TRUE)) # sample sizes

X1_prob <- 0.4 # probability of being female
X3_prob <- 0.7 # mechvent prob

b0 <- 0.4 # creatinine increase
bW <- -0.1 # average treatment effect
b1 <- -0.05 # prognostic - female
b2 <- 2 # prognostic - APACHE ish
b3 <- 0.3
b4 <- 0.1

s2 <- 1 # var prognostic
s4 <- 1
s_err <- 0.5 # var error term
s <- s_err + s2 # total variation (from apriori knowledge only?)

# function for generating the data


generate_dataset <- function(n) {
  W <- rbinom(n, 1, 0.5)
  X1 <- rbinom(n, 1, X1_prob)
  X2 <- rnorm(n, 0, s2)
  X3 <- rbinom(n, 1, X3_prob)
  X4 <- rnorm(n, 0, s4)
  
  err <- rnorm(n, 0, s_err)
  
  Y <- b0 + b1*X1 + b2*X2 + W*(bW + b3*X3 + b4*exp(-abs(X4))) + err

  p0 <- b0 + b1*X1 + b2*X2
  p1 <- b0 + b1*X1 + b2*X2 + (bW + b4*exp(abs(X4)))
  tau <- p1 - p0
  
  # add a bunch of variables with no relation to outcome or treatment
  X01 <- rnorm(n, 0, 1)
  X02 <- rnorm(n, 0, 1)
  X03 <- rnorm(n, 0, 1)
  cats <- sample(c("A", "B", "C"), size = n, replace = TRUE, prob = c(0.45, 0.3, 0.25))
  X04 <- as.integer(cats == "A")
  X05 <- as.integer(cats == "B")
  
  dataset <- as.data.frame(cbind(Y, W, X1, X2, X4, X01, X02, X03, X04, X05))
  truth <- as.data.frame(cbind(p0, p1, tau))
  
  return(list(dataset = dataset, truth = truth))
}

# generating the data ----

for (size in sizes) {
  # make sure bW is right size for power
  diff <- power.t.test(n = size/2, delta = NULL, sd = s, power = 0.75)$delta
  bW <- round(-diff, digits = 2)
  
  dataset <- lapply(1:sims, function(i) generate_dataset(size))
  saveRDS(dataset, file = paste0("live/data/continuous/scenario_10_", size, ".RDS"))
  
  # save the true DGM function for the oracle DR learner
  fmla <- "b0 + b1*X1 + b2*X2 + W*(bW + b3*X3 + b4*exp(-abs(X4)))"
  fmla <- gsub("\\b(X\\d+)\\b", "X$\\1", fmla)
  oracle_list <- list(fmla = fmla, b0 = b0, b1 = b1, b2 = b2, b4 = b4, bW = bW)
  saveRDS(oracle_list, file = paste0("live/data/continuous/scenario_10_", size, "_oracle.RDS"))
}
