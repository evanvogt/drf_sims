###############
# title: simulating data for scenario 5 - two-way interaction with HTE
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

sizes <- as.numeric(commandArgs(trailingOnly = TRUE)) # sample sizes

X1_prob <- 0.4 # probability of being female
X3_prob <- 0.7 # mechvent prob

b0 <- -0.4 # baseline log odds risk
b1 <- 0.5 # prognostic - female
b2 <- 0.5 # prognostic - APACHE ish
b34 <- -0.5 # two way interaction

# function for generating the data

generate_dataset <- function(n) {
  W <- rbinom(n, 1, 0.5)
  X1 <- rbinom(n, 1, X1_prob)
  X2 <- rnorm(n, 0, 1)
  X3 <- rbinom(n, 1, X3_prob)
  X4 <- rnorm(n, 0, 1)
  
  lp <- b0 + b1*X1 + b2*X2 + W*(bW + b34*X3*X4)
  prob <- plogis(lp)
  Y <- rbinom(n, 1, prob)
  
  p0 <- plogis(b0 + b1*X1 + b2*X2)
  p1 <- plogis(b0 + b1*X1 + b2*X2 + (bW + b34*X3*X4))
  tau <- p1 - p0
  
  # add a bunch of variables with no relation to outcome or treatment
  X01 <- rnorm(n, 0, 1)
  X02 <- rnorm(n, 0, 1)
  X03 <- rnorm(n, 0, 1)
  cats <- sample(c("A", "B", "C"), size = n, replace = TRUE, prob = c(0.45, 0.3, 0.25))
  X04 <- as.integer(cats == "A")
  X05 <- as.integer(cats == "B")
  
  dataset <- as.data.frame(cbind(Y, W, X1, X2, X3, X4, X01, X02, X03, X04, X05))
  truth <- as.data.frame(cbind(p0, p1, tau))
  
  return(list(dataset = dataset, truth = truth))
}

# generating the data ----

for (size in sizes) {
  cat(size)
  # adequately powered bW
  p1 <- plogis(b0)
  p2 <- power.prop.test(size/2, p2 = p1, power = 0.75)$p1 # makes sure that p2 is less than p1
  bW <- round(qlogis(p2) - b0, digits = 2)
  
  dataset <- lapply(1:sims, function(i) generate_dataset(size))
  cat(paste0("saving new dataset scenario 5 ", size))
  saveRDS(dataset, file = paste0("live/data/binary/scenario_5_", size, ".RDS"))
  
  # save the true DGM function for the oracle DR learner
  fmla <- "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b34*X$X3*X$X4)"
  oracle_list <- list(fmla = fmla, b0 = b0, b1 = b1, b2 = b2, b34 = b34, bW = bW)
  saveRDS(oracle_list, file = paste0("live/data/binary/scenario_5_", size, "_oracle.RDS"))
}
