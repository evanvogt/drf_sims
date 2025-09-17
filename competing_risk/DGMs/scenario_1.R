###############
# title: simulating data for secnario 1 - no HTE
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
cat(sizes)


X1_prob <- 0.4 # probability of being female

#discharge parameters
b0 <- 4 # baseline discharge time
bW <- -1 # treatment shortens ICU stay
b1 <- -0.05 # prognostic - female
b2 <- 0.7 # prognostic - APACHE ish

#mortality parameters
c0 <- 8 # baseline mortality time
cW <- 1 # treatment prolongs mortality
c1 <- 0.05 # prognostic - female
c2 <- -1 # prognostic

generate_dataset <- function(n) {

  W <- rbinom(n, 1, 0.5)
  X1 <- rbinom(n, 1, X1_prob)
  X2 <- rnorm(n, 0, 1)

  lp <- b0 + b1*X1 + b2*X2 + W*bW
  prob <- plogis(lp)
  Y <- rbinom(n, 1, prob)
  
  p0 <- plogis(b0 + b1*X1 + b2*X2)
  p1 <- plogis(b0 + b1*X1 + b2*X2 + bW)
  tau <- p1 - p0
  
  # add a bunch of variables with no relation to outcome or treatment
  X01 <- rnorm(n, 0, 1)
  X02 <- rnorm(n, 0, 1)
  X03 <- rnorm(n, 0, 1)
  cats <- sample(c("A", "B", "C"), size = n, replace = TRUE, prob = c(0.45, 0.3, 0.25))
  X04 <- as.integer(cats == "A")
  X05 <- as.integer(cats == "B")
  
  dataset <- as.data.frame(cbind(Y, W, X1, X2, X01, X02, X03, X04, X05))
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
  cat(paste0("saving new dataset scenario 1 ", size))
  saveRDS(dataset, file = paste0("live/data/binary/scenario_1_", size, ".RDS"))
  
  # save the true DGM function for the oracle DR learner with the correctly powered bW
  fmla <- "b0 + b1*X$X1 + b2*X$X2 + W*bW"
  oracle_list <- list(fmla = fmla, b0 = b0, b1 = b1, b2 = b2, bW = bW)
  saveRDS(oracle_list, file = paste0("live/data/binary/scenario_1_", size, "_oracle.RDS"))
}