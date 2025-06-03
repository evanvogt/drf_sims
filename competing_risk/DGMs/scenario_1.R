###############
# title: simulating data for cts outcome and no HTE
# date started: 07/01/2025
# date finished:
# author: Ellie Van Vogt
###############
rm(list = ls(all = T))
set.seed(1998)
# libraries ----
library(survival)
library(survRM2)
library(dplyr)

# paths ----
path <- "/rds/general/user/evanvogt/projects/nihr_drf_simulations"
setwd(path)

# parameters -----
sims <- 1000

sizes <- c(250, 500, 1000, 5000) # sample sizes



X1_prob <- 0.4 # probability of being female

#discharge parameters
b0 <- 4 # baseline discharge time
bW <- -1 # treatment shortens ICU stay
b1 <- -0.05 # prognostic - female
b2 <- 0.7 # prognostic - APACHE ish

#mortality parameters
c0 <- 2 # baseline mortality time
cW <- 1 # treatment prolongs mortality
c1 <- 0.05 # prognostic - female
c2 <- -1



# function for generating the data


generate_dataset <- function(n) {
  W <- rbinom(n, 1, 0.5)
  X1 <- rbinom(n, 1, X1_prob)
  X2 <- rnorm(n, 0, 1)
  
  
  scale_discharge <- exp(b0 + b1*X1 + b2*X2 + W*bW)
  scale_mortality <- exp(c0 + c1*X1 + c2*X2 + W*cW)
  
  shape_discharge <- 2
  shape_mortality <- 0.75
  
  time_discharge <- rweibull(n, shape_discharge, scale_discharge)
  time_mortality <- rweibull(n, shape_mortality, scale_mortality)

  summary(time_discharge)
  summary(time_mortality)
  
  hist(time_discharge)
  hist(time_mortality)
  
  time_event <- pmin(time_discharge, time_mortality)
  
  event_type <- ifelse(time_discharge < time_mortality, 1, 2)
  
  
  censoring_time <- runif(n, 0, 90)  
  censored <- time_event > censoring_time
  time_event[censored] <- censoring_time[censored]
  event_type[censored] <- 0 
  
  tau <- rmst2(time_event, event_type == 1, W, tau = 90)
  
  tau <- p1 - p0
  
  dataset <- as.data.frame(cbind(Y, W, X1, X2))
  truth <- as.data.frame(cbind(p0, p1, tau))
  
  return(list(dataset = dataset, truth = truth))
}


# generating the data ----

for (size in sizes) {
  dataset <- lapply(1:sims, function(i) generate_dataset(size))
  saveRDS(dataset, file = paste0("live/data/competing_risk//scenario_1_", size, ".RDS"))
}

# save the true DGM function for the oracle DR learner
fmla <- "b0 + b1*X$X1 + b2*X$X2 + W*bW"
oracle_list <- list(fmla = fmla, b0 = b0, b1 = b1, b2 = b2, bW = bW)
saveRDS(oracle_list, file = paste0("live/data/competing_risk//scenario_1_oracle.RDS"))