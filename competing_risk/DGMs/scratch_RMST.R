# figuring out the ATEs to specify for each given sample size
library(cmprsk)
library(survRM2)
library(mets)
library(dplyr)
# function to generate datasets with competing risks
# event horizon == censoring time here
simulate_cr_data <- function(n, shape1, scale1, ate1, shape2, scale2, ate2, studyend) {
  W <- rbinom(n, 1, 0.5)
  
  log_scale1 <- log(scale1) + ate1*W
  log_scale2 <- log(scale2) + ate2*W
  
  time1 <- rweibull(n, shape1, exp(log_scale1))
  time2 <- rweibull(n, shape2, exp(log_scale2))
  
  
  status <- ifelse(time1 < time2, 1, 2)
  time <- pmin(time1, time2)
  
  # administrative censoring + uniform censoring
  censtime <- runif(n, 0, 90)
  status <- ifelse(time < studyend & time < censtime, status, 0)
  
  df <- data.frame(cbind(W, time, status))
  return(df)
}

# analysis of RMST - subdistribution
analyse_sub_dist <- function(df) {
  # keeping those who have CE in the risk set
  df <- df %>% mutate(
    #time = ifelse(status == 2, studyend + 1, time), # people who have CE kept until after end of study
    event = factor(status, 0:2, c("censored", "EOI", "CE"))
  )

  fgfit <- finegray(Surv(time, event) ~ ., data = df, etype = "EOI")
  
  crfit <- crr(df$time, df$status, df$W)
  crr_times <- predict(crfit, df$W)
  
  cufit <- cuminc(df$time, df$status, df$W)
}

analyse_cause_spec <- function(df) {
  # keeping those who have CE in the risk set
  df <- df %>% mutate(
    event = factor(status, 0:2, c("censored", "EOI", "CE"))
  )
  
  ajfit <- survfit(Surv(time, event) ~ 1, df)
  res <- print(sub_model)
}

sim_trials <- lapply(seq(1, 1000), function(i) {
  df <- simulate_cr_data(500, 1.5, 4, -1, 0.5, 8, 1, 90)
  
  sd_res <- analyse_sub_dist(df)
  
  cs_res <- analyse_cause_spec(df)
  list(sub_dist = sd_res, cause_spec = cs_res)
})


df <- sim_trials[[1]]
