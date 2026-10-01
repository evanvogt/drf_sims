##########
# title: generating competing risk data using joint hazards (Beyersmann 2009)
##########
library(future.apply)
library(dplyr)
# correlated_covariates(): the Gaussian copula the missing-data and correlated
# sample-size sets draw from (see CORRELATED COVARIATES in that file)
source(here::here("R", "dgm_scenarios.R"))
# generating data from weibull distributions

# create dataframe of survival parameters
survival_scenario_params <- data.frame(
  scenario = 1:7,
  description = c(
    "ATE on E1 only",
    "ATE on E2 only",
    "HTE on E1, no ATE on E2",
    "HTE on E1, ATE on E2",
    "HTE on E2, no ATE on E1",
    "HTE on E2, ATE on E1",
    "HTE on both events"
  ),
  # Base parameters
  X1_prob = 0.4,
  X3_prob = 0.7,
  
  # Admin censoring (end of study - 6 months)
  admin_time = 180,
  
  # Event horizon (for truth calc)
  event_horizon = 28,
  
  # Baseline Weibull parameters for event 1
  # (25 rather than 15 since 2026-10-01: at 15 nearly everyone had an event by
  # the horizon, and the competition between the events swamped the X3 effects
  # on the RMTL scale - see ADEMP.md)
  shape1 = 2,
  scale1_base = 25,
  
  # Baseline Weibull parameters for event 2
  shape2 = 1.1,
  scale2_base = 45,
  
  # keep all effects on the log-scale
  # Prognostic effects
  b1_1 = -0.1,
  b2_1 = 0.1,
  b1_2 = -0.1,
  b2_2 = 0.1,
  
  # Treatment and predictive effects, written as log-HR / shape: a log-scale
  # coefficient b multiplies the hazard by exp(-shape * b), so dividing by the
  # shape gives both events the same log-HRs (+-0.35 ATE, +-0.7 more if X3).
  # E1 HR 1.42, 2.86 if X3; E2 HR 0.70, 0.35 if X3.
  # Treatment effects
  bW_1 = c(-0.35, 0, -0.35, -0.35, 0, -0.35, -0.35) / 2,
  bW_2 = c(0, 0.35, 0, 0.35, 0.35, 0.35, 0.35) / 1.1,

  # predictive effects
  b3_1 = c(0, 0, -0.7, -0.7, 0, 0, -0.7) / 2,
  b3_2 = c(0, 0, 0, 0, 0.7, 0.7, 0.7) / 1.1,
  
  stringsAsFactors = FALSE
)

# helper functions
cum_haz_all <- function(t, shape1, shape2, scale1, scale2) {
  (t/scale1)^shape1 + (t/scale2)^shape2
}

find_time <- function(u, shape1, shape2, scale1, scale2) {
  f <- function(t) cum_haz_all(t, shape1, shape2, scale1, scale2) + log(u)
  uniroot(f, c(0, 200))$root
}

hazard <- function(t, shape, scale) {
  shape * (1/scale)^shape * t^(shape - 1)
}

truth_individual <- function(shape1, shape2, scale1_1, scale1_0, scale2_1, scale2_0, horizon) {
  # hazards under each treatment
  h1_1 <- function(t) hazard(t, shape1, scale1_1)
  h1_0 <- function(t) hazard(t, shape1, scale1_0)
  
  h2_1 <- function(t) hazard(t, shape2, scale2_1)
  h2_0 <- function(t) hazard(t, shape2, scale2_0)
  
  # all-cause cumulative hazards
  cum_haz_all_1 <- function(t) cum_haz_all(t, shape1, shape2, scale1_1, scale2_1)
  cum_haz_all_0 <- function(t) cum_haz_all(t, shape1, shape2, scale1_0, scale2_0)
  
  # all cause survival functions
  S_all_1 <- function(t) exp(-cum_haz_all_1(t))
  S_all_0 <- function(t) exp(-cum_haz_all_0(t))
  
  # terms to integrate for getting CIFs
  f1_1 <- function(t) S_all_1(t) * h1_1(t)
  f1_0 <- function(t) S_all_0(t) * h1_0(t)
  
  f2_1 <- function(t) S_all_1(t) * h2_1(t)
  f2_0 <- function(t) S_all_0(t) * h2_0(t)
  
  # Cumulative incidence functions for each cause
  CIF1_1 <- Vectorize(function(t) integrate(f1_1, 0, t)$value)
  CIF1_0 <- Vectorize(function(t) integrate(f1_0, 0, t)$value)
  
  CIF2_1 <- Vectorize(function(t) integrate(f2_1, 0, t)$value)
  CIF2_0 <- Vectorize(function(t) integrate(f2_0, 0, t)$value)
  
  # cause-specific hazard estimands (actually sub-distribution but everyone calls it cause-specific UGH)
  # RMTL - expected time from event i to the event horizon
  RMTL1_1 <- integrate(function(u) CIF1_1(u), 0, horizon)$value
  RMTL1_0 <- integrate(function(u) CIF1_0(u), 0, horizon)$value
  RMTL2_1 <- integrate(function(u) CIF2_1(u), 0, horizon)$value
  RMTL2_0 <- integrate(function(u) CIF2_0(u), 0, horizon)$value
  
  # RMSTc - expected time free from any event (essentially a composite outcome)
  RMSTc_1 <- integrate(function(u) 1 - CIF1_1(u) - CIF2_1(u), 0, horizon)$value
  RMSTc_0 <- integrate(function(u) 1 - CIF1_0(u) - CIF2_0(u), 0, horizon)$value
  
  # RMSTi - expected time to event i, irrespective of competing events 
  RMST1_1 <- horizon - RMTL1_1
  RMST1_0 <- horizon - RMTL1_0
  RMST2_1 <- horizon - RMTL2_1
  RMST2_0 <- horizon - RMTL2_0
  
  # Cause-specific (net) RMST — ignores competing events, uses cause-specific survival only
  # Target estimand for IPW and csf_cs approaches
  S1_cs_1 <- function(t) exp(-(t/scale1_1)^shape1)
  S1_cs_0 <- function(t) exp(-(t/scale1_0)^shape1)
  S2_cs_1 <- function(t) exp(-(t/scale2_1)^shape2)
  S2_cs_0 <- function(t) exp(-(t/scale2_0)^shape2)

  RMST1_cs_1 <- integrate(function(u) S1_cs_1(u), 0, horizon)$value
  RMST1_cs_0 <- integrate(function(u) S1_cs_0(u), 0, horizon)$value
  RMST2_cs_1 <- integrate(function(u) S2_cs_1(u), 0, horizon)$value
  RMST2_cs_0 <- integrate(function(u) S2_cs_0(u), 0, horizon)$value

  list(RMTL1_1 = RMTL1_1, RMTL1_0 = RMTL1_0, RMTL2_1 = RMTL2_1, RMTL2_0 = RMTL2_0,
       RMSTc_1 = RMSTc_1, RMSTc_0 = RMSTc_0, RMST1_1 = RMST1_1, RMST1_0 = RMST1_0,
       RMST2_1 = RMST2_1, RMST2_0 = RMST2_0,
       RMST1_cs_1 = RMST1_cs_1, RMST1_cs_0 = RMST1_cs_0,
       RMST2_cs_1 = RMST2_cs_1, RMST2_cs_0 = RMST2_cs_0)
}

#' Generate Competing Risks Survival Data
#'
#' @param scenario Integer 1-7 specifying the data generation scenario
#' Covariates come from the copula (correlated_covariates(), R/dgm_scenarios.R)
#' at every rho, 0 included, so a run at rho = 0 and at rho = 0.5 is paired.
#' The copula draws X1-X5 and X01-X03; this study keeps X1, X2, X3 and X01-X03
#' (exchangeable correlation rho among those six) and drops X4/X5. X04/X05, the
#' categorical pair, are independent of everything. Marginals are unchanged:
#' X1 ~ Bernoulli(0.4), X2 ~ N(0, 1), X3 ~ Bernoulli(0.7), X01-X03 ~ N(0, 1).
#'
#' DRAW ORDER (part of the contract - runs are reproduced by index):
#'     W, Z-block (n x 8), cats, U, cause, [C]
#' cats comes before U so that both censoring settings of a run share every
#' covariate; only C is extra under censoring = TRUE, and it is drawn last.
#'
#' @param scenario Integer 1-7 specifying the data generation scenario
#' @param n Sample size
#' @param rho Exchangeable latent correlation of the copula (0 = independent)
#' @param return_truth Logical, whether to calculate true treatment effects
#' @param censoring Logical, whether to add uniform uninformative censoring
#' @return List containing dataset and optionally truth
generate_surv_data <- function(scenario, n, rho = 0, return_truth = TRUE,
                               censoring = FALSE) {
  # scenario check
  if (!scenario %in% 1:7) {
    stop("Scenario must be between 1 and 7")
  }

  # Get parameters for this scenario
  params <- survival_scenario_params[survival_scenario_params$scenario == scenario, ]

  # Generate treatment and covariates
  W <- rbinom(n, 1, 0.5)
  cv <- correlated_covariates(n, list(X1_prob = params$X1_prob,
                                      X3_prob = params$X3_prob, rho = rho,
                                      s2 = 1, s4 = 1, s5 = 1))
  X1 <- cv$X1
  X2 <- cv$X2
  X3 <- cv$X3
  # noise: X01-X03 correlated with X1-X3 (proxies, never in the hazards),
  # X04/X05 indicators of an independent 3-level factor
  X01 <- cv$X01
  X02 <- cv$X02
  X03 <- cv$X03
  cats <- sample(c("A", "B", "C"), size = n, replace = TRUE, prob = c(0.45, 0.3, 0.25))
  X04 <- as.integer(cats == "A")
  X05 <- as.integer(cats == "B")

  # Shape values (fixed)
  shape1 <- params$shape1
  shape2 <- params$shape2
  
  # Covariate & treatment effects on scale parameters
  log_scale1 <- log(params$scale1_base) + params$b1_1 * X1 + params$b2_1 * X2 + W * (params$bW_1 + params$b3_1 * X3)
  
  log_scale2 <- log(params$scale2_base) + params$b1_2 * X1 + params$b2_2 * X2 + W * (params$bW_2 + params$b3_2 * X3)
  
  scale1 <- exp(log_scale1)
  scale2 <- exp(log_scale2)
  
  # Solve for survival time (joint hazard)
  u <- runif(n)
  Y <- future_mapply(find_time, u, shape1, shape2, scale1, scale2)
  
  # Generate cause
  h1 <- hazard(Y, shape1, scale1)
  h2 <- hazard(Y, shape2, scale2)
  
  prob <- h1/(h1 + h2)
  
  D <- 2 - rbinom(n, 1, prob)
  
  # Administrative censoring (end of study)
  admin <- params$admin_time
  D <- ifelse(Y > admin, 0, D)
  Y <- ifelse(Y > admin, admin, Y)
  
  # Uniform uniformative censoring (if required)
  if (censoring) {
    censor_time <- runif(n, 1, admin)
    D <- ifelse(Y > censor_time, 0, D)
    Y <- pmin(Y, censor_time)
  }

  # Dataset for analysis
  dataset <- data.frame(
    Y = Y,
    D = as.integer(D),
    W = W,
    X1 = X1,
    X2 = X2,
    X3 = X3,
    X01 = X01,
    X02 = X02,
    X03 = X03,
    X04 = X04,
    X05 = X05
  )
  
  result <- list(dataset = dataset)
  
  # truth calculation for RMST
  if (return_truth) {
    result$truth <- surv_truth(params, X1, X2, X3)
  }
  return(result)
}

#' True potential-outcome restricted means and CATEs at given covariates
#'
#' Per unit, by numerical integration (truth_individual()), so it is exact
#' whatever the covariates' joint distribution. generate_surv_data() calls it
#' at the sampled units; surv_dgm_check.R at a covariate grid.
#'
#' @param params one row of survival_scenario_params
#' @param X1,X2,X3 covariate vectors of equal length
#' @return tibble, one row per unit: the RMTL / RMST levels under each arm and
#'   the tau_* contrasts
surv_truth <- function(params, X1, X2, X3) {
  # scales under each treatment
  log_scale1_0 <- log(params$scale1_base) + params$b1_1 * X1 + params$b2_1 * X2
  log_scale1_1 <- log_scale1_0 + (params$bW_1 + params$b3_1 * X3)

  log_scale2_0 <- log(params$scale2_base) + params$b1_2 * X1 + params$b2_2 * X2
  log_scale2_1 <- log_scale2_0 + (params$bW_2 + params$b3_2 * X3)

  scale1_0 <- exp(log_scale1_0)
  scale1_1 <- exp(log_scale1_1)

  scale2_0 <- exp(log_scale2_0)
  scale2_1 <- exp(log_scale2_1)

  horizon <- params$event_horizon
  truth <- future_mapply(truth_individual, params$shape1, params$shape2, scale1_1, scale1_0, scale2_1, scale2_0, horizon, SIMPLIFY = FALSE)
  truth <- bind_rows(truth)

  truth %>%
    mutate(
      tau_RMTL1 = RMTL1_1 - RMTL1_0,
      tau_RMTL2 = RMTL2_1 - RMTL2_0,
      tau_RMSTc = RMSTc_1 - RMSTc_0,
      tau_RMST1 = RMST1_1 - RMST1_0,
      tau_RMST2 = RMST2_1 - RMST2_0,
      tau_RMST1_cs = RMST1_cs_1 - RMST1_cs_0,
      tau_RMST2_cs = RMST2_cs_1 - RMST2_cs_0
      )
}
