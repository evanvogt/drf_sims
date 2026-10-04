##########
# title: single-event survival data - competing_risk's event 1 without event 2
##########
# The event time is event 1 of competing_risk/surv_dgm.R with the competing
# event removed: Weibull, shape 2, log-scale log 25 - 0.1 X1 + 0.1 X2 +
# W (bW + b3 X3). Every constant is read from survival_scenario_params, so the
# two studies cannot drift apart. The estimand is that study's net RMST1
# (tau_RMST1_cs).
#
# PAIRING WITH competing_risk. DRAW ORDER (part of the contract - runs are
# reproduced by index):
#     W, Z-block (n x 8), cats, U, [C]
# which is generate_surv_data()'s order without its `cause` draw. Under the
# same setup_rng_stream(run), W, the covariates and U are therefore those of
# the competing-risk run with the same run index at rho = 0, and T here solves
# Lambda1(T) = -log U where there it solves Lambda1(T) + Lambda2(T) = -log U, so
# T_single >= T_CR unit by unit. C is not shared: the parent draws `cause`
# before C. se_dgm_check.R asserts the pairing.

library(dplyr)
# survival_scenario_params, and through R/dgm_scenarios.R correlated_covariates()
source(here::here("competing_risk", "surv_dgm.R"))

# Three scenarios, built from the competing-risk parameters:
#   1 null           no treatment effect
#   2 constant       competing_risk scenario 1's event-1 effect: HR 1.42
#   3 heterogeneous  competing_risk scenario 3's event-1 effect: HR 1.42 / 2.86
#                    at X3 = 0 / 1
# "Constant" is on the hazard scale. On the RMST scale the CATE still varies a
# little with X1 and X2 (SD about 0.10 against a mean of -2.16; see ADEMP.md).
# Scenario 1 being the null is what makes cate_metrics()'s hard-coded "scenario
# 1 has no heterogeneity" convention correct here.
se_scenario_params <- with(survival_scenario_params, data.frame(
  scenario = 1:3,
  description = c(
    "No treatment effect",
    "Constant treatment effect (hazard scale)",
    "Heterogeneous treatment effect"
  ),
  X1_prob = X1_prob[1],
  X3_prob = X3_prob[1],
  admin_time = admin_time[1],
  event_horizon = event_horizon[1],
  shape = shape1[1],
  scale_base = scale1_base[1],
  b1 = b1_1[1],
  b2 = b2_1[1],
  # log-scale coefficients, as in the parent (log-HR / shape)
  bW = c(0, bW_1[1], bW_1[3]),
  b3 = c(0, 0, b3_1[3]),
  stringsAsFactors = FALSE
))

#' Restricted mean of a Weibull(shape, scale) to the horizon, in closed form
#'
#' integral_0^h exp(-(t / scale)^shape) dt
#'   = scale * Gamma(1 + 1/shape) * P(1/shape, (h / scale)^shape),
#' P the regularised lower incomplete gamma. For shape 2 this is
#' scale * sqrt(pi) / 2 * erf(h / scale).
rmst_weibull <- function(scale, shape, horizon) {
  scale * gamma(1 + 1 / shape) * pgamma((horizon / scale)^shape, 1 / shape)
}

#' Log-scale of the event time under each arm
se_log_scale <- function(params, X1, X2, X3, W) {
  log(params$scale_base) + params$b1 * X1 + params$b2 * X2 +
    W * (params$bW + params$b3 * X3)
}

#' True potential-outcome RMSTs and the CATE at given covariates
#'
#' @param params one row of se_scenario_params
#' @return tibble, one row per unit: RMST_0, RMST_1, tau_RMST
se_truth <- function(params, X1, X2, X3) {
  scale0 <- exp(se_log_scale(params, X1, X2, X3, 0))
  scale1 <- exp(se_log_scale(params, X1, X2, X3, 1))
  h <- params$event_horizon
  tibble(
    RMST_0 = rmst_weibull(scale0, params$shape, h),
    RMST_1 = rmst_weibull(scale1, params$shape, h)
  ) %>%
    mutate(tau_RMST = RMST_1 - RMST_0)
}

#' Generate single-event survival data
#'
#' @param scenario Integer 1-3
#' @param n Sample size
#' @param return_truth Logical, whether to calculate the true CATEs
#' @param censoring Logical, whether to add uniform uninformative censoring on
#'   top of the administrative censoring at 180 (the study end)
#' @return list(dataset, truth)
generate_se_data <- function(scenario, n, return_truth = TRUE,
                             censoring = FALSE) {
  if (!scenario %in% se_scenario_params$scenario) {
    stop("Scenario must be between 1 and 3")
  }
  params <- se_scenario_params[se_scenario_params$scenario == scenario, ]

  # treatment and covariates, in the parent's draw order (see the header)
  W <- rbinom(n, 1, 0.5)
  cv <- correlated_covariates(n, list(X1_prob = params$X1_prob,
                                      X3_prob = params$X3_prob, rho = 0,
                                      s2 = 1, s4 = 1, s5 = 1))
  X1 <- cv$X1
  X2 <- cv$X2
  X3 <- cv$X3
  # noise: X01-X03 from the copula, X04/X05 indicators of a 3-level factor
  X01 <- cv$X01
  X02 <- cv$X02
  X03 <- cv$X03
  cats <- sample(c("A", "B", "C"), size = n, replace = TRUE, prob = c(0.45, 0.3, 0.25))
  X04 <- as.integer(cats == "A")
  X05 <- as.integer(cats == "B")

  # event time by inversion: S(t) = exp(-(t / scale)^shape) = U
  scale <- exp(se_log_scale(params, X1, X2, X3, W))
  u <- runif(n)
  Y <- scale * (-log(u))^(1 / params$shape)
  D <- rep(1L, n)

  # administrative censoring (end of study)
  admin <- params$admin_time
  D <- ifelse(Y > admin, 0L, D)
  Y <- pmin(Y, admin)

  # uniform uninformative censoring (if required)
  if (censoring) {
    censor_time <- runif(n, 1, admin)
    D <- ifelse(Y > censor_time, 0L, D)
    Y <- pmin(Y, censor_time)
  }

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
  if (return_truth) {
    result$truth <- se_truth(params, X1, X2, X3)
  }
  result
}
