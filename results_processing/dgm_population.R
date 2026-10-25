##########
# title: the DGMs' true CATE over a trial population - appendix figures and tables
##########
# The DGMs' true CATE over a trial population, shared by thesis_figures/dgm.R
# and thesis_tables/dgm_tables.R so a figure and a table use the same draws.
# Reads no simulation output. rho = 0.5 only: continuous_corr_0.5 and
# binary_corr_0.5 scenarios 1-4, and competing risks' 7 scenarios (RMTL CATE,
# both events).
#
# Continuous / binary tau is reported as tau(x) - ATE, which does not depend on
# n (calibrate_bW() moves only the location); the ATE per n is in
# cate_summary(). Competing-risk tau is not calibrated, so it is reported as it
# is, in days. One copula draw serves all three outcomes; competing risks use
# its first DGM_N_POP_SURV rows.
#
# l_moments() and surv_moment_grid() give the SNR of the competing-risk check
# and rel_eff_snr.R (SD(tau) / SD(noise)).
#
# Source R/figures.R first (it loads dplyr and tidyr).

source(here::here("R", "dgm_scenarios.R"))
source(here::here("competing_risk", "surv_dgm.R"))

DGM_RHO <- 0.5
DGM_SCENARIOS <- 1:4
DGM_N_LEVELS <- c(100, 250, 500, 1000)
DGM_OUTCOMES <- c("continuous", "binary")
# copula draws for the population summaries. SD(tau) is exact (te_moments());
# the correlations and shares are Monte Carlo, to about 0.002 at this size
DGM_N_POP <- 5e5
# the competing-risk truth is interpolated per draw, 7 scenarios x 2 events
DGM_N_POP_SURV <- 2e5
DGM_SEED <- 2026

#' One scenario's parameters from a rho = DGM_RHO set
dgm_params <- function(outcome, scenario) {
  tbl <- resolve_set(corr_set(outcome, DGM_RHO))
  tbl[tbl$scenario == scenario, ]
}

#' The covariate population: one copula draw at DGM_RHO
#'
#' @return tibble of X1-X5, X01-X03, with a row id
dgm_covariates <- function(n_pop = DGM_N_POP, seed = DGM_SEED) {
  set.seed(seed)
  cv <- correlated_covariates(n_pop, dgm_params("continuous", 1))
  tibble(id = seq_len(n_pop), !!!cv)
}

#' Continuous / binary truth at each draw, scenarios DGM_SCENARIOS
#'
#' @param cv output of dgm_covariates()
#' @return long tibble: outcome, scenario, id, m0 (the control mean - the
#'   control risk for binary), tau_c = tau(x) - ATE
cate_population <- function(cv) {
  bind_rows(lapply(DGM_OUTCOMES, function(outcome) {
    binary <- outcome == "binary"
    bind_rows(lapply(DGM_SCENARIOS, function(s) {
      p <- dgm_params(outcome, s)
      g <- truth_at(p, 0, binary, cv$X1, cv$X2, cv$X3, cv$X4, cv$X5)$tau
      tibble(
        outcome = outcome,
        scenario = s,
        id = cv$id,
        m0 = control_mean(p, cv$X1, cv$X2, binary),
        tau_c = g - te_moments(p)$mean
      )
    }))
  }))
}

#' Exact tau(x) - ATE on a covariate grid, one scenario
#'
#' @param grid data frame with X3, X4, X5 (X1 and X2 do not enter tau)
#' @return the grid with tau_c alongside
cate_grid <- function(outcome, scenario, grid) {
  p <- dgm_params(outcome, scenario)
  g <- truth_at(p, 0, outcome == "binary", rep(0, nrow(grid)), rep(0, nrow(grid)),
                grid$X3, grid$X4, grid$X5)$tau
  mutate(grid, tau_c = g - te_moments(p)$mean)
}

#' Per outcome x scenario summaries of the continuous / binary truth
#'
#' @param pop output of cate_population()
#' @param cv output of dgm_covariates(), the same draw
#' @return tibble, one row per outcome x scenario:
#'   sd_tau            SD(tau), exact (te_moments())
#'   q01, q99          the 1st and 99th percentiles of tau - ATE over the
#'                     draws. Not te_range(): over the support (X4, X5 to
#'                     +-6 SD) scenario 3's X4 * X5 term reaches +-18
#'   cor_m0_tau        cor(m0, tau) over the draws
#'   r2_linear         R^2 of tau on X1-X5, X01-X03, linear: how much of the
#'                     heterogeneity a linear CATE model could capture
#'   ate_<n>           the ATE, bW + E[g], at each n in DGM_N_LEVELS
#'   opp_<n>           the share of the population whose tau has the opposite
#'                     sign to the ATE at that n
#'   cor_m0_tau and r2_linear are NA where tau is constant (scenario 1)
cate_summary <- function(pop, cv) {
  covs <- as.matrix(cv[, COPULA_VARS])
  bind_rows(lapply(DGM_OUTCOMES, function(outcome) {
    set <- corr_set(outcome, DGM_RHO)
    bind_rows(lapply(DGM_SCENARIOS, function(s) {
      p <- dgm_params(outcome, s)
      d <- pop[pop$outcome == outcome & pop$scenario == s, ]
      stopifnot(identical(d$id, cv$id))
      eg <- te_moments(p)$mean
      varies <- sd(d$tau_c) > 0
      ate <- sapply(DGM_N_LEVELS, function(n) calibrate_bW(p, n, calibration_for(set)) + eg)
      opp <- sapply(ate, function(a) mean((a + d$tau_c) * a < 0))
      tibble(
        outcome = outcome,
        scenario = s,
        sd_tau = sqrt(te_moments(p)$var),
        q01 = unname(quantile(d$tau_c, 0.01)),
        q99 = unname(quantile(d$tau_c, 0.99)),
        cor_m0_tau = if (varies) cor(d$m0, d$tau_c) else NA_real_,
        r2_linear = if (varies) summary(lm(d$tau_c ~ covs))$r.squared else NA_real_,
        !!!setNames(as.list(ate), paste0("ate_", DGM_N_LEVELS)),
        !!!setNames(as.list(opp), paste0("opp_", DGM_N_LEVELS))
      )
    }))
  }))
}

#' Competing-risk RMTL CATE of each event at each draw, all 7 scenarios
#'
#' surv_truth_grid() per scenario (numerical integration, in parallel), then
#' interpolated to the draws in X2 within each (X1, X3) cell
#' (surv_truth_at_draws()), as competing_risk/surv_dgm_check.R does.
#'
#' @param cv output of dgm_covariates(); its first DGM_N_POP_SURV rows are used
#' @param workers future workers for the grid truth
#' @return long tibble: scenario, event ("Event 1", "Event 2"), id, X3, tau
#'   (days) and tau_x3 = tau(X1, X2, X3 = 1) - tau(X1, X2, X3 = 0), the effect of
#'   switching X3 with X1 and X2 held
surv_population <- function(cv, workers = min(6, max(1, parallelly::availableCores() - 1))) {
  oplan <- future::plan(future::multisession, workers = workers)
  on.exit(future::plan(oplan), add = TRUE)

  cv <- cv[seq_len(min(DGM_N_POP_SURV, nrow(cv))), ]
  events <- c(tau_RMTL1 = "Event 1", tau_RMTL2 = "Event 2")
  one <- rep(1, nrow(cv))
  bind_rows(lapply(survival_scenario_params$scenario, function(s) {
    gt <- surv_truth_grid(survival_scenario_params[survival_scenario_params$scenario == s, ])
    bind_rows(lapply(names(events), function(col) {
      tibble(
        scenario = s,
        event = unname(events[col]),
        id = cv$id,
        X3 = cv$X3,
        tau = surv_truth_at_draws(gt, col, cv$X1, cv$X2, cv$X3),
        tau_x3 = surv_truth_at_draws(gt, col, cv$X1, cv$X2, one) -
          surv_truth_at_draws(gt, col, cv$X1, cv$X2, 0 * one)
      )
    }))
  }))
}

#' Per scenario x event summaries of the competing-risk truth
#'
#' @param spop output of surv_population()
#' @return tibble: mean and SD of tau, its mean within X3 = 0 and X3 = 1, and
#'   the standardised X3 effect E[tau_x3] (surv_dgm_check.R section 3)
surv_summary <- function(spop) {
  spop %>%
    group_by(scenario, event) %>%
    summarise(
      mean_tau = mean(tau),
      sd_tau = sd(tau),
      mean_x3_0 = mean(tau[X3 == 0]),
      mean_x3_1 = mean(tau[X3 == 1]),
      x3_effect = mean(tau_x3),
      .groups = "drop"
    )
}

#' The first two moments of L_k = (h - T)+ 1(D = k) given (x, w)
#'
#' E[L_k] = int_0^h (h - t) f_k(t) dt and E[L_k^2] = int_0^h (h - t)^2 f_k(t) dt,
#' f_k = S_all h_k the subdensity of event k (surv_dgm.R's hazards). E[L_k] is
#' the RMTL of event k, which truth_individual() gets by the nested integral of
#' the CIF.
#'
#' @param params one row of survival_scenario_params
#' @param scale1,scale2 the two events' Weibull scales at (x, w)
#' @return c(m1_1, var_1, m1_2, var_2), the mean and variance for each event
l_moments <- function(params, scale1, scale2) {
  h <- params$event_horizon
  shapes <- c(params$shape1, params$shape2)
  scales <- c(scale1, scale2)
  out <- numeric(0)
  for (k in 1:2) {
    f <- function(t) {
      exp(-cum_haz_all(t, params$shape1, params$shape2, scale1, scale2)) *
        hazard(t, shapes[k], scales[k])
    }
    m1 <- integrate(function(t) (h - t) * f(t), 0, h)$value
    m2 <- integrate(function(t) (h - t)^2 * f(t), 0, h)$value
    out[paste0(c("m1_", "var_"), k)] <- c(m1, m2 - m1^2)
  }
  out
}

#' l_moments() in both arms on surv_truth_grid()'s X2 x X1 x X3 grid
#'
#' @return the grid, with per event k: tau_RMTLk (the CATE) and noise_k, the
#'   variance of L_k averaged over the two arms
surv_moment_grid <- function(params, X2_grid = seq(-4.5, 4.5, length.out = 41)) {
  grid <- expand.grid(X2 = X2_grid, X1 = 0:1, X3 = 0:1)
  # the scales as surv_truth() builds them
  log_s1_0 <- log(params$scale1_base) + params$b1_1 * grid$X1 + params$b2_1 * grid$X2
  log_s2_0 <- log(params$scale2_base) + params$b1_2 * grid$X1 + params$b2_2 * grid$X2
  log_s1_1 <- log_s1_0 + params$bW_1 + params$b3_1 * grid$X3
  log_s2_1 <- log_s2_0 + params$bW_2 + params$b3_2 * grid$X3
  arm <- function(ls1, ls2) {
    t(mapply(function(a, b) l_moments(params, exp(a), exp(b)), ls1, ls2))
  }
  m0 <- arm(log_s1_0, log_s2_0)
  m1 <- arm(log_s1_1, log_s2_1)

  # E[L_k] is the RMTL: check it against truth_individual() at a few cells
  for (i in c(1, 62, 123, 164)) {
    ti <- truth_individual(params$shape1, params$shape2,
                           exp(log_s1_1[i]), exp(log_s1_0[i]),
                           exp(log_s2_1[i]), exp(log_s2_0[i]),
                           params$event_horizon)
    stopifnot(
      abs(m1[i, "m1_1"] - ti$RMTL1_1) < 1e-4, abs(m0[i, "m1_1"] - ti$RMTL1_0) < 1e-4,
      abs(m1[i, "m1_2"] - ti$RMTL2_1) < 1e-4, abs(m0[i, "m1_2"] - ti$RMTL2_0) < 1e-4
    )
  }

  mutate(
    grid,
    tau_RMTL1 = m1[, "m1_1"] - m0[, "m1_1"],
    tau_RMTL2 = m1[, "m1_2"] - m0[, "m1_2"],
    noise_1 = 0.5 * (m0[, "var_1"] + m1[, "var_1"]),
    noise_2 = 0.5 * (m0[, "var_2"] + m1[, "var_2"])
  )
}
