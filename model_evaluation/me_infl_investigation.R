##########
# title: investigation of Alaa & Van Der Schaar 2019 influence function corrected PEHE metric
##########
# Population-scale check of the proxy scores on this study's own DGM
# (scenarios 1-4, continuous_corr_0.5 - me_dgms.R) with ORACLE nuisances:
# mu0/mu1 are the true p0/p1 and pi is the true 0.5, so any gap between a
# score and true_pehe is the score's own bias, not nuisance error. N_BIG is
# large enough that Monte Carlo noise is negligible; bW is calibrated at
# CALIB_N, a study n, so the true ATE is the one the study sees.
#
# The scores are me_metrics.R's calc_infl_score()/calc_dr_risk() and the
# pseudo-outcome is me_nuisance.R's calculate_pseudos() - the same code that
# scores the candidates.

library(here)
source(here("model_evaluation", "me_dgms.R"))
source(here("model_evaluation", "me_nuisance.R")) # calculate_pseudos()
source(here("model_evaluation", "me_metrics.R")) # calc_infl_score(), calc_dr_risk(), study

N_BIG <- 2e6
CALIB_N <- 500 # the grid also has 250 and 1000
PI_TRUE <- 0.5

# one-step corrected PEHE against an arbitrary (possibly biased) surrogate
onestep <- function(th, sur, phi) {
  mean((sur - th)^2 + 2 * (sur - th) * (phi - sur))
}

# corrected version of the infl score:
calc_infl_score_corrected <- function(tau_hat, tau, pi, Y, W) {
  A <- W - pi
  C <- pi * (1 - pi)
  B <- 2 * W * A / C
  diff_tau <- tau - tau_hat
  
  pehe_infl <- (1 - B) * tau^2 + B * Y * diff_tau + 2 * A * diff_tau^2 - 2 * tau_hat^2
  
  mean(pehe_infl)
}

investigate_scenario <- function(scenario) {
  gen <- generate_me_scenario_data(scenario, N_BIG, seed = 1, calib_n = CALIB_N)
  Y <- gen$dataset$Y
  W <- gen$dataset$W
  X2 <- gen$dataset$X2
  tau <- gen$truth$tau

  oracle <- data.frame(
    mu0_DR = gen$truth$p0,
    mu1_DR = gen$truth$p1,
    pi = PI_TRUE
  )
  phi <- calculate_pseudos(oracle, Y, W)$phi

  cands <- list(
    truth = tau,
    ate = rep(mean(tau), N_BIG),
    half = mean(tau) + 0.5 * (tau - mean(tau)),
    mu1 = gen$truth$p1,
    tau_plus_X2 = tau + X2
  )

  message(sprintf(
    "scenario %d: bW = %.3f, sd(tau) = %.3f",
    scenario,
    gen$bW,
    sd(tau)
  ))

  bind_rows(lapply(names(cands), function(nm) {
    th <- cands[[nm]]
    tibble::tibble(
      scenario = scenario,
      candidate = nm,
      true_pehe = mean((th - tau)^2),
      dr = calc_dr_risk(th, phi),
      infl_tau = calc_infl_score(th, tau, PI_TRUE, Y, W),
      infl_phi = calc_infl_score(th, phi, PI_TRUE, Y, W),
      infl_v2_tau = calc_infl_score_corrected(th, tau, PI_TRUE, Y, W),
      infl_v2_phi = calc_infl_score_corrected(th, tau, PI_TRUE, Y, W),
      onestep_biased_sur = onestep(th, tau + 0.3 * X2, phi)
    )
  }))
}

res <- bind_rows(lapply(
  sort(unique(study$grid$scenario)),
  investigate_scenario
))

# results
print(mutate(res, across(where(is.double), ~ round(.x, 3))), n = Inf)

# table sorted by true PEHE
print(arrange(res, scenario, true_pehe), n = Inf)
