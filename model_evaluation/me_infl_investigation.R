##########
# title: investigation of Alaa & Van Der Schaar 2019 influence function corrected PEHE metric
##########
#   Rscript me_infl_investigation.R
#
# The proxy scores on this study's own DGM (scenarios 1-4,
# continuous_corr_0.5 - me_dgms.R) with ORACLE nuisances: mu0/mu1 are the true
# p0/p1 and pi is the true 0.5, so any gap between a score and true_pehe is
# the score's own bias, not nuisance error. N_BIG is large enough that Monte
# Carlo noise is negligible; bW is calibrated at CALIB_N, a study n, so the
# true ATE is the one the study sees. With the oracle nuisances, infl_tau is
# the IF score with the T-learner surrogate (mu1 - mu0 = tau). population.RDS
# is also thesis_tables/me_tables.R's me_if_population table.
#
# Saves population.RDS to INVESTIGATION_DIR, results/correlated/
# model_evaluation/investigation/.
#
# The same comparison on the simulated runs needs no script of its own:
# me_metrics.R already scores every candidate with the IF score under both
# surrogates (inflt, infl) and DR risk, and me_infl_investigation.qmd reads
# them from me_strat_metrics.RDS.
#
# me_metrics.R is deliberately NOT sourced: when <prefix>_all.RDS exists (on
# the cluster) sourcing it recomputes and overwrites <prefix>_metrics.RDS.
# calc_infl_score()/calc_dr_risk() are copied from it below instead.

library(here)
library(dplyr)
source(here("model_evaluation", "me_config.R")) # study
source(here("model_evaluation", "me_dgms.R"))
source(here("model_evaluation", "me_nuisance.R")) # calculate_pseudos()

N_BIG <- 2e6
CALIB_N <- 500 # the grid also has 250 and 1000
PI_TRUE <- 0.5

INVESTIGATION_DIR <- file.path(study$res_path, "investigation")

save_investigation <- function(x, file) {
  dir.create(INVESTIGATION_DIR, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(INVESTIGATION_DIR, file)
  saveRDS(x, path)
  message("saved ", path)
}

# ---- scores -----------------------------------------------------------------

# copied verbatim from me_metrics.R - see the header for why it isn't sourced
calc_infl_score <- function(tau_hat, tau, pi, Y, W) {
  A <- W - pi
  C <- pi * (1 - pi)
  B <- 2 * W * A / C
  diff_tau <- tau - tau_hat

  pehe_infl <- (1 - B) * tau^2 + B * Y * diff_tau - A * diff_tau^2 + tau_hat^2

  mean(pehe_infl)
}

# copied verbatim from me_metrics.R
calc_dr_risk <- function(tau_hat, tau_dr) {
  mean((tau_hat - tau_dr)^2)
}

# one-step corrected PEHE against an arbitrary (possibly biased) surrogate
onestep <- function(th, sur, phi) {
  mean((sur - th)^2 + 2 * (sur - th) * (phi - sur))
}

# ---- population check ---------------------------------------------------------

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
      n_pop = N_BIG,
      true_pehe = mean((th - tau)^2),
      # what the IF score tracks at pi = 0.5, up to a constant (me_metrics.R's
      # me_per_model(), mu1_mse)
      mu1_mse = mean((th - gen$truth$p1)^2),
      dr = calc_dr_risk(th, phi),
      infl_tau = calc_infl_score(th, tau, PI_TRUE, Y, W),
      infl_phi = calc_infl_score(th, phi, PI_TRUE, Y, W),
      onestep_biased_sur = onestep(th, tau + 0.3 * X2, phi)
    )
  }))
}

run_population <- function() {
  res <- bind_rows(lapply(
    sort(unique(study$grid$scenario)),
    investigate_scenario
  ))
  save_investigation(res, "population.RDS")

  # results
  print(mutate(res, across(where(is.double), ~ round(.x, 3))), n = Inf)

  # table sorted by true PEHE
  print(arrange(res, scenario, true_pehe), n = Inf)

  invisible(res)
}

# ---- run ----------------------------------------------------------------------

run_population()
