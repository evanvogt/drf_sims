##########
# title: what the single-event DGM implies, and that it is what it claims
##########
# Produces the numbers in ADEMP.md "What the DGM implies", after three checks
# that stop() if they fail:
#
#   A. rmst_weibull()'s closed form matches numerical integration;
#   B. se_truth()$tau_RMST is competing_risk's tau_RMST1_cs at the same
#      covariates (scenarios 1 / 2 / 3 here are parent scenarios 2 / 1 / 3 on
#      event 1: parent scenario 2 has no event-1 effect);
#   C. the pairing in se_dgm.R's header: under the same setup_rng_stream(),
#      W and X equal generate_surv_data(rho = 0)'s, and T_single >= T_CR.
#
# Then
#   1. event mix by the horizon per arm - generate_se_data() at n = 20,000 per
#      (scenario, censoring) cell, seeded once (illustrative, not the study's
#      streams);
#   2. population truths - se_truth() over N_POP covariate draws: mean and SD
#      of tau_RMST, the X3 = 1 minus X3 = 0 difference, control-arm RMST.
#
# Local, about a minute:  Rscript se_dgm_check.R  (from competing_risk/single_event/)

library(here)
library(dplyr)
library(future)
library(knitr)
source(here("R", "utils.R"))
source(here("competing_risk", "single_event", "se_dgm.R"))

plan(sequential)

N_MIX <- 20000   # units per event-mix cell
N_POP <- 200000  # covariate draws for the population truths
horizon <- se_scenario_params$event_horizon[1]
scenarios <- se_scenario_params$scenario

show <- function(df, digits = 2) {
  print(kable(df, format = "pipe", digits = digits))
  cat("\n")
}

# ---- A. closed form against integrate() -------------------------------------
scales <- exp(seq(log(5), log(80), length.out = 25))
closed <- rmst_weibull(scales, 2, horizon)
numeric_ <- vapply(scales, function(s) {
  integrate(function(t) exp(-(t / s)^2), 0, horizon, rel.tol = 1e-12)$value
}, numeric(1))
err_a <- max(abs(closed - numeric_))
cat(sprintf("A. rmst_weibull vs integrate: max |diff| = %.2e\n", err_a))
stopifnot(err_a < 1e-8)

# ---- B. truth equals the parent's net RMST1 ---------------------------------
set.seed(20261004)
n_b <- 300
X1 <- rbinom(n_b, 1, 0.4)
X2 <- rnorm(n_b)
X3 <- rbinom(n_b, 1, 0.7)
parent_of <- c(`1` = 2, `2` = 1, `3` = 3)
err_b <- vapply(scenarios, function(s) {
  ours <- se_truth(se_scenario_params[s, ], X1, X2, X3)$tau_RMST
  theirs <- surv_truth(survival_scenario_params[parent_of[[as.character(s)]], ],
                       X1, X2, X3)$tau_RMST1_cs
  max(abs(ours - theirs))
}, numeric(1))
cat(sprintf("B. tau_RMST vs parent tau_RMST1_cs, scenario %d: max |diff| = %.2e\n",
            scenarios, err_b), sep = "")
stopifnot(all(err_b < 1e-6))

# ---- C. pairing with competing_risk at rho = 0 --------------------------------
x_cols <- c("W", "X1", "X2", "X3", "X01", "X02", "X03", "X04", "X05")
for (run in c(1, 17)) {
  for (s in scenarios) {
    setup_rng_stream(run)
    se <- generate_se_data(s, 500, return_truth = FALSE)$dataset
    setup_rng_stream(run)
    cr <- generate_surv_data(parent_of[[as.character(s)]], 500, rho = 0,
                             return_truth = FALSE)$dataset
    stopifnot(identical(se[x_cols], cr[x_cols]), all(se$Y >= cr$Y - 1e-3))
  }
}
cat("C. pairing: W and X identical to competing_risk rho = 0, T_single >= T_CR\n\n")

# ---- 1. event mix -------------------------------------------------------------
cat("=== 1. Event mix by the horizon (n =", N_MIX, "per cell) ===\n\n")
set.seed(20261004)
mix <- bind_rows(lapply(scenarios, function(s) {
  bind_rows(lapply(c(FALSE, TRUE), function(cens) {
    d <- generate_se_data(s, N_MIX, return_truth = FALSE, censoring = cens)$dataset
    d %>%
      group_by(arm = ifelse(W == 1, "treated", "control")) %>%
      summarise(event_by_h = mean(D == 1 & Y <= horizon),
                censored_before_h = mean(D == 0 & Y < horizon),
                observed_past_h = mean(Y > horizon),
                .groups = "drop") %>%
      mutate(scenario = s, censoring = cens, .before = 1)
  }))
}))
show(mix, 3)

# ---- 2. population truths -----------------------------------------------------
cat("=== 2. Population truths (", N_POP, "covariate draws) ===\n\n")
set.seed(20261004)
cv <- correlated_covariates(N_POP, list(X1_prob = 0.4, X3_prob = 0.7, rho = 0,
                                        s2 = 1, s4 = 1, s5 = 1))
truths <- bind_rows(lapply(scenarios, function(s) {
  tr <- se_truth(se_scenario_params[s, ], cv$X1, cv$X2, cv$X3)
  tibble(
    scenario = s,
    description = se_scenario_params$description[s],
    RMST_0 = mean(tr$RMST_0),
    RMST_1 = mean(tr$RMST_1),
    tau_mean = mean(tr$tau_RMST),
    tau_sd = sd(tr$tau_RMST),
    tau_X3_0 = mean(tr$tau_RMST[cv$X3 == 0]),
    tau_X3_1 = mean(tr$tau_RMST[cv$X3 == 1]),
    tau_min = min(tr$tau_RMST),
    tau_max = max(tr$tau_RMST)
  )
}))
show(truths, 3)
