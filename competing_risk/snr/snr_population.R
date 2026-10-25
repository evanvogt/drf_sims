##########
# title: population SNR per scenario and hte_mult - competing-risk SNR check
##########
# The check's x axis: SD(tau) / SD(noise) for each event, in each scenario at
# each hte_mult, from the same definition and draws as
# results_processing/thesis_figures/rel_eff_snr.R (dgm_covariates(), first
# DGM_N_POP_SURV rows, rho = 0.5; surv_moment_grid()). hte_mult scales b3_1 and
# b3_2 as generate_surv_data() does; scenarios 1-2 have neither, so they are at
# hte_mult 1 only. Reads no simulation output and needs no cluster (about a
# minute per scenario).
#
# Writes to the results folder:
#   snr_population.RDS  per (scenario, hte_mult, target): hazard ratios at X3 = 0 / 1,
#                       SD(tau), SD(noise) and SNR
#
#   cd competing_risk/snr; Rscript snr_population.R [results dir]
#   default: study$res_path (results/competing_risk_snr)

library(here)
source(here("R", "figures.R"))
source(here("results_processing", "dgm_population.R"))

args <- commandArgs(trailingOnly = TRUE)
res_path <- if (length(args) >= 1) args[1] else
  file.path(dirname(here()), "results", "competing_risk_snr")
dir.create(res_path, showWarnings = FALSE, recursive = TRUE)

SCENARIOS <- 1:7
HTE_LEVELS <- c(1, 1.5, 2, 3)
# the grid's cells (snr_config.R): scenarios 1-2 have no X3 effect to scale
hte_levels_for <- function(scenario) if (scenario %in% 1:2) 1 else HTE_LEVELS
EVENTS <- c(tau_RMTL1 = "Event 1", tau_RMTL2 = "Event 2")

cv <- dgm_covariates()
cv <- cv[seq_len(min(DGM_N_POP_SURV, nrow(cv))), ]

cells <- do.call(rbind, lapply(SCENARIOS, function(sc) {
  data.frame(scenario = sc, hte_mult = hte_levels_for(sc))
}))

population <- bind_rows(lapply(seq_len(nrow(cells)), function(i) {
  sc <- cells$scenario[i]
  m <- cells$hte_mult[i]
  p <- survival_scenario_params[survival_scenario_params$scenario == sc, ]
  p$b3_1 <- p$b3_1 * m
  p$b3_2 <- p$b3_2 * m
  gt <- surv_moment_grid(p)
  bind_rows(lapply(seq_along(EVENTS), function(k) {
    tau <- surv_truth_at_draws(gt, names(EVENTS)[k], cv$X1, cv$X2, cv$X3)
    noise <- surv_truth_at_draws(gt, paste0("noise_", k), cv$X1, cv$X2, cv$X3)
    tibble(
      scenario = sc, hte_mult = m, target = unname(EVENTS[k]),
      # Weibull scale effect b on the log scale -> hazard ratio exp(-shape * b)
      hr_x3_0 = exp(-p[[paste0("shape", k)]] * p[[paste0("bW_", k)]]),
      hr_x3_1 = exp(-p[[paste0("shape", k)]] *
                      (p[[paste0("bW_", k)]] + p[[paste0("b3_", k)]])),
      sd_tau = sd(tau), sd_noise = sqrt(mean(noise))
    )
  }))
})) %>%
  mutate(snr = sd_tau / sd_noise)

message("Population SNR per cell:")
print(as.data.frame(population), digits = 3)
saveRDS(population, file.path(res_path, "snr_population.RDS"))
