##########
# title: metrics - correlated binary SNR check
##########
# As ../binary/bin_corr_metrics.R for the SNR_MODELS arms, plus each run's RMSE
# ratio to causal_forest (below 1 = better), as in competing_risk/snr/, and the
# RMSE against the truth of the nuisance estimates each arm is built on.
#
# The SNR on x is SD(tau) / sqrt(E[p(1 - p)]), as in rel_eff_snr.R, but from the
# clamped truth over one seeded draw of POP_N covariates (te_moments() knows
# nothing of the clamp). bW is calibrated at SNR_CALIB_N, so n does not enter.
#
# Writes to the results folder:
#   bin_snr_population.RDS  per (scenario, hte_mult): bW, ATE, SD(tau),
#                           SD(noise), SNR and the share clamped
#   bin_snr_metrics.RDS     one row per (scenario, n, hte_mult, run, model),
#                           with the nuisance RMSEs (nuisance_metrics())
#   bin_snr_ratios.RDS      RMSE ratios per (scenario, n, hte_mult, model):
#                           median, IQR, geometric mean, MCSE of the mean log
#   bin_snr_nuisance_hist.RDS  per (scenario, n, hte_mult, nuisance, source,
#                           bin): unit counts of the nuisance estimates and
#                           truths, pooled over the runs (nuisance_bins())
#
#   Rscript bin_snr_metrics.R [results dir]   # default: study$res_path

library(here)
source(here("sample_size", "snr", "bin_snr_config.R"))
source(here("R", "dgm_scenarios.R"))
source(here("R", "metrics.R"))
source(here("R", "cate_models.R"))

args <- commandArgs(trailingOnly = TRUE)
res_path <- if (length(args) >= 1) args[1] else study$res_path

CF_ARM <- "causal_forest"
METRIC_MODELS <- c(SNR_MODELS, "t_random_forest")
POP_N <- 1e6
POP_SEED <- 20261015

# ---- population SNR per cell ---------------------------------------------------

#' One (scenario, hte_mult) cell's population truth summaries
pop_cell <- function(scenario, hte_mult, cv) {
  tbl <- resolve_set(SNR_SET)
  p <- scale_hte(tbl[tbl$scenario == scenario, ], hte_mult)
  bW <- calibrate_bW(p, SNR_CALIB_N, calibration_for(SNR_SET))
  tr <- truth_at(p, bW, TRUE, cv$X1, cv$X2, cv$X3, cv$X4, cv$X5)
  sd_noise <- sqrt(0.5 * mean(tr$p0 * (1 - tr$p0)) + 0.5 * mean(tr$p1 * (1 - tr$p1)))
  tibble(
    scenario = scenario, hte_mult = hte_mult, bW = bW,
    ate = mean(tr$tau), sd_tau = sd(tr$tau), sd_noise = sd_noise,
    snr = sd(tr$tau) / sd_noise,
    share_clamped = mean(tr$p1 <= p$risk_lo | tr$p1 >= p$risk_hi)
  )
}

cells <- unique(study$grid[, c("scenario", "hte_mult")])
set.seed(POP_SEED)
pop_cv <- correlated_covariates(POP_N, resolve_set(SNR_SET)[1, ])
population <- bind_rows(lapply(seq_len(nrow(cells)), function(i) {
  pop_cell(cells$scenario[i], cells$hte_mult[i], pop_cv)
})) %>%
  arrange(scenario, hte_mult)
rm(pop_cv)

message("Population SNR per cell:")
print(as.data.frame(population), digits = 3)

# ---- nuisances ------------------------------------------------------------------

#' A run's per-unit nuisance estimates and their truths, as
#' list(nuisance = list(source = vector)): the sources are the arms the fits
#' belong to, and "truth"
#'
#' The true propensity is 0.5 (W ~ rbinom(n, 1, 0.5)) and the true marginal
#' mean E[Y | X] is (p0 + p1) / 2. causal_forest has grf's own propensity and
#' marginal mean; dr_random_forest has nuisance_rf()'s propensity (trimmed) and
#' per-arm outcome forests, with Y1.hat recovered through t_learner_tau(). A
#' unit's own-arm prediction is OOB, the other arm's comes from a forest that
#' never saw it, and both enter the pseudo-outcome. No other arm is listed:
#' t_random_forest shares dr_random_forest's fits, the oracles' are the truth,
#' and dr_semi_oracle's forests are nuisance_rf()'s estimator again.
nuisance_values <- function(sim_res) {
  tr <- sim_res$truth
  rf <- sim_res$nuisances_rf
  cf <- sim_res$causal_forest
  list(
    ps  = list(causal_forest = cf$W.hat, dr_random_forest = rf$W.hat, truth = 0.5),
    m   = list(causal_forest = cf$Y.hat.cf, truth = (tr$p0 + tr$p1) / 2),
    mu0 = list(dr_random_forest = rf$Y0.hat, truth = tr$p0),
    mu1 = list(dr_random_forest = rf$Y0.hat +
                 t_learner_tau(rf, sim_res$data$Y, sim_res$data$W),
               truth = tr$p1)
  )
}

#' RMSE against the truth of the nuisance estimates one arm is built on, over
#' all n units; NA where nuisance_values() has no such fit for the arm
nuisance_metrics <- function(model_res, model, sim_res) {
  # the oracle's outcome model is the truth, so a mismatch means truth and data are misaligned
  if (model == "dr_oracle") stopifnot(isTRUE(all.equal(model_res$Y0.hat, sim_res$truth$p0)))
  v <- nuisance_values(sim_res)
  rmse <- function(nu) {
    est <- v[[nu]][[model]]
    if (is.null(est)) NA_real_ else sqrt(mean((est - v[[nu]]$truth)^2))
  }
  tibble(ps_rmse = rmse("ps"), mu0_rmse = rmse("mu0"), mu1_rmse = rmse("mu1"),
         m_rmse = rmse("m"))
}

NUIS_BREAKS <- seq(0, 1, by = 0.01)

#' A run's nuisance estimates and truths binned on NUIS_BREAKS, one row per
#' (nuisance, source, non-empty bin). The propensity's truth, a constant 0.5,
#' is left out
nuisance_bins <- function(sim_res) {
  v <- nuisance_values(sim_res)
  v$ps$truth <- NULL
  bind_rows(lapply(names(v), function(nu) {
    tibble(nuisance = nu, source = rep(names(v[[nu]]), lengths(v[[nu]])),
           value = unlist(v[[nu]], use.names = FALSE))
  })) %>%
    mutate(bin = findInterval(value, NUIS_BREAKS, all.inside = TRUE)) %>%
    count(nuisance, source, bin, name = "count")
}

# ---- per-run metrics ------------------------------------------------------------

message("Reading bin_snr_all.RDS...")
all_results_df <- readRDS(file.path(res_path, "bin_snr_all.RDS"))

metrics <- compute_metrics(
  study, all_results_df, models = METRIC_MODELS,
  augment = add_t_learners,
  per_model = function(model_res, true_tau, model, sim_res, keys) {
    bind_cols(
      cate_metrics(model_res$tau, true_tau, keys$scenario),
      hte_test_metrics(model_res, sim_res, model),
      nuisance_metrics(model_res, model, sim_res)
    )
  }
)

# every unit's nuisance estimates and truths, pooled over a cell's runs
nuisance_hist <- compute_run_metrics(study, all_results_df, nuisance_bins) %>%
  group_by(scenario, n, hte_mult, nuisance, source, bin) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  mutate(mid = (NUIS_BREAKS[bin] + NUIS_BREAKS[bin + 1]) / 2)
rm(all_results_df)

# ---- paired per-run ratio to the causal forest ----------------------------------

ratios <- metrics %>%
  filter(model != CF_ARM) %>%
  select(scenario, n, hte_mult, run, model, rmse) %>%
  inner_join(
    metrics %>%
      filter(model == CF_ARM) %>%
      select(scenario, n, hte_mult, run, rmse_cf = rmse),
    by = c("scenario", "n", "hte_mult", "run")
  ) %>%
  mutate(ratio = rmse / rmse_cf) %>%
  filter(is.finite(ratio), ratio > 0) %>%
  group_by(scenario, n, hte_mult, model) %>%
  summarise(
    runs       = n(),
    median     = median(ratio),
    q25        = unname(quantile(ratio, 0.25)),
    q75        = unname(quantile(ratio, 0.75)),
    geo_mean   = exp(mean(log(ratio))),
    mcse_log   = sd(log(ratio)) / sqrt(n()),
    prop_better = mean(ratio < 1),
    .groups = "drop"
  ) %>%
  left_join(select(population, scenario, hte_mult, snr), by = c("scenario", "hte_mult")) %>%
  arrange(model, scenario, n, hte_mult)

print(as.data.frame(ratios), digits = 3)

saveRDS(population, file.path(res_path, "bin_snr_population.RDS"))
saveRDS(metrics, file.path(res_path, "bin_snr_metrics.RDS"))
saveRDS(ratios, file.path(res_path, "bin_snr_ratios.RDS"))
saveRDS(nuisance_hist, file.path(res_path, "bin_snr_nuisance_hist.RDS"))
print("metrics calculated!")
