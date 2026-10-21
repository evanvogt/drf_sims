##########
# title: metrics - correlated binary SNR check
##########
# As ../binary/bin_corr_metrics.R for the SNR_MODELS arms, plus each run's RMSE
# ratio to causal_forest (below 1 = better), as in competing_risk/snr/.
#
# The SNR on x is SD(tau) / sqrt(E[p(1 - p)]), as in rel_eff_snr.R, but from the
# CLAMPED truth over one seeded draw of POP_N covariates (te_moments() knows
# nothing of the clamp). bW is calibrated at SNR_CALIB_N, so n does not enter.
#
# Writes to the results folder:
#   bin_snr_population.RDS  per (scenario, hte_mult): bW, ATE, SD(tau),
#                           SD(noise), SNR and the share clamped
#   bin_snr_metrics.RDS     one row per (scenario, n, hte_mult, run, model)
#   bin_snr_ratios.RDS      RMSE ratios per (scenario, n, hte_mult, model):
#                           median, IQR, geometric mean, MCSE of the mean log
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

# ---- per-run metrics ------------------------------------------------------------

message("Reading bin_snr_all.RDS...")
all_results_df <- readRDS(file.path(res_path, "bin_snr_all.RDS"))

metrics <- compute_metrics(
  study, all_results_df, models = METRIC_MODELS,
  augment = add_t_learners,
  per_model = function(model_res, true_tau, model, sim_res, keys) {
    bind_cols(
      cate_metrics(model_res$tau, true_tau, keys$scenario),
      hte_test_metrics(model_res, sim_res, model)
    )
  }
)
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
print("metrics calculated!")
