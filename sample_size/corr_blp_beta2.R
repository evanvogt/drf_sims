##########
# title: correlated sample-size sets - extract the saved BLP beta.2 coefficients
##########
# The metrics scripts reduce each estimator's saved BLP_whole to p-values
# (BLP_p, BLP_p_os - R/metrics.R::hte_test_metrics()), dropping the beta.2
# estimate itself. This pulls the beta.2 row back out of every run as saved at
# estimation time (R/cate_models.R::run_blp_whole(): homoskedastic SEs,
# two-sided p), for the estimators that save one - CATE_MODELS. Nothing is
# refitted, so there are no T-learner, HC3 or true-CATE rows: those BLPs only
# ever existed at metrics time.
#
# Reads the collected <prefix>_all.RDS, as the metrics scripts do, so it runs
# after *_corr_collect.R and covers exactly the runs that collected.
#
# One row per (rho, scenario, n, run, model):
#   beta2_est, beta2_se, beta2_t, beta2_p   the beta.2 row of BLP_whole
#   df                                      residual df of the BLP regression
# A model whose BLP_whole is NULL (constant tau - see run_blp_whole()) gets a
# row of NAs, as BLP_p is NA there. Results saved before 2026-09-27 kept only
# the Estimate and Pr(>|t|) columns, so se, t and df would be NA for those;
# every correlated run is later than that.
#
# Writes <res_path>/<prefix>_blp_beta2.RDS for each outcome. Run from
# sample_size/correlated/:
#   Rscript corr_blp_beta2.R                 # both outcomes
#   Rscript corr_blp_beta2.R continuous      # or just one

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tibble)
})
source(here("R", "metrics.R"))  # CATE_MODELS, compute_metrics()

outcomes <- commandArgs(trailingOnly = TRUE)
if (length(outcomes) == 0) outcomes <- c("continuous", "binary")
stopifnot(all(outcomes %in% c("continuous", "binary")))

config_file <- c(continuous = "continuous/cts_corr_config.R",
                 binary     = "binary/bin_corr_config.R")

#' The beta.2 row of one model's saved BLP_whole, as a one-row tibble
beta2_row <- function(blp) {
  na <- tibble(beta2_est = NA_real_, beta2_se = NA_real_, beta2_t = NA_real_,
               beta2_p = NA_real_, df = NA_real_)
  if (is.null(blp) || !"beta.2" %in% rownames(blp)) return(na)
  b <- blp["beta.2", ]
  full <- length(b) >= 4  # the pre-2026-09-27 shape is Estimate, Pr(>|t|) only
  df <- attr(blp, "df")
  tibble(
    beta2_est = unname(b[1]),
    beta2_se  = if (full) unname(b[2]) else NA_real_,
    beta2_t   = if (full) unname(b[3]) else NA_real_,
    beta2_p   = unname(b[length(b)]),
    df        = if (is.null(df)) NA_real_ else as.numeric(df)
  )
}

for (outcome in outcomes) {
  source(here("sample_size", "correlated", config_file[[outcome]]))  # -> study

  all_results_df <- readRDS(file.path(study$res_path,
                                      paste0(study$prefix, "_all.RDS")))

  beta2 <- compute_metrics(
    study, all_results_df, models = CATE_MODELS,
    per_model = function(model_res, true_tau, model, sim_res, keys) {
      beta2_row(model_res$BLP_whole)
    }
  )
  rm(all_results_df)
  gc()

  out <- file.path(study$res_path, paste0(study$prefix, "_blp_beta2.RDS"))
  saveRDS(beta2, out)
  message(outcome, ": ", nrow(distinct(beta2, across(c(all_of(study$path_cols), run)))),
          " runs, ", nrow(beta2), " rows -> ", out)
}
