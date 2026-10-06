##########
# title: metrics for competing risk outcome
##########
# This study does not use compute_metrics(): its results nest as
# framework x target rather than one entry per model, and each target is scored
# against a different truth column. The point metrics themselves come from
# cate_metrics() so the bias convention matches the rest of the repo.

library(here)
source(here("competing_risk/surv_config.R"))
source(here("R", "metrics.R"))

# One future::multisession worker per task below; keep the jobscript's PBS
# ncpus/ompthreads in step with this (jobscripts/surv_metrics.sh). There are 28
# (rho, scenario, n, censoring) combos to share out; workers <= 1 runs
# sequentially instead.
workers <- 2

# The pseudo-value frameworks come in arms crossing two factors - see
# surv_models.R and the README:
#   whole_oob  whole-sample pseudo-values, whole-sample OOB fit (production parity)
#   whole_scf  whole-sample pseudo-values, single crossfit (the control)
#   cvps_scf   leave-one-fold-out pseudo-values, single crossfit
# The SuperLearner arms have no OOB analogue, so they are scf throughout and vary
# the pseudo-values only (sl_*_whole vs sl_*_cvps). The random survival forest
# DR-learner fits (Y, D) rather than pseudo-values, so it varies the fitting only
# (rsf_dr_oob vs rsf_dr_scf) - but it targets the same RMTL truth.
#
# NOTE: frameworks_run below is intersect(names(sim_res), frameworks), so a
# framework missing from ANY of these three lists is dropped silently rather than
# erroring. Add to all three together.
pseudo_arms <- c("whole_oob", "whole_scf", "cvps_scf")

frameworks <- c(
  "ipw", "csf_cs", "csf_sh",
  paste0("pseudo_cf_", pseudo_arms),
  paste0("pseudo_dr_", pseudo_arms),
  "sl_t_whole", "sl_t_cvps",
  "sl_dr_whole", "sl_dr_cvps",
  "sl_t_split",
  "rsf_dr_oob", "rsf_dr_scf"
)

# every pseudo-value framework shares the same targets and truth columns
pseudo_frameworks <- setdiff(frameworks, c("ipw", "csf_cs", "csf_sh"))
pseudo_targets    <- c("RMTL1", "RMTL2", "RMSTc")
pseudo_truth      <- c(RMTL1 = "tau_RMTL1", RMTL2 = "tau_RMTL2", RMSTc = "tau_RMSTc")

# which targets are valid per framework
framework_targets <- c(
  list(
    ipw    = c("RMST1", "RMST2", "RMSTc"),
    csf_cs = c("RMST1", "RMST2", "RMSTc"),
    csf_sh = c("RMST1", "RMST2")
  ),
  setNames(rep(list(pseudo_targets), length(pseudo_frameworks)), pseudo_frameworks)
)

# framework-specific truth column mapping.
# ipw and csf_cs remove competing events, so they target the cause-specific (net)
# RMST = integral of S*(t). csf_sh keeps competing events in the risk set
# (Fine-Gray) so it targets the subdistribution RMST = horizon - RMTL.
framework_truth_map <- c(
  list(
    ipw    = c(RMST1 = "tau_RMST1_cs", RMST2 = "tau_RMST2_cs", RMSTc = "tau_RMSTc"),
    csf_cs = c(RMST1 = "tau_RMST1_cs", RMST2 = "tau_RMST2_cs", RMSTc = "tau_RMSTc"),
    csf_sh = c(RMST1 = "tau_RMST1",    RMST2 = "tau_RMST2")
  ),
  setNames(rep(list(pseudo_truth), length(pseudo_frameworks)), pseudo_frameworks)
)

# C-statistic = (Kendall tau_b + 1) / 2, equivalent to Harrell's C for a
# continuous outcome.
c_statistic <- function(est, true) {
  (cor(true, est, method = "kendall", use = "pairwise.complete.obs") + 1) / 2
}

# No scenario here is null. cate_metrics() hard-codes scenario 1 as the
# no-heterogeneity scenario (the sample_size/ convention) and stores Pearson and
# Spearman there as 0, but this study's scenario 1 CATE varies with X1 and X2
# (ADEMP.md, "Performance measures"). So both are recomputed for every scenario.
# A constant truth (tau_RMST1_cs in scenarios 2 and 5, tau_RMST2_cs in 1 and 3)
# still gives NA, with cor()'s "standard deviation is zero" warning.
association <- function(est, true) {
  tibble(
    corr     = cor(true, est, use = "pairwise.complete.obs"),
    spearman = cor(true, est, method = "spearman", use = "pairwise.complete.obs"),
    c_stat   = c_statistic(est, true)
  )
}

#' One run's metrics: a row per (framework, target) the run carries
score_run <- function(sim_res, scenario) {

  truth          <- sim_res$truth
  frameworks_run <- intersect(names(sim_res), frameworks)

  map_dfr(frameworks_run, function(framework) {

    fw_data     <- sim_res[[framework]]
    targets_run <- intersect(names(fw_data), framework_targets[[framework]])

    map_dfr(targets_run, function(target) {

      model_tau <- fw_data[[target]]
      true_tau  <- truth[[framework_truth_map[[framework]][[target]]]]

      bind_cols(
        tibble(framework = framework, target = target),
        select(cate_metrics(model_tau, true_tau, scenario),
               -c(corr, spearman)),
        association(model_tau, true_tau)
      )
    })
  })
}

#' One (rho, scenario, n, censoring) combo's runs -> its metric rows
score_combo <- function(args) {
  map_dfr(args$results, function(entry) {
    mutate(score_run(entry$result, args$scenario),
           rho = args$rho, scenario = args$scenario, n = args$n,
           censoring = args$censoring, run = entry$run, .before = 1)
  })
}

message("Reading surv_all.RDS...")
all_results_df <- readRDS(file.path(study$res_path, "surv_all.RDS"))

# Bundle each combo's own slice into its own list element BEFORE handing it to
# future_map(), which ships .x[[i]] to the worker scoring it - mapping over row
# indices instead would make every worker export the whole of all_results_df as
# a captured global (see surv_nuisance_extract.R). Each run is also cut down to
# the truth and the framework estimates: its data and nuisances are what make
# surv_all.RDS large, and nothing here reads them.
combo_args <- pmap(
  list(all_results_df$results, all_results_df$rho, all_results_df$scenario,
       all_results_df$n, all_results_df$censoring),
  function(results, rho, scenario, n, censoring) {
    results <- map(results, function(entry) {
      keep <- intersect(names(entry$result), c("truth", frameworks))
      list(run = entry$run, result = entry$result[keep])
    })
    list(results = results, rho = rho, scenario = scenario, n = n,
         censoring = censoring)
  }
)
rm(all_results_df)
invisible(gc())

message("Scoring ", length(combo_args), " combos, ",
        sum(map_int(combo_args, ~ length(.x$results))), " runs total, ",
        workers, " worker(s)...")

if (workers > 1) {
  require(future)
  require(furrr)
  future::plan(future::multisession, workers = workers)
  on.exit(future::plan(future::sequential), add = TRUE)
  scored <- furrr::future_map(
    combo_args, score_combo,
    .options = furrr::furrr_options(packages = c("dplyr", "purrr", "tibble"))
  )
} else {
  scored <- map(combo_args, score_combo)
}

metrics <- bind_rows(scored) %>%
  # RMST and RMTL are inverses of each other, so label by event rather than scale
  mutate(target = case_when(target %in% c("RMST1", "RMTL1") ~ "Event 1",
                            target %in% c("RMST2", "RMTL2") ~ "Event 2",
                            target == "RMSTc" ~ "Combined"))

saveRDS(metrics, file.path(study$res_path, "surv_metrics.RDS"))
print("metrics calculated!")
