##########
# title: metrics - competing-risk SNR check
##########
# As competing_risk/surv_metrics.R, for the arms in SNR_MODELS: every arm is a
# pseudo-value-scale estimator, scored against tau_RMTL1 / tau_RMTL2 with
# cate_metrics(), so the bias convention matches the rest of the repo. Pearson
# and Spearman are recomputed here as in surv_metrics.R (association()).
#
# Then the comparison the study exists for: each run's RMSE of every DR arm
# (semi-oracle twins included) divided by the same run's RMSE of
# pseudo_cf_whole_oob, rel_eff_snr.R's ratio. Below 1, DR did better.
#
# The oracle DR-learner (snr_oracle_config.R) is its own study with its own
# collected file. When snr_oracle_all.RDS exists its arm is merged into each
# run here, matched on scenario, n, hte_mult and run, and scored with the rest;
# a run with no oracle file just lacks that arm.
#
# Writes to the results folder:
#   snr_metrics.RDS  one row per (scenario, n, hte_mult, run, framework, target)
#   snr_ratios.RDS   per (scenario, n, hte_mult, target, framework): the run-level RMSE
#                    ratios summarised - median and IQR, and the geometric mean
#                    exp(mean log ratio) with the MCSE of the mean log ratio

library(here)
source(here("competing_risk", "snr", "snr_config.R"))
source(here("R", "metrics.R"))

# its own environment: the oracle config reassigns `study`
oracle <- new.env()
source(here("competing_risk", "snr", "snr_oracle_config.R"), local = oracle)

CF_ARM <- "pseudo_cf_whole_oob"
TRUTH  <- c(RMTL1 = "tau_RMTL1", RMTL2 = "tau_RMTL2")
EVENT_LABELS <- c(RMTL1 = "Event 1", RMTL2 = "Event 2")

# C-statistic = (Kendall tau_b + 1) / 2, as surv_metrics.R
c_statistic <- function(est, true) {
  (cor(true, est, method = "kendall", use = "pairwise.complete.obs") + 1) / 2
}

association <- function(est, true) {
  tibble(
    corr     = cor(true, est, use = "pairwise.complete.obs"),
    spearman = cor(true, est, method = "spearman", use = "pairwise.complete.obs"),
    c_stat   = c_statistic(est, true)
  )
}

#' One run's metrics: a row per (framework, target) the run carries
score_run <- function(sim_res, scenario) {
  map_dfr(intersect(c(SNR_MODELS, oracle$ORACLE_ARM), names(sim_res)), function(framework) {
    fw <- sim_res[[framework]]
    map_dfr(intersect(names(TRUTH), names(fw)), function(target) {
      est  <- fw[[target]]
      true <- sim_res$truth[[TRUTH[[target]]]]
      bind_cols(
        tibble(framework = framework, target = target),
        select(cate_metrics(est, true, scenario), -c(corr, spearman)),
        association(est, true)
      )
    })
  })
}

message("Reading snr_all.RDS...")
all_results_df <- readRDS(file.path(study$res_path, "snr_all.RDS"))

oracle_file <- file.path(oracle$study$res_path, "snr_oracle_all.RDS")
if (file.exists(oracle_file)) {
  message("Merging ", oracle$ORACLE_ARM, " from snr_oracle_all.RDS...")
  oracle_df <- readRDS(oracle_file)
  cell_key <- function(df) paste(df$scenario, df$n, df$hte_mult, sep = "_")
  oracle_cells <- setNames(oracle_df$results, cell_key(oracle_df))
  all_results_df$results <- map2(
    all_results_df$results, cell_key(all_results_df),
    function(results, key) {
      cell <- oracle_cells[[key]]
      if (is.null(cell)) return(results)
      by_run <- setNames(cell, vapply(cell, function(x) as.character(x$run), ""))
      map(results, function(entry) {
        found <- by_run[[as.character(entry$run)]]
        if (!is.null(found)) {
          entry$result[[oracle$ORACLE_ARM]] <- found$result[[oracle$ORACLE_ARM]]
        }
        entry
      })
    }
  )
  rm(oracle_df, oracle_cells)
} else {
  message("No snr_oracle_all.RDS - scoring without ", oracle$ORACLE_ARM)
}

metrics <- bind_rows(pmap(
  list(all_results_df$results, all_results_df$scenario, all_results_df$n,
       all_results_df$hte_mult),
  function(results, scenario, n, hte_mult) {
    map_dfr(results, function(entry) {
      mutate(score_run(entry$result, scenario), scenario = scenario, n = n,
             hte_mult = hte_mult, run = entry$run, .before = 1)
    })
  }
)) %>%
  mutate(target = unname(EVENT_LABELS[target]))
rm(all_results_df)

# paired per-run ratio to the causal forest, then summarised per cell
ratios <- metrics %>%
  filter(framework != CF_ARM) %>%
  select(scenario, n, hte_mult, run, target, framework, rmse) %>%
  inner_join(
    metrics %>%
      filter(framework == CF_ARM) %>%
      select(scenario, n, hte_mult, run, target, rmse_cf = rmse),
    by = c("scenario", "n", "hte_mult", "run", "target")
  ) %>%
  mutate(ratio = rmse / rmse_cf) %>%
  group_by(scenario, n, hte_mult, target, framework) %>%
  summarise(
    runs       = n(),
    median     = median(ratio),
    q25        = unname(quantile(ratio, 0.25)),
    q75        = unname(quantile(ratio, 0.75)),
    geo_mean   = exp(mean(log(ratio))),
    mcse_log   = sd(log(ratio)) / sqrt(n()),
    prop_dr_better = mean(ratio < 1),
    .groups = "drop"
  ) %>%
  arrange(scenario, target, framework, n, hte_mult)

print(as.data.frame(ratios), digits = 3)

saveRDS(metrics, file.path(study$res_path, "snr_metrics.RDS"))
saveRDS(ratios, file.path(study$res_path, "snr_ratios.RDS"))
print("metrics calculated!")
