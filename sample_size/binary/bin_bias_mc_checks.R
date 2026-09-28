##########
# title: Monte Carlo checks on the binary study's ATE bias
##########
# At n = 500 every model, including dr_oracle, showed an ATE bias of about
# +0.007 in all of scenarios 1-4 (~4 MCSEs). These checks, run 2026-09-28,
# traced it to the outcome draws in the particular streams used, not to the
# estimators, the code or the seeding. See README.md, "Monte Carlo checks on
# the ATE bias", and the thesis appendix (sc:bin_bias_mc).
#
#   1. dr_oracle ATE bias by run block (1-100 from bin_1.sh, 101-500 from
#      bin_extra.sh), with a per-scenario MCSE
#   2. how tightly the scenarios move together within a run - they share
#      random numbers, because setup_rng_stream() seeds on the run alone
#   3. every saved scenario 1 dataset regenerated from its seed with the
#      current code, plus a model-free bias: difference in means minus the
#      sample ATE
#   4. the same model-free bias in fresh streams beyond the study's (no
#      models fitted), to test the DGM and the seeding themselves
#   5. whether the reuse of each run's stream across n correlates the n cells:
#      model-free errors over the fresh streams, then the saved per-run bias
#      and MSE of the fitted models
#   6. each model's ATE bias relative to dr_oracle's on the same dataset, with
#      the paired MCSE - the comparison that removes the shared noise
#
# Result on 2026-09-28: every dataset reproduced; the difference in means
# carried the same bias; fresh streams 501-5000 were unbiased at every n
# (n = 500: +0.0002, z = 0.4); cross-n correlations were |r| <= 0.019
# (model-free) and |mean Spearman r| < 0.09 (bias, MSE).
#
# MCSEs are per scenario, never pooled across scenarios: scenarios 1-4 are
# ~0.9 correlated within a run, so pooling them overstates z roughly 2x.
#
# Fits no models. Needs bin_all.RDS and bin_metrics.RDS (after bin_collect.R
# and bin_metrics.R). From sample_size/binary/:
#   Rscript bin_bias_mc_checks.R [n_streams]
# n_streams (default 5000, must exceed the study's 500 runs) is how many
# streams check 4 generates; streams 1-500 reproduce the study's datasets and
# act as the control. Writes bin_bias_mc_checks.RDS to study$res_path.

library(dplyr)
library(tidyr)
library(purrr)
library(here)

source(here("sample_size/binary/bin_config.R"))
source(here("R", "utils.R"))
source(here("sample_size", "binary", "bin_dgms.R"))

args <- as.numeric(commandArgs(trailingOnly = TRUE))
n_streams <- if (length(args) >= 1 && !is.na(args[1])) args[1] else 5000

study_n <- sort(unique(study$grid$n))
study_runs <- max(study$grid$run)
stopifnot(n_streams > study_runs)

run_block <- function(run) if_else(run <= 100, "1-100", paste0("101-", study_runs))

# difference in means minus the sample ATE: the ATE error with no model at all
dim_error <- function(data, tau) {
  mean(data$Y[data$W == 1]) - mean(data$Y[data$W == 0]) - mean(tau)
}

metrics <- readRDS(file.path(study$res_path, "bin_metrics.RDS")) %>%
  filter(scenario %in% 1:4)

# ---- 1. dr_oracle ATE bias by run block ----------------------------------------

oracle_blocks <- metrics %>%
  filter(model == "dr_oracle") %>%
  group_by(scenario, n, block = run_block(run)) %>%
  summarise(mcse = sd(ate_bias) / sqrt(n()), bias = mean(ate_bias),
            z = bias / mcse, .groups = "drop")

print(oracle_blocks, n = Inf)

# ---- 2. correlation of dr_oracle's per-run ATE bias across scenarios -------------

scenario_cor <- map(set_names(study_n), function(size) {
  metrics %>%
    filter(model == "dr_oracle", n == size) %>%
    select(scenario, run, ate_bias) %>%
    pivot_wider(names_from = scenario, values_from = ate_bias) %>%
    select(-run) %>%
    cor(use = "pairwise.complete.obs")
})

print(scenario_cor)

# ---- 3. reproduce the saved scenario 1 datasets; model-free bias ----------------

all_results <- readRDS(file.path(study$res_path, "bin_all.RDS"))

reproduce <- all_results %>%
  filter(scenario == 1) %>%
  unnest_longer(results) %>%
  mutate(
    run = map_int(results, "run"),
    # as bin_analysis.R draws it
    same = map2_lgl(results, n, function(r, size) {
      setup_rng_stream(r$run)
      gen <- generate_binary_scenario_data(1, size)
      isTRUE(all.equal(gen$dataset, r$result$data))
    }),
    dim_err = map_dbl(results, ~ dim_error(.x$result$data, .x$result$truth$tau))
  ) %>%
  select(n, run, same, dim_err)

rm(all_results)
gc()

reproduce_summary <- reproduce %>%
  group_by(n, block = run_block(run)) %>%
  summarise(mismatch = sum(!same), mcse = sd(dim_err) / sqrt(n()),
            dim_err = mean(dim_err), .groups = "drop")

print(reproduce_summary)
print(filter(reproduce, !same) %>% arrange(n, run))

# ---- 4. model-free bias in fresh streams ------------------------------------------

# walks the streams setup_rng_stream() would give runs 1..n_streams, without
# its restart from the base seed for every run
RNGkind("L'Ecuyer-CMRG")
set.seed(formals(setup_rng_stream)$seed)
stream <- .Random.seed

fresh <- vector("list", n_streams)
for (r in seq_len(n_streams)) {
  stream <- parallel::nextRNGStream(stream)
  fresh[[r]] <- bind_rows(lapply(study_n, function(size) {
    assign(".Random.seed", stream, envir = globalenv())
    gen <- generate_binary_scenario_data(1, size)
    tibble(run = r, n = size, dim_err = dim_error(gen$dataset, gen$truth$tau))
  }))
}
fresh <- bind_rows(fresh)

fresh_summary <- fresh %>%
  mutate(block = cut(run, c(0, 100, study_runs, n_streams),
                     labels = c("1-100", paste0("101-", study_runs),
                                paste0(study_runs + 1, "-", n_streams)))) %>%
  group_by(n, block) %>%
  summarise(mcse = sd(dim_err) / sqrt(n()), bias = mean(dim_err),
            z = bias / mcse, .groups = "drop")

print(fresh_summary, n = Inf)

# ---- 5. correlation across n -------------------------------------------------

# model-free, over every stream; with n_streams = 5000, |r| < ~0.028 is noise
cross_n_dim <- fresh %>%
  pivot_wider(names_from = n, values_from = dim_err) %>%
  select(-run) %>%
  cor()

print(cross_n_dim)

# the fitted models' per-run bias and MSE, Spearman because per-run MSE is
# right-skewed; with 500 runs, |r| < ~0.09 is noise. Read mean_r: the cells in
# a row share datasets, so they clear the threshold together or not at all.
n_pairs <- combn(study_n, 2, simplify = FALSE)

cross_n_fits <- metrics %>%
  select(scenario, model, n, run, bias, mse) %>%
  pivot_longer(c(bias, mse), names_to = "metric") %>%
  group_by(scenario, model, metric) %>%
  group_modify(function(d, key) {
    w <- pivot_wider(d, names_from = n, values_from = value)
    map_dfr(n_pairs, function(p) {
      a <- w[[as.character(p[1])]]
      b <- w[[as.character(p[2])]]
      ok <- !is.na(a) & !is.na(b)
      tibble(pair = paste(p, collapse = " vs "), n_runs = sum(ok),
             r = cor(a[ok], b[ok], method = "spearman"))
    })
  }) %>%
  ungroup()

cross_n_fits_summary <- cross_n_fits %>%
  group_by(metric, pair) %>%
  summarise(mean_r = mean(r), max_abs_r = max(abs(r)),
            share_beyond = mean(abs(r) > 1.96 / sqrt(n_runs)), .groups = "drop")

print(cross_n_fits_summary)

# ---- 6. ATE bias relative to dr_oracle ------------------------------------------

oracle_ate_bias <- metrics %>%
  filter(model == "dr_oracle") %>%
  select(scenario, n, run, oracle = ate_bias)

vs_oracle <- metrics %>%
  filter(model != "dr_oracle") %>%
  left_join(oracle_ate_bias, by = c("scenario", "n", "run")) %>%
  mutate(rel = ate_bias - oracle) %>%
  group_by(scenario, n, model) %>%
  summarise(mcse = sd(rel, na.rm = TRUE) / sqrt(sum(!is.na(rel))),
            rel_bias = mean(rel, na.rm = TRUE), z = rel_bias / mcse,
            .groups = "drop") %>%
  arrange(desc(abs(z)))

print(vs_oracle, n = 20)

saveRDS(
  list(oracle_blocks = oracle_blocks, scenario_cor = scenario_cor,
       reproduce = reproduce, reproduce_summary = reproduce_summary,
       fresh = fresh, fresh_summary = fresh_summary, cross_n_dim = cross_n_dim,
       cross_n_fits = cross_n_fits, cross_n_fits_summary = cross_n_fits_summary,
       vs_oracle = vs_oracle),
  file.path(study$res_path, "bin_bias_mc_checks.RDS")
)

print("Monte Carlo checks complete!")
