##########
# title: metrics for the interim-analysis validation study
##########
# get_results() (R/pipeline.R) already attaches the grouping columns
# (scenario, n, interim_prop, run) to every collected run, replacing the
# hand-rolled four-level nested loops this used to take.
#
# This does NOT reuse R/metrics.R::compute_metrics() - that function's
# contract (per_model(model_res, true_tau, model, sim_res, keys)) assumes one
# flat set of model results and a sim_res$truth$tau per run. A cts_val run
# instead saves results1/results2/validations (two fitted chunks, already
# compared against each other by cts_val_analysis.R), so the flattening below
# is purpose-built - it reuses only the grouping-column plumbing from the
# study config.

library(here)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
source(here("validation/cts_val_config.R"))

all_results_df <- readRDS(file.path(study$res_path, "cts_val_all.RDS"))

#' Flatten the collected results into one row per (combination, run, model)
#'
#' @param extractor function(validations, key_row) returning a tibble of rows
#'   for one run's validations object, `key_row` giving the grouping columns
flatten_validations <- function(all_results_df, extractor) {
  group_cols <- study$path_cols

  df <- all_results_df %>%
    unnest_longer(results) %>%
    mutate(run = map_int(results, ~ .x$run),
           sim_res = map(results, ~ .x$result)) %>%
    select(-results)

  map_dfr(seq_len(nrow(df)), function(i) {
    key_row <- df[i, c(group_cols, "run"), drop = FALSE]
    extractor(df$sim_res[[i]]$validations, key_row)
  })
}

# ---- subgroup persistence: p-values for the top/bottom-10% responder groups
subgroups_metrics <- flatten_validations(all_results_df, function(val, key) {
  map_dfr(names(val$subgroups), function(model) {
    vals <- val$subgroups[[model]]
    bind_cols(key, tibble(model = model,
                          top_pval = unname(vals["top"]),
                          bottom_pval = unname(vals["bottom"])))
  })
})

# ---- CATE variance stability between chunks
variances_metrics <- flatten_validations(all_results_df, function(val, key) {
  map_dfr(names(val$variances), function(model) {
    vals <- val$variances[[model]]
    bind_cols(key, tibble(model = model,
                          vt1 = unname(vals["vt1"]),
                          vt2 = unname(vals["vt2"]),
                          var_change = unname(vals["vt2"]) - unname(vals["vt1"])))
  })
})

# ---- variable-importance rank stability between chunks
var_imp_metrics <- flatten_validations(all_results_df, function(val, key) {
  map_dfr(names(val$var_imps), function(model) {
    vi_df <- val$var_imps[[model]]
    bind_cols(key[rep(1, nrow(vi_df)), , drop = FALSE], tibble(model = model), vi_df)
  })
})

metrics <- list(subgroups = subgroups_metrics,
                variances = variances_metrics,
                var_imps = var_imp_metrics)

saveRDS(metrics, file.path(study$res_path, "cts_val_metrics.RDS"))
print("metrics calculated!")
