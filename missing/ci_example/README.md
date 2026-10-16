# Confidence intervals after multiple imputation

A single cell of the missing-data design, asking a narrower question: once the
covariates have been multiply imputed, **how should the bootstrap intervals from
each imputation be combined?**

| | |
|---|---|
| scenarios | 1–4, the main scenarios |
| n | 500, MAR, `prop = 0.3`, `type = both` |
| method | `multiple_imputation` only |
| runs | 100 |
| array | **400 jobs** |
| results | `../results/missing/ci_example/scenario_<k>/500/both/0.3/MAR/multiple_imputation/` |
| figures | `cts_miss_ci_results.R` / `.qmd` — every metric, to `results/all_figures/`; the diagnostic counterpart to the chapter script `results_processing/thesis_figures/miss_cts_ci.R`. Both show all four scenarios. |

The grid is ordered by `(scenario, run)` rather than left in `expand.grid`
order. The array index is a row number, so reordering renumbers every job.

## The three pooling strategies

Each of 50 imputed datasets gets its own half-sample bootstrap; `combine_mi_ci()`
in `R/bootstrap_ci.R` then pools them three ways so they can be compared:

| strategy | how |
|---|---|
| `pooled` | empirical quantiles of the bootstrap replicates stacked across imputations |
| `mib` | Rubin's rules — within + between variance, critical value averaged over the per-imputation maxima |
| `hybrid` | one variance and one critical value from the stacked draws |

**Imputation within arm, with the outcome.** The 50 imputations are drawn
separately in each treatment arm with Y as a predictor (`impute_by_arm()`,
`R/missingness.R`), as in the main design - see `missing/ADEMP.md`.

**Coverage on complete and incomplete units.** The truth is τ at the
unamputed covariates, which an incomplete unit's interval cannot pin down from
the observed values, so marginal coverage and the point metrics are also
reported on the complete units (`_cu`) and the incomplete units (`_iu`)
separately (`interval_metrics_split()`, `R/metrics.R`). Each run saves its
missingness mask for this.

## To add: grid-based simultaneous coverage

`simultaneous_coverage` here is **per-unit** — it asks whether the band covers
every unit of the sample the model was fitted on. `sample_size/confidence_intervals/`
also evaluates coverage over a fixed covariate grid, which is the more
meaningful target for a band and does not move with the realised sample. The
machinery already exists in `R/dgm_scenarios.R`; this study just does not thread
it through.

Three call sites are needed, copying
`sample_size/confidence_intervals/continuous/cts_corr_ci_analysis.R`:

1. `cts_miss_ci_analysis.R` — build the query points and their truth:
   ```r
   Z_query     <- build_query_grid(scenario, set = "continuous_missing", ...)
   grid_truth  <- build_query_grid_truth(scenario, set = "continuous_missing",
                                         gen$bW, Z_query)
   ```
   and store `grid_truth` on the results object.
2. Pass `Z_query` into the model fits, so each returns `grid_lb` / `grid_ub`
   alongside its per-unit interval.
3. `cts_miss_ci_metrics.R` — emit the grid interval as its own row, labelled
   `"<model>_grid"`, the pattern `cts_corr_ci_metrics.R` already uses:
   ```r
   if (!is.null(model_res$grid_lb)) {
     bind_cols(tibble(model = paste0(model, "_grid")),
               interval_metrics(model_res$grid_lb, model_res$grid_ub,
                                sim_res$grid_truth$tau))
   }
   ```

Until that lands, the simultaneous-coverage numbers here and in
`sample_size/confidence_intervals/` are measuring different things and should
not be compared directly. `cts_miss_ci_results.qmd` says so on the relevant
section.

## Gotchas

**The `mib` and `hybrid` critical values are the `1 - alpha` quantile, not
`1 - alpha/2`.** They are quantiles of a maximum of *absolute* roots, so both
tails are already in it, as in `simultaneous_band()`. `pooled` uses `alpha/2`
and `1 - alpha/2`, correctly — it takes quantiles of signed replicates.

**The `hybrid` margin is `sqrt(lambda_hat) * S_star`**, as in the other two
strategies and `simultaneous_band()`. `sqrt(lambda_hat * S_star)` takes the
square root of the critical value and gives bands about half as wide. The
pooled bands are saved without their draws, so a wrong band cannot be patched
after the fact: it needs a re-run.

**`alpha` is an explicit argument of `combine_mi()`**, not read from the
global environment.
