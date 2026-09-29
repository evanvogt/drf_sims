# Confidence intervals after multiple imputation

A single cell of the missing-data design, asking a narrower question: once the
covariates have been multiply imputed, **how should the bootstrap intervals from
each imputation be combined?**

| | |
|---|---|
| scenarios | 1, 3, 4, 5, 6 (of the missing-data set, which uses the main study's numbers; 1–5 before 2026-09-26) |
| n | 500, MAR, `prop = 0.3`, `type = both` |
| method | `multiple_imputation` only |
| runs | 100 |
| array | **500 jobs** |
| results | `../results/missing/ci_example/scenario_<k>/500/both/0.3/MAR/multiple_imputation/` |
| figures | `cts_miss_ci_results.R` / `.qmd` — every metric, to `results/all_figures/`; the diagnostic counterpart to the chapter script `results_processing/thesis_figures/miss_cts_ci.R`. Keeps all five scenarios, where the chapter script shows four. |

The grid is ordered by `(scenario, run)` rather than left in `expand.grid` order.
That was already true and is preserved — the array index is a row number, so
reordering renumbers every job.

## The three pooling strategies

Each of 50 imputed datasets gets its own half-sample bootstrap; `combine_mi_ci()`
in `R/bootstrap_ci.R` then pools them three ways so they can be compared:

| strategy | how |
|---|---|
| `pooled` | empirical quantiles of the bootstrap replicates stacked across imputations |
| `mib` | Rubin's rules — within + between variance, critical value averaged over the per-imputation maxima |
| `hybrid` | one variance and one critical value from the stacked draws |

## To add: grid-based simultaneous coverage

`simultaneous_coverage` here is **per-unit** — it asks whether the band covers
every unit of the sample the model was fitted on. `confidence_intervals/` now
also evaluates coverage over a fixed covariate grid, which is the more
meaningful target for a band and does not move with the realised sample. The
machinery already exists in `R/dgm_scenarios.R`; this study just does not thread
it through.

Three call sites are needed, copying
`sample_size/confidence_intervals/continuous/cts_ci_analysis.R`:

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
   `"<model>_grid"`, the pattern `cts_ci_metrics.R` already uses:
   ```r
   if (!is.null(model_res$grid_lb)) {
     bind_cols(tibble(model = paste0(model, "_grid")),
               interval_metrics(model_res$grid_lb, model_res$grid_ub,
                                sim_res$grid_truth$tau))
   }
   ```

Until that lands, the simultaneous-coverage numbers here and in
`confidence_intervals/` are measuring different things and should not be
compared directly. `cts_miss_ci_results.qmd` says so on the relevant section.

## Gotchas

**`alpha` used to be a free variable in `combine_mi()`; it is now an explicit argument.**

**The `mib` and `hybrid` critical values are now the `1 - alpha` quantile (2026-09-28), not `1 - alpha/2`.** They are quantiles of a maximum of *absolute* roots, so both tails are already in it; `1 - alpha/2` gave 97.5% bands at `alpha = 0.05`. Same fix as `simultaneous_band()`, so each imputation's own band changes too. `pooled` still uses `alpha/2` and `1 - alpha/2`, correctly — it takes quantiles of signed replicates. Results from before this date are 97.5% bands.

**The hybrid margin looks like a typo.** It computes
`sqrt(lambda_hat * S_star)` where the other two strategies and
`simultaneous_band()` use `sqrt(lambda_hat) * S_star` — i.e. it takes the square
root of the critical value. Preserved as written, but flagged: if the hybrid
intervals look oddly narrow, this is why.

**The old `AUX` / `AUX-Y` mechanism spellings are gone** (2026-09-28): the names are now MAR / MNAR-Y0 / MNAR-tau and anything else is an error. This study runs MAR only. **Its covariates changed the same day:** `continuous_missing` now draws correlated covariates with X01–X03 as auxiliaries (`missing/ADEMP.md`), so every dataset here changes too.

**The header claimed imputations were reduced "from 50 to 20" to keep generation tractable; the code actually uses 50 and always has.**

## Status

**Needs a re-run** — the crossfitting strategy change to `R/cate_models.R`
(see root README Methods/Status) affects the estimators used here even
though there's no SuperLearner arm. Bug F does not apply (no SuperLearner
arm), but bug O does: it changed the `continuous_missing` baseline and `bW`
this study generates from (`missing/continuous/README.md`). Archive any old
tree first with `R/archive_old_results.R` (root `README.md`, Status, step 0) -
it would use the pre-2026-09-26 scenario numbers (1-5, now 1 and 3-6).
The 2026-09-28 correlated covariates are a further reason: every dataset
changes.
