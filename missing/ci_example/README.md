# Confidence intervals after multiple imputation

A single cell of the missing-data design, asking a narrower question: once the
covariates have been multiply imputed, **how should the bootstrap intervals from
each imputation be combined?**

| | |
|---|---|
| scenarios | 1–5 (the reduced set) |
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
`confidence_intervals/continuous/cts_ci_analysis.R`:

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

**`alpha` used to be a free variable.** `combine_mi()` read `alpha` from the
global environment at five places — it was never one of its arguments — and only
worked because `cts_miss_ci_analysis.R` happened to define `alpha <- 0.05` at top
level. It is now an explicit argument of `combine_mi_ci()`. Any other caller
would have silently used the wrong level or failed outright.

**The hybrid margin looks like a typo.** It computes
`sqrt(lambda_hat * S_star)` where the other two strategies and
`simultaneous_band()` use `sqrt(lambda_hat) * S_star` — i.e. it takes the square
root of the critical value. Preserved as written, but flagged: if the hybrid
intervals look oddly narrow, this is why.

**This study spells the mechanisms `AUX` / `AUX-Y`** where the others say
`MNAR` / `MNAR-Y`. Both are accepted and normalised in `R/missingness.R`. Its
grid only ever runs MAR, so the unobserved-`U` branches never fire here — which
is also why its divergent `p1` definition (it did not subtract the `U` term) made
no difference to anything it produced.

**The header used to claim imputations were reduced "from 50 to 20"** to keep
generation tractable. The code set `n_imp <- 50` and the completion message
printed `"<n_imp> 50 imputed datasets"`. 50 is what ran, and 50 is kept.

## Status

**Needs a re-run** — the crossfitting strategy change to `R/cate_models.R`
(see root README Methods/Status) affects the estimators used here even
though there's no SuperLearner arm. Unaffected by the *bug ledger* itself:
bug F does not apply (no SuperLearner arm), and the DGM is the continuous
one, so neither of those is a separate reason to re-run.
