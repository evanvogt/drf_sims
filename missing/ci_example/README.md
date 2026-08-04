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

Unaffected by the bug ledger — no SuperLearner arm, so bug F does not apply, and
the DGM is the continuous one. No re-run needed.
