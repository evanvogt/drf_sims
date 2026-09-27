# ADEMP — confidence intervals for the CATE

Three studies: `continuous/`, `binary/` and `optimal_sf/`.

## Aims

- **continuous / binary:** assess the coverage and width of half-sample
  bootstrap simultaneous confidence bands for the CATE, and how they depend on
  the bootstrap forests' `sample.fraction` (`CI_sf`).
- **optimal_sf:** choose `CI_sf` by a data-driven calibration that uses the
  point estimate as a plug-in truth.

## Data-generating mechanisms

The `continuous/` or `binary/` DGM (see their `ADEMP.md`), scenarios 1–10.

| | continuous / binary | optimal_sf |
|---|---|---|
| scenarios | 1–10 | 1–10 |
| n | 500, 1000 | 500, 1000 |
| `CI_sf` | 0.05 to 0.5 by 0.05 (design factor) | same values, as calibration candidates |
| repetitions | 100; 500 for scenarios 1–4 | 100 |

Full factorial over the rows above, per outcome type.

## Estimands

- Unit-level CATE for every unit in the simulated sample.
- CATE at a fixed covariate query grid (`build_query_grid()`), same points in every run (continuous / binary).
- optimal_sf: the estimated CATE `tau.hat` (plug-in target).

## Methods

Estimators (as `continuous/`, no SuperLearner, no HTE tests):
`causal_forest`, `dr_random_forest`, `dr_oracle`, `dr_semi_oracle`.

Intervals:

- **Half-sample bootstrap band** (`R/bootstrap_ci.R`): `CI_boot = 200` draws;
  refit the second stage on half the sample with nuisances held fixed; roots
  `tau_full − tau_half` standardised by their bootstrap SD; critical value =
  `1 − α/2` quantile of the per-draw maximum over units; α = 0.05.
- **`causal_forest_inbuilt`:** pointwise normal interval from grf's variance estimate.

optimal_sf (`find_optimal_sf()`): resample pseudo-outcome residuals around
`tau.hat` within folds, re-estimate, build the band at each candidate `CI_sf`,
and pick the value whose coverage of `tau.hat` is closest to 0.95.

## Performance measures

- marginal coverage (proportion of units covered by their own interval)
- simultaneous coverage (band covers every unit; the target the method controls)
- mean interval length
- query grid only: bias-eliminated marginal and simultaneous coverage
  (coverage of the across-run mean estimate)
- optimal_sf: mean coverage of `tau.hat` per candidate `CI_sf`, and the selected `CI_sf`

Nominal coverage 0.95.
