# Confidence intervals for CATE estimates

A point estimate of a heterogeneous treatment effect is not much use without an
interval. These studies build **simultaneous** intervals — bands that cover every
unit at once, not each unit separately — using a half-sample bootstrap, and ask
how well they cover.

| folder | |
|---|---|
| `continuous/` | continuous outcome |
| `binary/` | binary outcome — **carries bug A**, see its README |
| `optimal_sf/` | calibrating the bootstrap's `sample.fraction` |

## The method

For each of `CI_boot` draws: take half of each fold, refit the second stage on
that half, and form the *half-sample root* `tau_full - tau_half`. Standardise
each root by its own bootstrap SD, take the maximum over units within each draw,
and use the `1 - alpha` quantile of that maximum as a **single** critical
value. That scalar is what makes the band simultaneous. The maximum is of
absolute roots, so both tails are already in it — `1 - alpha`, not
`1 - alpha/2`, gives a two-sided `1 - alpha` band. (Before 2026-09-28 the code
used `1 - alpha/2`, i.e. 97.5% bands at `alpha = 0.05` — any results from
before then, archived ones included, are 97.5% bands.)

Implemented in `R/bootstrap_ci.R` — `cf_half_boot` for the causal forest,
`rf_half_boot` for the DR-learner second stage, `simultaneous_band` for the
shared final step.

The causal forest also reports its own variance estimates, so the metrics score
a normal-approximation interval alongside the bootstrap one, labelled
`causal_forest_inbuilt`.

## Design

| | |
|---|---|
| scenarios | 1–10 |
| n | 500, 1000 |
| `CI_sf` | 0.05 to 0.5 in steps of 0.05 — the sweep `optimal_sf/` consumes |
| runs | 100 |
| array | **20,000 jobs** per outcome, split across two jobscripts |

`CI_sf` is the `sample.fraction` passed to the half-sample forests. It is swept
rather than fixed because the right value is not known a priori — too small and
the forests are noisy, too large and the half-samples stop being independent
enough.

## What these studies do *not* do

No SuperLearner arm, and no BLP or independence tests (`profile = "ci"`). The
question is interval coverage, not heterogeneity detection, and the bootstrap
already dominates the runtime.

## Metrics

`marginal_coverage` (proportion of units their own interval covers),
`simultaneous_coverage` (0/1 — did the band cover everything at once), and
`mean_ci_length`. Nominal is 0.95 for both coverages; simultaneous coverage is
the one the method is constructed to control.

`continuous/` and `binary/` additionally score a `"<model>_grid"` row against
a fixed covariate query grid (`R/dgm_scenarios.R::build_query_grid()`) rather
than the per-run sampled units. Those grid rows carry two more columns,
`be_marginal_coverage` and `be_simultaneous_coverage` — coverage of the same
interval against the across-run *mean* point estimate at each grid point, in
place of the true tau ("bias-eliminated coverage"; see
https://joonho112.github.io/simsum-mini-course/06-metrics-inference.html#sec-becoverage).
This isolates whether interval *width* is correctly calibrated, independent
of point-estimate bias. It's grid-only because BE-coverage needs a fixed
estimand replicated identically across runs — the query grid points are, but
the per-run sampled units are not (a fresh sample is drawn every run) — so
`be_marginal_coverage`/`be_simultaneous_coverage` are `NA` on every non-`_grid`
row (`optimal_sf/` never builds a query grid, so it has none of this).

## Correlated covariates

Built 2026-10-04 as `sample_size/correlated/confidence_intervals/`. All three
studies here (`continuous/`, `binary/`, `optimal_sf/`) are rerun on
`correlated/`'s scenarios 1–4 at ρ = 0 and 0.5, with this folder's design
otherwise: n ∈ {500, 1000}, 100 runs, and the full `CI_sf` sweep. That is
16,000 jobs per CI outcome and 1,600 per optimal_sf outcome. Its ρ = 0 arm has
this study's distribution but is not paired with these runs. See that
folder's README.

## Status

**Archive the old results first** - `R/archive_old_results.R` (root `README.md`, Status, step 0). They predate the current DGM and use the pre-2026-09-26 scenario numbers, so running into that tree would mix old and new results.

All three studies are also moved by the per-arm (T-learner) outcome models
the DR-learners now use (2026-09-27): `dr_random_forest` and `dr_semi_oracle`
fit one outcome forest per arm (`R/cate_models.R::t_learner_rf`). No
SuperLearner arm runs here (`sl_lib = NULL`), so the library change does not
apply.

`continuous/` — unaffected by the *bug ledger*, but **needs a re-run** for
the crossfitting strategy change to `R/cate_models.R` (see root README
Methods/Status).

`binary/` — **re-runs entirely.** See `binary/README.md`.

`optimal_sf/` — **both variants re-run.** `cts` for the crossfitting change
alone; `bin` for that plus bug A/the DGM issue. See
`optimal_sf/README.md`.
