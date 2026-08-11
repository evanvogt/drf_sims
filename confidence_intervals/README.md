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
and use the `1 - alpha/2` quantile of that maximum as a **single** critical
value. That scalar is what makes the band simultaneous.

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

## Status

`continuous/` — unaffected by the *bug ledger*, but **needs a re-run** for
the crossfitting strategy change to `R/cate_models.R` (see root README
Methods/Status).

`binary/` — **re-runs entirely.** See `binary/README.md`.

`optimal_sf/` — **both variants re-run.** `cts` for the crossfitting change
alone; `bin` for that plus bug A/the DGM issue. See
`optimal_sf/README.md`.
