# Binary outcome — sample size study

The continuous study, on a logit scale. Same estimators, same crossfitting, same
sample-size sweep; the outcome is `rbinom(n, 1, plogis(lp))` and the estimand is
a **risk difference**.

## Design

| | |
|---|---|
| scenarios | **1, 3, 8, 9** — a subset, not all ten |
| n | 100, 250, 500, 1000 |
| runs | 100 |
| array | **1,600 jobs** |
| results | `../results/binary/scenario_<k>/<n>/res_sim_<run>.RDS` |

The four scenarios are the ones the chapter reports: null, simple, complex and
non-linear (`SS_SCENARIO_LABELS` in `R/figures.R`).

The coefficient table differs from the continuous study — `b0 = -0.4`,
`b1 = 0.5`, `b2 = 0.5` rather than the continuous values — because the same
numbers on a logit scale would saturate `plogis`. Scenario 10's treatment effect
also differs (`exp(X4)` rather than `exp(-abs(X4))`), and `bW` is calibrated with
`power.prop.test` rather than `power.t.test`. Those are the only differences from
`continuous/`.

## The grid was declared three ways

Before the configs existed, this study's grid appeared in three places and they
did not agree:

| | |
|---|---|
| `bin_analysis.R` | `scenario = c(1:10)` |
| `bin_check.R` / `bin_collect.R` | `scenario = c(1, 3, 8, 9)` |
| `jobscripts/bin_1.sh` | `#PBS -J 1-1600` |

4 scenarios × 4 sample sizes × 100 runs is exactly 1,600, so `c(1, 3, 8, 9)` is
the design and the analysis script's `c(1:10)` was the stale line.

**The results on disk are the intended ones** — four scenarios, 100 runs each.
The stale `c(1:10)` did not corrupt them; it was simply out of step with the
grid the runs were actually launched from.

`bin_config.R` now declares `c(1, 3, 8, 9)` once and every script reads it from
there, so the three-way drift cannot recur.

## Files

Same layout as `continuous/`. `bin_models.R` sets `family = binomial()` and
`oracle_link = "logit"` — the oracle formula from `bin_dgms.R` is a linear
predictor, so the model code applies `plogis`. (`missing/binary/` uses the
opposite convention; see `R/README.md`.)

## Running it

```bash
qsub binary/jobscripts/bin_1.sh     # 1-1600
Rscript binary/bin_check.R
qsub binary/jobscripts/bin_collect.sh
```

## Status

**Re-runs required for bug F only**, which changes `dr_superlearner` as in
`continuous/`. Only that arm should move; the harness can confirm the other four
are unchanged.

The DGM is unaffected by the bug ledger — bug A is the *confidence-interval*
binary study, not this one. `bias` also flips sign when the metrics are
regenerated (bug G), but that needs no cluster time.
