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

## ⚠ The grid disagreed three ways

Before the configs existed, this study's grid was declared in three places and
they did not match:

| | |
|---|---|
| `bin_analysis.R` | `scenario = c(1:10)` → 4,000 rows |
| `bin_check.R` / `bin_collect.R` | `scenario = c(1, 3, 8, 9)` → 1,600 rows |
| `jobscripts/bin_1.sh` | `#PBS -J 1-1600` |

4 scenarios × 4 sample sizes × 100 runs is exactly 1,600, so `c(1, 3, 8, 9)` was
the intent and the analysis script's `c(1:10)` was stale.

**Consequence for the results currently on disk.** `expand.grid` varies its first
column fastest, so submitting indices 1–1600 against the 4,000-row grid ran
*runs 1–40 of all ten scenarios*, not runs 1–100 of four. `bin_collect.R` then
looked for scenarios 1, 3, 8 and 9 and found 40 runs in each. **The study
silently has 40 replicates per cell, not 100** — Monte Carlo error is about 1.6×
larger than intended.

`bin_config.R` now declares `c(1, 3, 8, 9)` once, so index `i` means the same
thing everywhere. Worth confirming the replicate count against
`../results/binary/` before trusting any pre-existing figures.

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

**Re-runs required**, for two independent reasons:

1. bug F changes `dr_superlearner` (as in `continuous/`)
2. the grid fix means the intended design — 4 scenarios × 100 runs — has never
   actually been run
