# Continuous outcome — sample size study

How well do doubly-robust and forest-based CATE estimators recover heterogeneous
treatment effects as sample size grows, for a continuous outcome?

This is the reference study: the binary, missing-data and confidence-interval
studies are all variations on it.

## Design

Ten scenarios varying the structure of the CATE, crossed with four sample sizes,
100 runs each — **4,000 array jobs** — plus runs 101–500 for scenarios 1, 3, 8
and 9, another 6,400 appended as grid rows 4001–10400.

| | |
|---|---|
| scenarios | 1–10 (see `R/dgm_scenarios.R`, `DESC_10`) |
| n | 100, 250, 500, 1000 |
| runs | 100; **500 for scenarios 1, 3, 8, 9** |
| folds | 4 at n=100, 5 at n=250, else 10 |
| results | `../results/continuous/scenario_<k>/<n>/res_sim_<run>.RDS` |

Folds are reduced at small n because the double-crossfitting procedure fits
nuisances over all `C(V,2)` fold pairs — 45 fits at V=10 — and the training
sets become too small otherwise.

The SuperLearner library also shrinks at n=100 (`SL.earth` and `SL.ranger` are
dropped), which is why `dr_superlearner` is sometimes filtered out of the n=100
figures.

### Scenarios

| # | CATE structure |
|---|---|
| 1 | no HTE (ATE only) |
| 2 | simple, binary variable (X3) |
| 3 | simple, continuous variable (X4) |
| 4 | two variables, additive |
| 5 | continuous × binary interaction |
| 6 | single effects + interaction |
| 7 | continuous × continuous interaction |
| 8 | single effects + a different interaction |
| 9 | cosine |
| 10 | exponential |

Each dataset also carries five deliberately unrelated covariates (`X01`–`X05`),
so the estimators have to find the signal rather than being handed it.

## Estimators

`causal_forest`, `dr_random_forest`, `dr_oracle`, `dr_semi_oracle`,
`dr_superlearner` — all defined in `R/cate_models.R`, all using double
crossfitting for the nuisances. The oracle uses the true outcome model and a
known propensity of 0.5; the semi-oracle knows only the propensity.

## Files

| file | role |
|---|---|
| `cts_config.R` | the parameter grid and results path — **the** definition |
| `cts_dgms.R` | names this study's slice of `R/dgm_scenarios.R` |
| `cts_models.R` | `family = gaussian()`, `profile = "base"` |
| `cts_analysis.R` | array entry point; one row of the grid per index |
| `cts_check.R` | finds missing runs, writes `jobscripts/failed_ids.txt`, and updates `-J` and the resource request in the rerun jobscript |
| `cts_collect.R` | gathers per-run files into `cts_all.RDS` |
| `cts_metrics.R` | computes `cts_metrics.RDS`, plus `cts_true_cate_tests.RDS` — see below |
| `results_cts.R`, `cts_results.Rmd` | summaries |

### True-CATE HTE test evaluation

`cts_true_cate_tests.RDS` reruns the BLP and independence tests
(`run_true_cate_tests()`, `R/cate_models.R`) against the *true* CATE and true
nuisances (`truth$tau`, `truth$p0`, `W.hat = 0.5`) instead of an estimator's
fitted ones — one `BLP_p`/`indep_cate` row per (scenario, n, run), with no
per-model dimension, since nothing here is estimated. This isolates the
tests' own size/power from any estimator's error: scenario 1 is the null
(no heterogeneity), scenarios 2-10 the alternative.

Scenario 1's true CATE is *exactly* constant (no estimation noise to give it
apparent variance), so `BLP_p` is `NA` for every scenario-1 run — `GenericML::BLP()`
cannot identify its interaction coefficient when tau has zero variance (the
same degenerate-tau guard as `hte_test_metrics()`'s `BLP_whole = NULL`
elsewhere) — while `indep_cate` still returns a real p-value there.

## Running it

```bash
qsub continuous/jobscripts/cts_1.sh     # 1-4000
qsub continuous/jobscripts/cts_extra.sh # 4001-10400: runs 101-500, scenarios 1, 3, 8, 9
Rscript continuous/cts_check.R          # writes failed_ids.txt if any are missing
qsub continuous/jobscripts/cts_collect.sh
qsub continuous/jobscripts/cts_metrics.sh
```

## Sizing the array job

`cts_1.sh` shipped with **unmeasured** `#PBS -l` lines: `ncpus=2:ompthreads=2` was
requested to match `cts_analysis.R`'s `workers <- 2`, but `num.threads` was never
set anywhere in `R/cate_models.R` (the shared estimator engine this study, and
six others, source), so grf always defaulted to using every visible core per
`multisession` worker regardless of what PBS allocated — the 2x2 pairing was a
guess, not a measurement. `R/cate_models.R`'s `cate_methods()` now accepts an
explicit `num.threads` argument (default `NULL`, so every other study using it
is unaffected), and `cts_analysis.R` forwards it from an optional trailing CLI
arg, the same pattern `crossfitting/cf_analysis.R` uses.

`cts_1.sh` now asks for one core and 2gb for an hour, and runs
`cts_analysis.R <i> 1 1` (one worker, one grf thread). Those figures are set by
hand. The `syrup` profiling sweep that was meant to measure them didn't work
for this study (see the root README's "Resource profiling (removed)"). When
changing them, size for the `n = 1000` cell, not a middle value: the array
spans `n ∈ {100, 250, 500, 1000}` under one `#PBS -l` line. Keep
`workers × grf_threads` within `ncpus`, and change the trailing args and the
`#PBS -l` line together.

## Status

**Re-runs required** — for the crossfitting strategy change to
`R/cate_models.R` (see root README Methods/Status), which moves all five
estimator arms, and separately for bug F, now fixed permanently, which changes
`dr_superlearner`: the second-stage SuperLearner library was pretested and
the result discarded, so failing algorithms were never dropped. Only the
`dr_superlearner` arm moves for bug F specifically — the harness can confirm
the other four are unchanged there.

Nothing else in this study was affected by the *bug ledger*. `bias` also
changes sign when the metrics are regenerated (bug G), but that needs no
cluster time.
