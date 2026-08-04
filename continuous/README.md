# Continuous outcome — sample size study

How well do doubly-robust and forest-based CATE estimators recover heterogeneous
treatment effects as sample size grows, for a continuous outcome?

This is the reference study: the binary, missing-data and confidence-interval
studies are all variations on it.

## Design

Ten scenarios varying the structure of the CATE, crossed with four sample sizes,
100 runs each — **4,000 array jobs**.

| | |
|---|---|
| scenarios | 1–10 (see `R/dgm_scenarios.R`, `DESC_10`) |
| n | 100, 250, 500, 1000 |
| runs | 100 |
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
| `cts_check.R` | finds missing runs, writes `jobscripts/failed_ids.txt` |
| `cts_collect.R` | gathers per-run files into `cts_all.RDS` |
| `cts_metrics.R` | computes `cts_metrics.RDS` |
| `results_cts.R`, `cts_results.Rmd` | summaries |

## Running it

```bash
qsub continuous/jobscripts/cts_1.sh     # 1-4000
Rscript continuous/cts_check.R          # writes failed_ids.txt if any are missing
qsub continuous/jobscripts/cts_collect.sh
qsub continuous/jobscripts/cts_metrics.sh
```

## Status

**Re-runs required.** Fixing bug F (`PRETEST_STAGE2`) changes `dr_superlearner`:
the second-stage SuperLearner library was pretested and the result discarded, so
failing algorithms were never dropped. Only that arm should move — the harness
can confirm the other four are unchanged.

Nothing else in this study was affected by the bug ledger. `bias` also changes
sign when the metrics are regenerated (bug G), but that needs no cluster time.
