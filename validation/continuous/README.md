# Interim-analysis validation — continuous outcome

The only validation study run today. See `validation/README.md` for what a
validation study checks and why.

## Design

| | |
|---|---|
| scenario | 3 (continuous × binary interaction — see `R/dgm_scenarios.R`, `DESC_10`) |
| n | 1000, split into two chunks of `n * interim_prop` and `n * (1 - interim_prop)` |
| interim_prop | 0.25, 0.5, 0.75 |
| runs | 100 — **300 array jobs** |
| folds | 5 below n=250 per chunk, else 10 |
| results | `../results/cts_val/scenario_3/1000/<interim_prop>/res_sim_<run>.RDS` |

Only scenario 3 and n = 1000 are run today; both are still columns of the grid
(`cts_val_config.R`) so extending to more scenarios or sizes is a one-line
change, not a rewrite.

## Estimators

Only `causal_forest` and `dr_random_forest`, both from `R/cate_models.R` — the
question here is whether one estimator's own findings persist across chunks,
not which estimator is most accurate, so the oracle/semi-oracle/SuperLearner
arms this repo's other studies carry are not run.

Each estimator's per-run object also carries `te_vims`: a treatment-effect
variable importance measure computed by refitting the second stage with each
covariate dropped in turn and comparing out-of-fold prediction error. This is
the one piece of estimation logic that is specific to this study — it isn't
duplicated anywhere else in the repo, so it stays in `cts_val_models.R` rather
than moving into `R/`.

## Files

| file | role |
|---|---|
| `cts_val_config.R` | the parameter grid and results path — **the** definition |
| `cts_val_dgms.R` | names this study's slice of `R/dgm_scenarios.R` |
| `cts_val_models.R` | wraps `R/cate_models.R`'s causal_forest/DR-RF + the TE-VIM helpers |
| `cts_val_analysis.R` | array entry point; fits both chunks, computes the three chunk comparisons |
| `cts_val_check.R` | finds missing runs, writes `jobscripts/failed_ids.txt`, and updates `-J` and the resource request in the rerun jobscript |
| `cts_val_collect.R` | gathers per-run files into `cts_val_all.RDS` |
| `cts_val_metrics.R` | flattens `validations` into tidy `cts_val_metrics.RDS` |
| `results_cts_val.R` | summary plots |

## Running it

```bash
qsub validation/continuous/jobscripts/cts_val_1.sh     # 1-300
Rscript validation/continuous/cts_val_check.R          # writes failed_ids.txt if any are missing
qsub validation/continuous/jobscripts/cts_val_collect.sh
qsub validation/continuous/jobscripts/cts_val_metrics.sh
```

## Status

**Needs a fresh run** — for two independent reasons. First, the crossfitting
strategy change to `R/cate_models.R` (see root README Methods/Status) affects
both estimators used here (`causal_forest`, `dr_random_forest`). Second, this
folder was rebuilt from a pre-restructure fork (standalone copies of the DGM
and the CATE estimators, grid/index arithmetic done by hand in
`validation.sh`, hardcoded absolute cluster paths) onto the shared `R/`
pattern every other study now uses. Three things changed on purpose, not by
accident, as part of that rebuild:

- **Causal-forest nuisance.** The old code built the causal forest's outcome
  nuisance from a regression of `Y` on `(W, X)` with the observed `W` plugged
  back in. `R/cate_models.R::run_causal_forest()` instead uses `Y.hat.cf`, a
  separate `X`-only regression — the input `causal_forest()` actually expects.
  This is very likely a bug fix, but it changes the causal-forest arm's
  numbers (and its TE-VIMs, which reuse the same nuisances). Any results from
  before this refactor are not comparable and should not be reused.
- **Results location.** Results now write to `../results/cts_val/...`
  (outside the repo, alongside every other study), not the old
  `live/results/validation/...` cluster path. There is no migration step —
  see the point above.
- **The collect/metrics handoff was broken and is now just fixed.**
  `collect_validation.R` wrote `validation_all_tidy.RDS`; `metrics_validation.R`
  read `validation_all.RDS` with an incompatible nested-list structure, and
  also carried a stale, wider scenario/size grid than anything actually run.
  `cts_val_collect.R` and `cts_val_metrics.R` now share one grid
  (`cts_val_config.R`) and one file (`cts_val_all.RDS`).

**Not implemented.** The "compare HTE tests between chunks" comparison was
stubbed in the original code and still is. Both estimators' per-run objects
now carry `BLP_whole`/`independence_cate`/`independence_po` in the same shape
`R/metrics.R::hte_test_metrics()` consumes, which is what a fourth chunk
comparison (alongside subgroups/variance/var-imps) would build on.

**Not in scope for this refactor.** `R/regression_check.R`'s old-vs-new
harness compares a behaviour-*preserving* refactor against a baseline; the
causal-forest nuisance change above means there is no valid "before" to check
against here. Worth adding once a fresh baseline exists.
