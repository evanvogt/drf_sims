# `R/` — the shared library

Every study sources from here. Nothing in this directory knows about a specific
study: the differences between them are arguments, not forks.

## Files

| file | holds |
|---|---|
| `utils.R` | `setup_rng_stream`, `collate_predictions`, `scatter_folds`, `trim_ps`, `timed` |
| `dgm_scenarios.R` | scenario tables + `generate_scenario_data()`, `get_oracle_info()` |
| `missingness.R` | `introduce_missingness`, `handle_missingness` (by-arm imputation via `impute_by_arm`), `generate_and_process_data`, `generate_reference_data`, `amputation_mask`, `u_term_by_completeness` |
| `cate_models.R` | the DR-learner family, causal forest, `cate_methods()`, `combine_mi` |
| `sl_library.R` | per-nuisance SuperLearner libraries (`sl_libraries()`), custom learner wrappers, `sl_fit_predict()`, `pretest_superlearner()` |
| `bootstrap_ci.R` | `cf_half_boot`, `rf_half_boot`, `combine_mi_ci`, `find_optimal_sf` |
| `metrics.R` | `cate_metrics`, `cate_metrics_split` (complete / incomplete units, missing-data studies), `interval_metrics`, `interval_metrics_split`, `compute_metrics` |
| `pipeline.R` | `study_config`, `get_results`, `check_failed`, `grid_indices` |
| `figures.R` | display labels, palette/theme, `summarise_metrics`, plot helpers |
| `regression_check.R` | the old-vs-new equivalence harness |
| `bin_verify_hte.R` | checks the binary risk-difference scenarios and re-derives `RD_SCALE`, `RD_SCALE_MISS`/`BU_MISS` and `RD_SCALE_CORR` (`Rscript R/bin_verify_hte.R`) |
| `calibration_report.R` | planned vs realised power and `bW` per scenario |
| `archive_old_results.R` | packs results trees into `results/_archive/<label>/*.tar` (`--trees`, `--label`) |

## `cate_methods()` — the axes

| argument | what it changes | who uses it |
|---|---|---|
| `family` | SuperLearner outcome model: `binomial()` adds `method.NNloglik` | binary studies |
| `ipw` | grf `sample.weights` / SuperLearner `obsWeights`. `NULL` is the unweighted path exactly | the missing-data IPW arm |
| `ci` | `list(boot=, sf=, alpha=)` turns on the half-sample bootstrap | the CI studies |
| `models` | `NULL` runs every arm; a subset of `CATE_MODELS` runs only those (e.g. no oracle-type arms) | applied analyses |
| `sl_plan` | a future plan (e.g. `tweak(multisession, workers = 5)`) set for the SuperLearner arm alone, then switched back: its learners are single-threaded and want more workers than the grf arms, whose memory is the limit. Results are unchanged (furrr's seeds do not depend on the worker count) | applied analyses |
| `X_ps` | propensity-only covariates (e.g. calendar time under a drifting allocation): added to every *estimated* propensity model - `nuisance_rf`'s, the causal forest's (which then takes `nuisance_rf`'s `W.hat` instead of its own), the SuperLearner arm's - and to nothing else. `all_cate_surv_models()` takes the same argument, along with `models` and `estimands` | applied analyses |

The last three default to the behaviour without them exactly.

Every oracle formula in `dgm_scenarios.R` returns the outcome mean
`E[Y | X, W]` — the linear predictor for a continuous outcome, the risk for a
binary one — so `run_dr_oracle()` applies no link.

### Orchestration profiles

Which post-estimation tests run is set by `profile` (`PROFILES` in
`cate_models.R`):

| | `base` | `ci` | `missing` |
|---|---|---|---|
| causal forest variance | no | yes | yes |
| BLP / independence tests, every model | yes | no | yes |
| SuperLearner, oracle and semi-oracle skipped when `X` has NAs | no | no | yes |

The CI profile drops the tests only to save runtime. The SuperLearner arm and
the bootstrap are set by `sl_lib` and `ci`, under any profile, so applied
analyses use `missing` with `ci` set.

`multiple_imputation` runs keep each imputation's tests unpooled, in `mi_tests`
(`mi_test_table()`). The pooling rule is still to be decided; see
`missing/README.md`.

## The parameter grid contract

Each study declares its grid **once**, in `<study>/<prefix>_config.R`. The PBS
array index is a **row number** of `study$grid`, so:

> Never filter or reorder the grid after construction.

Doing so renumbers every job, so index `i` would mean different things in the
analysis and check scripts. To run a subset, select indices and leave the
numbering alone:

```r
idx <- grid_indices(study, method = "complete_data")
```

## DGM draw order

`generate_scenario_data()` draws in this order, and the order is part of the
contract — every study reproduces runs by index through `setup_rng_stream()`:

```
W, X1, X2, X3, X4, X5, [U], [err], X01, X02, X03, cats
```

`U` only for the MNAR mechanisms (MNAR-Y0, MNAR-tau), `err` only for
continuous outcomes. The correlated sets (`*_corr_<ρ>`, `continuous_missing`,
`binary_missing`) draw `W, Z-block (X1–X5, X01–X03), [U], [err], cats` instead
(`correlated_covariates()`). Every scenario draws and returns `X1`–`X5`,
whether or not its treatment effect uses them. `regression_check.R`
fingerprints the generated dataset, not just the estimates, so a change here
fails loudly.

Both outcomes are one model, `E[Y | x, W] = m0(x) + W·τ(x)` (`control_mean()`
gives `m0`). For a binary outcome `m0` is a logistic scaled into
`[p0_lo, p0_hi]`, so the treatment effect adds on the risk-difference scale;
`sample_size/README.md` has the design and why.

## Verifying a change

```bash
Rscript tests/run_tests.R               # quick invariants, no results needed
Rscript R/regression_check.R baseline   # before touching anything
Rscript R/regression_check.R verify     # after each step - every study identical
Rscript crossfitting/cf_testing.R       # independent check of the estimators
```

The harness runs each study in its own subprocess, because
`run_all_cate_methods` is defined in each study's models file and sourcing two
studies into one session would silently give the second one's definition to
both.

Scope: it proves the changed code reproduces the **current** code on **this
machine**. It is not a claim about the cluster's numbers — R 4.5.3 here versus
4.3.2 there. Re-capture the baseline if you change machine.
