# `R/` — the shared library

Every study sources from here. Nothing in this directory knows about a specific
study: the differences between them are arguments, not forks.

This folder exists because they used to be forks. The same CATE estimation code
lived in seven files, the same DGM in four, and the same collect/check
boilerplate in eight — `sample_size/continuous/cts_models.R` and `sample_size/binary/bin_models.R`
differed in **two** places out of 438 lines. Consolidating removed roughly 3,000
lines and, more usefully, removed the possibility of the copies drifting apart,
which is where most of the bugs in the ledger below came from.

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

The repo-root `utils.R` is a two-line shim onto `R/utils.R`, so existing
`source(here("utils.R"))` calls keep working.

## `cate_methods()` — the axes

The seven forked model files differed on exactly four things. Three are now
arguments:

| argument | what it changes | who uses it |
|---|---|---|
| `family` | SuperLearner outcome model: `binomial()` adds `method.NNloglik` | binary studies |
| `ipw` | grf `sample.weights` / SuperLearner `obsWeights`. `NULL` is the unweighted path exactly | the missing-data IPW arm |
| `ci` | `list(boot=, sf=, alpha=)` turns on the half-sample bootstrap | the CI studies |

Two more were added for applied analyses (2026-10-06). Both default to the old
behaviour exactly:

| argument | what it changes |
|---|---|
| `models` | `NULL` runs every arm; a subset of `CATE_MODELS` runs only those (e.g. no oracle-type arms) |
| `X_ps` | propensity-only covariates (e.g. calendar time under a drifting allocation): added to every *estimated* propensity model - `nuisance_rf`'s, the causal forest's (which then takes `nuisance_rf`'s `W.hat` instead of its own), the SuperLearner arm's - and to nothing else. `all_cate_surv_models()` takes the same argument, along with `models` and `estimands` |

The fourth was where the oracle arm's inverse link lives, and it is gone. Every
oracle formula in `dgm_scenarios.R` returns the outcome mean `E[Y | X, W]` —
the linear predictor for a continuous outcome, the risk for a binary one — so
`run_dr_oracle()` applies no link. While the binary DGM was on the logit scale
its oracle formulas were linear predictors, and an `oracle_link` argument said
whether to apply `plogis`. That argument was not implied by `family`, and
getting it wrong was bug M; with one convention there is nothing left to get
wrong.

### Orchestration profiles

The variants also disagreed about which post-estimation tests run. Those
disagreements look like drift rather than design, so they are reproduced exactly
rather than harmonised — changing them would move published numbers.

| | `base` | `ci` | `missing` | `ci_mi` | `full` |
|---|---|---|---|---|---|
| causal forest variance | no | yes | yes | yes | yes |
| causal forest BLP/independence | yes | no | yes | no | yes |
| `dr_random_forest` BLP/independence | yes | no | yes | no | yes |
| oracle / semi-oracle tests | yes | no | yes | no | yes |
| SuperLearner arm | yes | no | yes (if `X` complete) | no | yes (if `X` complete) |
| half-sample bootstrap | no | yes | no | yes | with `ci` |
| nuisance row means | yes | no | yes | yes | |

`full` (2026-10-06) is `missing` under a name for applied use: tests, causal
forest variance and, given `ci`, the bootstrap intervals all from one run. The
"SuperLearner arm" and "bootstrap" rows are really set by `sl_lib` and `ci`,
which work under any profile; the CI profiles only switch the tests off.

The `dr_random_forest` row used to be the odd one: `missing` alone set
`dr_rf_tests = FALSE`, so `BLP_p` was `NA` for exactly one model in those two
studies. Nothing marked it as deliberate, and the decision has been taken —
every model carries the tests where possible.

The results made before that were back-filled in place by a one-off patch
rather than re-run. Those results are now archived, the re-run carries the
tests natively, and the patch scripts were removed after commit `e7b1d59` (see
`missing/binary/README.md` for the history).

`multiple_imputation` runs keep each imputation's tests unpooled, in `mi_tests`
(`mi_test_table()`). The pooling rule is still to be decided; see
`missing/README.md`.

The CI profiles keep the tests off; that one **is** deliberate.

## The parameter grid contract

Each study declares its grid **once**, in `<study>/<prefix>_config.R`. The PBS
array index is a **row number** of `study$grid`, so:

> Never filter or reorder the grid after construction.

Doing so renumbers every job, which is what made index `i` mean different things
in the analysis and check scripts (bug D). To run a subset, select indices and
leave the numbering alone:

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
continuous outcomes. The missing-data sets (`continuous_missing`,
`binary_missing`) have correlated covariates since 2026-09-28 and draw
`W, Z-block (X1–X5, X01–X03), [U], [err], cats` instead
(`correlated_covariates()`); every other set is unchanged.
`regression_check.R` fingerprints the generated dataset, not just the
estimates, so a change here fails loudly. Every scenario draws and returns
`X1`–`X5` since 2026-09-27, whether or not its treatment effect uses them;
before then `X3`/`X4`/`X5` were drawn only where needed, so earlier runs of
scenarios 1, 2, 4–8 and 10 are not reproducible from the current code.

Both outcomes are one model, `E[Y | x, W] = m0(x) + W·τ(x)` (`control_mean()`
gives `m0`). For a binary outcome `m0` is a logistic scaled into
`[p0_lo, p0_hi]`, so the treatment effect adds on the risk-difference scale;
`sample_size/binary/README.md` has the design and why.

## Fixed bugs

Legacy flags for bugs A, F, and the missing/binary fork have been deleted; the code now always uses the fixed behaviour.

- **bug K** — `pretest_superlearner()` could return an empty SuperLearner library, crashing downstream calls; fixed.
- **bug L** — `run_blp_whole()` had no `tryCatch`, crashed on degenerate CATEs; fixed.

## Verifying a change

```bash
Rscript R/regression_check.R baseline   # before touching anything
Rscript R/regression_check.R verify     # after each step - must be 8/8
Rscript crossfitting/cf_testing.R       # independent check of the estimators
```

The harness runs each study in its own subprocess, because `run_all_cate_methods`
is defined in seven shim files and sourcing two studies into one session would
silently give the second one's definition to both.

Scope: it proves the refactored code reproduces the **current** code on **this
machine**. It is not a claim about the cluster's numbers — R 4.5.3 here versus
4.3.2 there. Re-capture the baseline if you change machine.
