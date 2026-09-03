# `R/` — the shared library

Every study sources from here. Nothing in this directory knows about a specific
study: the differences between them are arguments, not forks.

This folder exists because they used to be forks. The same CATE estimation code
lived in seven files, the same DGM in four, and the same collect/check
boilerplate in eight — `continuous/cts_models.R` and `binary/bin_models.R`
differed in **two** places out of 438 lines. Consolidating removed roughly 3,000
lines and, more usefully, removed the possibility of the copies drifting apart,
which is where most of the bugs in the ledger below came from.

## Files

| file | holds |
|---|---|
| `utils.R` | `setup_rng_stream`, `collate_predictions`, `scatter_folds`, `trim_ps`, `timed` |
| `dgm_scenarios.R` | scenario tables + `generate_scenario_data()`, `get_oracle_info()` |
| `missingness.R` | `introduce_missingness`, `handle_missingness`, `generate_and_process_data` |
| `cate_models.R` | the DR-learner family, causal forest, `cate_methods()`, `combine_mi` |
| `bootstrap_ci.R` | `cf_half_boot`, `rf_half_boot`, `combine_mi_ci`, `find_optimal_sf` |
| `metrics.R` | `cate_metrics`, `interval_metrics`, `compute_metrics` |
| `pipeline.R` | `study_config`, `get_results`, `check_failed`, `grid_indices` |
| `figures.R` | display labels, palette/theme, `summarise_metrics`, plot helpers |
| `regression_check.R` | the old-vs-new equivalence harness |

The repo-root `utils.R` is a two-line shim onto `R/utils.R`, so existing
`source(here("utils.R"))` calls keep working.

## `cate_methods()` — the four axes

The seven forked model files differed on exactly four things, which are now
arguments:

| argument | what it changes | who uses it |
|---|---|---|
| `family` | SuperLearner outcome model: `binomial()` adds `method.NNloglik` | binary studies |
| `oracle_link` | whether `plogis` is applied to the oracle formula **here** | see below |
| `ipw` | grf `sample.weights` / SuperLearner `obsWeights`. `NULL` is the unweighted path exactly | the missing-data IPW arm |
| `ci` | `list(boot=, sf=, alpha=)` turns on the half-sample bootstrap | the CI studies |

**`oracle_link` is not the same as `family`, and conflating them is a bug.** The
repo has two conventions for where the inverse link lives:

- `binary/` and `confidence_intervals/binary/` — `get_binary_oracle_info()`
  returns a **linear predictor**, and the model code applies `plogis`.
  `oracle_link = "logit"`.
- `missing/binary/` — the oracle formula string already contains `plogis(...)`.
  `oracle_link = "identity"`.

Both give the same answer. Applying `plogis` to a formula that already contains
it does not.

### Orchestration profiles

The variants also disagreed about which post-estimation tests run. Those
disagreements look like drift rather than design, so they are reproduced exactly
rather than harmonised — changing them would move published numbers.

| | `base` | `ci` | `missing` | `ci_mi` |
|---|---|---|---|---|
| causal forest variance | no | yes | yes | yes |
| causal forest BLP/independence | yes | no | yes | no |
| `dr_random_forest` BLP/independence | yes | no | yes | no |
| oracle / semi-oracle tests | yes | no | yes | no |
| SuperLearner arm | yes | no | yes (if `X` complete) | no |
| half-sample bootstrap | no | yes | no | yes |
| nuisance row means | yes | no | yes | yes |

The `dr_random_forest` row used to be the odd one: `missing` alone set
`dr_rf_tests = FALSE`, so `BLP_p` was `NA` for exactly one model in those two
studies. Nothing marked it as deliberate, and the decision has been taken —
every model carries the tests where possible.

The 19,800 results already on disk were **not** re-run for it. Both tests are
deterministic (`GenericML::BLP` is an OLS with a sandwich vcov;
`coin::independence_test` is asymptotic under `teststat = "quadratic"`), and the
saved files retain `nuisances_rf`, `data` and `tau`, so `R/patch_hte_tests.R`
recomputed the three fields in place. `missing/patch_hte_verify.R` is the proof:
it runs the same row twice from one seed, with the flag off and on, and checks
that every arm's `tau` is byte-identical and that the patch reproduces the
re-run's test values exactly. One field is not recoverable —
`dr_random_forest$variance`, which the old inline branch discarded — so patched
files lack it and newly run ones have it. Nothing reads it.

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
W, X1, X2, [X3], [X4], [X5], [U], [err], X01, X02, X03, cats
```

`X3`/`X4`/`X5` only when the scenario needs them, `U` only for the MNAR
mechanisms, `err` only for continuous outcomes. `regression_check.R`
fingerprints the generated dataset, not just the estimates, so a change here
fails loudly.

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
