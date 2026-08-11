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
| `dr_random_forest` BLP/independence | **yes** | no | **no** | no |
| oracle / semi-oracle tests | yes | no | yes | no |
| SuperLearner arm | yes | no | yes (if `X` complete) | no |
| half-sample bootstrap | no | yes | no | yes |
| nuisance row means | yes | no | yes | yes |

The `dr_random_forest` row is the odd one: the base studies get BLP and
independence tests for that model and the missing-data studies do not, so
`BLP_p` is `NA` for exactly one model in those studies. Nothing marks it as
deliberate. **Worth a decision.**

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

## Legacy flags — to be deleted once the re-runs land

These reproduce known bugs. They currently default to the **fixed** behaviour;
flipping one back reproduces exactly what produced the results now on the
cluster, which is useful while the corrected runs are in flight.

| flag | file | reproduces |
|---|---|---|
| `LEGACY_BIN_CI_PARAMS` | `dgm_scenarios.R` | bug A — the binary CI study ran on the continuous coefficient table |
| `LEGACY_BIN_MISS` | `dgm_scenarios.R` | the three-defect `missing/binary` fork (coefficients, `power.t.test` calibration, log-odds truth) |
| `PRETEST_STAGE2` | `cate_models.R` | bug F — `stage_2_sl` discarding the pretested SuperLearner library |
| `BIAS_SIGN` | `metrics.R` | bug G — `bias` as `true - est` |

## Known but unfixed

Not gated by a legacy flag - there's no fixed behaviour to fall back to yet,
just an open crash.

**`pretest_superlearner()` can return an empty library, and its callers don't
guard against that.** `pretest_superlearner()` (line 366) drops any candidate
algorithm that errors, warns, or returns all-`NA` on a 2-fold inner CV; if
every candidate fails on a given fold it returns `character(0)`.
`nuisance_sl()` (line 306) and `stage_2_sl()` (line 408) both pass that result
straight into `SuperLearner(..., SL.library = <possibly empty>)`, which then
crashes (`arguments imply differing number of rows: 0, 1`, inside
`future_map`'s per-fold apply). Confirmed reproducing on unmodified
`missing/binary/bin_miss_analysis.R` at `scenario = 1, mechanism = MAR, run =
1` for both `complete_cases` (n = 331) and `mean_imputation` (n = 500) - see
"Known issue found while profiling" in `missing/binary/README.md` for the full
repro and a walk-through of the root cause. Not yet checked in the other
studies that call `pretest_superlearner()` (`binary/`, `continuous/`,
`crossfitting/`, `case_study/`). Distinct from bug F, which is about which
library variant is used once pretesting leaves at least one survivor, not
about what happens when it leaves zero - candidate "bug H" once it's fixed
(bug G is already taken).

**When the corrected results are in and written up, remove:**

1. all four flags and their `if` branches
2. `SCENARIO_SETS$binary_missing` (the legacy table) and the
   `binary_missing_fixed` indirection in `resolve_set()`
3. the `LEGACY_BIN_MISS` branch in `calibration_for()`
4. `stage_2_sl`'s dead vector branch — every caller passes a matrix, so the
   `single` branch has never executed in any study
5. the `bias_sign` argument of `cate_metrics()`

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
