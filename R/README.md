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

## Fixed since the flags above were written

Not gated by a legacy flag - both were straight crashes with no "fixed
behaviour to fall back to" (nothing to reproduce), so nothing to flip back.
See the root `README.md` bug ledger for the full write-up (letters K, L).

- **bug K** — `pretest_superlearner()` could return an empty library
  (`character(0)`) when every candidate warned or errored on a fold; its
  callers (`nuisance_sl()`, `stage_2_sl()`, and `crossfitting/cf_models.R`'s
  own `sl_nuisance_fit()`/`stage2_crossfit_sl()`) then passed that straight
  into a live `SuperLearner(..., SL.library = character(0))`, which crashed
  building a 0-column predictions `data.frame()`. Fixed by giving
  `pretest_superlearner()` a floor: falls back to `"SL.mean"` (asserted
  directly, not re-run through the pretest loop) if nothing else survives.
  Confirmed via `missing/binary/bin_miss_analysis.R` (see its README) and
  `continuous/cts_analysis.R` (23 of the 24 array IDs that failed on the
  cluster for this reason). Distinct from bug F, which is about *which*
  library variant is used once pretesting leaves at least one survivor, not
  about what happens when it leaves zero.
- **bug L** — `run_blp_whole()` called `GenericML::BLP()` with no
  `tryCatch`, unlike `run_independence_test_whole()` next to it. `BLP()`
  crashes (`subscript out of bounds`) whenever the fitted CATE (`tau`) is
  exactly constant: the `beta.2` regressor it builds
  (`(W - W.hat) * (tau - mean(tau))`) becomes identically zero, `lm()` marks
  it aliased, and `sandwich::vcovHC()` drops that coefficient's row/column
  entirely rather than keeping it as `NA`, so `GenericML`'s internal
  name-based indexing throws. Fixed by wrapping the `BLP()` call in
  `tryCatch`, returning `NULL` on error — `R/metrics.R`'s
  `hte_test_metrics()` already maps `BLP_whole = NULL` to `BLP_p = NA`, and
  `NA` is the statistically correct answer here (there's no fitted
  coefficient to attach a p-value to when `tau` has zero variance).
  Confirmed via `continuous/` array id 3529 (scenario 9, n=100, run=89).

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
