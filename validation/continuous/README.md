# Interim-analysis validation — continuous outcome

The only validation study run today. See `validation/README.md` for what a
validation study checks and why.

## Design

| | |
|---|---|
| scenario | 3 (continuous × binary interaction — see `R/dgm_scenarios.R`, `DESC_10`) |
| n | 1000, split into two chunks of `n * interim_prop` and `n * (1 - interim_prop)` |
| interim_prop | 0.25 to 0.75 in steps of 0.05 — 11 interim points |
| runs | 100 — **1100 array jobs**, 1h walltime |
| folds | 5 below n=250 per chunk, else 10 |
| results | `../results/cts_val/scenario_3/1000/<interim_prop>/res_sim_<run>.RDS` |

`cts_val_config.R` rounds the `interim_prop` sequence to 2 dp, which is
load-bearing rather than cosmetic: `interim_prop` is a `path_cols` entry, so its
`as.character()` form becomes a results directory name, and an unrounded
`seq(by = 0.05)` yields values like `0.30000000000000004` that `get_results()`
can no longer match back to their grid row.

Only scenario 3 and n = 1000 are run today; both are still columns of the grid
(`cts_val_config.R`) so extending to more scenarios or sizes is a one-line
change, not a rewrite.

## Estimators

Only `causal_forest` and `dr_random_forest`, both from `R/cate_models.R` — the
question here is whether one estimator's own findings persist across chunks,
not which estimator is most accurate, so the oracle/semi-oracle/SuperLearner
arms this repo's other studies carry are not run.

## Variable importance

Each estimator's per-run object carries **two** importance measures, one column
per covariate each. They ask different questions and need not agree, which is
why both are kept — each gets its own chunk-1-vs-chunk-2 rank comparison, and
each nominates its own top covariate for the interaction test below.

- `te_vims` — refit the second stage with each covariate dropped in turn and
  compare out-of-bag prediction error. How much worse is the CATE predicted
  without this covariate?
- `shap_vims` — fit a tuned xgboost surrogate to the estimated CATEs and take
  exact TreeSHAP of it, averaging |SHAP| over units. How much of each unit's
  estimated CATE is attributable to this covariate? This is "Strategy 3" /
  indirect SHAP from Svensson et al.'s SHAP_CATE, adapted from that paper's
  `cvboost3()` with a reduced tuning grid (3 learning rates × 3 depths, 5-fold,
  ≤2000 rounds) — the paper's 24-combination search runs four times per array
  job across 1100 jobs. Because it reads only `(X, tau)` it applies unchanged to
  both estimators, unlike the TE-VIMs which need an estimator-specific refit.

Both are larger-is-more-important, so `rank()` means the same thing for each.

**New dependencies**: `xgboost` and `SHAPforxgboost`, used nowhere else in this
repo. `renv` is not wired up here (see `.claude/CLAUDE.md`), so they must be
installed into the cluster's `sim-env` R library by hand — `cts_val_testing.R`
check 1 is there to catch their absence before 1100 jobs go out.

## What is compared across chunks

`cts_val_analysis.R` saves four comparisons per run, under `validations`:

| | |
|---|---|
| `subgroups` | fit a tree on chunk 1's top/bottom-10%-CATE responder groups, predict them into chunk 2, test the `W x subgroup` interaction there |
| `variances` | `var(tau)` in each chunk, against the true `Var(tau) = 1` |
| `var_imps` | each measure's covariate ranking in each chunk, and the rank change |
| `top_var_tests` | chunk 1's most important covariate, interaction-tested on the remaining participants — both as a continuous `W x X_top` interaction and after a median split |

The estimation logic behind these (both importance measures and
`interaction_pval()`) is specific to this study — it isn't duplicated anywhere
else in the repo, so it stays in `cts_val_models.R` rather than moving into `R/`.

## Files

| file | role |
|---|---|
| `cts_val_config.R` | the parameter grid and results path — **the** definition |
| `cts_val_dgms.R` | names this study's slice of `R/dgm_scenarios.R` |
| `cts_val_models.R` | wraps `R/cate_models.R`'s causal_forest/DR-RF + the TE-VIM, TreeSHAP and interaction-test helpers |
| `cts_val_analysis.R` | array entry point; fits both chunks, computes the four chunk comparisons |
| `cts_val_run.R` | runs the whole grid in one RStudio session, 8 rows at a time — the no-queue alternative to `cts_val_1.sh` |
| `cts_val_testing.R` | pre-submission verification — dependencies, grid, and the helpers above |
| `cts_val_profile.R` | resource sweep over `(interim_prop, workers, grf_threads)` — run in an RStudio session, not a job |
| `cts_val_profile_summary.R` | turns the sweep into PBS directives and writes them into the jobscripts |
| `cts_val_check.R` | finds missing runs, writes `jobscripts/failed_ids.txt`, and updates `-J` and the resource request in the rerun jobscript |
| `cts_val_collect.R` | gathers per-run files into `cts_val_all.RDS` |
| `cts_val_metrics.R` | flattens `validations` into tidy `cts_val_metrics.RDS` |
| `results_cts_val.R` | summary plots |
| `cts_val_results.qmd` | the written-up report |

## Running it

```bash
Rscript validation/continuous/cts_val_testing.R full   # do this first - see below
qsub validation/continuous/jobscripts/cts_val_1.sh     # 1-1100
Rscript validation/continuous/cts_val_check.R          # writes failed_ids.txt if any are missing
qsub validation/continuous/jobscripts/cts_val_collect.sh
qsub validation/continuous/jobscripts/cts_val_metrics.sh
```

`cts_val_testing.R` runs six cheap structural checks; `full` adds one real
replicate end to end and reports how long it took, which is the number to read
against the jobscript's 1h walltime before submitting. It exits non-zero on any
failure. Drop `full` for a quick local smoke test.

### Without the queue

When the queue is slower than the work, `cts_val_run.R` runs the same 1100 rows
inside one interactive session instead. Request an RStudio session with 8 cores
and 64gb — the same session `cts_val_profile.R` asks for — then:

```r
source(here::here("validation", "continuous", "cts_val_run.R"))
```

It runs 8 rows at a time, each as its own `Rscript cts_val_analysis.R <i> 1 1`
subprocess — the identical command `cts_val_1.sh` gives one array index, so a
row produced here and a row produced by the array are the same calculation and
land in the same place. Budget roughly 7h for a full grid.

Rows that already have a results file are skipped (the same missing-run logic
`cts_val_check.R` uses), so an interrupted session — or one that hits its
walltime — is resumed by sourcing the file again, and so is a row that failed.
Edit `ids` at the top to run a subset, e.g.
`grid_indices(study, interim_prop = 0.25)`; trailing arguments do the same
non-interactively (`Rscript cts_val_run.R 1 2 3`). Per-row logs go to
`<results>/validation/continuous/session_logs/`.

It never writes `jobscripts/failed_ids.txt` and never edits `cts_val_rerun.sh` —
those belong to the queue path. Afterwards, `cts_val_check.R`, collect and
metrics run exactly as they would after an array run.

### Sizing the job

`cts_val_1.sh`'s `ncpus=5` only ever mirrored a hardcoded `workers <- 5`; it was
never measured, and until recently no `num.threads` reached grf at all, so each
of those five workers spawned forests on every visible core whatever PBS had
allocated. `cts_val_profile.R` measures the alternatives. Unlike
`continuous/cts_profile.R` it is **not** an array job — there is no
`cts_val_profile.sh`. Ask for an interactive RStudio session with 8 cores and
64gb, then:

```r
source(here::here("validation", "continuous", "cts_val_profile.R"))
source(here::here("validation", "continuous", "cts_val_profile_summary.R"))
```

The sweep is 28 cells over `(interim_prop, run, workers, grf_threads)` and writes
one `prof_<i>.RDS` per cell as it goes, skipping any already on disk — an
interrupted session is resumed by sourcing it again. Set `cells <- 1` at the top
for a single trial cell first, which also gives a per-cell time to budget the
rest against. The summary prints the observed-against-allocated CPU table (the
direct read on whether `ncpus=5` was ever honoured) and a cost-against-`ncpus`
table, then **edits `cts_val_1.sh` and `cts_val_rerun.sh` in place** — check
`git diff` before submitting.

`workers` and `grf_threads` reach `cts_val_analysis.R` as trailing arguments on
the `Rscript` line, so the PBS resource request and the R-level parallelism
cannot drift apart. Called without them it falls back to what it did before
(`workers = 5`, grf unthrottled).

## Status

**Needs a fresh run** — for several independent reasons. The most recent set:

- **`bottom_pval` was never a subgroup test.** The analysis script pulled the
  interaction p-value correctly for the top-10% group (`pvals_top[4]`) but read
  `pvals_bottom[1]` — the *intercept* — for the bottom-10% one. Every
  `bottom_pval` in results generated before this fix is meaningless, and the
  report excluded the column entirely. `interaction_pval()` now indexes the
  `W:v` coefficient by name and returns `NA` when the subgroup carries no
  contrast, so both arms are real tests and neither can silently read the wrong
  row again.
- **The interim grid is finer.** `interim_prop` was three points
  (0.25/0.5/0.75); it is now eleven, 0.25 to 0.75 in 0.05 steps, so replication
  rates read as a curve in interim size. 300 array jobs became 1100.
- **A second importance measure.** Surrogate TreeSHAP now runs alongside the
  TE-VIMs (see Variable importance above), and every var-imp row carries a
  `measure` column that older results do not have.
- **A fourth chunk comparison.** `top_var_tests` — see What is compared across
  chunks.

Before that, the crossfitting strategy change to `R/cate_models.R` (see root
README Methods/Status) affected both estimators used here, and this folder was
rebuilt from a pre-restructure fork (standalone copies of the DGM and the CATE
estimators, grid/index arithmetic done by hand in `validation.sh`, hardcoded
absolute cluster paths) onto the shared `R/` pattern every other study now
uses. Three things changed on purpose, not by accident, as part of that
rebuild:

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
`R/metrics.R::hte_test_metrics()` consumes, which is what a fifth chunk
comparison (alongside subgroups/variance/var-imps/top-var) would build on.

**Not in scope for this refactor.** `R/regression_check.R`'s old-vs-new
harness compares a behaviour-*preserving* refactor against a baseline; the
causal-forest nuisance change above means there is no valid "before" to check
against here. Worth adding once a fresh baseline exists.
