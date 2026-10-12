# Interim-analysis validation — continuous outcome

See `validation/README.md` for what a validation study checks and why. The
binary arm (`../binary/`) runs the same design, and the code both arms share
lives in `../val_common.R`.

## Design

| | |
|---|---|
| covariates | correlated, rho = 0.5 — the `continuous_corr_0.5` set (`R/dgm_scenarios.R`, CORRELATED COVARIATES), the same set `sample_size/correlated/` runs as its primary arm |
| scenarios | 1, 2, 3, 4 (`R/dgm_scenarios.R`, `DESC_10`): 1 the null, no HTE (since 2026-10-12); 2 simple HTE, `b4 * X4`; 3 complex, `b3 * X3 + b4 * X4 + b45 * X4 * X5`; 4 non-linear, `b4 * cos(X4)`. Scenario 2 was scenario 3 before the 2026-09-26 renumbering |
| n | 1000, **one** trial split into its first `n * interim_prop` participants and the rest |
| interim_prop | 0.25 to 0.75 in steps of 0.05 — 11 interim points |
| runs | 100 per scenario — **4400 array jobs**, 1 h walltime (see Sizing the job) |
| folds | DR SuperLearner only: 5 for a chunk below 500 rows, else 10 |
| results | `../results/validation/continuous/rho_0.5/scenario_<k>/1000/<interim_prop>/res_sim_<run>.RDS` |

Each run draws one trial of `n` and splits it. Because the seed is the run alone,
a run's chunk 1 at 0.25 is nested in its chunk 1 at 0.30, so every curve across
`interim_prop` is paired within run. Before the 2026-10-08 changes the two chunks
were two separate draws, and since `generate_scenario_data()` calibrates `bW`
to 80% power at the `n` it is given, each chunk was its own trial with its own
ATE.

`cts_val_config.R` rounds the `interim_prop` sequence to 2 dp, which is
load-bearing rather than cosmetic: `interim_prop` is a `path_cols` entry, so its
`as.character()` form becomes a results directory name, and an unrounded
`seq(by = 0.05)` yields values like `0.30000000000000004` that `get_results()`
can no longer match back to their grid row.

Only rho = 0.5, scenarios 1–4 and n = 1000 are run today; all three are
columns of the grid (`cts_val_config.R`) so extending to rho = 0, more scenarios
or sizes is a one-line change, not a rewrite. (The grid briefly carried rho = 0
too, on 2026-10-09. Its `rho_0/` results are outside the grid now, so check
and collect ignore them.) The same run number draws the same participants in
every scenario; only the CATE changes.

## Estimators

`causal_forest`, `dr_random_forest` and `dr_superlearner`, all from
`R/cate_models.R`. The question here is whether one estimator's own findings
persist across chunks, not which estimator is most accurate, so the oracle and
semi-oracle arms this repo's other studies carry are not run.

The DR SuperLearner is fit as in the sample-size studies: `nuisance_sl()` and
`run_dr_superlearner()`, single leave-one-fold-out crossfitting with the same
folds for nuisances and stage 2, and `sl_libraries()` sized to each chunk. Every
chunk is ≥ 250 rows, so it gets the full libraries apart from `SL.earth` in the
outcome model below 500. Its pseudo-outcome comes from its own SuperLearner
nuisances, and its TE-VIMs are scored against that rather than against
`nuisance_rf()`'s.

## Variable importance

Each estimator's per-run object carries **two** importance measures, one column
per covariate each. They ask different questions and need not agree, which is
why both are kept — each gets its own chunk-1-vs-chunk-2 rank comparison, and
each nominates its own top covariate for the interaction test below.

- `te_vims` — refit the second stage with each covariate dropped in turn and
  compare out-of-sample prediction error. How much worse is the CATE predicted
  without this covariate? Out-of-bag for the two forests. For the DR
  SuperLearner it is a crossfit refit of its stage-2 SuperLearner on the same
  folds (`get_te_vims_superlearner()`), one per covariate × fold. Each fold
  reuses the library the full fit's pretest left there instead of pretesting
  again. This is by far the most expensive step in a row.
- `shap_vims` — fit a tuned xgboost surrogate to the estimated CATEs and take
  exact TreeSHAP of it, averaging |SHAP| over units. How much of each unit's
  estimated CATE is attributable to this covariate? This is "Strategy 3" /
  indirect SHAP from Svensson et al.'s SHAP_CATE, adapted from that paper's
  `cvboost3()` with a reduced tuning grid (3 learning rates × 3 depths, 5-fold,
  ≤2000 rounds) — the paper's 24-combination search runs four times per array
  job across 4400 jobs. Because it reads only `(X, tau)` it applies unchanged to
  every estimator, unlike the TE-VIMs which need an estimator-specific refit.

Both are larger-is-more-important, so `rank()` means the same thing for each.

**New dependencies**: `xgboost` and `SHAPforxgboost`, used nowhere else in this
repo. They are in `environment.yml`; a fresh install of `sim-env` must include them, and
they must otherwise be installed into the cluster's `sim-env` R library by hand — `cts_val_testing.R`
check 1 is there to catch their absence before 4400 jobs go out.

## What is compared across chunks

`cts_val_analysis.R` saves four comparisons per run, under `validations`:

| | |
|---|---|
| `subgroups` | fit a tree on chunk 1's top/bottom-10%-CATE responder groups, predict them into chunk 2, test the `W x subgroup` interaction there |
| `variances` | `var(tau)` in each chunk, against the true `Var(tau)` (1 in scenario 2; the report computes it per scenario) |
| `var_imps` | each measure's covariate ranking in each chunk, and the rank change |
| `top_var_tests` | chunk 1's most important covariate, interaction-tested on the remaining participants — as a continuous `W x X_top` interaction (`p_cts`), the same adjusted for every other covariate (`p_cts_adj`), after a median split (`p_split`), after a split learned on chunk 1 (`p_dd`, with `dd_cuts` and `dd_prop`), and as a spline projection learned on chunk 1 (`p_proj`) |

Since 2026-10-11 every interaction test (the subgroup tests and all five
top-covariate forms) is saved twice, with classical and with HC3 standard
errors: `subgroups[[model]]` is a two-row `data.frame(se, top, bottom)` and
`top_var_tests[[model]]` leads with an `se` column (`chunk_validations_se()`).
`variances` and `var_imps` do not depend on the standard errors and are saved
once.

**Estimates and the direction of a replication** (since 2026-10-12). A
p-value alone cannot say whether chunk 2 found the effect chunk 1 claimed, or
one pointing the other way. Each test therefore saves four more numbers next to
its p-value, named `<stat>_<test>`. The tests are `top`/`bottom` for the
subgroups and `cts`, `cts_adj`, `split`, `dd`, `proj` for the top-covariate
forms. They come after the existing columns, whose names and values are
unchanged:

| column | what it is |
|---|---|
| `est_<test>` | the chunk-2 `W:v` estimate (the same under both kinds of standard error) |
| `se_<test>` | its standard error (classical or HC3, per row) |
| `ref1_<test>` | the same contrast in chunk 1's *estimated* tau: what chunk 1 claimed |
| `true2_<test>` | the same contrast in chunk 2's *true* tau: what the test estimates; 0 in scenario 1 |

W is randomised, so `W:v` estimates the slope of the linear projection of tau
on v, and `ref1`/`true2` are that projection with tau-hat or true tau as the
response (`tau_slope()`, `tau_slope_adj()`). For a subgroup, that is the
difference in mean tau inside vs outside. `ref1` uses the form's v built in
chunk 1 with the same helper (the tree's own groups, chunk 1's median split,
the tree and spline evaluated at chunk 1's `X_top`). So `ref1_top > 0`,
`ref1_bottom < 0` and `ref1_dd > 0` by construction. A replication "in the
right direction" has `sign(est) == sign(ref1)`, and also `== sign(true2)` if
it is right in truth. The report's "Direction of the replications" section
summarises both. `true2` needs `results2$truth`, which the analysis scripts
attach; `chunk_validations()` called without it gives NA.

**Scenario 1, the null** (since 2026-10-12). tau = bW for everyone, so every
significant test is a type I error and every replication rate is a size. That
is the reference the other scenarios' rates are read against.

`p_cts_adj` exists because of the correlation. With rho = 0.5, E[tau | X5]
moves with X5 because X5 is correlated with X4, so a marginal `W x X5` test
finds a real interaction even though X5 modifies nothing. A run that nominates
a proxy of X4 then "replicates" on `p_cts`. `p_cts_adj` is the `W:X_top`
coefficient of `Y ~ W * (all covariates)` and should not flag it.

**The data-driven forms** let chunk 1 decide the shape of the test; chunk 2's
outcomes are seen only by the test, so both stay honest.

- `p_dd`: a depth-2 regression tree of chunk 1's estimated CATE on `X_top`
  alone (`cp = 0`, leaves ≥ 10% of the chunk), its leaves grouped into above /
  below chunk 1's mean CATE. The split can be one threshold, a binary
  covariate's 0/1, or middle vs. tails. `dd_cuts` records the cut points,
  `dd_prop` chunk 2's share on the high side. A flat chunk 1 CATE gives no
  split and an NA.
- `p_proj`: a smoothing spline (4 df) of chunk 1's estimated CATE on `X_top`,
  evaluated at chunk 2's `X_top`, tested as a 1-df `W x g(X_top)` interaction. A
  covariate with ≤ 4 distinct values (X1, X3) is scored as itself, which is then
  the `p_cts` test.

**Scenario 4 is symmetric in X4.** `cos` is even and X4 is centred on 0, so the
best linear fit of tau on X4 has slope 0 and E[tau | X4 > median] = E[tau | X4 <
median]. `p_cts`, `p_cts_adj` and `p_split` therefore have no power there even
when X4 is correctly nominated; their rate should sit near 5%. This is kept on
purpose and reported as a finding about those tests. `p_dd` and `p_proj` are the
forms that can replicate it.

**`p_split` on a binary covariate** (fixed 2026-10-11). The split was
`xt > median(xt)`, which is all 0 for X3 ~ Bernoulli(0.7) (median 1), so
`p_split` was NA whenever X3 was nominated — the likely winner in scenario 3,
where X3 is the largest single source of heterogeneity. `median_split()` now uses a two-valued
covariate as its own split. Nothing changes for continuous covariates, or for
X1 (P = 0.4, median 0).

The comparisons are `chunk_validations()`, and the estimation logic behind them
(both importance measures and `interaction_pval()`), in `../val_common.R`,
shared with the binary arm. It is specific to these studies, so it stays in
`validation/` rather than moving into `R/`.

**HC3 interaction tests.** This arm calls `chunk_validations_se()`, so every
interaction test (the subgroup tests and the five top-covariate forms) is
saved as a t-test on classical standard errors and on
`sandwich::vcovHC(type = "HC3")` ones, as in the binary arm. HC3 is the
primary kind; the report compares the two. The noise is homoskedastic, but `Y ~ W * v` leaves out the
prognostic X1 and X2 and, in the treated arm, the spread of tau within each
level of `v`. At rho = 0.5 both of those vary with the subgroup, so the
residual variance differs across the four `W x v` cells. The W:v estimate is a
contrast of the four cell means, its variance is dominated by the small cells
(a 10% responder subgroup is ~12–40 per cell in chunk 2), and the classical
pooled variance comes almost entirely from the large ones. HC3 is the cell-wise
(Welch-type) variance there, slightly conservative. For `p_split`, `p_cts` and
`p_cts_adj` the two give nearly the same answer, so the switch costs little
there and keeps the two arms' tests the same. Classical standard errors until
2026-10-08, HC3 only until 2026-10-11.

## Files

| file | role |
|---|---|
| `cts_val_config.R` | the parameter grid and results path — **the** definition |
| `cts_val_dgms.R` | names this study's slice of `R/dgm_scenarios.R` |
| `cts_val_models.R` | `run_all_cate_methods()`: `../val_common.R`'s `fit_val_methods()` with `family = gaussian()` |
| `cts_val_analysis.R` | array entry point; splits the trial, fits both chunks, runs `chunk_validations_se()` |
| `../val_common.R` | shared with `binary/`: the split, the estimators and importance measures, the interaction tests (including the data-driven ones), the four comparisons |
| `cts_val_run.R` | runs the whole grid in one RStudio session, 8 rows at a time — the no-queue alternative to `cts_val_1.sh` |
| `cts_val_testing.R` | pre-submission verification — dependencies, grid, and the helpers above |
| `cts_val_check.R` | finds missing runs, writes `jobscripts/failed_ids.txt`, and updates `-J` and the resource request in the rerun jobscript |
| `cts_val_collect.R` | gathers per-run files into `cts_val_all.RDS` |
| `cts_val_metrics.R` | flattens `validations` into tidy `cts_val_metrics.RDS` |
| `results_cts_val.R` | summary plots, one scenario at a time (`Rscript results_cts_val.R <scenario>`), HC3 |
| `cts_val_results.qmd` | the written-up report, one render per scenario (Quarto `scenario` param) |

## Running it

```bash
Rscript validation/continuous/cts_val_testing.R full   # do this first - see below
qsub validation/continuous/jobscripts/cts_val_1.sh     # 1-4400
Rscript validation/continuous/cts_val_check.R          # writes failed_ids.txt if any are missing
qsub validation/continuous/jobscripts/cts_val_collect.sh
qsub validation/continuous/jobscripts/cts_val_metrics.sh
cd validation/continuous
quarto render cts_val_results.qmd -P scenario:1 --output cts_val_results_sc1.html   # and 2, 3, 4
```

`cts_val_testing.R` runs eight cheap checks; `full` adds one real
replicate end to end and reports how long it took, which is the number to read
against the jobscript's walltime before submitting. It exits non-zero on any
failure. Drop `full` for a quick local smoke test.

**Before the 2026-10-12 run, archive every old result, both arms.** Nothing
saved so far was made by the current code. Files from before 2026-10-11 carry
HC3 p-values only and no `p_dd`/`p_proj`. None carry the interaction estimates.
And adding scenario 1 reordered the grid. `check_failed()` would count old
files as done, and `get_results()` would collect them next to the new ones. So
archive both arms' whole trees: every `rho_*/scenario_*/`, the old `rho_0/`,
the collected `*_val_all.RDS` / `*_val_metrics.RDS` and `session_logs/`.
`archive_old_results.R --trees` writes one `.tar` per tree under
`results/_archive/<label>/`, and deletes a tree only once its tar lists exactly
its files. From the repo root, dry run then `--apply`:

```bash
Rscript R/archive_old_results.R --label validation_pre_2026-10-12 \
  --trees validation/continuous validation/binary
Rscript R/archive_old_results.R --label validation_pre_2026-10-12 \
  --trees validation/continuous validation/binary --apply
```

or, if the trees are big, as a job. The jobscript adds `--apply` itself, and
`ARCHIVE_ARGS` must stay comma-free:

```bash
cd R/jobscripts
qsub -v ARCHIVE_ARGS="--label validation_pre_2026-10-12 --trees validation/continuous validation/binary" archive_old_results.sh
```

To restore a tree, run `tar -xf _archive/validation_pre_2026-10-12/<name>.tar`
from `results/`.

### Without the queue

When the queue is slower than the work, `cts_val_run.R` runs the same 4400 rows
inside one interactive session instead. Request an RStudio session with 8 cores
and 64gb, then:

```r
source(here::here("validation", "continuous", "cts_val_run.R"))
```

It runs 8 rows at a time, each as its own `Rscript cts_val_analysis.R <i> 1 1`
subprocess — the identical command `cts_val_1.sh` gives one array index, so a
row produced here and a row produced by the array are the same calculation and
land in the same place. A full grid took roughly 7h before the DR SuperLearner
was added; with it, budget 4400 × (the replicate time `cts_val_testing.R full`
reports) / 8.

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
allocated. `cts_val_1.sh` now runs one worker with one grf thread
(`cts_val_analysis.R <i> 1 1`) on `ncpus=1`. That is set by hand, as is the
walltime. A row took roughly 3 min with the two forest estimators. The DR
SuperLearner's TE-VIM refits (10 covariates × 5–10 folds of a 9-learner
SuperLearner, per chunk) make it longer; the walltime was raised to 20 min,
then to 1 h on 2026-10-09.
`cts_val_testing.R full` passes the
same `1 1`, so the replicate time it reports is the array job's. Check 6 also
prints a step-by-step timing of one 250-row chunk. If a row is too slow on one
core, raise `ncpus`, `ompthreads` and the `workers` argument together
(`cts_val_analysis.R <i> k 1` on `ncpus=k:ompthreads=k`). The SuperLearner
folds and the TE-VIM refits all run through `future_map()`, so they spread
across the extra workers. The `syrup`
resource sweep meant to measure the alternatives didn't work for this study
(see `docs/DEVLOG.md`'s "Resource profiling (removed)"). Read the replicate time
that `cts_val_testing.R full` reports against the jobscript's walltime instead.

`workers` and `grf_threads` reach `cts_val_analysis.R` as trailing arguments on
the `Rscript` line, so the PBS resource request and the R-level parallelism
cannot drift apart. Called without them it falls back to what it did before
(`workers = 5`, grf unthrottled).

## Status

**Needs a fresh run** of all 4400 rows. Archive both arms' whole results trees
first (Running it, above).

Changed 2026-10-12:

- **Scenario 1, the null**, added to the grid: 4 x 11 x 100 = 4400 array jobs
  (was 3300). The grid's row order changed with it, and row 1 is now scenario
  1 at 0.25, run 1. Every significant test there is a type I error.
- **Interaction estimates and reference contrasts**: `est_`, `se_`, `ref1_`,
  `true2_<test>` for every subgroup and top-covariate test — see What is
  compared across chunks. The p-values are unchanged and keep their names
  (checked against the previous `chunk_validations()` on the same fits). The
  report gains a "Direction of the replications" section, and handles a
  scenario without modifiers.

Changed 2026-10-11:

- **Scenarios 2–4**, rho = 0.5 only. The grid had carried rho = 0 as well
  since 2026-10-09; it is back to the primary arm, with scenarios 3 (complex)
  and 4 (non-linear) added. 3300 array jobs.
- **Classical and HC3 standard errors**, both saved
  (`chunk_validations_se()`) — see What is compared across chunks. The
  subgroup and top-covariate tables of a results file change shape (an `se`
  column), and `cts_val_metrics.R` reads the new shape only.
- **Two data-driven top-covariate tests**, `p_dd` and `p_proj` — see What is
  compared across chunks. Scenario 4's symmetric `cos(X4)` is where they
  matter.
- **`p_split` on a binary covariate** — see What is compared across chunks.
- **The report renders per scenario** (`-P scenario:<k>`), with the
  `Var(tau)` reference line and the true-modifier set derived from the
  scenario rather than hard-coded to scenario 2, and a classical-vs-HC3
  section.

Before that, the run needed for several independent reasons: bug O changed the
continuous DGM, and the scenario is now numbered 2, not 3 (root `README.md`,
Status, step 0).

Changed 2026-10-08, ahead of the rerun:

- **Correlated covariates, rho = 0.5 only.** The study now draws from
  `continuous_corr_0.5` rather than the independent `continuous` set. `rho` is
  a grid and path column (`rho_0.5/`), and `cts_val_analysis.R` writes through
  `combo_dir()` rather than a hand-built path.
- **One trial, split.** See Design. The chunks used to be two separate draws,
  each with `bW` calibrated to its own size.
- **`p_cts_adj`** in `top_var_tests` — see What is compared across chunks.
- **DR SuperLearner added** as a third estimator, with crossfit TE-VIMs — see
  Estimators and Variable importance. The TE-VIM scoring shared by all three
  is now one helper, `te_vim_scores()`; the two forest TE-VIMs compute exactly
  what they did before.
- **Shared code moved to `../val_common.R`** when the binary arm was added:
  the split, `fit_val_methods()`, the helpers and `chunk_validations()` (the
  four comparisons, formerly inline in `cts_val_analysis.R`). The comparisons
  were checked identical to the inline code at four interim points.
- **HC3 interaction tests** (`robust = TRUE`), after the move above — see What
  is compared across chunks. Every p-value in `subgroups` and `top_var_tests`
  changes; the other two comparisons do not.
- **`cts_val_rerun.sh` reset** to what `check_failed()` computes from
  `cts_val_1.sh`. It still asked for `ncpus=6:ompthreads=5:mem=12gb` from the
  old 5-worker setup, and `check_failed()` only ever raises a rerun's resources.

Before that:

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
stubbed in the original code and still is. All three estimators' per-run objects
now carry `BLP_whole`/`independence_cate`/`independence_po` in the same shape
`R/metrics.R::hte_test_metrics()` consumes, which is what a fifth chunk
comparison (alongside subgroups/variance/var-imps/top-var) would build on.

**Not in scope for this refactor.** `R/regression_check.R`'s old-vs-new
harness compares a behaviour-*preserving* refactor against a baseline; the
causal-forest nuisance change above means there is no valid "before" to check
against here. Worth adding once a fresh baseline exists.
