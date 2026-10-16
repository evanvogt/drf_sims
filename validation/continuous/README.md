# Interim-analysis validation — continuous outcome

See `validation/README.md` for what a validation study checks and why. The
binary arm (`../binary/`) runs the same design, and the code both arms share
lives in `../val_common.R`.

## Design

| | |
|---|---|
| covariates | correlated, rho = 0.5 — the `continuous_corr_0.5` set (`R/dgm_scenarios.R`, CORRELATED COVARIATES), the same set `sample_size/` runs as its primary arm |
| scenarios | 1, 2, 3, 4 (`R/dgm_scenarios.R`, `DESC_10`): 1 the null, no HTE; 2 simple HTE, `b4 * X4`; 3 complex, `b3 * X3 + b4 * X4 + b45 * X4 * X5`; 4 non-linear, `b4 * cos(X4)` |
| n | 1000, **one** trial split into its first `n * interim_prop` participants and the rest |
| interim_prop | 0.25 to 0.75 in steps of 0.05 — 11 interim points |
| runs | 100 per scenario — **4400 array jobs**, 1 h walltime (see Sizing the job) |
| folds | DR SuperLearner only: 5 for a chunk below 500 rows, else 10 |
| results | `../results/validation/continuous/rho_0.5/scenario_<k>/1000/<interim_prop>/res_sim_<run>.RDS` |

Each run draws one trial of `n` and splits it. Because the seed is the run alone,
a run's chunk 1 at 0.25 is nested in its chunk 1 at 0.30, so every curve across
`interim_prop` is paired within run. Two separate draws would not do: since
`generate_scenario_data()` calibrates `bW` to 80% power at the `n` it is given,
each chunk would be its own trial with its own ATE.

`cts_val_config.R` rounds the `interim_prop` sequence to 2 dp, which is
load-bearing rather than cosmetic: `interim_prop` is a `path_cols` entry, so its
`as.character()` form becomes a results directory name, and an unrounded
`seq(by = 0.05)` yields values like `0.30000000000000004` that `get_results()`
can no longer match back to their grid row.

Only rho = 0.5, scenarios 1–4 and n = 1000 are run; all three are columns of
the grid (`cts_val_config.R`), so extending to rho = 0, more scenarios or sizes
is a one-line change. The same run number draws the same participants in every
scenario; only the CATE changes.

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

`xgboost` and `SHAPforxgboost` are used only by these studies. They are in
`environment.yml`; `cts_val_testing.R` check 1 catches their absence before
4400 jobs go out.

## What is compared across chunks

`cts_val_analysis.R` saves four comparisons per run, under `validations`:

| | |
|---|---|
| `subgroups` | fit a tree on chunk 1's top/bottom-10%-CATE responder groups, predict them into chunk 2, test the `W x subgroup` interaction there |
| `variances` | `var(tau)` in each chunk, against the true `Var(tau)` (1 in scenario 2; the report computes it per scenario) |
| `var_imps` | each measure's covariate ranking in each chunk, and the rank change |
| `top_var_tests` | chunk 1's most important covariate, interaction-tested on the remaining participants — as a continuous `W x X_top` interaction (`p_cts`), the same adjusted for every other covariate (`p_cts_adj`), after a median split (`p_split`), after a split learned on chunk 1 (`p_dd`, with `dd_cuts` and `dd_prop`), and as a spline projection learned on chunk 1 (`p_proj`) |

Every interaction test (the subgroup tests and all five
top-covariate forms) is saved twice, with classical and with HC3 standard
errors: `subgroups[[model]]` is a two-row `data.frame(se, top, bottom)` and
`top_var_tests[[model]]` leads with an `se` column (`chunk_validations_se()`).
`variances` and `var_imps` do not depend on the standard errors and are saved
once.

**Estimates and the direction of a replication.** A p-value alone cannot say
whether chunk 2 found the effect chunk 1 claimed, or one pointing the other
way. Each test therefore saves four more numbers next to its p-value, named
`<stat>_<test>`. The tests are `top`/`bottom` for the subgroups and `cts`,
`cts_adj`, `split`, `dd`, `proj` for the top-covariate forms:

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

**Scenario 1, the null.** tau = bW for everyone, so every
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

**`p_split` on a binary covariate.** `median_split()` uses a two-valued
covariate as its own split: `xt > median(xt)` is all 0 for X3 ~ Bernoulli(0.7)
(median 1), and X3 is the likely nominee in scenario 3.

The comparisons are `chunk_validations()`, and the estimation logic behind them
(both importance measures and `interaction_pval()`), in `../val_common.R`,
shared with the binary arm.

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
there and keeps the two arms' tests the same.

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
| `study.R check` | finds missing runs, writes `jobscripts/failed_ids.txt`, and updates `-J` and the resource request in the rerun jobscript |
| `study.R collect` | gathers per-run files into `cts_val_all.RDS` |
| `cts_val_metrics.R` | flattens `validations` into tidy `cts_val_metrics.RDS` |
| `results_cts_val.R` | summary plots, one scenario at a time (`Rscript results_cts_val.R <scenario>`), HC3 |
| `../val_results.qmd` | the written-up report, shared with `binary/`, one render per scenario (Quarto `outcome` and `scenario` params) |

## Running it

```bash
Rscript validation/continuous/cts_val_testing.R full   # do this first - see below
qsub validation/continuous/jobscripts/cts_val_1.sh     # 1-4400
Rscript study.R check validation/continuous/cts_val_config.R          # writes failed_ids.txt if any are missing
qsub validation/continuous/jobscripts/cts_val_collect.sh
qsub validation/continuous/jobscripts/cts_val_metrics.sh
cd validation
quarto render val_results.qmd -P outcome:continuous -P scenario:1 --output cts_val_results_sc1.html   # and 2, 3, 4
```

`cts_val_testing.R` runs eight cheap checks; `full` adds one real
replicate end to end and reports how long it took, which is the number to read
against the jobscript's walltime before submitting. It exits non-zero on any
failure. Drop `full` for a quick local smoke test.

`check_failed()` counts any results file on the grid's paths as done, so if
the code changes in a way that makes old results incomparable, archive the
trees first with `R/archive_old_results.R --label <label> --trees
validation/continuous validation/binary` (dry run, then `--apply`).

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
land in the same place. Budget 4400 × (the replicate time
`cts_val_testing.R full` reports) / 8.

Rows that already have a results file are skipped (the same missing-run logic
`study.R check` uses), so an interrupted session — or one that hits its
walltime — is resumed by sourcing the file again, and so is a row that failed.
Edit `ids` at the top to run a subset, e.g.
`grid_indices(study, interim_prop = 0.25)`; trailing arguments do the same
non-interactively (`Rscript cts_val_run.R 1 2 3`). Per-row logs go to
`<results>/validation/continuous/session_logs/`.

It never writes `jobscripts/failed_ids.txt` and never edits `cts_val_rerun.sh` —
those belong to the queue path. Afterwards, check, collect and
metrics run exactly as they would after an array run.

### Sizing the job

`cts_val_1.sh` runs one worker with one grf thread
(`cts_val_analysis.R <i> 1 1`) on `ncpus=1`, with a 1 h walltime, both set by
hand. The DR SuperLearner's TE-VIM refits (10 covariates × 5–10 folds of a
9-learner SuperLearner, per chunk) dominate a row's time.
`cts_val_testing.R full` passes the
same `1 1`, so the replicate time it reports is the array job's. Check 6 also
prints a step-by-step timing of one 250-row chunk. If a row is too slow on one
core, raise `ncpus`, `ompthreads` and the `workers` argument together
(`cts_val_analysis.R <i> k 1` on `ncpus=k:ompthreads=k`). The SuperLearner
folds and the TE-VIM refits all run through `future_map()`, so they spread
across the extra workers.

`workers` and `grf_threads` reach `cts_val_analysis.R` as trailing arguments on
the `Rscript` line, so the PBS resource request and the R-level parallelism
cannot drift apart. Called without them it uses `workers = 5`, grf
unthrottled.

## Notes

- **Causal-forest nuisance.** `R/cate_models.R::run_causal_forest()` uses
  `Y.hat.cf`, a separate `X`-only regression — the input `causal_forest()`
  expects — not a regression of `Y` on `(W, X)` with the observed `W` plugged
  back in. The causal forest's TE-VIMs reuse the same nuisances.
- **Not implemented: comparing HTE tests between chunks.** All three
  estimators' per-run objects carry `BLP_whole`/`independence_cate`/
  `independence_po` in the shape `R/metrics.R::hte_test_metrics()` consumes,
  which is what a fifth chunk comparison would build on.
- **No `R/regression_check.R` entry yet.** Worth adding once a baseline exists.
