# Missing covariates — continuous outcome

The main missing-data design. See `missing/README.md` for the mechanisms,
handling methods and the shared bug fixes (B, C, D — all of which originated
here).

| | |
|---|---|
| array | **12,600 jobs**: `cts_miss_1.sh` (1–9900, scenarios 1, 3, 4, 5) and `cts_miss_extra.sh` (9901–12600, scenario 2) |
| results | `../results/missing/continuous/scenario_<k>/<n>/<type>/<prop>/<mechanism>/<method>/` |
| metrics | `cts_miss_metrics.RDS`, including `rel_efficiency` against the `complete_data` arm; plus `cts_miss_true_cate_tests.RDS`, the true-CATE HTE test evaluation |

## Files

| file | role |
|---|---|
| `cts_miss_config.R` | the grid — **the** definition, never filtered |
| `cts_miss_dgms.R` | names the `continuous_missing` scenario set + the missingness machinery |
| `cts_miss_models.R` | `family = gaussian()`, `profile = "missing"`, `ipw` threaded through |
| `cts_miss_analysis.R` | array entry point |
| `cts_miss_results.R` / `.qmd` | every metric, to `results/all_figures/` — the diagnostic counterpart to the chapter script `results_processing/thesis_figures/miss_cts.R` |
| `mi_scratch.R` | exploratory, unmaintained |
| `results_cts_miss.{R,Rmd,html}` | **superseded** by `cts_miss_results.R`/`.qmd` — kept for reference, renamed to the repo's legacy convention (cf. `continuous/results_cts.R`). Points at the old HPC path `/rds/general/...` and at `results/new_format/metrics_cts_miss_df.RDS`, neither of which exists any more; scenarios are `"scenario_5"` strings, and only bias and MSE are plotted. It was called `cts_miss_results.*`, so its `.html` was being clobbered by the new `.qmd`'s render output. |

The `ipw` argument is the only thing separating this study's estimators from
`continuous/`: it becomes grf `sample.weights` and SuperLearner `obsWeights`.
`ipw = NULL` reproduces the unweighted path exactly.

## Gotchas

**Scenario numbers are the main study's** (since 2026-09-26; see
`missing/README.md` for the old numbering). Scenario 2 was appended later, as a
second block in `cts_miss_config.R`, so rows 1–9900 kept their meaning.

**`multiple_imputation` returns a list of 50 datasets, not one.** The analysis
script fits each and Rubin-combines with `combine_mi()`; only
`causal_forest`, `dr_random_forest` and `dr_semi_oracle` are combined. Each of
those arms also saves `mi_tests`, the 50 imputations' HTE tests, unpooled: the
pooling rule is still to be decided (`missing/README.md`), so the arm's
`BLP_p`/`indep_*` metrics stay `NA` until it is.

**The row-dropping methods drop rows from the truth too.** `complete_cases` and
`IPW` return `retained_indices`, and `generate_and_process_data()` subsets
`truth` to match — otherwise estimates and truth would be misaligned.

**`cts_miss_true_cate_tests.RDS` runs the BLP and independence tests on the
true CATE and true nuisances instead of an estimator's** (`truth$tau`,
`truth$p0`, `W.hat = 0.5` — see `sample_size/continuous/README.md`'s "True-CATE HTE test
evaluation"; its "HTE tests" describes the test columns, including the
one-sided, HC3 `BLP_p_os`), one row per (scenario, n, type, prop, mechanism,
method, run), no per-model dimension. Every scenario-1 true-CATE test is `NA`:
the true CATE is constant there. `method == "multiple_imputation"` rows are `NA`/`NA`
there, the same pooling question as the estimated-CATE tests above: `data` is a
list of 50 imputed data.frames, with no single covariate matrix to test against.

**`dr_random_forest` used to carry no BLP or independence test in this study**,
unlike `continuous/`, so `BLP_p` was `NA` for that one model. Copy-paste drift
rather than a decision, and the decision has been taken: every model carries the
tests where possible. `PROFILES$missing` (`R/cate_models.R`) now sets
`dr_rf_tests = TRUE`, so the re-run writes them natively. The results made
before that were back-filled in place by a one-off patch rather than re-run,
and are now archived; that history is written up once, in
`missing/binary/README.md`, since `profile = "missing"` covers both studies.
The re-run has no patch step.

## Running it

```bash
qsub missing/continuous/jobscripts/cts_miss_1.sh       # 1-9900
qsub missing/continuous/jobscripts/cts_miss_extra.sh   # 9901-12600, scenario 2
Rscript missing/continuous/cts_miss_check.R
qsub missing/continuous/jobscripts/cts_miss_collect.sh
qsub missing/continuous/jobscripts/cts_miss_metrics.sh
```

Then, for the figures:

```bash
Rscript missing/continuous/cts_miss_results.R           # every metric, to results/all_figures/
quarto render missing/continuous/cts_miss_results.qmd   # the same, as a browsable report
```

To run only the `complete_data` reference arm, take
`grid_indices(study, method = "complete_data")` (1,100 indices) and submit those
— do **not** filter the grid.

## Status

**All 12,600 rows re-run** (`cts_miss_1.sh`, then `cts_miss_extra.sh` for
scenario 2), after archiving the old tree with `R/archive_old_results.R` (root
`README.md`, Status, step 0). Bug O, below, supersedes rows 1–9900. The old
tree uses the pre-2026-09-26 numbers, so running into it would put the new
scenario 2 on top of the old scenario 2 (now 5). Then collect/metrics - no
patch step, see `missing/README.md` Status.

**Re-runs required** — for the crossfitting strategy change to
`R/cate_models.R` (see root README Methods/Status), which moves all five
estimator arms, and separately for bug F (`dr_superlearner` only). Bugs B and C
were collection/metrics problems, so re-running `cts_miss_collect.R` and
`cts_miss_metrics.R` over the existing per-run files recovers the mechanisms
that were being missed and populates `rel_efficiency` — no cluster time needed
for those two.

**Also re-run for bug O** — the `continuous_missing` table now shares the main
study's baseline (`b0 = 0.4`, `b1 = −0.5`, `b2 = 1`; scenario 4 alone had
`b2 = 1` before) and its `bW` calibration: every trial is planned for 80% power
under homogeneity, so every scenario has a true ATE of −0.29 at n = 500. See
`sample_size/continuous/README.md`'s "Outcome model and `bW` calibration". With complete
data, realised power is 0.81 in scenario 1 and 0.63–0.76 elsewhere. MNAR-Y's `U`
term is heterogeneity the plan knows nothing about, like the rest, so it is
left out of the calibration too. That keeps one `bW` and one truth per scenario
across every mechanism, and lowers realised power under MNAR-Y to 0.54–0.67.
With `PROFILES$missing` already set, re-run results carry the
`dr_random_forest` tests, so the back-fill patch only ever mattered for results
made before bug O.
