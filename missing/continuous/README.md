# Missing covariates — continuous outcome

The main missing-data design. See `missing/README.md` for the mechanisms and
handling methods.

| | |
|---|---|
| array | **12,600 jobs**: `cts_miss_1.sh` (1–9900, scenarios 1, 3, 4, 5) and `cts_miss_extra.sh` (9901–12600, scenario 2) |
| results | `../results/missing/continuous/scenario_<k>/<n>/<type>/<prop>/<mechanism>/<method>/` |
| metrics | `cts_miss_metrics.RDS`, including `rel_efficiency` against the `complete_data` arm; plus `cts_miss_true_cate_tests.RDS`, the true-CATE HTE test evaluation |

## Calibration

The `continuous_missing` table shares `sample_size/`'s baseline (`b0 = 0.4`,
`b1 = −0.5`, `b2 = 1`) and its `bW` calibration (`sample_size/README.md`,
"Calibration of bW"): every trial is planned for 80% power under homogeneity.
With correlated covariates the planned SD carries Cov(X1, X2), which moves the
planned effect to −0.27, and scenario 3's `bW` moves further because
E[X4·X5] = 0.5. The MNAR `U` term is heterogeneity the plan knows nothing
about, like the rest, so it is left out of the calibration too; that keeps one
`bW` and one truth per scenario across every mechanism. At n = 500, from
`Rscript R/calibration_report.R`:

| scenario | `bW` | realised power | MNAR-Y0 | MNAR-tau |
|---|---|---|---|---|
| 1 | −0.27 | 0.81 | 0.54 | — |
| 2 | −0.27 | 0.66 | 0.46 | 0.54 |
| 3 | −1.42 | 0.56 | 0.41 | 0.47 |
| 4 | −0.87 | 0.76 | 0.51 | 0.61 |
| 5 | −1.67 | 0.68 | 0.47 | 0.56 |

MNAR-Y0's `U` (`bU = 1`, as prognostic as X2) adds variance 1 to *both* arms,
the plan knows nothing about it, and so it costs the most power.

## Files

| file | role |
|---|---|
| `cts_miss_config.R` | the grid — **the** definition, never filtered |
| `cts_miss_dgms.R` | names the `continuous_missing` scenario set + the missingness machinery |
| `cts_miss_models.R` | `family = gaussian()`, `profile = "missing"`, `ipw` threaded through |
| `cts_miss_analysis.R` | array entry point |
| `../miss_results.R` / `.qmd` | shared with `binary/`: every metric, to `results/all_figures/` — the diagnostic counterpart to the chapter script `results_processing/thesis_figures/miss_cts.R` |

The `ipw` argument is the only thing separating this study's estimators from
`sample_size/continuous/`: it becomes grf `sample.weights` and SuperLearner
`obsWeights`. `ipw = NULL` reproduces the unweighted path exactly.

## Gotchas

**Scenario 2 is a second block** in `cts_miss_config.R`, appended after the
others, so rows 1–9900 keep their meaning.

**`multiple_imputation` returns a list of 50 datasets, not one.** The analysis
script fits each and Rubin-combines with `combine_mi()`; only
`causal_forest`, `dr_random_forest` and `dr_semi_oracle` are combined. Each of
those arms also keeps each imputation's CATE estimates (`tau_mi`, n × 50) and
the Rubin variance's two components (`variance_within`, `variance_between`),
and saves `mi_tests`, the 50 imputations' HTE tests, unpooled: the
pooling rule is still to be decided (`missing/README.md`), so the arm's
`BLP_p`/`indep_*` metrics stay `NA` until it is.

**The row-dropping methods drop rows from the truth too.** `complete_cases` and
`IPW` return `retained_indices`, and `generate_and_process_data()` subsets
`truth` to match — otherwise estimates and truth would be misaligned.

**`cts_miss_true_cate_tests.RDS` runs the BLP and independence tests on the
true CATE and true nuisances instead of an estimator's** (`truth$tau`,
`truth$p0`, `W.hat = 0.5` — see `sample_size/README.md`'s "True-CATE
tests"; its "HTE tests" describes the test columns, including the
one-sided, HC3 `BLP_p_os`), one row per (scenario, n, type, prop, mechanism,
method, run), no per-model dimension. Every scenario-1 true-CATE test is `NA`:
the true CATE is constant there. `method == "multiple_imputation"` rows are `NA`/`NA`
there, the same pooling question as the estimated-CATE tests above: `data` is a
list of 50 imputed data.frames, with no single covariate matrix to test against.

**Every model carries the HTE tests**, `dr_random_forest` included
(`PROFILES$missing` sets `dr_rf_tests = TRUE` in `R/cate_models.R`).

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
Rscript missing/miss_results.R continuous         # every metric, to results/all_figures/
quarto render missing/miss_results.qmd -P outcome:continuous --output cts_miss_results.html
```

To run only some methods, take e.g. `grid_indices(study, method = "complete_data")`
(1,100 indices) or `grid_indices(study, method = c("IPW", "regression"))` (2,800)
and submit those — do **not** filter the grid.
