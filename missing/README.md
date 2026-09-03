# Missing covariates

What happens to CATE estimation when covariates are missing, and which handling
method costs least?

Three studies:

| folder | question |
|---|---|
| `continuous/` | the main design — continuous outcome, 8 handling methods × 3 mechanisms |
| `binary/` | the same design on a logit scale |
| `ci_example/` | how to build a *confidence interval* after multiple imputation |

## The design (continuous and binary)

| | |
|---|---|
| scenarios | 1, 2, 4, 5 — a **renumbered** five-scenario set, not a subset of the main ten |
| n | 500 |
| type | `both` (prognostic and predictive covariates amputated) |
| prop | 0.3 |
| mechanism | MAR, MNAR, MNAR-Y |
| method | 8 handling methods + `complete_data` reference |
| runs | 100 |
| array | **9,900 jobs** (scenario 1 has no MNAR-Y, so 100 rows are dropped) |

Scenario numbering here is **not** the main study's. Scenario `k` corresponds to
main-study scenario 1, 2, 4, 8, 9 respectively — so "scenario 4" means something
different in this folder than in `continuous/`.

## Results files

Each study has a `*_results.R` / `*_results.qmd` pair showing **every** metric
its `*_metrics.R` computes — bias, ATE bias, MSE, RMSE, MAE, relative
efficiency, correlation, Spearman, sign accuracy, the three HTE-test p-values,
and the `n_na` diagnostic. Figures go to `results/all_figures/missing/<study>/`.

These are the diagnostic view. `results_processing/thesis_figures/miss_*.R` is
the chapter view: the two or three panels that get printed, written to
`results/thesis_figures/`. The two never overwrite each other.

Labels, palette, the mean ± MCSE summary and the point-range panel are shared
via `R/figures.R` — `MISS_SCENARIO_LABELS`, `METHOD_LABELS`,
`MECHANISM_LABELS`, `STRATEGY_LABELS` and `point_range_plot()` in particular.
Rename an estimator or a handling method there and every figure follows.

### Mechanisms

- **MAR** — missingness depends on the observed covariates
- **MNAR** — driven entirely by an unobserved `U`
- **MNAR-Y** — `U` also enters the treatment effect, so missingness is related to
  the outcome. Not defined for scenario 1, which has no treatment effect
  heterogeneity to relate to.

`missing/ci_example/` still calls these `AUX` / `AUX-Y`. Both spellings are
accepted and normalised in `R/missingness.R`.

### Handling methods

`complete_cases`, `mean_imputation`, `missforest`, `regression`,
`missing_indicator`, `IPW`, `multiple_imputation`, `none` (let the estimator
handle it), plus `complete_data` — a reference arm with **no missingness
introduced at all**, so the others can be scored against complete-data
performance via `rel_efficiency`.

## Open gap: the `multiple_imputation` arm has no HTE tests, for any model

Not a defect in one estimator — a hole in the arm. `multiple_imputation` runs
fit each of the 50 imputed datasets and Rubin-combine with `combine_mi()`
(`R/cate_models.R`), which returns only `tau` and `variance`. The analysis
scripts then build `results` from those combined objects alone, so
`nuisances_rf` is never saved. Consequently `BLP_p`, `indep_cate` and
`indep_po` are `NA` for **every** model on all 1,100 MI runs per study.

Unlike the `dr_random_forest` gap, this one **cannot** be patched: there is
nothing on disk to recompute a BLP from. Closing it needs two things:

1. **A methodological decision.** What is "the" heterogeneity test across 50
   imputations? Rubin-combine the BLP coefficients and their standard errors,
   combine the per-imputation p-values (Fisher, Stouffer), or something else.
   Each answers a slightly different question and none is the obvious default.
2. **A re-run** of the MI arm — 1,100 runs per study, the most expensive method
   in the grid.

Until then the `NA`s are honest and should be read as "not computed", not as
"the test failed". `R/patch_hte_tests.R` detects these runs and refuses them,
which is why `check_all.R` reports `patchable_jobs` of 8,800 rather than 9,900.

The same gap carries over to the true-CATE HTE test evaluation
(`*_true_cate_tests.RDS`, `true_cate_test_row()` in `R/cate_models.R` — see
`continuous/README.md`): `multiple_imputation` rows are `NA`/`NA` there too,
for the same reason (`data` is a list of 50 imputed data.frames, not one),
even though that evaluation needs no nuisances from `nuisances_rf` at all —
it is the multiple-imputation *pooling* question in point 1 above, not a
missing-nuisance problem, that is unresolved for those rows.

## Bugs fixed here

**Bug B** — `cts_miss_collect.R` looked for `AUX`/`AUX-Y` where the DGM said `MNAR`/`MNAR-Y`; two mechanisms were silently collected as empty. Fixed.

**Bug C** — `rel_efficiency` was `NA` everywhere; the metrics scripts built the reference from `method == "complete_data"`, which the collect grid did not contain. Fixed.

**Bug D** — the array index meant two different things: both analysis scripts filtered `method == "complete_data"` *after* `expand.grid`, renumbering every row. A `failed_ids.txt` would have resubmitted the wrong parameters. The grid now lives in `<prefix>_config.R` and is never filtered. Fixed.

## Status

`missing/continuous` — re-runs for the crossfitting strategy change (all five
arms) and separately for bug F (`dr_superlearner` only).

`missing/binary` — re-runs entirely; see its own README, it had three further
defects plus a `pretest_superlearner()` crash (bug K), now fixed.

**Both** additionally owe the `dr_random_forest` HTE back-fill — a one-off
in-place repair of finished results, not a re-run. Track it in
`check_all_studies.md`'s `patch_status` column; a study reads `complete` on the
run counts while still owing this, which is exactly why that column exists. See
"Patched: every model now carries the HTE tests" in `missing/binary/README.md`.

`missing/ci_example` — see its own README.
