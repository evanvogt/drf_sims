# Missing covariates

What happens to CATE estimation when covariates are missing, and which handling
method costs least?

Three studies:

| folder | question |
|---|---|
| `continuous/` | the main design — continuous outcome, 8 handling methods × 3 mechanisms |
| `binary/` | the same design with a binary outcome, the treatment effect on the risk-difference scale |
| `ci_example/` | how to build a *confidence interval* after multiple imputation |

## The design (continuous and binary)

| | |
|---|---|
| scenarios | 1–5 — the main study's scenarios 1–5, same numbers |
| n | 500 |
| type | `both` (prognostic and predictive covariates amputated) |
| prop | 0.3 |
| mechanism | MAR, MNAR, MNAR-Y |
| method | 8 handling methods + `complete_data` reference |
| runs | 100 |
| array | **12,600 jobs** — rows 1–9,900 are scenarios 1, 3, 4, 5 (scenario 1 has no MNAR-Y, so 900 rows are dropped); rows 9,901–12,600 are scenario 2 |

Scenario numbering is the main study's: the missing-data set is the main set's
scenarios 1–6, so scenario `k` here is scenario `k` in `continuous/` and
`binary/`. The scenarios of interest, 1–4, are all run. Scenario 2 (simple HTE
on the continuous `X4`) was added after the others had been run, so it is the
second block of the grid, rows 9,901–12,600.

**Before 2026-09-26 this folder numbered them differently**, and results saved
before then use the old numbers:

| now | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| before 2026-09-26 | 1 | 6 | 4 | 5 | 2 | 3 |

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

## Open decision: pooling the `multiple_imputation` arm's HTE tests

`multiple_imputation` runs fit each of the 50 imputed datasets and
Rubin-combine with `combine_mi()` (`R/cate_models.R`), which pools only `tau`
and `variance`. Each imputation's fit runs the BLP and independence tests, but
until 2026-09-27 the analysis scripts threw them away, so `BLP_p`, `indep_cate`
and `indep_po` were `NA` for **every** model on all 1,400 MI runs per study.
Those runs kept no nuisances either, so the gap could not be patched.

**From the re-run on, the per-imputation tests are saved.** Each MI arm
(`causal_forest`, `dr_random_forest`, `dr_semi_oracle`) carries `mi_tests`, a
50-row table from `mi_test_table()` with the ingredients any candidate pooling
rule needs:

| columns | for |
|---|---|
| `blp_estimate`, `blp_se`, `blp_df` (β₂ and its residual df) | Rubin's rules on the BLP coefficient |
| `indep_{cate,po}_stat`, `indep_{cate,po}_df` (chi-square) | a statistic-pooling rule such as D2 |
| `blp_p`, `indep_cate_p`, `indep_po_p` | a p-value combination rule (Fisher, Stouffer, median p) |

A `NA` BLP row is bug L's degenerate-tau fallback. A failed independence test
reads `p = 1`, `stat = 0`, `df = NA`.

**No pooling rule is applied yet.** What "the" heterogeneity test across 50
imputations should be is still a methodological decision. The options answer
slightly different questions and none is the obvious default. It will be made
at metrics time, reading `mi_tests` from the collected results, with no
re-run. Until then `hte_test_metrics()` sees no `BLP_whole` /
`independence_*` on MI arms and reports `NA`, which means "not pooled yet",
not "the test failed".

The same gap carries over to the true-CATE HTE test evaluation
(`*_true_cate_tests.RDS`, `true_cate_test_row()` in `R/cate_models.R` — see
`continuous/README.md`): `multiple_imputation` rows are `NA`/`NA` there too,
for the same reason (`data` is a list of 50 imputed data.frames, not one),
even though that evaluation needs no nuisances from `nuisances_rf` at all.
It is the same pooling question as above, not a missing-nuisance problem, that
is unresolved for those rows. Its BLP half does not need X, and Y and W are
never imputed, so that half could be filled in from the saved results once
the question is settled.

## Bugs fixed here

**Bug B** — `cts_miss_collect.R` looked for `AUX`/`AUX-Y` where the DGM said `MNAR`/`MNAR-Y`; two mechanisms were silently collected as empty. Fixed.

**Bug C** — `rel_efficiency` was `NA` everywhere; the metrics scripts built the reference from `method == "complete_data"`, which the collect grid did not contain. Fixed.

**Bug D** — the array index meant two different things: both analysis scripts filtered `method == "complete_data"` *after* `expand.grid`, renumbering every row. A `failed_ids.txt` would have resubmitted the wrong parameters. The grid now lives in `<prefix>_config.R` and is never filtered. Fixed.

**Bug M** — `missing/binary`'s `dr_oracle` was handed log-odds as outcome predictions: `6b06db3` dropped the `plogis` from its oracle formulas but kept `oracle_link = "identity"`. Fixed; no finished result was affected — see `missing/binary/README.md`. Since the risk-difference DGM every oracle formula returns the outcome mean, and the `oracle_link` argument is gone, so this mismatch can no longer be written.

**Bug N** — `missing/binary`'s MNAR-Y truth was evaluated at U = 0 rather than averaged over U. The two agree on the identity scale, not the logit one, so about 0.02–0.03 of every binary MNAR-Y arm's bias was the truth's. Fixed in the generator, and repaired from the saved `p0`/`p1` when `bin_miss_metrics.R` ran. **Moot since the risk-difference DGM:** U now enters as `bU·tanh(U)`, which has mean zero, so the U-free truth *is* the average over U, on the binary scale as on the continuous one. The quadrature and the metrics-time repair are gone. See `missing/binary/README.md`.

## Status

**Both re-run all 12,600 rows**, after archiving their old trees
(`R/archive_old_results.R`; root `README.md`, Status, step 0):

- `missing/continuous` — bug O changed `continuous_missing`'s baseline and `bW`
  (see `continuous/README.md`), which moves every dataset, on top of the
  crossfitting strategy change and bug F. Its finished rows 1–9,900 are
  superseded, not "untouched".
- `missing/binary` — bug P and then the risk-difference DGM
  (`binary/README.md`), on top of the three defects and bug K in its own
  README. Submit on the risk-difference code.

Archiving first matters here in particular: the old trees use the
pre-2026-09-26 numbers, so the new scenario 2 would land on the paths of the
old scenario 2 (now scenario 5), and the check script would count the old files
as the new scenario being done.

```bash
qsub missing/continuous/jobscripts/cts_miss_1.sh       # rows 1-9,900
qsub missing/continuous/jobscripts/cts_miss_extra.sh   # rows 9,901-12,600, scenario 2
qsub missing/binary/jobscripts/bin_miss_1.sh
qsub missing/binary/jobscripts/bin_miss_extra.sh
```

Then check, collect and metrics as usual. **There is no patch step.** New runs
carry every model's HTE tests natively (`PROFILES$missing`), and the MI arm
saves its per-imputation tests (above). The back-fill applied to the archived
results, and the scripts that did it (removed after commit `e7b1d59`), are
summarised under "Every model carries the HTE tests" in
`missing/binary/README.md`.

`missing/ci_example` — see its own README.
