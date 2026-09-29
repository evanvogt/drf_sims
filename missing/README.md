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
| covariates | correlated (copula, ρ = 0.5), X01–X03 auxiliary — see below |
| mechanism | MAR, MNAR-Y0, MNAR-tau |
| method | 8 handling methods + `complete_data` reference |
| runs | 100 |
| array | **12,600 jobs** — rows 1–9,900 are scenarios 1, 3, 4, 5 (scenario 1 has no MNAR-tau, so 900 rows are dropped); rows 9,901–12,600 are scenario 2 |

Scenario numbering is the main study's: the missing-data set is the main set's
scenarios 1–6, so scenario `k` here is scenario `k` in `continuous/` and
`binary/`. The scenarios of interest, 1–4, are all run. Scenario 2 (simple HTE
on the continuous `X4`) was added after the others had been run, so it is the
second block of the grid, rows 9,901–12,600.

**Runs stay at 100 for every scenario**, unlike `continuous/` and `binary/`,
which take scenarios 1–4 to 500. More runs cost far more here: one run is one
job per mechanism × method, so each extra run of scenarios 1–4 adds 99 jobs
(scenario 1 has no MNAR-tau), against 16 in the main studies (4 scenarios × 4
sample sizes). Going to 500 runs for scenarios 1–4 would add 39,600 jobs per study (52,200
in all, about 4× the current 12,600), so 79,200 across `continuous/` and
`binary/`. That is four more array scripts per study under the 10,000-subjob
limit. The `multiple_imputation` jobs are the costliest, since each fits 50
imputed datasets. If more runs are ever needed, append runs 101+ as a new block
at the end of the grid, as the main studies do, so existing row numbers keep
their meaning.

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

### Covariates are correlated here (since 2026-09-28)

Unlike the main studies, the missing-data sets draw X1–X5 and X01–X03 from a
Gaussian copula with exchangeable latent correlation ρ = 0.5 (X1 and X3 by
thresholding). X01–X03 are **auxiliaries**: correlated with X1–X5, in neither
the outcome model nor the CATE, and never amputated, so they are what
`regression` imputation and the `IPW` model have to work with. X04/X05 stay pure
noise. Before this every covariate was independent, which left the imputation
methods nothing to impute from and made MAR indistinguishable from MCAR as far
as the CATE is concerned. See `ADEMP.md`.

### Mechanisms

- **MAR** — missingness depends on the observed covariates. With correlated
  covariates the missing values differ from the observed ones, so this is the
  arm that tests the imputation methods.
- **MNAR-Y0** — missingness is driven by an unobserved `U` (independent of X)
  that also shifts the **control outcome**: missingness is related to
  prognosis.
- **MNAR-tau** — `U` enters the **treatment effect** instead: missingness is
  related to benefit. Not defined for scenario 1, which has no treatment effect
  heterogeneity to relate to. Displayed as "MNAR-τ".

These replaced MAR / MNAR / MNAR-Y on 2026-09-28. The old "MNAR" had `U` in
neither X nor Y, so it was MCAR in effect, and was dropped; the old "MNAR-Y" is
MNAR-tau. The old names (and `ci_example`'s older `AUX` / `AUX-Y`) are now
rejected with an error rather than mapped.

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
reads `p = 1`, `stat = 0`, `df = NA`. From the `is_constant()` guard on, a
constant tau gives `NA` for both tests (`p`, `stat` and `df` all `NA` for the
independence tests); rows saved before it can instead carry an independence
p ≈ 0 for a constant tau, which coin returns with only a warning (see
`sample_size/continuous/README.md`, "True-CATE HTE test evaluation").

**No pooling rule is applied yet.** What "the" heterogeneity test across 50
imputations should be is still a methodological decision. The options answer
slightly different questions and none is the obvious default. It will be made
at metrics time, reading `mi_tests` from the collected results, with no
re-run. Until then `hte_test_metrics()` sees no `BLP_whole` /
`independence_*` on MI arms and reports `NA`, which means "not pooled yet",
not "the test failed".

The same gap carries over to the true-CATE HTE test evaluation
(`*_true_cate_tests.RDS`, `true_cate_test_row()` in `R/cate_models.R` — see
`sample_size/continuous/README.md`): `multiple_imputation` rows are `NA`/`NA` there too,
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
  (see `sample_size/continuous/README.md`), which moves every dataset, on top of the
  crossfitting strategy change and bug F. Its finished rows 1–9,900 are
  superseded, not "untouched".
- `missing/binary` — bug P and then the risk-difference DGM
  (`sample_size/binary/README.md`), on top of the three defects and bug K in its own
  README. Submit on the risk-difference code.
- Both, and `ci_example` (which runs on `continuous_missing`), also for the
  2026-09-28 redesign: correlated covariates with X01–X03 as auxiliaries, and
  the mechanisms MAR / MNAR-Y0 / MNAR-tau. Every dataset changes, including the
  `complete_data` arm's.

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
