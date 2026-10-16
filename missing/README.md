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

Scenario `k` here is scenario `k` in `sample_size/`. The scenarios of interest,
1–4, are all run. Scenario 2 (simple HTE on the continuous `X4`) was added
after the others, so it is appended as the second block of the grid, rows
9,901–12,600, and the earlier rows keep their meaning.

**Runs stay at 100 for every scenario**, unlike `sample_size/`, which takes
scenarios 1–4 to 500. One run here is one job per mechanism × method, so each
extra run of scenarios 1–4 adds 99 jobs (scenario 1 has no MNAR-tau), against
16 in `sample_size/` (4 scenarios × 4 sample sizes): 500 runs would add 39,600
jobs per study. The `multiple_imputation` jobs are the costliest, since each
fits 50 imputed datasets. If more runs are ever needed, append runs 101+ as a
new block at the end of the grid, so existing row numbers keep their meaning.

Results saved before 2026-09-26 use older scenario numbers; the map is in
`docs/DEVLOG.md`.

### Covariates are correlated

The missing-data sets draw X1–X5 and X01–X03 from a Gaussian copula with
exchangeable latent correlation ρ = 0.5 (X1 and X3 by thresholding). X01–X03
are **auxiliaries**: correlated with X1–X5, in neither the outcome model nor
the CATE, and never amputated, so they are what `regression` imputation and the
`IPW` model have to work with. X04/X05 stay pure noise. With independent
covariates the imputation methods would have nothing to impute from, and MAR
would be indistinguishable from MCAR as far as the CATE is concerned. See
`ADEMP.md`.

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

Any other mechanism name is rejected with an error.

**The mechanisms are paired within a run.** `U` is drawn under MAR too, and not
used, so a run's three mechanisms share W, every covariate, the continuous
error and X04/X05; MAR and MNAR outcomes differ by `U`'s term alone. The
amputation still differs by design.

**Binary `U`.** `missing/binary` has its own HTE scale (`RD_SCALE_MISS`,
`R/dgm_scenarios.R`), derived at n = 500 alone, which leaves room for
`bU = 0.12`, `sU = 2`; complete cases' ATE bias under MNAR-tau is about a fifth
of the RD. The binary mechanisms stay far weaker than the continuous ones - a
Bernoulli outcome's own variance dominates - so compare the two outcomes'
mechanism effects with that in mind (strength table in `ADEMP.md`;
`miss_dgm_checks.R`).

### How the methods are scored

The truth is τ at the *unamputed* covariates, which no method can recover for
a unit whose effect modifiers are missing. Scored on all analysed units, the
imputation methods carry an error floor on those units that `complete_cases`
and `IPW` - scored on their complete units only - never face. So every run
saves its missingness mask (`miss_mask`, 500 × 5) and `retained_indices`, and
every point metric comes on all units, the complete units (`_cu`, the subset
every method has - the comparison across methods) and the incomplete units
(`_iu`). The headlines are `rel_efficiency_cu` and `bias_diff_complete_cu`.
Under MNAR-tau a secondary truth that adds E[U-term | complete or incomplete]
is scored too (`rmse_r` etc.). Details in `ADEMP.md`, "Estimands" and
"Performance measures".

`miss_dgm_checks.R` (run from the repo root, `Rscript
missing/miss_dgm_checks.R`) checks the DGM's claims at n = 20,000: the
copula, the amputation, how strongly each mechanism selects, the pairing, and
each handling method's signature on complete and incomplete units.

### Handling methods

`complete_cases`, `mean_imputation`, `missforest`, `regression`,
`missing_indicator`, `IPW`, `multiple_imputation`, `none` (let the estimator
handle it), plus `complete_data` — a reference arm with **no missingness
introduced at all**, so the others can be scored against complete-data
performance via `rel_efficiency_cu`. The `complete_data` arm still saves the
mask the other arms' amputation leaves, so it can be split the same way.

`missforest`, `multiple_imputation` and `regression` impute **within each
treatment arm, with Y as a predictor** (`impute_by_arm()` in
`R/missingness.R`): leaving Y and W out of the imputation model flattens the
heterogeneity the study is trying to estimate. `regression` stays single and
deterministic, predicting each incomplete covariate from X01–X05 and Y.

`IPW`'s model for P(complete) is `X01–X05 + W * Y`, and its weights are
stabilised by the marginal proportion complete. Without Y and W the weights are
near constant under MNAR and IPW reproduces complete cases; with `W * Y` it
halves complete cases' MNAR-tau ATE bias (continuous).

The aim is that every handling method is one a researcher could reasonably use
in a secondary HTE analysis of an RCT. One consequence: under MNAR-Y0 /
MNAR-tau, `U` drives both the missingness and Y. Any imputation or weighting
model that conditions on Y (`regression`, `missforest`, `multiple_imputation`,
`IPW`) therefore opens a collider path, X → Y ← U → R. Only the methods that
ignore Y (complete cases, mean imputation, the indicator method) keep the
MCAR-like view of X under MNAR-Y0 - see `ADEMP.md`.

## Open decision: pooling the `multiple_imputation` arm's HTE tests

`multiple_imputation` runs fit each of the 50 imputed datasets and
Rubin-combine with `combine_mi()` (`R/cate_models.R`), which pools only `tau`
and `variance`. It also keeps what it pooled: `tau_mi`, each imputation's CATE
estimates, n × 50, and `variance_within` / `variance_between`, W̄ and B — for
Barnard–Rubin df, the fraction of missing information, or another pooling rule.

The per-imputation tests are saved. Each MI arm (`causal_forest`,
`dr_random_forest`, `dr_semi_oracle`) carries `mi_tests`, a 50-row table from
`mi_test_table()` with the ingredients any candidate pooling rule needs:

| columns | for |
|---|---|
| `blp_estimate`, `blp_se`, `blp_df` (β₂ and its residual df) | Rubin's rules on the BLP coefficient |
| `blp_se_hc3` (β₂'s HC3 SE; same estimate and df) | Rubin's rules for the `BLP_p_os` version of the test |
| `indep_{cate,po}_stat`, `indep_{cate,po}_df` (chi-square) | a statistic-pooling rule such as D2 |
| `blp_p`, `blp_p_os`, `indep_cate_p`, `indep_po_p` | a p-value combination rule (Fisher, Stouffer, median p) |

`blp_se_hc3` and `blp_p_os` match the headline `BLP_p_os` metric: one-sided,
HC3 SEs. For every other arm that is recomputed at metrics time from the run's
saved nuisances (`blp_inputs()`), but the MI arm keeps no per-imputation
nuisances, so `mi_test_table()` computes it at estimation time instead.

A `NA` BLP row is `run_blp_whole()`'s degenerate-tau fallback. A failed
independence test reads `p = 1`, `stat = 0`, `df = NA`. A constant tau gives
`NA` for both tests (`p`, `stat` and `df` all `NA` for the independence tests;
see `sample_size/README.md`, "True-CATE tests" under "HTE tests").

**No pooling rule is applied yet.** What "the" heterogeneity test across 50
imputations should be is still a methodological decision. The options answer
slightly different questions and none is the obvious default. It will be made
at metrics time, reading `mi_tests` from the collected results, with no
re-run. Until then `hte_test_metrics()` sees no `BLP_whole` /
`independence_*` on MI arms and reports `NA`, which means "not pooled yet",
not "the test failed".

The same gap carries over to the true-CATE HTE test evaluation
(`*_true_cate_tests.RDS`, `true_cate_test_row()` in `R/cate_models.R` — see
`sample_size/README.md`): `multiple_imputation` rows are `NA`/`NA` there too,
because `data` is a list of 50 imputed data.frames, not one. Its BLP half does
not need X, and Y and W are never imputed, so that half could be filled in from
the saved results once the question is settled.

## Running

```bash
qsub missing/continuous/jobscripts/cts_miss_1.sh       # rows 1-9,900
qsub missing/continuous/jobscripts/cts_miss_extra.sh   # rows 9,901-12,600, scenario 2
qsub missing/binary/jobscripts/bin_miss_1.sh
qsub missing/binary/jobscripts/bin_miss_extra.sh
```

Then check, collect and metrics as usual (see each study's README). New runs
carry every model's HTE tests natively (`PROFILES$missing`), and the MI arm
saves its per-imputation tests (above).

## Results files

`miss_results.R` / `miss_results.qmd` (one copy for both outcomes: an
`outcome` argument / Quarto param) show **every** metric a study's
`*_metrics.R` computes — bias, ATE bias, MSE, RMSE, MAE, relative
efficiency, correlation, Spearman, sign accuracy, the three HTE-test p-values,
and the `n_na` diagnostic. Figures go to `results/all_figures/missing/<study>/`.

These are the diagnostic view. `results_processing/thesis_figures/miss_*.R` is
the chapter view: the two or three panels that get printed, written to
`results/thesis_figures/`. The two never overwrite each other.

Labels, palette, the mean ± MCSE summary and the point-range panel are shared
via `R/figures.R` — `MISS_SCENARIO_LABELS`, `METHOD_LABELS`,
`MECHANISM_LABELS`, `STRATEGY_LABELS` and `point_range_plot()` in particular.
Rename an estimator or a handling method there and every figure follows.
