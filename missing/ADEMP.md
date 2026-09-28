# ADEMP — missing covariates

Three studies: `continuous/` and `binary/` (main design) and `ci_example/`
(intervals after multiple imputation).

## Aims

- **continuous / binary:** assess how covariate missingness, under different
  mechanisms, affects CATE estimation, and compare missing-data handling methods.
- **ci_example:** compare ways of pooling bootstrap confidence intervals for the
  CATE across multiply imputed datasets.

## Data-generating mechanisms

Complete data: the `sample_size/` continuous or binary DGM (see `sample_size/ADEMP.md`),
scenarios numbered as there.

| | continuous / binary | ci_example (continuous) |
|---|---|---|
| scenarios | 1–5 (no HTE; X4; X3 + X4 + X4·X5; cos(X4); X3) | 1, 3, 4, 5, 6 (6: X3 + X4) |
| n | 500 | 500 |
| amputed covariates | `both`: X1–X5 in every scenario, effect modifiers or not (X01–X05 never) | `both` |
| proportion incomplete | 0.3 (`mice::ampute`) | 0.3 |
| mechanism | MAR, MNAR, MNAR-Y | MAR |
| repetitions | 100 | 100 |

Mechanisms:

- **MAR:** missingness depends on the observed covariates.
- **MNAR:** missingness driven only by an unobserved `U ~ N(0, 1)`.
- **MNAR-Y:** as MNAR, and `U` also enters the treatment effect
  (continuous: `+ U`; binary: `+ 0.08·tanh(U)`). Not run for scenario 1.

`bW` calibrated as in the complete-data study, ignoring `U`.

## Estimands

- Unit-level CATE given the observed-data covariates (averaged over `U`), for
  every unit in the analysed sample (continuous: mean difference; binary: risk difference).
- Presence of heterogeneity, for the tests (continuous / binary).

## Methods

**continuous / binary:** handling method × estimator.

Handling methods: `complete_cases`, `mean_imputation`, `missforest`,
`regression` imputation, `missing_indicator`, `IPW` (complete cases,
inverse-probability weighted), `multiple_imputation` (50 imputations, Rubin's
rules), `none` (estimator handles NAs), and `complete_data` (no missingness;
reference).

Estimators: `causal_forest`, `dr_random_forest`, `dr_superlearner`,
`dr_oracle`, `dr_semi_oracle` (as `continuous/`, V = 10). Arms needing a
complete covariate matrix are skipped under `none`; under
`multiple_imputation` only `causal_forest`, `dr_random_forest` and
`dr_semi_oracle` are pooled.

**ci_example:** `multiple_imputation` only (50 imputations). Per imputation,
half-sample bootstrap simultaneous bands (`CI_boot = 200`, `CI_sf = 0.5`,
α = 0.05) for `causal_forest`, `dr_random_forest`, `dr_oracle`,
`dr_semi_oracle`, pooled three ways:

| strategy | pooling |
|---|---|
| `pooled` | empirical quantiles of bootstrap replicates stacked across imputations |
| `mib` | Rubin's rules variance; critical value averaged over imputations |
| `hybrid` | one variance and one critical value from the stacked draws |

## Performance measures

**continuous / binary:**

- bias, ATE bias, relative biases
- MSE, RMSE, MAE
- Pearson and Spearman correlation, sign accuracy
- relative efficiency (MSE / MSE of `complete_data`) and bias relative to `complete_data`
- HTE test rejection rates (`BLP_p`, `indep_cate`, `indep_po`); not yet pooled
  for `multiple_imputation`
- `n_na` (units with no estimate)

**ci_example:** marginal coverage, simultaneous coverage (per unit; nominal
0.95), mean interval length, plus the point metrics above.
