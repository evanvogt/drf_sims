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

`bW` calibrated as in the complete-data study, ignoring `U`. `U` is drawn
(after X5, before the continuous error) only under MNAR and MNAR-Y, so within
a run the complete data under MAR differ from those under the MNAR mechanisms.

### Amputation

Implementation: `R/missingness.R::introduce_missingness()`, one call to
`mice::ampute` per dataset, after the complete data are generated.

**Variables.** Y, W and the noise covariates X01–X05 are set aside before
`ampute` is called: they are never made missing and, under MAR, do not drive
the missingness either. Under `type = "both"` all of X1–X5 can be missing.
Under MNAR / MNAR-Y, `U` is appended as an extra, always-observed column for
the call and dropped afterwards; it is never seen by the handling methods or
estimators.

**Patterns.** Every observed/missing combination of X1–X5 except all-observed:

| mechanism | patterns | note |
|---|---|---|
| MAR | 30 | all-missing also dropped: a pattern must leave something observed to weight |
| MNAR, MNAR-Y | 31 | all-missing kept, since `U` stays observed in it |

Each covariate is missing in about half the patterns, so each is missing for
about 15% of units, and 30% of units have at least one missing covariate
(complete-case n ≈ 350).

**`ampute` settings** (package defaults unless stated):

- `prop = 0.3`, `bycases = TRUE`: 30% of *units* are incomplete.
- `freq`: equal, each unit is first allocated to one pattern with probability
  1/30 (MAR) or 1/31 (MNAR).
- `std = TRUE`: variables standardised before scoring.
- `cont = TRUE`, `type = "RIGHT"`: within its pattern, a unit becomes
  incomplete with probability logistic in its weighted sum score, shifted so
  the overall proportion is 0.3; high scores are the most likely to be missing.
- `mech = "MAR"` in every call; the mechanism is set by the weights:

| mechanism | weighted sum score |
|---|---|
| MAR | `ampute` default: weight 1 on each of X1–X5 *observed* in the pattern, 0 on those made missing |
| MNAR, MNAR-Y | user-supplied: weight 1 on `U`, 0 on X1–X5, in every pattern |

**What each mechanism implies** (checked over 20 draws of scenario 3):

- **MAR:** units with higher X1–X5 are more often incomplete (correlation of
  incompleteness with each covariate ≈ 0.11–0.14). Since X1 and X2 are
  prognostic, incompleteness is also correlated with Y (≈ 0.12).
- **MNAR:** `U` enters neither the covariates nor the outcome, so
  incompleteness is independent of X, W and Y (all correlations ≈ 0): in
  Rubin's taxonomy this arm is MCAR, with `U` only a device for generating it.
- **MNAR-Y:** incompleteness is independent of X but, through `U`, correlated
  with the treated arm's outcome (≈ 0.10 with Y overall): units with the
  largest individual treatment effects are the most likely to be incomplete.
  Not correlated with the estimand, which averages over `U`.

## Estimands

- Unit-level CATE given the observed-data covariates (averaged over `U`), for
  every unit in the analysed sample (continuous: mean difference; binary: risk difference).
- Presence of heterogeneity, for the tests (continuous / binary).

## Methods

**continuous / binary:** handling method × estimator.

### Missing-data handling

Implementation: `R/missingness.R::handle_missingness()`, applied to the
amputed dataset before any estimator sees it. Every method except
`complete_data` sees the same amputed dataset within a (scenario, mechanism,
run), since each run is seeded by its run index alone.

| method | what is done | analysed n |
|---|---|---|
| `complete_data` | no amputation; the complete-data reference | 500 |
| `complete_cases` | drop units with any missing covariate | ≈ 350 |
| `IPW` | complete cases, weighted by 1 / P̂(complete) | ≈ 350 |
| `mean_imputation` | each missing value replaced by its column's observed mean | 500 |
| `missing_indicator` | mean imputation, plus a 0/1 indicator per amputed covariate | 500 |
| `regression` | single deterministic regression imputation | 500 |
| `missforest` | single random-forest imputation (`missForest`) | 500 |
| `multiple_imputation` | 50 imputations by `mice` random forests, analysed separately and pooled | 500 |
| `none` | NAs passed to the estimator | 500 |

Per method:

- **`complete_cases`:** `complete.cases()` over all columns; the truth is
  subset to the same units, so the CATE is scored on the retained units only.
- **`IPW`:** a logistic regression of the complete-case indicator on the
  *fully observed* covariates, fit on all 500 units; complete cases are then
  weighted by the inverse of their fitted probability. Under `both` the fully
  observed covariates are X01–X05 alone, which are unrelated to missingness
  under every mechanism, so the weights are close to constant (≈ 1 / 0.7).
  The weights reach every fit that accepts them: grf `sample.weights` in every
  forest (per-arm outcome forests, propensity forest, causal forest, stage-2
  forests) and SuperLearner `obsWeights` in every SuperLearner fit. Truth
  subset as `complete_cases`.
- **`mean_imputation`:** column means over observed values, so the binary X1
  and X3 are imputed with a fraction (their observed prevalence).
- **`missing_indicator`:** as `mean_imputation`, plus `X1_missing`–
  `X5_missing` appended as covariates, so every estimator sees 15 covariates
  instead of 10.
- **`regression`:** for each incomplete covariate, `VIM::regressionImp` fits
  a linear model (`lm`; X1 and X3 are numeric, so linear too) of that
  covariate on the *fully observed* covariates among the observed units and
  fills in its predictions, with no residual noise. Y and W are excluded from
  the imputation model. Under `both` the predictors are X01–X05 alone, which
  are independent of X1–X5, so the imputations are close to the observed means.
- **`missforest`:** `missForest` on the ten covariates (Y and W excluded from
  the imputation model), binary covariates as factors so they are imputed by
  classification forests and stay 0/1. Package defaults: 100 trees, `mtry`
  ⌊√10⌋ = 3, up to 10 iterations, stopping at the first iteration whose change
  in the imputed values is larger than the previous one's. One completed
  dataset.
- **`multiple_imputation`:** `mice(m = 50, method = "rf")`, i.e.
  `mice.impute.rf` for every incomplete covariate (10 trees; each imputation
  is an observed donor value drawn from the matching leaves, so X1 and X3 stay
  0/1), default 5 iterations. The predictor matrix uses every other covariate
  to impute each one, and sets the Y and W columns to 0, so **the imputation
  model omits the outcome and the treatment**. Each of the 50 completed
  datasets is analysed separately (no `dr_oracle`, no `dr_superlearner`); per
  unit, the pooled point estimate is the mean of the 50 CATE estimates and the
  pooled variance is Rubin's `W̄ + (1 + 1/50)·B`, with W̄ the mean of grf's
  per-imputation variance estimates (`R/cate_models.R::combine_mi`). HTE tests
  are run per imputation and saved unpooled (`mi_test_table`).
- **`none`:** the amputed data go straight to the estimators. Only the grf
  estimators accept NAs (grf splits on missingness natively, missing
  incorporated in attributes), so only `causal_forest` and `dr_random_forest`
  are run.

### Estimators

`causal_forest`, `dr_random_forest`, `dr_superlearner`, `dr_oracle`,
`dr_semi_oracle`, fit exactly as in `sample_size/` (see
`sample_size/ADEMP.md`), with the n ≥ 500 SuperLearner libraries and V = 10.
Covariates are every column after Y and W of the handled dataset. The
T-learners are not derived here.

`dr_oracle` evaluates the true outcome-mean formula at the covariates *as
handled*: the imputed values for the imputation methods, so it isolates the
nuisance-model error but not the imputation error.

| method | `causal_forest` | `dr_random_forest` | `dr_superlearner` | `dr_oracle` | `dr_semi_oracle` |
|---|---|---|---|---|---|
| `complete_data`, `complete_cases`, `IPW`, single imputations | ✓ | ✓ | ✓ | ✓ | ✓ |
| `multiple_imputation` | ✓ (pooled) | ✓ (pooled) | — | — | ✓ (pooled) |
| `none` | ✓ | ✓ | — | — | — |

### ci_example

`multiple_imputation` only, MAR only: amputation and imputation exactly as
above (50 `mice` random-forest imputations, Y and W excluded from the
imputation model). Unlike the main design, `dr_oracle` is fit on each
imputation (at the imputed covariates). Per imputation, half-sample bootstrap simultaneous bands (`CI_boot = 200`, `CI_sf = 0.5`,
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
