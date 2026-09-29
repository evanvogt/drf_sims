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
scenarios numbered as there - the same outcome model, CATE and baseline - **except
that the covariates are correlated**.

| | continuous / binary | ci_example (continuous) |
|---|---|---|
| scenarios | 1–5 (no HTE; X4; X3 + X4 + X4·X5; cos(X4); X3) | 1, 3, 4, 5, 6 (6: X3 + X4) |
| n | 500 | 500 |
| covariates | copula, exchangeable latent ρ = 0.5; X01–X03 auxiliary | as main |
| amputed covariates | `both`: X1–X5 in every scenario, effect modifiers or not (X01–X05 never) | `both` |
| proportion incomplete | 0.3 (`mice::ampute`) | 0.3 |
| mechanism | MAR, MNAR-Y0, MNAR-τ | MAR |
| repetitions | 100 | 100 |

**Covariates.** A latent `Z ~ N(0, R)` over (X1, X2, X3, X4, X5, X01, X02,
X03), R exchangeable with ρ = 0.5 (`rho` in the scenario table). X1 and X3
are thresholds of their latent columns at the main study's prevalences (0.4,
0.7); X2, X4, X5 are scaled latent columns (SD 1); X01–X03 are the latent
columns themselves. X04/X05 (the categorical pair) are drawn independently as
in the main study. So:

- the prognostic X1, X2 and the effect modifiers X3–X5 are correlated with one
  another: on the latent scale each is 44% predictable (R²) from the other
  seven;
- X01–X03 are **auxiliary**: in neither the outcome model nor the CATE, never
  amputated, but correlated with X1–X5 (a latent R² ≈ 0.375 from the three), so
  they carry information for `regression` imputation and the `IPW` model;
- the truth is unchanged as a function: τ(X3, X4, X5) and m0(X1, X2) as in the
  main study, and X01–X03 add nothing to either given X1–X5.

Correlation changes E[g] only where the CATE multiplies two modifiers
(scenario 3's X4·X5), and the planned outcome SD / control event rate through
Cov(X1, X2); `calibrate_bW()` accounts for both (`te_grid()` and
`baseline_grid()` integrate over the correlated latent), so the true ATE is
still the planned effect. Before 2026-09-28 every covariate was independent:
imputation then had nothing to impute from, and MAR selected on X only, which
leaves the CATE unbiased - see "What each mechanism implies".

Mechanisms:

- **MAR:** missingness depends on the observed covariates.
- **MNAR-Y0:** missingness driven only by an unobserved `U ~ N(0, 1)`,
  independent of X, which also shifts the control outcome
  (continuous: `+ U`; binary: `+ 0.08·tanh(U)`, in both arms).
- **MNAR-τ** (`MNAR-tau` in the grid and paths): as MNAR-Y0, but `U` enters
  the treatment effect instead (same terms, treated arm only). Not run for
  scenario 1.

`bW` calibrated as in the complete-data study, ignoring `U`. `U` is drawn
(after the covariate block, before the continuous error) only under the MNAR
mechanisms, so within a run the complete data under MAR differ from those
under the MNAR mechanisms. Draw order: `W, Z-block, [U], [err], cats`.

These mechanisms replaced MAR / MNAR / MNAR-Y on 2026-09-28. The old MNAR had
`U` in neither X nor Y, so it was MCAR; the old MNAR-Y is MNAR-τ.

### Amputation

Implementation: `R/missingness.R::introduce_missingness()`, one call to
`mice::ampute` per dataset, after the complete data are generated.

**Variables.** Y, W and X01–X05 (the auxiliaries X01–X03 and the noise pair
X04/X05) are set aside before
`ampute` is called: they are never made missing and, under MAR, do not drive
the missingness either. Under `type = "both"` all of X1–X5 can be missing.
Under MNAR-Y0 / MNAR-τ, `U` is appended as an extra, always-observed column for
the call and dropped afterwards; it is never seen by the handling methods or
estimators.

**Patterns.** Every observed/missing combination of X1–X5 except all-observed:

| mechanism | patterns | note |
|---|---|---|
| MAR | 30 | all-missing also dropped: a pattern must leave something observed to weight |
| MNAR-Y0, MNAR-τ | 31 | all-missing kept, since `U` stays observed in it |

Each covariate is missing in about half the patterns, so each is missing for
about 15% of units, and 30% of units have at least one missing covariate
(complete-case n ≈ 350).

**`ampute` settings** (package defaults unless stated):

- `prop = 0.3`, `bycases = TRUE`: 30% of *units* are incomplete.
- `freq`: equal, each unit is first allocated to one pattern with probability
  1/30 (MAR) or 1/31 (MNAR-Y0, MNAR-τ).
- `std = TRUE`: variables standardised before scoring.
- `cont = TRUE`, `type = "RIGHT"`: within its pattern, a unit becomes
  incomplete with probability logistic in its weighted sum score, shifted so
  the overall proportion is 0.3; high scores are the most likely to be missing.
- `mech = "MAR"` in every call; the mechanism is set by the weights:

| mechanism | weighted sum score |
|---|---|
| MAR | `ampute` default: weight 1 on each of X1–X5 *observed* in the pattern, 0 on those made missing |
| MNAR-Y0, MNAR-τ | user-supplied: weight 1 on `U`, 0 on X1–X5, in every pattern |

**What each mechanism implies** (one draw of scenario 2 at n = 20,000, both
outcomes; "CC" is a complete-case `lm(Y ~ W·X4 + X1 + X2 + X3 + X5)`, which is
correctly specified for the continuous CATE):

- **MAR:** units with higher X1–X5 are more often incomplete (correlation of
  incompleteness with each covariate ≈ 0.27). With correlated covariates the
  *missing values themselves* are shifted: the missing X4 values average
  0.39 SD above the observed ones, so mean imputation is biased and the
  imputation models have something to learn. Selection is on X only, so CC
  stays unbiased for the CATE (continuous: W and W·X4 within 1 SE).
- **MNAR-Y0:** incompleteness is independent of X (|correlation| ≤ 0.01) but,
  through `U`, correlated with the outcome in *both* arms (continuous ≈ 0.26;
  binary ≈ 0.02–0.04). The shift is the same in both arms (W is randomised),
  so it differences out of the CATE: CC unbiased (continuous: W −0.013 ± 0.018,
  W·X4 +0.010 ± 0.018), and since X is independent of missingness the
  imputation methods face MCAR-like X. **Predicted benign for the CATE** -
  missingness tied to prognosis costs power (U is unplanned outcome variance)
  and makes the missingness indicators prognostic, not biased.
- **MNAR-τ:** incompleteness is independent of X but correlated with the
  treated arm's outcome only (continuous ≈ 0.27, control ≈ 0): the units with
  the largest individual effects are the most likely to be incomplete. The
  estimand averages over `U`, so CC is biased by a constant shift -
  continuous W −0.27 ± 0.01 against a true ATE of −0.27, binary −0.020 against
  an RD of −0.118. No handling method that uses only X can remove it.

Before 2026-09-28 (independent covariates, and an "MNAR" arm with `U` in
neither X nor Y) MAR's missing values were marginally distributed, the old
MNAR was MCAR, and only MNAR-Y (now MNAR-τ) could bias the CATE.

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
  observed covariates are X01–X05 alone. Under MAR the auxiliaries X01–X03 are
  correlated with the X1–X5 values that drive missingness, so the model is
  informative but misspecified (ampute's pattern-wise score is not a logistic
  function of X01–X05). Under MNAR-Y0 / MNAR-τ missingness depends only on
  `U`, independent of every covariate, so the weights are close to constant
  (≈ 1 / 0.7). (Before the 2026-09-28 correlated covariates they were close to
  constant under every mechanism.)
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
  the imputation model. Under `both` the predictors are X01–X05 alone: the
  auxiliaries X01–X03 explain part of each amputed covariate (latent
  R² ≈ 0.375), X04/X05 nothing. (Before the 2026-09-28 correlated covariates
  none of them did, and the imputations were close to the observed means.)
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
