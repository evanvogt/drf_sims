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
that the covariates are correlated**, and that binary scenario 4's HTE is
smaller (scale 0.179 rather than 0.204, to make room for `U` - see "Binary `U`
has its own scale" below).

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
- **MNAR-Y0:** missingness driven only by an unobserved `U`, independent of X,
  which also shifts the control outcome, in both arms (continuous: `+ U`,
  `U ~ N(0, 1)`; binary: `+ 0.12·tanh(U)`, `U ~ N(0, 2²)`).
- **MNAR-τ** (`MNAR-tau` in the grid and paths): as MNAR-Y0, but `U` enters
  the treatment effect instead (same terms, treated arm only). Not run for
  scenario 1.

`bW` calibrated as in the complete-data study, ignoring `U`. Since 2026-09-29
`U` is drawn under **every** mechanism (after the covariate block, before the
continuous error) and used only under MNAR. So within a run the three
mechanisms share W, every covariate, the continuous error and X04/X05; the
MAR and MNAR outcomes differ by `U`'s term alone, and the comparison across
mechanisms is paired. Draw order: `W, Z-block, U, [err], cats`. (Before then
`U` was drawn only under MNAR, which shifted the error and X04/X05, so MAR and
MNAR were paired on W and X1–X5 / X01–X03 only.)

**Binary `U` has its own scale (since 2026-09-29).** A risk-difference effect
must keep every risk inside [0.01, 0.99], so `bU` is capped by the same
bounds as the HTE. `missing/binary` now has its own `RD_SCALE_MISS` rather
than the main study's `RD_SCALE`, derived at n = 500 alone. That frees the
floor room the main studies need at n = 100, and `U` spends it: `bU = 0.12`
(was 0.08) is the largest value that leaves scenarios 2, 3, 5 and 6 their
`sample_size/binary` HTE, and only scenario 4, where the ceiling binds, is
smaller (scale 0.179 against 0.204). `sU = 2` pushes `tanh(U)` towards ±1,
which strengthens the selection at no cost to the bounds (the amputation
standardises `U`). Both are re-derived by `sample_size/binary/bin_verify_hte.R`.
The binary mechanisms remain much weaker than the continuous ones, because a
Bernoulli outcome's own variance dominates - see the strength table below.

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

**What each mechanism implies.** From `missing/miss_dgm_checks.R`: one draw of
scenario 2 at n = 20,000 per outcome, bW calibrated at the studies' n = 500 so
the true ATE is the study's (continuous −0.27, binary RD −0.118). "CC" is
complete cases analysed with a correctly specified
`lm(Y ~ W * covariates)`, whose CATE is scored against the truth.

| | continuous | binary |
|---|---|---|
| AUC of incompleteness from X1–X5 (MAR / MNAR) | 0.71 / 0.51 | 0.71 / 0.51 |
| AUC of incompleteness from X01–X05, what IPW sees (MAR) | 0.66 | 0.66 |
| missing X4 minus observed X4, SD (MAR) | +0.30 | +0.36 |
| AUC of incompleteness from Y (MAR / MNAR-Y0 / MNAR-τ) | 0.54 / 0.66 / 0.58 | 0.52 / 0.54 / 0.52 |
| `U`'s share of Var(Y) (MNAR-Y0) | 45% | 4% |
| E[U-term \| complete], E[U-term \| incomplete] | −0.25, +0.60 | −0.02, +0.05 |
| CC ATE bias under MNAR-τ, % of the ATE | −89% | about −20% |

- **MAR:** units with higher X1–X5 are more often incomplete. With correlated
  covariates the *missing values themselves* are shifted, so mean imputation
  is biased and the imputation models have something to learn. Selection is
  on X only, so CC stays unbiased for the CATE.
- **MNAR-Y0:** incompleteness is independent of X but, through `U`, tied to
  the outcome in *both* arms. The shift is the same in both arms (W is
  randomised), so it differences out of the CATE: CC unbiased, and the
  imputation methods face MCAR-like X. **Predicted benign for the CATE** -
  missingness tied to prognosis costs power (U is unplanned outcome variance)
  and makes the missingness indicators prognostic, not biased. For the binary
  outcome it is close to MCAR.
- **MNAR-τ:** incompleteness is independent of X but tied to the treated
  arm's outcome only: the units with the largest individual effects are the
  most likely to be incomplete. Against the primary truth (averaged over `U`)
  CC is biased by a constant shift, E[U-term | complete]. No handling method
  that uses only X can remove it; one that uses the missingness itself can
  learn it (see "Estimands").

**Signatures: complete vs incomplete units.** The same check fits the
correctly specified lm to each handled dataset and scores its CATE on the
complete and the incomplete units separately. In continuous scenario 2 under
MAR the complete-unit RMSE is 0.02–0.12 for every method, but on the
incomplete units it is 0.50–0.66 for every imputation method - the floor from
scoring against τ at covariates they never saw - while CC has no incomplete
units to be scored on. Hence the split scoring under "Performance measures".
By-arm imputation with the outcome brings the slope of the CATE on the truth
to 0.90 (`multiple_imputation`) and 0.95 (`missforest`), against 0.86 for
mean imputation and 0.84 for the pre-2026-09-29 MI without Y and W.

Before 2026-09-28 (independent covariates, and an "MNAR" arm with `U` in
neither X nor Y) MAR's missing values were marginally distributed, the old
MNAR was MCAR, and only MNAR-Y (now MNAR-τ) could bias the CATE.

## Estimands

- **Unit-level CATE τ(X) at the true, unamputed covariates, averaged over
  `U`**, for every unit in the analysed sample (continuous: mean difference;
  binary: risk difference). This is the same target under every mechanism and
  every handling method: the CATE as a function of X.
- Presence of heterogeneity, for the tests (continuous / binary).

**What that target implies for incomplete units.** No method sees an
incomplete unit's missing effect modifiers, so none can recover its τ(X) -
at best E[τ(X) | what is observed]. On those units every imputation method
therefore carries an error floor that more data does not remove (about 0.5-0.7
RMSE in continuous scenario 2, against about 0.02 on complete units -
`missing/miss_dgm_checks.R`), while `complete_cases` and `IPW` are scored on
their complete units alone. That is why the performance measures are split by
completeness (below). Until 2026-09-29 this section described the target as
"the CATE given the observed-data covariates", which is not what the truth
was.

**MNAR-τ: two candidate targets.** Under MNAR-τ, `U` drives the missingness
and sits in the treatment effect, so an incomplete unit's expected effect is
τ(X) + E[U-term | incomplete] and a complete unit's τ(X) + E[U-term |
complete] (continuous: about +0.59 and −0.25; binary: about +0.05 and −0.02).
- **Primary (current): τ(X).** The population CATE as a function of X, the
  same target as under MAR and MNAR-Y0. Methods that can see whether a unit is
  complete (`missing_indicator`, `none`, and forests on mean-imputed point
  masses) are scored against a target that ignores information they
  legitimately use, so they can look worse for learning something real.
- **Secondary: τ_R = τ(X) + E[U-term | complete or incomplete].** The CATE
  given the covariates *and* the unit's completeness. Scored as `bias_r`,
  `ate_bias_r`, `rmse_r`, `rmse_r_cu`, `rmse_r_iu` on MNAR-τ rows. Under it,
  complete cases' "bias" is their failure to adjust for selection, and the
  methods that learn from missingness are rewarded; but the target then
  depends on the mechanism, so absolute errors no longer compare across
  mechanisms.
The constants come from `u_term_by_completeness()` (`R/missingness.R`), a
large simulation of the DGM. Both truths are computed at metrics time, so
which one the chapter leads with can be decided with no re-run.

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
| `missforest` | single random-forest imputation (`missForest`), within each arm, Y as a predictor | 500 |
| `multiple_imputation` | 50 imputations by `mice` random forests within each arm, Y as a predictor, analysed separately and pooled | 500 |
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
  This is deliberately the weak single-imputation baseline: no outcome, no
  other amputed covariate, no arm-specific model.
- **Imputing within arm, with the outcome (since 2026-09-29).** `missforest`
  and `multiple_imputation` impute **separately in each treatment arm**, with
  Y among the predictors (`impute_by_arm()`, `R/missingness.R`). An imputation
  model that leaves Y out imputes X independently of the outcome, which pulls
  every X-Y association, and so the heterogeneity, towards zero; a pooled
  model with Y but main effects only would still flatten the W×X
  interactions that are the CATE. Imputing by arm keeps both, and is the
  standard advice for imputing effect modifiers. W is constant within an arm,
  so it is not a predictor. Until 2026-09-29 both methods excluded Y and W:
  in continuous scenario 2 at large n that gave a slope of the CATE on the
  truth of 0.84 for `multiple_imputation`, against 0.90 now
  (`missing/miss_dgm_checks.R`).
- **`missforest`:** `missForest` on each arm's covariates plus Y, binary
  columns (X1, X3, and Y for a binary outcome) as factors so they are imputed
  by classification forests and stay 0/1. Package defaults: 100 trees, `mtry`
  ⌊√11⌋ = 3, up to 10 iterations, stopping at the first iteration whose
  change in the imputed values is larger than the previous one's. One
  completed dataset.
- **`multiple_imputation`:** `mice(m = 50, method = "rf")` in each arm, i.e.
  `mice.impute.rf` for every incomplete covariate (10 trees; each imputation
  is an observed donor value drawn from the matching leaves, so X1 and X3 stay
  0/1), default 5 iterations. The predictor matrix uses every other column of
  the arm, Y included, to impute each one. Imputation *i* of the two arms is
  recombined into the original row order. Each of the 50 completed
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
above (50 `mice` random-forest imputations within each arm, Y as a
predictor). Unlike the main design, `dr_oracle` is fit on each
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

- bias, ATE bias, relative ATE bias
- MSE, RMSE, MAE
- Pearson and Spearman correlation, sign accuracy
- HTE test rejection rates (`BLP_p`, `indep_cate`, `indep_po`); not yet pooled
  for `multiple_imputation`
- `n_na` (units with no estimate)

**Every point metric three ways** (`cate_metrics_split()`, `R/metrics.R`;
since 2026-09-29), using the missingness mask each run now saves:

| suffix | units | use |
|---|---|---|
| none | all analysed units (~350 for `complete_cases` / `IPW`, 500 otherwise) | not comparable across those two groups |
| `_cu` | complete units only | **the comparison across methods** - every method has them |
| `_iu` | incomplete units only (`NA` for `complete_cases` / `IPW`) | the error floor of scoring against τ at unseen covariates |

**Comparisons with `complete_data`** (same scenario, mechanism, run, model):

- `rel_efficiency_cu` = MSE_cu / MSE_cu of `complete_data` - the headline
  efficiency measure. `rel_efficiency` (all units) is `NA` for
  `complete_cases` and `IPW`, which analyse fewer units; `rel_efficiency_iu`
  on the incomplete units.
- `bias_diff_complete_cu` = bias_cu − bias_cu of `complete_data`: the bias the
  missingness and the handling add. A difference, not a ratio: the old
  `rel_bias_complete` divided by a complete-data bias that is often near zero.

**Comparing mechanisms.** The complete / incomplete split is the
amputation's, so it is a different set of units under each mechanism (MAR's
complete units have lower X1–X5, so a different spread of τ); and MNAR-Y0's
`U` adds outcome variance to the `complete_data` arm itself (continuous power
81% → 54%). So mechanisms are compared only through the within-mechanism
measures above, never through absolute `_cu` / `_iu` errors.

`rel_bias_cate`, the per-unit `(est − true) / true`, is still computed but not
reported: the true CATE crosses zero in scenarios 2–5.

**MNAR-τ secondary truth:** `bias_r`, `ate_bias_r`, `rmse_r`, `rmse_r_cu`,
`rmse_r_iu` - see "Estimands".

**ci_example:** marginal coverage (also on complete and incomplete units,
`_cu` / `_iu`), simultaneous coverage (per unit; nominal 0.95), mean interval
length, plus the point metrics above.
