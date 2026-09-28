# ADEMP — sample size studies

Five studies, sharing one set of data-generating mechanisms
(`R/dgm_scenarios.R`) and one set of estimators (`R/cate_models.R`):

| study | folder | question |
|---|---|---|
| continuous | `continuous/` | CATE estimation and HTE tests, continuous outcome |
| binary | `binary/` | the same, binary outcome (risk-difference scale) |
| CI, continuous | `confidence_intervals/continuous/` | bootstrap confidence bands for the CATE |
| CI, binary | `confidence_intervals/binary/` | the same, binary outcome |
| optimal_sf | `confidence_intervals/optimal_sf/` | data-driven choice of the bootstrap `sample.fraction`, both outcomes |

## Aims

- **continuous / binary:** assess how well DR-learner and causal forest
  estimators recover the CATE in a two-arm RCT, across CATE structures and
  sample sizes, and the size and power of post-estimation heterogeneity tests.
- **CI continuous / binary:** assess the coverage and width of half-sample
  bootstrap simultaneous confidence bands for the CATE, and how they depend on
  the bootstrap forests' `sample.fraction` (`CI_sf`).
- **optimal_sf:** choose `CI_sf` per dataset by a calibration that uses the
  point estimate as a plug-in truth.

## Data-generating mechanisms

### Common to both outcomes

- `W ~ Bernoulli(0.5)`
- Prognostic: `X1 ~ Bernoulli(0.4)`, `X2 ~ N(0, 1)`
- Potential effect modifiers (drawn in every scenario, whether or not τ uses
  them): `X3 ~ Bernoulli(0.7)`, `X4, X5 ~ N(0, 1)`
- Noise: `X01`–`X03 ~ N(0, 1)`; `X04`, `X05` indicators of a 3-level factor
  (0.45 / 0.30 / 0.25)
- `E[Y | x, W] = m0(x) + W·τ(x)`, `τ(x) = bW + g(x)`
- `bW` set so the true ATE gives 80% power in an unadjusted test at that n
  (`calibrate_bW()`)

### Continuous

- `Y = 0.4 − 0.5·X1 + X2 + W·τ(x) + ε`, `ε ~ N(0, 0.5²)`
- Two-sample t-test calibration: ATE −0.65, −0.41, −0.29, −0.20 at
  n = 100, 250, 500, 1000

### Binary

- `Y ~ Bernoulli(m0(x) + W·τ(x))`
- Control risk `m0(x) = 0.34 + 0.36·plogis(−1.72 − 0.5·X1 + X2)`; E[m0] = 0.40
- `g` is the continuous modifier with signs reversed, scaled by `RD_SCALE`,
  with X4 and X5 entering through `tanh` (t4 = tanh(X4), t5 = tanh(X5))
- Two-proportion test calibration: RD −0.248, −0.164, −0.118, −0.085 at
  n = 100, 250, 500, 1000
- Every treated risk stays within [0.01, 0.99]

### Scenarios

| scenario | description | continuous g(x) | binary g(x) |
|---|---|---|---|
| 1 | null | 0 | 0 |
| 2 | simple, continuous | −X4 | 0.082·t4 |
| 3 | complex | 2·X3 + 0.5·X4 − 0.5·X4·X5 | −0.102·X3 − 0.026·t4 + 0.026·t4·t5 |
| 4 | non-linear | cos(X4) | −0.204·cos(X4) |
| 5 | simple, binary | 2·X3 | −0.274·X3 |
| 6 | two variables | 0.3·X3 − X4 | −0.023·X3 + 0.075·t4 |
| 7 | binary × continuous | X3·X4 | −0.082·X3·t4 |
| 8 | single effects + interaction | 2·X3 + 0.5·X4 − 0.5·X3·X4 | −0.274·X3 − 0.069·t4 + 0.069·X3·t4 |
| 9 | continuous × continuous | −0.5·X4·X5 | 0.082·t4·t5 |
| 10 | exponential | 0.3·X3 + 0.1·exp(−\|X4\|) | −0.179·X3 − 0.060·exp(−\|X4\|) |

Scenarios 1–4 are the ones reported. Results saved before 2026-09-26 use the
old numbering (see the header of `R/dgm_scenarios.R`).

### Design

Full factorial over the rows below, per outcome type:

| | continuous / binary | CI continuous / binary | optimal_sf |
|---|---|---|---|
| scenarios | 1–10 | 1–10 | 1–10 |
| n | 100, 250, 500, 1000 | 500, 1000 | 500, 1000 |
| `CI_sf` | — | 0.05 to 0.5 by 0.05 (design factor) | chosen per run from the same values |
| repetitions | 100; 500 for scenarios 1–4 | 100 | 100 |

Each run is seeded by its run index alone (`setup_rng_stream(run)`), so a given
(scenario, n, run) is the same dataset in every study and every `CI_sf` cell.

## Estimands

- Unit-level CATE `τ(x_i)` for every unit in the simulated sample: a mean
  difference (continuous) or risk difference
  `P(Y=1 | x, W=1) − P(Y=1 | x, W=0)` (binary).
- ATE, as the sample mean of the CATE.
- Presence of heterogeneity (H0: τ constant), for the tests (continuous /
  binary).
- CI and optimal_sf studies: the CATE at a fixed covariate query grid
  (`build_query_grid()`) over the scenario's effect modifiers, all other
  covariates at 0; the same points in every run.
- optimal_sf calibration only: the estimated CATE `tau.hat` (plug-in target).

## Methods

### Estimators

| method | description | used in |
|---|---|---|
| `causal_forest` | grf `causal_forest`, internal cross-fitting | all but optimal_sf |
| `dr_random_forest` | DR-learner; outcome model fit per arm (T-learner RF, own-arm predictions OOB), RF propensity, OOB RF second stage | all |
| `dr_superlearner` | DR-learner; SuperLearner per-arm outcome models and propensity, single leave-one-fold-out crossfit (V = 4, 5, 10 at n = 100, 250, ≥ 500), SuperLearner second stage | continuous / binary |
| `dr_oracle` | DR-learner with the true outcome model (binary: true risk) and propensity 0.5 | all but optimal_sf |
| `dr_semi_oracle` | DR-learner with the known propensity 0.5 and per-arm outcome forests specified as in `dr_random_forest` (refitted) | all but optimal_sf |
| `t_random_forest` | T-learner $\hat\mu_1 - \hat\mu_0$ from `dr_random_forest`'s per-arm outcome forests (own-arm predictions OOB); derived at metrics time, not refitted | continuous / binary |
| `t_superlearner` | T-learner $\hat\mu_1 - \hat\mu_0$ from `dr_superlearner`'s per-arm SuperLearners (out-of-fold); derived at metrics time, not refitted | continuous / binary |

Estimated propensities in the DR-learners (`dr_random_forest`,
`dr_superlearner`) are trimmed to [0.05, 0.95]; the causal forest's internal
propensity is not.

### Model fitting

Implementation: `R/cate_models.R::cate_methods()`, called through each study's
`*_models.R` shim. The binary and continuous studies differ only in the
SuperLearner `family` (see below).

**Inputs.** Every learner receives the 10 covariates X1–X5, X01–X05 as they
are generated (no scaling or feature engineering, except the interaction lasso
below). W enters the causal forest as the treatment; it is never a covariate
of an outcome model, which is fit per arm.

**DR pseudo-outcome** (`dr_pseudo`), shared by every DR-learner:

$$\hat\phi_i = \hat\mu_1(x_i) - \hat\mu_0(x_i) + \frac{(Y_i - \hat\mu_{W_i}(x_i))(W_i - \hat e(x_i))}{\hat e(x_i)(1 - \hat e(x_i))}$$

The CATE is the second-stage regression of $\hat\phi$ on X.

**Forests.** All forests are grf `regression_forest` / `causal_forest` at
grf's defaults: 2000 trees, `sample.fraction` 0.5, honest splitting
(`honesty.fraction` 0.5, `honesty.prune.leaves` TRUE), `min.node.size` 5,
`mtry` = min(⌈√p⌉ + 20, p) = all 10 covariates, `alpha` 0.05,
`ci.group.size` 2, no parameter tuning. Only the CI studies' bootstrap forests
change a setting (`sample.fraction = CI_sf`).

**Per estimator:**

- `causal_forest`: `causal_forest(X, Y, W)` with no supplied nuisances, so grf
  fits `Y.hat` and `W.hat` itself as OOB predictions from regression forests of
  500 trees (max(50, 2000/4)); τ̂ is the forest's OOB prediction. Splits use
  grf's default `stabilize.splits = TRUE`.
- `dr_random_forest` (whole sample, no sample splitting):
  - $\hat\mu_0, \hat\mu_1$: one regression forest of Y on X per arm, fit on that
    arm's rows. A unit's own-arm prediction is OOB; its other-arm prediction is
    from a forest that never saw it.
  - $\hat e$: regression forest of W on X, OOB predictions, trimmed.
  - Stage 2: regression forest of $\hat\phi$ on X; τ̂ is its OOB prediction,
    with grf's OOB variance estimate saved.
- `dr_oracle`: $\mu_w(x)$ is the true outcome mean evaluated at W = w (the
  DGM's oracle formula, `get_*_oracle_info()`); e = 0.5. Stage 2 as
  `dr_random_forest`.
- `dr_semi_oracle`: per-arm outcome forests as in `dr_random_forest` (a separate
  fit); e = 0.5. Stage 2 as `dr_random_forest`.
- `dr_superlearner` (single cross-fit; V = 4, 5, 10 folds at n = 100, 250,
  ≥ 500; rows assigned to folds in contiguous blocks):
  - Stage 1, per fold k: on the rows outside k, one SuperLearner of Y on X per
    arm and one of W on X; predict fold k's rows; trim $\hat e$; form
    $\hat\phi$ for fold k.
  - Stage 2, per fold k: SuperLearner of $\hat\phi$ on X fit on the rows outside
    k, predicting fold k. Stage 1 and stage 2 use the same folds.
- `t_random_forest`, `t_superlearner`: $\hat\mu_1 - \hat\mu_0$, recovered at
  metrics time from the saved DR nuisances by removing the residual term from
  $\hat\phi$ (`R/metrics.R::add_t_learners`).

**SuperLearner settings** (`R/sl_library.R::sl_fit_predict`):

| model | family | meta-learner |
|---|---|---|
| propensity | binomial | `method.NNloglik` |
| outcome, continuous | gaussian | `method.NNLS` |
| outcome, binary | binomial | `method.NNloglik` |
| CATE (stage 2) | gaussian | `method.NNLS` |

- Ensemble weights from SuperLearner's internal 10-fold CV (package default),
  within each cross-fitting training set.
- Before each fit, each candidate is fit alone with 2-fold CV
  (`pretest_superlearner`); learners that error or give non-finite predictions
  are dropped and recorded in `sl_dropped`. If all fail, the library is
  `SL.mean` alone.
- If a SuperLearner fit itself errors, every prediction is the mean of its
  training outcome, with a warning, recorded in `sl_dropped` as `(whole fit)`.
- If an outcome or propensity fit returns all-zero predictions, they are
  replaced by the training mean (for the outcome model, that arm's mean), with
  a warning.

SuperLearner libraries, one per nuisance (`R/sl_library.R::sl_libraries`),
chosen without reference to the scenarios' HTE shapes. The outcome model is fit
in each arm separately, so its library is sized to about half the training
rows (~37 per arm at n = 100):

| | n = 100 | n = 250 | n ≥ 500 |
|---|---|---|---|
| propensity | mean, glm | same | same |
| outcome (per arm) | mean, lasso, ranger (min.node.size 25) | + glm, gam | + earth |
| CATE (stage 2) | mean, glm, lasso (lambda.min and lambda.1se), gam, ranger (min.node.size 25) | + lasso on pairwise interactions (both tunings), earth | same |

`family = binomial()` for the binary outcome model. Each stage-2 lasso is
included at both tunings and SuperLearner's CV weighs them. Learners that error
or give non-finite predictions on a fold are dropped before fitting
(`pretest_superlearner`) and saved as `sl_dropped`.

Candidate learners (SuperLearner wrapper defaults unless stated; p = 10):

| learner | definition |
|---|---|
| `SL.mean` | training mean |
| `SL.glm` | `glm` on all main effects |
| `SL.glmnet` | lasso (`alpha` 1) on main effects, λ chosen by 10-fold `cv.glmnet` at `lambda.min` |
| `SL.glmnet.1se` | as `SL.glmnet`, at `lambda.1se` (custom wrapper) |
| `SL.glmnet.int` | lasso on all main effects and pairwise interactions (`model.matrix(~ .^2)`), `lambda.min` (custom wrapper) |
| `SL.glmnet.int.1se` | as `SL.glmnet.int`, at `lambda.1se` (custom wrapper) |
| `SL.gam` | `gam` package; smoothing spline with 2 df for covariates with more than 4 unique values, linear otherwise |
| `SL.earth` | MARS: `degree` 2, `penalty` 3, `nk` = max(21, 2p + 1), backward pruning |
| `SL.ranger.ns25` | ranger, 500 trees, `mtry` ⌊√p⌋ = 3, `min.node.size` 25 (custom wrapper; the default is 5); probability forest for a binomial outcome |

### Heterogeneity tests (continuous / binary)

Per method, on the full sample (`run_blp_whole`,
`run_independence_test_whole`). Also run on the true CATE and nuisances
(`cts_true_cate_tests.RDS`, `bin_true_cate_tests.RDS`).

- **BLP** (GenericML `BLP`, Chernozhukov et al.): weighted least squares, weights
  $1/(\hat e(1-\hat e))$, of Y on an intercept, the baseline proxy
  $\hat\mu_0(x)$, $(W - \hat e)$ and $(W - \hat e)(\hat\tau(x) - \bar{\hat\tau})$.
  The test is on the last coefficient, β2. $\hat e$ and $\hat\mu_0$ are the
  estimator's own stage-1 nuisances: `nuisances_rf` for `causal_forest`,
  `dr_random_forest` and `t_random_forest`; `nuisances_sl` for
  `dr_superlearner` and `t_superlearner`; e = 0.5 and the arm's own
  $\mu_0$ for the oracles.
  - `BLP_p`: two-sided, homoskedastic OLS standard errors (saved at
    estimation time).
  - `BLP_p_os`: one-sided (H1: β2 > 0), HC3 standard errors (recomputed at
    metrics time from the saved nuisances).
- **Independence test** (coin `independence_test(τ ~ X, teststat =
  "quadratic")`, asymptotic χ² reference): of the estimated CATE against the
  covariates (`indep_cate`), and of the pseudo-outcome against the covariates
  (`indep_po`; reported once per pseudo-outcome, so not for `causal_forest` or
  the T-learners).
- A CATE that is constant (up to rounding) gives NA for both tests.

### Intervals (CI studies)

No SuperLearner arm and no HTE tests. The point estimators are fit exactly as
above; since `family` only reaches SuperLearner, the continuous and binary CI
studies fit the same models. The causal forest also saves grf's variance
estimate (for `causal_forest_inbuilt`).

- **Half-sample bootstrap band** (`R/bootstrap_ci.R`), B = 200, α = 0.05: per
  draw, refit the second stage (DR-learners: regression forest on the fixed
  pseudo-outcomes; causal forest: with `Y.hat`, `W.hat` fixed) on an
  unstratified half sample, with `sample.fraction = CI_sf`. Roots
  `tau_full − tau_half` (in-half units OOB, the rest out of sample) are
  standardised by their bootstrap SD; critical value = `1 − α/2` quantile of
  the per-draw maximum over units. Built over the sample's units and,
  separately, over the query grid.
- **`causal_forest_inbuilt`:** pointwise normal interval from grf's variance
  estimate, α = 0.05.

### Calibration (optimal_sf)

`dr_random_forest` only. `find_optimal_sf()`: for each candidate `CI_sf`, 50
times, resample the pseudo-outcome residuals around `tau.hat` (whole sample,
with replacement), re-estimate the CATE, build the band (B = 100) and record
its coverage of `tau.hat`. Pick the candidate whose mean coverage is closest
to `1 − α`, then build the final band on the observed data at that value
(B = 200). α = 0.10 here, so the target is 0.90.

## Performance measures

### continuous / binary

Per run, over the units in the sample (`R/metrics.R::cate_metrics`), then
averaged over runs:

- bias (mean of estimate − truth), ATE bias, relative ATE bias, relative CATE bias
- MSE, RMSE, MAE
- Pearson and Spearman correlation with the true CATE (set to 0 in scenario 1)
- sign accuracy
- HTE tests: p-values `BLP_p`, `BLP_p_os`, `indep_cate`, `indep_po` → rejection rates at
  0.05 (type I error in scenario 1, power in 2–10)

### CI continuous / binary

Per run (`R/metrics.R::interval_metrics`), then averaged over runs. Nominal
coverage 0.95.

- marginal coverage (proportion of units covered by their own interval)
- simultaneous coverage (band covers every unit; the target the method controls)
- mean interval length
- query grid only: bias-eliminated marginal and simultaneous coverage
  (coverage of the across-run mean estimate)

### optimal_sf

No metrics script yet. Each run saves the selected `CI_sf`, the calibration
curve (mean coverage of `tau.hat` and mean width per candidate) and the final
band over units and query grid alongside the truth, so coverage of the true
CATE (nominal 0.90) can be scored as in the CI studies.
