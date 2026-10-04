# ADEMP — correlated covariates

Six studies, all on one data-generating mechanism: scenarios 1–4 of the
sample-size studies, with the covariates drawn from a Gaussian copula at
latent correlation ρ ∈ {0, 0.5}. Each study is a parent sample-size study rerun
on this DGM. Its estimators, intervals, tests and metrics are the parent's
code, unchanged. Only the covariates' joint distribution differs.

| study | folder | prefix | parent | question |
|---|---|---|---|---|
| correlated continuous | `continuous/` | `cts_corr_*` | `sample_size/continuous/` | CATE estimation and HTE tests, continuous outcome |
| correlated binary | `binary/` | `bin_corr_*` | `sample_size/binary/` | the same, binary outcome (risk-difference scale) |
| correlated CI, continuous | `confidence_intervals/continuous/` | `cts_corr_ci_*` | `sample_size/confidence_intervals/continuous/` | bootstrap simultaneous bands for the CATE across `CI_sf` |
| correlated CI, binary | `confidence_intervals/binary/` | `bin_corr_ci_*` | `sample_size/confidence_intervals/binary/` | the same, binary outcome |
| correlated optimal_sf, continuous | `confidence_intervals/optimal_sf/` | `cts_corr_ci_sf_*` | `sample_size/confidence_intervals/optimal_sf/` | data-driven choice of `CI_sf`, scored against the true CATE |
| correlated optimal_sf, binary | `confidence_intervals/optimal_sf/` | `bin_corr_ci_sf_*` | the same | the same, binary outcome |

Folders are relative to `sample_size/correlated/`. The cts/bin studies were
first run on 2026-10-01, and the CI and optimal_sf studies on 2026-10-04
(`R/study_registry.R`).

## Aims

In the parent studies every covariate is independent, so prognosis (m0) and
effect modification (τ) never move together. Here they do. The aim of each
study is to measure, within each ρ, the same performance its parent measures,
and then how that performance changes from ρ = 0 to ρ = 0.5.

- **Correlated continuous / binary:**
  - how well DR-learner, T-learner and causal forest estimators recover the
    CATE in a two-arm RCT, across scenarios 1–4 and n = 100 to 1000;
  - the size and power of the post-estimation heterogeneity tests, on both the
    estimated and the true CATE;
  - whether separating prognosis from effect modification gets harder when the
    two are correlated.
- **Correlated CI continuous / binary:** the coverage and width of half-sample
  bootstrap simultaneous confidence bands for the CATE, across the bootstrap
  forests' `sample.fraction` (`CI_sf`), and whether this changes with ρ.
- **Correlated optimal_sf continuous / binary:** whether a per-dataset choice of
  `CI_sf` still reaches nominal 0.90 coverage of the *true* CATE when m0
  tracks τ. The choice is made by a calibration that uses the point estimate as
  a plug-in truth, so it cannot see that estimate's bias, and at ρ = 0.5
  that bias is where any growth is expected.

## Data-generating mechanisms

The implementation is in `R/dgm_scenarios.R` ("CORRELATED SAMPLE-SIZE SETS").
The scenario sets are `continuous_corr_<ρ>` and `binary_corr_<ρ>`
(`corr_set()`), one per ρ in `CORR_RHOS`.

### Outcome model

Both outcomes use

`E[Y | x, W] = m0(x) + W·τ(x)`, `τ(x) = bW + g(x)`, `W ~ Bernoulli(0.5)`.

X1 and X2 are purely prognostic. X3, X4 and X5 are the potential effect
modifiers, and each scenario's g uses some of them. All ten covariates
X1–X5 and X01–X05 are drawn and returned in every scenario, and every
estimator sees all ten.

**Continuous.**

- `Y = 0.4 − 0.5·X1 + X2 + W·τ(x) + ε`, `ε ~ N(0, 0.5²)`.

**Binary.** The effect is on the risk-difference scale.

- `Y ~ Bernoulli(m0(x) + W·τ(x))`.
- Control risk: `m0(x) = 0.34 + 0.36·plogis(−1.72 − 0.5·X1 + X2)`.
  - E[m0] = 0.40, and m0 lies in [0.34, 0.70].
  - X1 and X2 are therefore only weakly prognostic.
- g is the continuous g with its signs reversed and scaled by `RD_SCALE`
  (0.082, 0.051, 0.204 for scenarios 2, 3, 4).
- X4 and X5 enter g through tanh (t4 = tanh(X4), t5 = tanh(X5)), which
  keeps the risks bounded without changing the covariates' distributions.
- The event is harmful, and treatment lowers its risk.

### Scenarios

| scenario | description | continuous g(x) | binary g(x) |
|---|---|---|---|
| 1 | null | 0 | 0 |
| 2 | simple, continuous | −X4 | 0.082·t4 |
| 3 | complex | 2·X3 + 0.5·X4 − 0.5·X4·X5 | −0.102·X3 − 0.026·t4 + 0.026·t4·t5 |
| 4 | non-linear | cos(X4) | −0.204·cos(X4) |

These are coefficient for coefficient the parent studies' scenarios 1–4
(`SCENARIO_SETS$continuous[1:4, ]` and `SCENARIO_SETS$binary[1:4, ]`, plus a
`rho` column). τ(x) and m0(x) are therefore the same functions as in the parent
studies.

### Covariates: the copula

- A latent `Z ~ N(0, R)` is drawn over (X1, X2, X3, X4, X5, X01, X02, X03),
  with R exchangeable at ρ (`correlated_covariates()`, `exch_cor()`). Z is one
  `rnorm()` block of n × 8 draws, multiplied by `chol(R)`.
- X1 and X3 are thresholds of their latent columns at the parent studies'
  prevalences, so X1 ~ Bernoulli(0.4) and X3 ~ Bernoulli(0.7) marginally.
- X2, X4 and X5 are the latent columns, so each is N(0, 1) marginally.
- X01–X03 are the latent columns themselves, N(0, 1) marginally:
  - noise for both the outcome and the CATE;
  - correlated with X1–X5, so they are partial proxies for the modifiers.
- X04 and X05 are indicators of a 3-level factor (0.45 / 0.30 / 0.25), drawn
  independently of everything else.

Every marginal is the same as in the parent studies; only the dependence
changes. At ρ = 0.5 the observed correlations are about 0.5 between continuous
covariates, and 0.3–0.4 for pairs involving X1 or X3.

**ρ = 0 is the paired independent arm.** The covariates are independent, as
in the parent studies, but drawn through the copula branch in the copula's
order, so each run is paired with the same run at ρ = 0.5 (see Seeding and
pairing below).

### Treatment-effect calibration

`bW` is set so that the true ATE equals the effect an unadjusted test would
have 80% power to detect at that n, planned as if the effect were homogeneous
(`calibrate_bW()`):

- **Continuous:**
  - The power calculation is a two-sample t-test with n/2 per arm.
  - The planned SD is
    `sqrt(b1²·p(1 − p) + b2²·s2² + 2·b1·b2·Cov(X1, X2) + s_err²)`.
  - The covariance term (`cov_X1X2()`) is zero at ρ = 0. Var(g) is left out
    of the plan on purpose.
  - bW is rounded to 2 dp.
- **Binary:**
  - The power calculation is a two-proportion test with n/2 per arm, at the
    control event rate E[m0].
  - bW is rounded to 3 dp.
- **Both:** `bW = δ − E[g]`, with E[g] (`te_moments()`) and E[m0] taken over
  the correlated latent. At ρ = 0, bW equals the parent studies' bW at every n.

At ρ = 0.5, bW and the ATE move through two channels:

- E[g], in scenario 3, where g contains the product X4·X5;
- the planned SD or control event rate, through Cov(X1, X2).

From `R/calibration_report.R`, plus a large-sample simulation for SD(τ) and
cor(m0, τ):

| | ρ = 0 | ρ = 0.5 |
|---|---|---|
| **continuous** | | |
| true ATE, n = 100 / 250 / 500 / 1000 | −0.65 / −0.41 / −0.29 / −0.20 | −0.60 / … / −0.19 |
| realised power, scenarios 2 / 3 / 4 | 0.65–0.67 / 0.61–0.63 / 0.76–0.77 | 0.65–0.66 / 0.55–0.56 / 0.75–0.76 |
| SD(τ), scenario 3 | 1.16 | 1.36 |
| cor(m0, τ), scenarios 2 / 3 / 4 | 0 / 0 / 0 | −0.44 / 0.39 / 0.01 |
| **binary** | | |
| true RD, n = 100 / 250 / 500 / 1000 | −0.248 / −0.164 / −0.118 / −0.085 | −0.247 / … / −0.085 |
| power | 0.80 | 0.80 |
| cor(m0, τ), scenarios 2 / 3 / 4 | 0 / 0 / 0 | 0.39 / −0.34 / 0.06 |
| lowest treated risk, scenario 3, n = 100 | 0.011 | **0.007** |

Realised power for the continuous outcome falls below 0.80 because the plan
ignores Var(g), as a trial planned under homogeneity would. For the binary
outcome, realised power equals planned power, because each arm's variance is
fixed by its marginal risk.

**Binary scenario 3 floor.**

- Every treated risk stays within [0.01, 0.99] (`RD_EPS` = 0.01), with one
  exception: at ρ = 0.5, E[g] rises in scenario 3 and bW falls, so the lowest
  treated risk over the covariate support at n = 100 is 0.007.
- That is inside [0, 1] but below `RD_EPS`. It is accepted on purpose, so that
  τ(x) stays the parent binary study's τ(x): a scale of 0.049 would restore
  the floor, at the cost of a smaller HTE.
- `sample_size/binary/bin_verify_hte.R` check 7 holds this one cell to > 0, and
  every other cell to `RD_EPS`.

### Design

Full factorial, per outcome:

| | correlated cts / bin | correlated CI cts / bin | correlated optimal_sf cts / bin |
|---|---|---|---|
| scenarios | 1–4 | 1–4 | 1–4 |
| n | 100, 250, 500, 1000 | 500, 1000 | 500, 1000 |
| ρ | 0, 0.5 | 0, 0.5 | 0, 0.5 |
| `CI_sf` | — | 0.05 to 0.5 by 0.05 (design factor) | chosen per run from the same 10 values |
| repetitions | 500 | 100 | 100 |
| grid rows per outcome | 16,000 | 16,000 | 1,600 |
| ρ boundary in the grid | rows 1–8000 / 8001–16000 | rows 1–8000 / 8001–16000 | rows 1–800 / 801–1600 |

ρ varies slowest in every grid. The jobscripts split at the PBS array cap
(1–10000, 10001–16000), not at the ρ boundary.

### Seeding and pairing

- Each run is seeded by its run index alone (`setup_rng_stream(run)`).
- The draw order is W, then the n × 8 latent block, then the continuous noise
  ε, then the categorical pair X04/X05.
- So run *r* at ρ = 0 and at ρ = 0.5 shares:
  - W;
  - the raw normals of the latent block, multiplied by a different Cholesky
    factor;
  - ε (continuous outcome);
  - X04 and X05.
- Binary Y is coupled across ρ only through `rbinom`'s shared uniforms. About
  91–99% of outcomes agree.
- In the CI studies, run *r* is also the same dataset at every `CI_sf`.
- Scenarios share seeds within a run, so their errors are correlated.
- The copula draws in a different order from the independent-covariate
  generator. The ρ = 0 runs are therefore **not** paired with the parent
  studies' runs: they come from the same distribution, but are different draws.

## Estimands

- **Unit-level CATE** `τ(x_i)` for every unit in the simulated sample:
  - continuous outcome: a mean difference;
  - binary outcome: a risk difference,
    `P(Y = 1 | x, W = 1) − P(Y = 1 | x, W = 0)`.

  The truth is computed per run at that run's correlated covariate draws, with
  the ρ-specific bW (`truth_at()`), and saved with the run.
- **ATE**, as the sample mean of the CATE.
- **Presence of heterogeneity** (H0: τ constant), for the tests (cts / bin).
  Scenario 1 is the null and scenarios 2–4 are alternatives.
- **CI and optimal_sf:** the CATE at a fixed covariate query grid
  (`build_query_grid()`), with the same points in every run and at both ρ.
  - The grid crosses whichever of X3 (0, 1), X4 and X5 (−2, −1, 0, 1, 2) the
    scenario's τ uses, with every other covariate at 0.
  - That gives 1 point in scenario 1, 5 in scenario 2, 50 in scenario 3 and 5
    in scenario 4.
  - The bands are also scored over the sample's own units.
- **optimal_sf calibration only:** the estimated CATE τ̂, used as a plug-in
  target.
- **The comparison of interest** in every study is the change in each
  performance measure from ρ = 0 to ρ = 0.5.

**Query-grid caveat.**

- The grid is the parent studies' grid, unchanged. τ at its points is still
  the truth at ρ = 0.5.
- But a point such as (X4 = 2, X5 = −2) is about 4 SD from the centre of the
  correlated latent at ρ = 0.5, against about 2 SD at ρ = 0, so the forests
  extrapolate there.
- A fall in grid coverage at ρ = 0.5 can therefore be extrapolation rather
  than miscalibration. **The per-unit band is the primary ρ comparison.**

## Methods

The estimators are fit exactly as in the parent studies. The cts / bin analysis
scripts source `sample_size/continuous/cts_models.R` and
`sample_size/binary/bin_models.R`. The CI scripts source the parent CI studies'
`*_ci_models.R`. Model code is in `R/cate_models.R` (`cate_methods()`), and
SuperLearner set-up is in `R/sl_library.R`.

### Estimators

| method | description | used in |
|---|---|---|
| `causal_forest` | grf `causal_forest`, internal cross-fitting | cts / bin, CI |
| `dr_random_forest` | DR-learner: per-arm outcome forests (own-arm predictions OOB), RF propensity, OOB RF second stage | all six |
| `dr_superlearner` | DR-learner: SuperLearner per-arm outcome models and propensity, single leave-one-fold-out cross-fit, SuperLearner second stage | cts / bin |
| `dr_oracle` | DR-learner with the true outcome mean (binary: the true risk) and propensity 0.5 | cts / bin, CI |
| `dr_semi_oracle` | DR-learner with the known propensity 0.5 and per-arm outcome forests specified as in `dr_random_forest` (refitted) | cts / bin, CI |
| `t_random_forest` | T-learner μ̂1 − μ̂0 from `dr_random_forest`'s per-arm outcome forests; derived at metrics time, not refitted | cts / bin |
| `t_superlearner` | T-learner μ̂1 − μ̂0 from `dr_superlearner`'s per-arm SuperLearners (out of fold); derived at metrics time, not refitted | cts / bin |
| `causal_forest_inbuilt` | `causal_forest`'s own pointwise interval from grf's variance estimate (an interval, not a separate fit) | CI |

The CI and optimal_sf studies have no SuperLearner arm. SuperLearner's
`family` is the only difference between the continuous and binary model code,
so the continuous and binary CI studies fit identical models. The oracles'
outcome formula is the ρ-specific one (`get_oracle_info()` with the run's bW).

### Model fitting

**Inputs.**

- Every learner receives the 10 covariates X1–X5 and X01–X05 as generated,
  with no scaling or feature engineering. The interaction lassos below are
  the one exception.
- W enters the causal forest as the treatment. It is never a covariate of an
  outcome model, which is fit per arm.

**DR pseudo-outcome** (`dr_pseudo`), shared by every DR-learner:

$$\hat\phi_i = \hat\mu_1(x_i) - \hat\mu_0(x_i) + \frac{(Y_i - \hat\mu_{W_i}(x_i))(W_i - \hat e(x_i))}{\hat e(x_i)(1 - \hat e(x_i))}$$

The CATE is the second-stage regression of φ̂ on X. Estimated propensities in
`dr_random_forest` and `dr_superlearner` are trimmed to [0.05, 0.95]. The
causal forest's internal propensity is not trimmed.

**Forests.**

- All forests are grf `regression_forest` / `causal_forest` at grf's defaults:
  - 2000 trees, `sample.fraction` 0.5;
  - honest splitting (`honesty.fraction` 0.5, `honesty.prune.leaves` TRUE);
  - `min.node.size` 5, `alpha` 0.05, `ci.group.size` 2;
  - `mtry` = min(⌈√p⌉ + 20, p), which is all 10 covariates;
  - no parameter tuning.
- Only the CI studies' bootstrap forests change a setting
  (`sample.fraction = CI_sf`).

**Per estimator.**

- `causal_forest`:
  - fit as `causal_forest(X, Y, W)` with no supplied nuisances;
  - grf fits `Y.hat` and `W.hat` itself, as OOB predictions from regression
    forests of 500 trees;
  - τ̂ is the forest's OOB prediction, with grf's default
    `stabilize.splits = TRUE`.
- `dr_random_forest` (whole sample, no sample splitting):
  - μ̂0 and μ̂1: one regression forest of Y on X per arm, fit on that arm's rows.
    A unit's own-arm prediction is OOB, and its other-arm prediction comes from
    a forest that never saw it.
  - ê: a regression forest of W on X, OOB predictions, trimmed.
  - Stage 2: a regression forest of φ̂ on X. τ̂ is its OOB prediction, and
    grf's OOB variance estimate is saved.
- `dr_oracle`:
  - μ_w(x) is the true outcome mean at W = w;
  - e = 0.5;
  - stage 2 as `dr_random_forest`.
- `dr_semi_oracle`:
  - per-arm outcome forests as in `dr_random_forest` (a separate fit);
  - e = 0.5;
  - stage 2 as `dr_random_forest`.
- `dr_superlearner`:
  - A single cross-fit with V = 4, 5, 10 folds at n = 100, 250, ≥ 500. Rows
    are assigned to folds in contiguous blocks.
  - Stage 1, per fold k: on the rows outside k, fit one SuperLearner of Y on X
    per arm and one of W on X. Predict fold k's rows, trim ê, and form φ̂ for
    fold k.
  - Stage 2, per fold k: fit a SuperLearner of φ̂ on X on the rows outside k,
    and predict fold k. Stage 1 and stage 2 use the same folds.
- `t_random_forest`, `t_superlearner`: μ̂1 − μ̂0, recovered at metrics time
  from the saved DR nuisances by removing the residual term from φ̂
  (`add_t_learners()` in `R/metrics.R`).

### SuperLearner (`dr_superlearner` only)

`sl_fit_predict()`:

| model | family | meta-learner |
|---|---|---|
| propensity | binomial | `method.NNloglik` |
| outcome, continuous | gaussian | `method.NNLS` |
| outcome, binary | binomial | `method.NNloglik` |
| CATE (stage 2) | gaussian | `method.NNLS` |

**Fitting rules.**

- Ensemble weights come from SuperLearner's internal 10-fold CV, within each
  cross-fitting training set.
- Before each fit, each candidate is fit alone with 2-fold CV
  (`pretest_superlearner`). Learners that error or give non-finite
  predictions are dropped and recorded in `sl_dropped`. If every candidate
  fails, the library is `SL.mean` alone.
- If a SuperLearner fit itself errors, every prediction is the mean of its
  training outcome, with a warning, recorded in `sl_dropped` as
  `(whole fit)`.
- If an outcome or propensity fit returns all-zero predictions, they are
  replaced by the training mean (for an outcome model, that arm's mean), with
  a warning.

**Libraries.** There is one library per nuisance (`sl_libraries()`), chosen
without reference to the scenarios' HTE shapes. The outcome model is fit in
each arm separately, so its library is sized to about half the training rows
(about 37 per arm at n = 100):

| | n = 100 | n = 250 | n ≥ 500 |
|---|---|---|---|
| propensity | mean, glm | same | same |
| outcome (per arm) | mean, lasso, ranger (min.node.size 25) | + glm, gam | + earth |
| CATE (stage 2) | mean, glm, lasso (lambda.min and lambda.1se), gam, ranger (min.node.size 25) | + lasso on pairwise interactions (both tunings), earth | same |

**Candidate learners.** These use SuperLearner's wrapper defaults unless
stated, with p = 10:

| learner | definition |
|---|---|
| `SL.mean` | training mean |
| `SL.glm` | `glm` on all main effects |
| `SL.glmnet` | lasso (`alpha` 1) on main effects, λ by 10-fold `cv.glmnet` at `lambda.min` |
| `SL.glmnet.1se` | as `SL.glmnet`, at `lambda.1se` (custom wrapper) |
| `SL.glmnet.int` | lasso on all main effects and pairwise interactions (`model.matrix(~ .^2)`), `lambda.min` (custom wrapper) |
| `SL.glmnet.int.1se` | as `SL.glmnet.int`, at `lambda.1se` (custom wrapper) |
| `SL.gam` | `gam` package; smoothing spline with 2 df for covariates with more than 4 unique values, linear otherwise |
| `SL.earth` | MARS: `degree` 2, `penalty` 3, `nk` = max(21, 2p + 1), backward pruning |
| `SL.ranger.ns25` | ranger, 500 trees, `mtry` ⌊√p⌋ = 3, `min.node.size` 25 (custom wrapper); a probability forest for a binomial outcome |

### Heterogeneity tests (correlated cts / bin)

These are run per method on the full sample (`run_blp_whole()`,
`run_independence_test_whole()`).

**BLP** (GenericML `BLP`, Chernozhukov et al.):

- A weighted least squares regression of Y, with weights 1/(ê(1 − ê)), on:
  - an intercept;
  - the baseline proxy μ̂0(x);
  - (W − ê);
  - (W − ê)(τ̂(x) − mean τ̂).
- The test is on the last coefficient, β2.
- ê and μ̂0 are the estimator's own stage-1 nuisances:
  - the forest nuisances for `causal_forest`, `dr_random_forest` and
    `t_random_forest`;
  - the SuperLearner nuisances for `dr_superlearner` and `t_superlearner`;
  - e = 0.5 and the arm's own μ0 for the oracles.
- Two p-values are reported:
  - `BLP_p`: two-sided, homoskedastic OLS standard errors, saved at estimation
    time.
  - `BLP_p_os`: one-sided (H1: β2 > 0), HC3 standard errors, recomputed at
    metrics time from the saved nuisances.

**Independence test** (coin `independence_test(τ ~ X, teststat =
"quadratic")`, with an asymptotic χ² reference):

- `indep_cate`: the estimated CATE against the covariates.
- `indep_po`: the DR pseudo-outcome against the covariates. It is reported
  once per pseudo-outcome, so it is not reported for `causal_forest` (which
  shares `dr_random_forest`'s) or for the T-learners.

A CATE that is constant up to rounding gives NA for both tests.

**True-CATE tests.** The BLP (both p-values) and `indep_cate` are also run on
the true CATE, with the true nuisances: tau = τ, μ̂0 = the true m0, ê = 0.5.
These use the same test functions and call shape (`run_true_cate_tests()`,
`true_cate_test_row()`), and the true CATE is the one computed for that run at
that ρ. They show how the tests behave with no estimation error. Scenario 1's
true CATE is constant, so these tests are NA there.

### Intervals (correlated CI)

**Half-sample bootstrap band** (`R/bootstrap_ci.R`), with B = 200 and
α = 0.05:

1. In each draw, refit the second stage on an unstratified half sample, with
   `sample.fraction = CI_sf`:
   - DR-learners: a regression forest on the fixed pseudo-outcomes;
   - causal forest: with `Y.hat` and `W.hat` fixed.
2. Form the roots `tau_full − tau_half`. In-half units are predicted OOB, and
   the rest out of sample.
3. Standardise the roots by their bootstrap SD.
4. Take as the critical value the 1 − α quantile of the per-draw maximum, over
   units, of the absolute standardised roots. The absolute value makes the
   band two-sided, so this gives a 95% simultaneous band.

The band is built over the sample's units and, separately, over the query
grid.

**`causal_forest_inbuilt`:** a pointwise normal interval from grf's variance
estimate, at α = 0.05.

### Calibration (correlated optimal_sf)

This is `dr_random_forest` only (`find_optimal_sf()`), with α = 0.10, so the
target is 0.90.

1. For each candidate `CI_sf` in 0.05, 0.10, …, 0.50, repeat 50 times:
   - resample the pseudo-outcome residuals around τ̂, whole sample, with
     replacement;
   - re-estimate the CATE from τ̂ plus the resampled residuals;
   - build the half-sample band (B = 100);
   - record its marginal coverage of τ̂ (the plug-in truth) and its mean width.
2. Pick the candidate whose mean coverage is closest to 0.90.
3. Build the final band on the observed data at that `CI_sf` (B = 200), over
   the sample's units and the query grid.

## Performance measures

### Correlated continuous / binary

Each measure is computed per run over the units in the sample
(`cate_metrics()`), then averaged over runs within each (ρ, scenario, n,
model) cell:

- **bias:** mean of (τ̂ − τ) over units.
- **ATE bias:** mean τ̂ − mean τ.
- **relative ATE bias:** ATE bias / true ATE.
- **relative CATE bias:** mean of (τ̂ − τ)/τ over units with τ ≠ 0.
- **MSE, RMSE, MAE** of τ̂ against τ.
- **Pearson and Spearman correlation** of τ̂ with τ. These are set to 0 in
  scenario 1, where τ is constant.
- **sign accuracy:** the proportion of units where sign(τ̂) = sign(τ).
- **`n_na`:** the number of units with no estimate.
- **HTE tests:** `BLP_p`, `BLP_p_os`, `indep_cate` and `indep_po` for each
  estimator, plus the true-CATE `BLP_p`, `BLP_p_os` and `indep_cate`.
  - Each becomes a rejection rate at 0.05: type I error in scenario 1, power
    in scenarios 2–4.
  - The number of runs with an NA p-value is reported alongside.
- **ATE bias relative to `dr_oracle`:** each estimator's ATE bias minus
  `dr_oracle`'s on the same run. This removes the noise the estimators share
  through a common dataset.

### Correlated CI continuous / binary

Per run (`interval_metrics()`), then averaged over runs within each
(ρ, scenario, n, `CI_sf`, model) cell. Nominal coverage is 0.95.

- **marginal coverage:** the proportion of units covered by their own interval.
- **simultaneous coverage:** whether the band covers every unit; this is the
  target the method controls.
- **mean interval length.**
- **query grid** (`<model>_grid` rows): the same three measures, plus
  bias-eliminated marginal and simultaneous coverage (`be_interval_metrics()`).
  - These are coverage of the across-run mean estimate at each grid point
    (`grid_be_reference()`), computed within each ρ.
  - Bias drops out, which separates miscalibration from extrapolation bias.

The per-unit band is the primary ρ comparison (see the query-grid caveat
under Estimands).

### Correlated optimal_sf continuous / binary

Each run's final band is scored against the true CATE with the CI measures
above, over units and over the grid, at nominal 0.90. Each row also records:

- `optimal_sf`: the chosen `CI_sf`;
- `plugin_coverage`: the calibration's own mean coverage of τ̂ at the chosen
  value;
- `plugin_ci_width`: the calibration's mean width at the chosen value.

The quantity of interest is the gap between `plugin_coverage`, which is
calibrated to 0.90, and the true-τ coverage. It shows how much of τ̂'s bias the
calibration misses, and whether that gap grows at ρ = 0.5.

### Monte Carlo error and comparisons

- **The effect of ρ (ρ = 0.5 − ρ = 0):** use paired per-run differences. The
  MCSE is sd(diff)/√runs.
- **Within one ρ:**
  - compare estimators against `dr_oracle` on the same run;
  - compute MCSEs per scenario, never pooled, because scenarios share seeds
    and so their errors are correlated;
  - rejection rates and coverages are per-run 0/1 indicators with a binomial
    MCSE.
- **ρ = 0 against the parent study:** these have the same distribution but are
  not paired, so use the unpaired MCSE. This serves as a sanity check that the
  copula branch reproduces the parent study.
- **Summary plots:** mean ± qnorm(0.975) × MCSE.

### Where they are reported

| report | contents |
|---|---|
| `corr_results.qmd` | both outcomes: levels by ρ, paired differences, HTE-test rejection rates, ATE bias relative to `dr_oracle`, and the ρ = 0 sanity check |
| `continuous/cts_corr_results.qmd`, `binary/bin_corr_results.qmd` | one outcome each: every metric by ρ, paired differences, the true-CATE tests, NA tables, the ρ = 0 sanity check, and a headline table |
| `confidence_intervals/<outcome>/*_corr_ci_results.qmd` | coverage and length by ρ across `CI_sf`, paired ρ differences, interval types, the best `CI_sf` by ρ, the optimal_sf pick and its true-τ coverage, the ρ = 0 sanity check, and a headline table |
