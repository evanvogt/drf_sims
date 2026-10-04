# ADEMP — single-event survival

One study, `competing_risk/single_event/` (`se_*` scripts). The DGM is
`competing_risk/`'s event 1 with event 2 removed (`se_dgm.R`). Most estimators
are that study's own arms, reused unchanged from `competing_risk/surv_models.R`;
`se_models.R` adds the pieces a single event needs.

## Aims

- Compare pseudo-value CATE estimators with the causal survival forest and the
  random survival forest (RSF) DR-learner, in a two-arm RCT with one
  time-to-event outcome and no competing event.
- Compare DR-learners with T-learners built on the same outcome models
  (random forest, SuperLearner and RSF).
- Give the competing-risk study a reference point. Its net RMST1 estimand is
  this study's estimand, so differences between the two studies come from the
  competing event rather than from survival estimation.

## Data-generating mechanisms

Implementation: `se_dgm.R::generate_se_data()`, parameters in
`se_scenario_params`, built from `competing_risk`'s `survival_scenario_params`
so the two cannot drift apart.

### Covariates and treatment

The same as `competing_risk/` at ρ = 0, and drawn by the same code
(`correlated_covariates()`, `R/dgm_scenarios.R`):

- `W ~ Bernoulli(0.5)`, independent of everything (an RCT)
- Prognostic: `X1 ~ Bernoulli(0.4)`, `X2 ~ N(0, 1)`
- Effect modifier: `X3 ~ Bernoulli(0.7)`
- Noise: `X01`–`X03 ~ N(0, 1)`; `X04`, `X05` indicators of a 3-level factor
  (0.45 / 0.30 / 0.25)

Every estimator sees the 8 covariates X1, X2, X3, X01–X05. There is no ρ
factor.

### Event time

Weibull, proportional hazards, with effects on the log scale:

$$S(t \mid x, w) = \exp\{-(t / s)^{a}\}, \qquad \log s = \log 25 - 0.1\,X_1 + 0.1\,X_2 + w\,(b_W + b_3 X_3), \qquad a = 2$$

T is drawn by inversion, T = s·(−log U)^{1/a}, `U ~ Uniform(0, 1)`.

**Censoring.** Administrative at 180 (the study end) in every run. With
`censoring = TRUE`, also `C ~ Uniform(1, 180)`, independent of X and W. Only
censoring before the horizon of 28 matters to the estimand.

### Scenarios

| scenario | description | b_W | b_3 | HR, X3 = 0 / 1 | as in `competing_risk/` |
|---|---|---|---|---|---|
| 1 | no treatment effect | 0 | 0 | 1 / 1 | scenario 2's event 1 |
| 2 | constant treatment effect | −0.175 | 0 | 1.42 / 1.42 | scenario 1's event 1 |
| 3 | heterogeneous treatment effect | −0.175 | −0.35 | 1.42 / 2.86 | scenario 3's event 1 |

The treatment raises the hazard, as in `competing_risk/`. "Constant" is on the
hazard scale. On the RMST scale the CATE still varies a little with X1 and X2,
because a restricted mean is not linear in them (SD 0.10 against a mean of
−2.16). An RMST CATE that is exactly constant is not possible with a Weibull
proportional-hazards model and prognostic covariates.

### Design

| | |
|---|---|
| scenarios | 1–3 |
| n | 500 |
| censoring | TRUE (uniform + administrative), FALSE (administrative only) |
| horizon | 28 |
| repetitions | 500 per (scenario, censoring) |
| array | 3,000 jobs (`se_config.R`, `jobscripts/se_1.sh`) |
| folds | V = 10, contiguous blocks of 50 rows |

**Seeding.** Each run is seeded by its run index alone (`setup_rng_stream(run)`).
Draw order: `W, Z (n × 8 copula normals), cats, U, [C]`. So:

- the three scenarios and the two censoring settings of a run share W, the
  covariates and U, and comparisons across them are paired;
- **the run is paired with `competing_risk/`'s ρ = 0 run of the same index.**
  The draw order is that study's without its `cause` draw, so W, X and U are
  identical, and T here solves Λ₁(T) = −log U where there it solves
  Λ₁(T) + Λ₂(T) = −log U. So T_single ≥ T_CR for every unit. C is not shared,
  because the parent draws `cause` before C.

`se_dgm_check.R` asserts both the pairing and the truth (below).

### What the DGM implies

From `Rscript se_dgm_check.R`. The event mix comes from `generate_se_data()`
at n = 20,000 per (scenario, censoring) cell, seeded once, so it is
illustrative. The truths are population values over 200,000 covariate draws.

**Event mix by the horizon.**

| scenario | arm | event by 28 | censored before 28 (censoring on) | observed past 28, censoring off / on |
|---|---|---|---|---|
| any | control | 0.74 | 0.10 | 0.26 / 0.22 |
| 1 | treated | 0.73 | 0.10 | 0.27 / 0.22 |
| 2 | treated | 0.86 | 0.09 | 0.14 / 0.13 |
| 3 | treated | 0.94 | 0.07 | 0.06 / 0.06 |

- Censoring on removes about 10% of each arm before the horizon. That is light
  censoring, so the censored and uncensored cells may not separate the methods
  much.
- **Without censoring, every pseudo-value is exactly min(Y, 28).** The
  pseudo-value arms are then regressions on the observed truncated time, and
  the pseudo-value construction does nothing.
- Someone is always observed past 28: about 55–65 controls per run of 500,
  and at least about 14 treated units even in scenario 3 with censoring on. The
  NA pseudo-values behind `competing_risk/`'s route C need nobody past the
  horizon, so they cannot occur in practice. `pseudo_rmst()` stops the run if
  one does.

**The true CATE.** Control-arm RMST 19.13 days.

| scenario | treated RMST | τ mean | τ SD | τ at X3 = 0 | τ at X3 = 1 | τ range |
|---|---|---|---|---|---|---|
| 1 | 19.13 | 0 | 0 | 0 | 0 | 0 |
| 2 | 16.96 | −2.16 | 0.10 | −2.16 | −2.16 | −2.27 to −1.42 |
| 3 | 13.89 | −5.24 | 2.02 | −2.16 | −6.56 | −6.65 to −1.42 |

These equal `competing_risk/`'s τ_RMST1_cs (to 1e-13 at the same covariates)
in its scenarios 2, 1 and 3.

## Estimands

The unit-level CATE on restricted mean survival time at horizon τ = 28,
evaluated at each sampled unit's (X1, X2, X3):

$$\tau(x) = \int_0^{28} S(t \mid x, 1)\,dt - \int_0^{28} S(t \mid x, 0)\,dt$$

For a Weibull the restricted mean has a closed form,
s·Γ(1 + 1/a)·P(1/a, (28/s)^a), with P the regularised lower incomplete gamma
(`rmst_weibull()`). For a = 2 this is s·(√π/2)·erf(28/s). The truth column is
`tau_RMST`. The ATE is the sample mean of the CATE. Point estimation only: no
intervals, no HTE tests.

## Methods

Implementation: `se_models.R::all_cate_se_models()`, called once per run by
`se_analysis.R` with `n_folds = 10`, `horizon = 28` and `sl_libraries(500)`.
Inputs are X (the 8 covariates as generated), W, the observed time Y and the
status D ∈ {0 censored, 1 event}.

Ten arms (`framework` in the results):

| arm | learner | outcome model | fitting | from |
|---|---|---|---|---|
| `csf` | causal survival forest | — | grf-internal | `csf_cs(event = 1)` |
| `pseudo_cf_whole_oob` | causal forest on pseudo-values | — | grf-internal | `pseudo_cf_whole_oob` |
| `pseudo_dr_whole_oob` | DR-learner | regression forest on pseudo-values, per arm | whole-sample OOB | `nuisance_pseudo_rf_oob` + `stage2_whole_rf` |
| `pseudo_t_whole_oob` | T-learner | 〃 | 〃 | `t_from_dr()` |
| `sl_dr_whole` | DR-learner | SuperLearner on pseudo-values, per arm | single crossfit, both stages | `nuisance_pseudo_sl` + `pseudo_dr_sl` |
| `sl_t_whole` | T-learner | 〃 | single crossfit | `t_from_dr()` |
| `rsf_dr_oob` | DR-learner | survival forest (randomForestSRC) on (Y, D), per arm | whole-sample OOB | `nuisance_rsf_se_oob` + `stage2_whole_rf` |
| `rsf_t_oob` | T-learner | 〃 | 〃 | `t_from_dr()` |
| `rsf_dr_scf` | DR-learner | 〃 | single crossfit, both stages | `nuisance_rsf_se_scf` + `stage_2_rf_scf` |
| `rsf_t_scf` | T-learner | 〃 | single crossfit | `t_from_dr()` |

Forests are grf at its defaults, apart from the RSF outcome models
(randomForestSRC at its defaults: 500 trees, `nodesize` 15, log-rank
splitting). Every estimated propensity in a DR-learner is trimmed to
[0.05, 0.95] (`trim_ps`).

- **Pseudo-values** (`pseudo_rmst()`): `pseudo::pseudomean(Y, D, 28)`,
  jackknife pseudo-values of the RMST from one Kaplan–Meier fit on all n. No
  covariates are needed because censoring is independent of X and W. Only
  whole-sample pseudo-values: the parent's whole/cvps comparison is not
  repeated here.
- **`csf`**: `causal_survival_forest(X, Y, W, D, target = "RMST",
  horizon = 28)`, OOB predictions. grf estimates the censoring itself.
- **DR-learners**: the DR pseudo-outcome with the pseudo-value θ in place of Y,
  as in `competing_risk/` (its ADEMP has the formula), then a regression of it
  on X. Per-arm outcome models (T-learner). The correction term always uses
  the whole-sample θ.
- **RSF outcome model**: a survival forest per arm on (Y, D) censored at 28,
  `ntime = 0`, μ̂ = ∫₀²⁸ Ŝ(t) dt from the step function, the last interval
  running from the last event time to 28 (`rsf_rmst()`). Own-arm predictions
  OOB in `_oob`, all from training folds in `_scf`. Pseudo-values enter only
  the correction term. ê and stage 2 are the grf ones of
  `pseudo_dr_whole_oob` / `stage_2_rf_scf`, so the outcome model is the only
  thing that differs from the pseudo-value RF DR-learner.
- **T-learners** (`t_from_dr()`): μ̂₁ − μ̂₀ from the DR-learner's own nuisances,
  recovered exactly as po − (θ − μ̂_W)(W − ê)/(ê(1 − ê)). No extra fitting, and
  the same honesty as the DR-learner (OOB or crossfit predictions). The
  SuperLearner T-learner is the same estimator as `competing_risk/`'s
  `sl_t_whole`.
- **SuperLearner** settings are the parent's (`sl_libraries(500)`, pretested,
  per nuisance; see `competing_risk/ADEMP.md`). In an RCT the propensity
  library's meta-learner often gives every learner zero weight, and the
  failsafe then uses mean(W). That is inherited behaviour and expected.

**Saved per run:** τ̂ per arm; the pseudo-values; the nuisances of the RF, SL
and RSF DR-learners (`po`, `pseudo.hat`, `pseudo0.hat`, `pseudo.hat.cf`,
`W.hat`); fold indices; the data and the truth.

## Performance measures

Per run, per arm, against `tau_RMST` (`se_metrics.R`), using `cate_metrics()`
so the conventions match the other studies:

- bias, ATE bias, relative ATE bias, relative CATE bias
- MSE, RMSE, MAE
- Pearson and Spearman correlation with the true CATE, sign accuracy
- C-statistic, (Kendall τ_b + 1) / 2
- `tau_sd`, the within-run SD of τ̂. In scenario 1 the truth is 0, so this is
  the spurious heterogeneity an arm reports.
- `n_na`

`se_results.qmd` averages each over runs, with Monte Carlo SE
sd / √(non-NA runs), plotted as mean ± 1.96 MCSE.

- **Scenario 1 is a true null.** `cate_metrics()` and `c_statistic()`
  hard-code scenario 1 as having no heterogeneity (correlations 0, C 0.5). In
  `competing_risk/` that convention is wrong for its scenario 1. Here it is
  right.
- **Scenario 2's ranking measures are near-meaningless.** Correlation, sign
  accuracy and the C-statistic are measured against a truth with SD 0.10. Read
  bias and RMSE there.
