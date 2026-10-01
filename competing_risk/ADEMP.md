# ADEMP — competing risks

One study, `competing_risk/` (`surv_*` scripts). The DGM is its own
(`surv_dgm.R`), not a `sample_size/` one, and most estimators are specific to
it (`surv_models.R`), borrowing the shared DR-learner and SuperLearner
machinery from `R/cate_models.R` and `R/sl_library.R`.

## Aims

- Assess how well forest- and SuperLearner-based estimators recover the
  unit-level CATE on restricted-mean time scales in a two-arm RCT where the
  event of interest (EOI) can be pre-empted by a competing event (CE), with and
  without right censoring.
- Compare three ways of handling the competing event: as censoring
  (cause-specific / net), kept in the risk set (subdistribution), or through
  pseudo-values of the cumulative incidence (restricted mean time lost, RMTL).
- Show what an analysis of one event reports when the treatment acts only on
  the other (scenarios 2, 5 and 6): an apparent effect on the EOI that is
  really a change in who survives long enough to have it.
- Secondary: how the pseudo-values should be built for a meta-learner -
  whole-sample against leave-one-fold-out pseudo-values, and whole-sample OOB
  against single crossfit fitting.

## Data-generating mechanisms

Implementation: `surv_dgm.R::generate_surv_data()`, parameters in
`survival_scenario_params`.

### Covariates and treatment

- `W ~ Bernoulli(0.5)`
- Prognostic for both events: `X1 ~ Bernoulli(0.4)`, `X2 ~ N(0, 1)`
- Effect modifier (drawn in every scenario, whether or not it is used):
  `X3 ~ Bernoulli(0.7)`
- Noise: `X01`–`X03 ~ N(0, 1)`; `X04`, `X05` indicators of a 3-level factor
  (0.45 / 0.30 / 0.25)

So every estimator sees 8 covariates (X1, X2, X3, X01–X05), all independent.
Unlike `sample_size/` there is no X4 / X5: the only effect modifier is the
binary X3.

### Event times

Cause-specific hazards are Weibull and proportional, cause k = 1 (EOI) or
2 (CE):

$$\lambda_k(t \mid x, w) = \frac{a_k}{s_k}\left(\frac{t}{s_k}\right)^{a_k - 1}, \qquad \log s_k = \log s_{k0} + b_{1k} X_1 + b_{2k} X_2 + w\,(b_{Wk} + b_{3k} X_3)$$

| | shape a_k | baseline scale s_k0 | b_1k (X1) | b_2k (X2) |
|---|---|---|---|---|
| EOI (k = 1) | 2 | 15 | −0.1 | 0.1 |
| CE (k = 2) | 1.1 | 45 | −0.1 | 0.1 |

Effects are on the log *scale*, so a shift b multiplies the hazard by
exp(−a_k·b): X1 raises the EOI hazard by 1.22 and the CE hazard by 1.12; one
SD of X2 lowers them by 0.82 and 0.90.

Generation follows Beyersmann et al. (2009): the event time T solves
Λ1(T) + Λ2(T) = −log U, `U ~ Uniform(0, 1)`, by `uniroot` on (0, 200)
(`find_time()`), and the cause is `Bernoulli(λ1(T) / (λ1(T) + λ2(T)))`.

**Censoring.** Administrative at 180 in every run. With `censoring = TRUE`,
also `C ~ Uniform(1, 180)`, independent of X and W. Only censoring before the
horizon of 28 matters to the estimands.

### Scenarios

Treatment effects on the log scale, and the hazard ratios they imply:

| scenario | description | b_W1 | b_31 | b_W2 | b_32 | EOI HR, X3 = 0 / 1 | CE HR, X3 = 0 / 1 |
|---|---|---|---|---|---|---|---|
| 1 | ATE on EOI only | −0.7 | 0 | 0 | 0 | 4.06 / 4.06 | 1 / 1 |
| 2 | ATE on CE only | 0 | 0 | 0.7 | 0 | 1 / 1 | 0.46 / 0.46 |
| 3 | HTE on EOI, no ATE on CE | −0.7 | −0.7 | 0 | 0 | 4.06 / 16.4 | 1 / 1 |
| 4 | HTE on EOI, ATE on CE | −0.7 | −0.7 | 0.7 | 0 | 4.06 / 16.4 | 0.46 / 0.46 |
| 5 | HTE on CE, no ATE on EOI | 0 | 0 | 0.7 | 0.7 | 1 / 1 | 0.46 / 0.21 |
| 6 | HTE on CE, ATE on EOI | −0.7 | 0 | 0.7 | 0.7 | 4.06 / 4.06 | 0.46 / 0.21 |
| 7 | HTE on both | −0.7 | −0.7 | 0.7 | 0.7 | 4.06 / 16.4 | 0.46 / 0.21 |

The treatment *raises* the EOI hazard and lowers the CE hazard. "ATE" and
"HTE" in the descriptions refer to the hazard scale, i.e. whether X3 modifies
the effect. On the time scales the estimands use, every effect is
heterogeneous. X1 and X2 shift the baseline hazards, and a restricted mean is
not linear in them, so the CATE varies with X1 and X2 even in scenarios 1 and
2 (see "What the DGM implies" below).

### Design

| | |
|---|---|
| scenarios | 1–7 |
| n | 500 |
| censoring | TRUE (uniform + administrative), FALSE (administrative only) |
| horizon | 28 |
| repetitions | 500 per (scenario, censoring) |
| array | 7,000 jobs (`surv_config.R`) |
| folds | V = 10 (`n_folds = 10` at n ≥ 300), contiguous blocks of 50 rows |

Runs 1–100 are the original 1,400-job array. Runs 101–500 were added later,
through `jobscripts/surv_extra.sh` or `surv_run.R` (see `README.md`). The
results report's text (`surv_results.qmd`) still says 100 runs.

**Seeding.** Each run is seeded by its run index alone
(`setup_rng_stream(run)`). Draw order: `W, X1, X2, X3, U, cause, [C], X01–X03,
cats`. So within a run all 14 (scenario, censoring) cells share W, X1–X3 and
the uniform the event time is solved from. Scenarios differ only by their
parameters, so comparisons across scenarios are paired. `censoring = TRUE`
consumes n more uniforms before the noise covariates, so X01–X05 differ
between the two censoring settings of a run.

### What the DGM implies

From `generate_surv_data()` at n = 20,000 per cell (event mix) and one draw at
n = 500 (truth). Seeded once, not with the study's streams, so these numbers
are illustrative.

**Event mix by the horizon.** The control arm is the same in every scenario:
EOI 0.76, CE 0.22, still event-free at 28 about 0.015. With censoring on that
becomes EOI 0.71–0.72, CE 0.21–0.22, censored before 28 about 0.055, and
about 0.014 observed past 28. Treated arm:

| scenario | EOI by 28 | CE by 28 | event-free at 28 | censored before 28 (censoring on) | observed past 28 (censoring on) |
|---|---|---|---|---|---|
| 1 | 0.89 | 0.11 | 0.000 | 0.028 | 0.000 |
| 2 | 0.87 | 0.11 | 0.023 | 0.061 | 0.019 |
| 3 | 0.93 | 0.07 | 0.000 | 0.014 | 0.000 |
| 4 | 0.97 | 0.03 | 0.000 | 0.019 | 0.000 |
| 5 | 0.90 | 0.07 | 0.026 | 0.064 | 0.020 |
| 6 | 0.97 | 0.03 | 0.000 | 0.032 | 0.000 |
| 7 | 0.97 | 0.03 | 0.000 | 0.015 | 0.000 |

- Almost everyone has an event before the horizon. The horizon, not the
  censoring, truncates the restricted means, and censoring is light: at most
  6% of a treated arm, 6% of controls.
- Wherever the treatment acts on the EOI (scenarios 1, 3, 4, 6, 7), no treated
  unit reaches 28, and about 1.5% of controls do: 3–4 units of 500. The data
  past the horizon are therefore thin. When there are none, the
  `pseudo::pseudoyl()` leave-one-out risk set empties at the horizon and
  returns NA pseudo-values (README "Known issues", route C).
- In the same scenarios the treated arm has few CEs (3–11%), so its RMTL2
  pseudo-values take only a handful of distinct values. That is the
  degenerate SuperLearner cell of route A'.

**The true CATEs.** Mean (SD) over the 500 units, in days. Control-arm levels:
RMTL1 12.5, RMTL2 4.5, RMSTc 11.0 (they sum to 28), net RMST1 12.8, net
RMST2 21.2.

| scenario | τ_RMTL1 | τ_RMTL2 | τ_RMSTc | τ_RMST1_cs | τ_RMST2_cs | cor(τ_RMTL1, τ_RMTL2) |
|---|---|---|---|---|---|---|
| 1 | 6.85 (0.37) | −1.79 (0.14) | −5.05 (0.51) | −6.40 (0.62) | 0 | +1.00 |
| 2 | 1.36 (0.18) | −2.28 (0.09) | 0.91 (0.08) | 0 | 3.32 (0.30) | −1.00 |
| 3 | 9.80 (1.93) | −2.72 (0.62) | −7.08 (1.46) | −8.67 (1.65) | 0 | −0.83 |
| 4 | 10.58 (1.76) | −3.65 (0.35) | −6.93 (1.56) | −8.62 (1.69) | 3.32 (0.30) | −0.64 |
| 5 | 1.87 (0.40) | −3.11 (0.55) | 1.23 (0.24) | 0 | 4.62 (0.93) | −0.92 |
| 6 | 8.41 (0.40) | −3.69 (0.36) | −4.73 (0.48) | −6.40 (0.63) | 4.58 (0.92) | −0.20 |
| 7 | 10.73 (1.93) | −3.89 (0.50) | −6.84 (1.57) | −8.56 (1.72) | 4.56 (0.96) | −0.79 |

τ_RMST1 and τ_RMST2 (subdistribution) are −τ_RMTL1 and −τ_RMTL2 exactly.

- **An effect on one cause moves the other's CIF.** In scenarios 2 and 5 the
  treatment does not touch the EOI hazard, yet τ_RMTL1 is +1.4 and +1.9: fewer
  CEs leave more units at risk of the EOI. The net truth τ_RMST1_cs is
  identically 0 there. Likewise τ_RMTL2 is −1.8 and −2.7 in scenarios 1 and 3,
  where τ_RMST2_cs is identically 0.
- **Scenario 5's RMTL1 heterogeneity is all induced** (cor −0.92 with
  τ_RMTL2), so an estimator that tracks the real, event-2, signal lands
  anti-correlated with τ_RMTL1. See README "Scenario 5 flips the sign".
- **"ATE" scenarios are not null on the time scale.** Scenario 1's τ_RMTL1
  has SD 0.37 and τ_RMST1_cs 0.62, all from X1 and X2. That is weak
  heterogeneity, but not zero (see "Performance measures").
- **Where X3 acts.** Mean τ_RMTL1 by X3 = 0 / 1 is 6.9 / 11.0 in scenario 3
  and 8.0 / 12.0 in scenario 7, but 8.0 / 8.6 in scenario 6: X3 modifying the
  CE barely reaches the EOI's RMTL (SD 0.40 against scenario 1's 0.37).

## Estimands

The unit-level CATE on four restricted-mean time scales, horizon τ = 28,
evaluated at each sampled unit's (X1, X2, X3) by numerical integration
(`truth_individual()`). The ATE is the sample mean of the CATE. With F_k the
cumulative incidence of cause k and S the all-cause survival:

| estimand | definition | truth columns | scored for |
|---|---|---|---|
| RMTL_k (restricted mean time lost to cause k) | ∫₀^τ F_k(t \| x, w) dt | `tau_RMTL1`, `tau_RMTL2` | every pseudo-value arm |
| subdistribution RMST_k | τ − RMTL_k | `tau_RMST1`, `tau_RMST2` | `csf_sh` |
| net (cause-specific) RMST_k | ∫₀^τ exp(−Λ_k(t \| x, w)) dt | `tau_RMST1_cs`, `tau_RMST2_cs` | `ipw`, `csf_cs` |
| RMSTc (event-free, composite) | ∫₀^τ S(t \| x, w) dt | `tau_RMSTc` | `ipw`, `csf_cs`, every pseudo-value arm |

Each CATE is the treated minus the control value at the same x.

- **Net RMST is a hypothetical estimand.** It is the RMST in a world where the
  other cause has been removed. The DGM's cause-specific-hazard construction
  makes it well defined in the simulation, but in real data it is identified
  only under independent latent event times, which cannot be checked.
- **Sign.** An RMTL CATE is positive when the event comes sooner; an RMST
  CATE is negative. Each arm is scored against its own truth, so the sign of a
  correlation means the same in every family, but the families' "Event 1" rows
  are different estimands with different spreads (scenario 1: SD 0.37 for
  τ_RMTL1, 0.62 for τ_RMST1_cs). Absolute errors therefore do not compare
  directly across families.
- No HTE tests and no intervals: the study is point estimation only.

## Methods

Implementation: `surv_models.R::all_cate_surv_models()`, called once per run
by `surv_analysis.R` with `n_folds = 10`, `horizon = 28` and
`sl_libraries(500)`. Inputs are X (the 8 covariates as generated, no scaling),
W, the observed time Y and the status D ∈ {0 censored, 1 EOI, 2 CE}.

Fourteen arms (`framework` in the results):

| arm | family | how the CE is handled | targets | fitting | pseudo-values |
|---|---|---|---|---|---|
| `ipw` | IPCW causal forest | censoring, inverse-probability weighted | net RMST1, RMST2; RMSTc | grf-internal (`cf_default`) | — |
| `csf_cs` | causal survival forest | censoring | net RMST1, RMST2; RMSTc | grf-internal | — |
| `csf_sh` | causal survival forest | kept in the risk set | subdistribution RMST1, RMST2 | grf-internal | — |
| `pseudo_cf_whole_oob` | causal forest on pseudo-values | Aalen–Johansen CIF (KM for RMSTc) | RMTL1, RMTL2, RMSTc | grf-internal | whole |
| `pseudo_cf_whole_scf` | 〃 | 〃 | 〃 | single crossfit | whole |
| `pseudo_cf_cvps_scf` | 〃 | 〃 | 〃 | single crossfit | leave-one-fold-out |
| `pseudo_dr_whole_oob` | RF DR-learner | 〃 | 〃 | whole-sample OOB (`oob_oob`) | whole |
| `pseudo_dr_whole_scf` | 〃 | 〃 | 〃 | single crossfit, both stages | whole |
| `pseudo_dr_cvps_scf` | 〃 | 〃 | 〃 | single crossfit, both stages | leave-one-fold-out (training only) |
| `sl_t_whole` | SuperLearner T-learner | 〃 | 〃 | single crossfit | whole |
| `sl_t_cvps` | 〃 | 〃 | 〃 | single crossfit | leave-one-fold-out |
| `sl_t_split` | 〃 | 〃 | 〃 | 3-way split | recomputed on V − 2 folds |
| `sl_dr_whole` | SuperLearner DR-learner | 〃 | 〃 | single crossfit, both stages (`scf_scf`) | whole |
| `sl_dr_cvps` | 〃 | 〃 | 〃 | single crossfit, both stages | leave-one-fold-out (training only) |

A fifteenth, `sl_dr_split` (DR-learner on split pseudo-observations, Cwiling
et al. 2025), is **disabled**: `pseudoyl()` returns NA for the max-time unit
of each split, and nothing guards it (README "Known issues").

All forests are grf at its defaults (as in `sample_size/ADEMP.md`: 2000 trees,
honest, `min.node.size` 5, `mtry` = all 8 covariates here). Every estimated
propensity in a DR-learner is trimmed to [0.05, 0.95] (`trim_ps`); the causal
forests' internal propensities are not.

### Direct survival estimators

- **`ipw`** (`cf_ipw`), following the grf tutorial:
  - Censoring weights: a `survival_forest` of time to censoring
    (D = 0 as the event) on (X, W); each unit's OOB probability of remaining
    uncensored at min(Y, τ), weight 1 / max(p, 0.001). With
    `censoring = FALSE` these are all about 1.
  - For event k, a second survival forest for the time to the *other* cause,
    weighted the same way, and the two weights multiplied. Units that were
    censored or had the competing event (at any time, including after τ) are
    dropped.
  - `causal_forest(X, min(Y, τ), W, sample.weights = w)` on the kept units:
    OOB τ̂ for them, a newdata prediction for the dropped ones.
  - Composite: censoring weights only, every uncensored unit kept.
- **`csf_cs`**: `causal_survival_forest(X, Y, W, 1{D = k}, target = "RMST",
  horizon = 28)` on the whole sample, the CE treated as censoring (grf
  estimates the covariate-dependent censoring itself). Composite: event =
  either cause. OOB predictions.
- **`csf_sh`**: event indicator 1{D = k}, with a CE's time moved to τ + 1, so
  that the unit stays in the risk set and never has the EOI before the
  horizon. That targets τ − RMTL_k. If any unit is censored before τ, every
  censored unit is dropped and the rest weighted by the censoring weights
  above; dropped units get a newdata prediction. Otherwise, whole sample, no
  weights.

### Pseudo-values

- **whole** (`pseudo_all`): `pseudo::pseudoyl(Y, D, 28)` gives jackknife
  pseudo-values of RMTL1 and RMTL2 from the Aalen–Johansen estimator;
  `pseudo::pseudomean(Y, 1{D ∈ {1, 2}}, 28)` gives RMSTc from the
  Kaplan–Meier one. One fit on all n. No covariates are needed because
  censoring is independent of X and W. There is no NA guard (route C).
- **cvps** (`pseudo_crossfit`): the same, for each fold k, on the rows outside
  k. An n × V matrix, NA on fold k's own rows: a jackknife pseudo-value for a
  unit needs that unit in the sample. NAs from `pseudoyl()`/`pseudomean()` are
  replaced by the whole-sample value. That leaks the held-out fold, so it is
  counted (`results$pseudos$cf_cv$n_na_fallback`).

`cvps` + OOB is not a cell: the cvps matrix exists only inside a fold loop.
`whole_scf` is the control that separates the pseudo-value factor from the
fitting factor.

### Pseudo-value causal forests

`causal_forest(X, θ, W)` with the pseudo-value θ as the outcome, grf's own
Y.hat and W.hat. `whole_oob`: whole sample, OOB τ̂. `*_scf`: per fold k, fit
on the rows outside k (θ = the whole vector or column k of the cvps matrix),
predict fold k.

### DR-learners

DR pseudo-outcome (`dr_pseudo`), with the pseudo-value θ in place of Y:

$$\hat\phi_i = \hat\mu_1(x_i) - \hat\mu_0(x_i) + \frac{(\theta_i - \hat\mu_{W_i}(x_i))(W_i - \hat e(x_i))}{\hat e(x_i)(1 - \hat e(x_i))}$$

The outcome model is fit separately in each arm (a T-learner), as in
`R/cate_models.R`.

- **`pseudo_dr_whole_oob`**: `t_learner_rf` (one regression forest per arm;
  own-arm predictions OOB, other-arm from a forest that never saw the unit),
  ê from an OOB regression forest of W on X, stage 2 an OOB regression forest
  of φ̂ on X (`stage2_whole_rf`). The production `dr_random_forest`'s
  strategy.
- **`pseudo_dr_*_scf`**: per fold k, per-arm forests and the W forest are fit
  on the rows outside k and predict fold k, giving φ̂ there. Stage 2
  (`stage_2_rf_scf`) is a regression forest on the rows outside k, predicting
  k.
- **`sl_dr_*`**: as `*_scf`, with SuperLearner for the per-arm outcome models,
  the propensity and stage 2 (`stage_2_sl`), on the same folds.
- **The correction term always uses the whole-sample θ.** cvps has no
  pseudo-value for the held-out rows, so in the DR arms "cvps" changes only the
  pseudo-values the nuisance regressions are *trained* on. Do not report the
  comparison as more than that (README "What the comparison actually
  measures").

### T-learners (SuperLearner)

- **`sl_t_whole`, `sl_t_cvps`**: per fold k, one SuperLearner of θ on X per
  arm, fit on the rows outside k; τ̂ = μ̂1 − μ̂0 on fold k.
- **`sl_t_split`** (after Cwiling et al. 2025, Algorithm 2): per fold k, fold
  k + 1 (wrapping) is set aside as the "KM set" and the other V − 2 folds
  train, with pseudo-values recomputed on those training folds alone; τ̂ on
  fold k. The T-learner never uses the KM set: split pseudo-observations only
  enter a DR correction term, i.e. the disabled arm. So in practice this arm is
  `sl_t_cvps` trained on 8 folds rather than 9. NA training pseudo-values are
  filled from the whole-sample vector and counted (`n_na_fallback`); rows
  still NA are dropped from the fit.

### SuperLearner settings

`sl_libraries(500)`, one library per nuisance (`R/sl_library.R`; learner
definitions in `sample_size/ADEMP.md`):

| nuisance | library | family / meta-learner | internal CV |
|---|---|---|---|
| propensity | mean, glm | binomial / `method.NNloglik` | V = 5 |
| outcome (θ, per arm) | mean, lasso, ranger (min.node.size 25), glm, gam, earth | gaussian / `method.NNLS` | V = 5 |
| CATE (stage 2) | mean, glm, lasso (lambda.min and lambda.1se), gam, ranger (min.node.size 25), lasso on pairwise interactions (both tunings), earth | gaussian / `method.NNLS` | V = 10 (package default) |

- Every fit is pretested (`pretest_superlearner`, 2-fold CV): learners that
  error or give non-finite predictions are dropped, falling back to `SL.mean`.
  On a constant outcome even that fallback fails (route A').
- Predictions come from the fit itself (`sl_fit_predict`); a fit that errors
  predicts its training mean. In the DR nuisances and `sl_t_split`, all-zero
  predictions (every learner weighted 0) are replaced by the training mean of
  that arm (or of W). `sl_t_whole`, `sl_t_cvps` and stage 2 have no such
  check.
- The outcome models are gaussian even though θ for RMTL is bounded on
  [0, 28]: pseudo-values can fall outside the range.

### Saved per run

τ̂ per arm and estimand; for the five DR nuisance arms, `po`, `pseudo.hat`,
`pseudo0.hat`, `pseudo.hat.cf` (a marginal regression of θ on X, saved but not
used by any estimator) and `W.hat`; both sets of pseudo-values with
`n_na_fallback`; fold indices; the data and the truth.

## Performance measures

Per run, per arm and target, against that arm's truth column
(`surv_metrics.R`, `framework_truth_map`), using `cate_metrics()` so the
conventions match the other studies:

- bias (mean of estimate − truth), ATE bias, relative ATE bias, relative CATE
  bias
- MSE, RMSE, MAE
- Pearson and Spearman correlation with the true CATE, sign accuracy
- C-statistic, (Kendall τ_b + 1) / 2: Harrell's C for a continuous outcome
- `n_na` (units with no estimate)

Targets are labelled by event rather than scale ("Event 1" = RMST1 or RMTL1,
"Event 2", "Combined"). Which truth each label means:

| arm | Event 1 | Event 2 | Combined |
|---|---|---|---|
| `ipw`, `csf_cs` | `tau_RMST1_cs` | `tau_RMST2_cs` | `tau_RMSTc` |
| `csf_sh` | `tau_RMST1` | `tau_RMST2` | — |
| pseudo-value arms | `tau_RMTL1` | `tau_RMTL2` | `tau_RMSTc` |

`surv_results.qmd` averages each measure over runs, with Monte Carlo SE
sd / √(non-NA runs), plotted as mean ± 1.96 MCSE. "Combined" is dropped from
the report.

**Cells with a constant truth.** `tau_RMST1_cs` ≡ 0 in scenarios 2 and 5, and
`tau_RMST2_cs` ≡ 0 in scenarios 1 and 3. In those `ipw` / `csf_cs` cells the
correlations and the C-statistic are NA (`surv_metrics.R` warns "the standard
deviation is zero"), relative biases are NA, and sign accuracy is about 0,
since `sign(true) = 0` matches no non-zero estimate. Drop them; they are not
failures. The README notes this for scenario 5's Event 1 only.

**Scenario 1 is treated as null, but is not.** `cate_metrics()` and
`c_statistic()` hard-code scenario 1 as the no-heterogeneity scenario (the
`sample_size/` convention): Pearson and Spearman are set to 0 and the
C-statistic to 0.5. Here scenario 1's CATE varies with X1 and X2 (SD 0.37 for
τ_RMTL1, 0.62 for τ_RMST1_cs), so those three columns are placeholders in
scenario 1, not measurements. Bias and the error measures are unaffected.
Fixing it needs only a metrics rerun, not a simulation rerun.

**Failed runs.** Runs that error produce no results and drop out of every
summary. `jobscripts/failed_ids.txt` lists 169 indices of 7,000 at the last
check. The deterministic failures are not random: route C fires only with
nobody observed past the horizon, which happens with censoring on and in
scenarios 1, 3, 4, 6 and 7. So the summaries for those cells are conditional
on someone surviving to 28. Check the coverage table in `surv_results.qmd`
before reading anything else.

**Diagnostics, not performance measures:** `n_na_fallback` (how much the cvps
and split arms leaned on whole-sample pseudo-values), and the nuisance
summaries from `surv_nuisance_extract.R` / `surv_nuisance_figs.R`: propensity
overlap and trimming rate, and the rate of extreme DR pseudo-outcomes.
