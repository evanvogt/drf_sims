# ADEMP — competing risks

One study, `competing_risk/` (`surv_*` scripts). The DGM is its own
(`surv_dgm.R`), not a `sample_size/` one, and most estimators are specific to
it (`surv_models.R`), borrowing the shared DR-learner and SuperLearner
machinery from `R/cate_models.R` and `R/sl_library.R`.

## Aims

- Assess how well forest- and SuperLearner-based estimators recover the
  unit-level CATE on restricted-mean time scales in a two-arm RCT where the
  event 1 (E1) can be pre-empted by a competing event 2 (E2), with and
  without right censoring.
- Compare three ways of handling the competing event: as censoring
  (cause-specific / net), kept in the risk set (subdistribution), or through
  pseudo-values of the cumulative incidence (restricted mean time lost, RMTL).
- Show what an analysis of one event reports when the treatment acts only on
  the other (scenarios 2, 5 and 6): an apparent effect on event 1 that is
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

Cause-specific hazards are Weibull and proportional, cause k = 1 (E1) or
2 (E2):

$$\lambda_k(t \mid x, w) = \frac{a_k}{s_k}\left(\frac{t}{s_k}\right)^{a_k - 1}, \qquad \log s_k = \log s_{k0} + b_{1k} X_1 + b_{2k} X_2 + w\,(b_{Wk} + b_{3k} X_3)$$

| | shape a_k | baseline scale s_k0 | b_1k (X1) | b_2k (X2) |
|---|---|---|---|---|
| E1 (k = 1) | 2 | 25 | −0.1 | 0.1 |
| E2 (k = 2) | 1.1 | 45 | −0.1 | 0.1 |

Effects are on the log *scale*, so a shift b multiplies the hazard by
exp(−a_k·b): X1 raises the E1 hazard by 1.22 and the E2 hazard by 1.12; one
SD of X2 lowers them by 0.82 and 0.90.

Generation follows Beyersmann et al. (2009): the event time T solves
Λ1(T) + Λ2(T) = −log U, `U ~ Uniform(0, 1)`, by `uniroot` on (0, 200)
(`find_time()`), and the cause is `Bernoulli(λ1(T) / (λ1(T) + λ2(T)))`.

**Censoring.** Administrative at 180 in every run. With `censoring = TRUE`,
also `C ~ Uniform(1, 180)`, independent of X and W. Only censoring before the
horizon of 28 matters to the estimands.

### Scenarios

Treatment effects are set as log hazard ratios, ±0.35 for the ATE and ±0.7
more when X3 = 1, and divided by the shape to get the log-scale coefficients
(`surv_dgm.R` writes them as `c(...) / 2` and `c(...) / 1.1`), so both events
get the same log-HRs:

| scenario | description | b_W1 | b_31 | b_W2 | b_32 | E1 HR, X3 = 0 / 1 | E2 HR, X3 = 0 / 1 |
|---|---|---|---|---|---|---|---|
| 1 | ATE on E1 only | −0.175 | 0 | 0 | 0 | 1.42 / 1.42 | 1 / 1 |
| 2 | ATE on E2 only | 0 | 0 | 0.318 | 0 | 1 / 1 | 0.70 / 0.70 |
| 3 | HTE on E1, no ATE on E2 | −0.175 | −0.35 | 0 | 0 | 1.42 / 2.86 | 1 / 1 |
| 4 | HTE on E1, ATE on E2 | −0.175 | −0.35 | 0.318 | 0 | 1.42 / 2.86 | 0.70 / 0.70 |
| 5 | HTE on E2, no ATE on E1 | 0 | 0 | 0.318 | 0.636 | 1 / 1 | 0.70 / 0.35 |
| 6 | HTE on E2, ATE on E1 | −0.175 | 0 | 0.318 | 0.636 | 1.42 / 1.42 | 0.70 / 0.35 |
| 7 | HTE on both | −0.175 | −0.35 | 0.318 | 0.636 | 1.42 / 2.86 | 0.70 / 0.35 |

The treatment *raises* the E1 hazard and lowers the E2 hazard. "ATE" and
"HTE" in the descriptions refer to the hazard scale, i.e. whether X3 modifies
the effect. On the time scales the estimands use, every effect is
heterogeneous. X1 and X2 shift the baseline hazards, and a restricted mean is
not linear in them, so the CATE varies with X1 and X2 even in scenarios 1 and
2 (see "What the DGM implies" below).

**Why these values (retuned 2026-10-01).** Until 2026-10-01 the E1 scale was
15 and every effect was ±0.7 on the log scale: E1 HR 4.06 / 16.4, E2 HR
0.46 / 0.21. That was correct on the hazard scale, but on the RMTL scale,
which 12 of the 14 arms target, the scenarios ran together. Two things caused
it. The same coefficient times shape 2 against 1.1 gave E1 a log-HR 1.8
times E2's. And 98% of controls had an event by the horizon, so the two
events competed heavily and an effect on either one moved both CIFs. Scenario 3
("no ATE on E2") then had more X3-driven RMTL2 heterogeneity than scenario 6
("HTE on E2"), and scenarios 4 and 6 could not be told apart. The fix was a
smaller, log-HR-matched set of effects and a lower E1 hazard (scale 25),
which cuts the competition.

How well the scenarios separate on the RMTL scale is measured as the smallest
X3 difference in the true CATE on an event the label says is heterogeneous,
divided by the largest *induced* X3 difference on an event the label says is
not (population values, X3 = 1 minus X3 = 0, in days):

| scenario | RMTL1 X3 diff, old | RMTL1, now | RMTL2 X3 diff, old | RMTL2, now |
|---|---|---|---|---|
| 3 | 4.15 | 3.97 | −1.29 (induced) | −0.85 (induced) |
| 4 | 3.68 | 4.10 | −0.63 (induced) | −0.65 (induced) |
| 5 | 0.70 (induced) | 0.68 (induced) | −1.17 | −1.99 |
| 6 | 0.57 (induced) | 0.80 (induced) | −0.69 | −1.87 |
| 7 | 4.01 | 5.05 | −0.99 | −2.22 |
| **separation** | 5.3 | 4.9 | **0.53** | **2.2** |

The induced differences can't be removed: in competing risks an effect on one
cause always moves the other cause's CIF. What the retune achieves is that every
labelled effect is now at least twice the size of any induced one.

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
results report's text (`surv_results.qmd`) still says 100 runs. **All results
produced before 2026-10-01 come from the old parameters** (see "Why these
values") and have to be regenerated, with the full 7,000-job array.

**Seeding.** Each run is seeded by its run index alone
(`setup_rng_stream(run)`). Draw order: `W, X1, X2, X3, U, cause, [C], X01–X03,
cats`. So within a run all 14 (scenario, censoring) cells share W, X1–X3 and
the uniform the event time is solved from. Scenarios differ only by their
parameters, so comparisons across scenarios are paired. `censoring = TRUE`
consumes n more uniforms before the noise covariates, so X01–X05 differ
between the two censoring settings of a run.

### What the DGM implies

For the parameters as of 2026-10-01. The event mix comes from
`generate_surv_data()` at n = 20,000 per (scenario, censoring) cell, seeded
once rather than with the study's streams, so it is illustrative. The truths
are population values: `truth_individual()` evaluated over the covariate
distribution (X1 and X3 exactly, X2 on 25 normal quantiles).

**Event mix by the horizon.** The control arm is the same in every scenario:
E1 0.52–0.54, E2 0.32–0.34, still event-free at 28 about 0.14. With censoring
on, that becomes E1 0.48–0.50, E2 0.30–0.31, censored before 28 about 0.08,
and about 0.12 observed past 28. Treated arm:

| scenario | E1 by 28 | E2 by 28 | event-free at 28 | censored before 28 (censoring on) | observed past 28 (censoring on) |
|---|---|---|---|---|---|
| 1 | 0.63 | 0.29 | 0.083 | 0.073 | 0.071 |
| 2 | 0.58 | 0.25 | 0.168 | 0.091 | 0.139 |
| 3 | 0.72 | 0.24 | 0.035 | 0.058 | 0.031 |
| 4 | 0.78 | 0.18 | 0.038 | 0.062 | 0.036 |
| 5 | 0.63 | 0.17 | 0.197 | 0.099 | 0.165 |
| 6 | 0.74 | 0.15 | 0.114 | 0.081 | 0.100 |
| 7 | 0.83 | 0.13 | 0.045 | 0.068 | 0.039 |

- The horizon still truncates most of the restricted means, but 4–20% of each
  arm is event-free at 28. Censoring now matters more than it did: 6–10% of
  an arm is censored before the horizon.
- Every cell has units observed past 28: about 30 controls per run and at
  least about 8 treated, even in scenario 3 with censoring on. The chance that
  nobody is past 28 in a run of 500 is below 1e-17 in every cell, so route C
  (README "Known issues", NA whole-sample pseudo-values) should not recur. It
  needed scenarios where no treated unit reached 28.
- The fewest treated E2s are in scenario 7, about 13% or roughly 32 units, against 3% before.
  So the near-constant treated-arm RMTL2 pseudo-values behind route A' should
  be much rarer.

**The true CATEs.** Population mean (SD), in days. Control-arm levels:
RMTL1 6.9, RMTL2 5.7, RMSTc 15.4 (they sum to 28), net RMST1 19.1, net
RMST2 21.2.

| scenario | τ_RMTL1 | τ_RMTL2 | τ_RMSTc | τ_RMST1_cs | τ_RMST2_cs | cor(τ_RMTL1, τ_RMTL2) |
|---|---|---|---|---|---|---|
| 1 | 1.82 (0.09) | −0.35 (0.07) | −1.47 (0.03) | −2.16 (0.10) | 0 | −0.96 |
| 2 | 0.51 (0.11) | −1.49 (0.08) | 0.97 (0.03) | 0 | 1.74 (0.14) | −0.99 |
| 3 | 4.59 (1.82) | −0.94 (0.42) | −3.65 (1.43) | −5.24 (2.02) | 0 | −0.94 |
| 4 | 5.29 (1.89) | −2.21 (0.37) | −3.08 (1.58) | −5.24 (2.02) | 1.74 (0.14) | −0.86 |
| 5 | 0.99 (0.38) | −2.88 (0.93) | 1.89 (0.60) | 0 | 3.41 (1.14) | −0.91 |
| 6 | 2.98 (0.48) | −3.06 (0.88) | 0.08 (0.50) | −2.16 (0.10) | 3.41 (1.14) | −0.90 |
| 7 | 5.96 (2.33) | −3.31 (1.05) | −2.65 (1.30) | −5.24 (2.02) | 3.41 (1.14) | −0.99 |

τ_RMST1 and τ_RMST2 (subdistribution) are −τ_RMTL1 and −τ_RMTL2 exactly.

- **An effect on one cause moves the other's CIF.** In scenarios 2 and 5 the
  treatment does not touch the E1 hazard, yet τ_RMTL1 is +0.5 and +1.0: fewer
  E2s leave more units at risk of E1. The net truth τ_RMST1_cs is
  identically 0 there. Likewise τ_RMTL2 is −0.35 and −0.94 in scenarios 1 and
  3, where τ_RMST2_cs is identically 0.
- **Scenario 5's RMTL1 heterogeneity is all induced** (cor −0.91 with
  τ_RMTL2), so an estimator that tracks the real, event-2, signal lands
  anti-correlated with τ_RMTL1. See README "Scenario 5 flips the sign".
- **Where X3 acts.** See the separation table under "Why these values": on
  RMTL1 the X3 difference is 4.0–5.1 in scenarios 3, 4 and 7 against an induced
  0.7–0.8 in 5 and 6. On RMTL2 it is 1.9–2.2 in scenarios 5, 6 and 7 against
  an induced 0.65–0.85 in 3 and 4. On the net scale, X3 acts only where the
  labels say.
- **"ATE" scenarios are close to null on the time scale, but not null.**
  Scenario 1's τ_RMTL1 has SD 0.09 and τ_RMST1_cs 0.10, all from X1 and X2
  (see "Performance measures").
- **Scenario 6's composite effect averages out.** The E1 and E2 effects
  offset on RMSTc, so τ_RMSTc has mean 0.08 but changes sign with X3: −0.67
  at X3 = 0 and +0.39 at X3 = 1.

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
  are different estimands with different spreads (scenario 3: SD 1.82 for
  τ_RMTL1, 2.02 for τ_RMST1_cs). Absolute errors therefore do not compare
  directly across families.
- No HTE tests and no intervals: the study is point estimation only.

## Methods

Implementation: `surv_models.R::all_cate_surv_models()`, called once per run
by `surv_analysis.R` with `n_folds = 10`, `horizon = 28` and
`sl_libraries(500)`. Inputs are X (the 8 covariates as generated, no scaling),
W, the observed time Y and the status D ∈ {0 censored, 1 E1, 2 E2}.

Fourteen arms (`framework` in the results):

| arm | family | how the competing event is handled | targets | fitting | pseudo-values |
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
  horizon = 28)` on the whole sample, the competing event treated as censoring (grf
  estimates the covariate-dependent censoring itself). Composite: event =
  either cause. OOB predictions.
- **`csf_sh`**: event indicator 1{D = k}, with a competing event's time moved to τ + 1, so
  that the unit stays in the risk set and never has event k before the
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
C-statistic to 0.5. Here scenario 1's CATE varies with X1 and X2 (SD 0.09 for
τ_RMTL1, 0.10 for τ_RMST1_cs), so those three columns are placeholders in
scenario 1, not measurements. Bias and the error measures are unaffected.
Fixing it needs only a metrics rerun, not a simulation rerun.

**Failed runs.** Runs that error produce no results and drop out of every
summary. Under the pre-2026-10-01 parameters, `jobscripts/failed_ids.txt`
listed 169 indices of 7,000. Those failures were not random: route C fires
only when nobody is observed past the horizon, which happened with censoring on
and in scenarios 1, 3, 4, 6 and 7, so those cells' summaries were conditional
on someone surviving to 28. Under the current parameters route C should not
fire (see "What the DGM implies"), but check the coverage table in
`surv_results.qmd` before reading anything else.

**Diagnostics, not performance measures:** `n_na_fallback` (how much the cvps
and split arms leaned on whole-sample pseudo-values), and the nuisance
summaries from `surv_nuisance_extract.R` / `surv_nuisance_figs.R`: propensity
overlap and trimming rate, and the rate of extreme DR pseudo-outcomes.
