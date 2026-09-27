# ADEMP — continuous outcome, sample size study

## Aims

Assess how well DR-learner and causal forest estimators recover the CATE in a
two-arm RCT with a continuous outcome, across CATE structures and sample sizes,
and the size and power of the post-estimation heterogeneity tests.

## Data-generating mechanisms

- `W ~ Bernoulli(0.5)`
- Prognostic: `X1 ~ Bernoulli(0.4)`, `X2 ~ N(0, 1)`
- Potential effect modifiers (drawn in every scenario, whether or not τ uses them): `X3 ~ Bernoulli(0.7)`, `X4, X5 ~ N(0, 1)`
- Noise: `X01`–`X03 ~ N(0, 1)`; `X04`, `X05` indicators of a 3-level factor (0.45 / 0.30 / 0.25)
- `Y = 0.4 − 0.5·X1 + X2 + W·τ(x) + ε`, `ε ~ N(0, 0.5²)`
- `τ(x) = bW + g(x)`; `bW` set so the true ATE gives 80% power in an
  unadjusted two-sample t-test (ATE −0.65, −0.41, −0.29, −0.20 at n = 100, 250, 500, 1000)

| scenario | g(x) |
|---|---|
| 1 | 0 (no HTE) |
| 2 | −X4 |
| 3 | 2·X3 + 0.5·X4 − 0.5·X4·X5 |
| 4 | cos(X4) |
| 5 | 2·X3 |
| 6 | 0.3·X3 − X4 |
| 7 | X3·X4 |
| 8 | 2·X3 + 0.5·X4 − 0.5·X3·X4 |
| 9 | −0.5·X4·X5 |
| 10 | 0.3·X3 + 0.1·exp(−\|X4\|) |

Full factorial: scenario (1–10) × n (100, 250, 500, 1000).
Repetitions: 100 per cell; 500 for scenarios 1–4 (the reported null, simple,
complex and non-linear scenarios).

## Estimands

- Unit-level CATE `τ(x_i)` for every unit in the simulated sample.
- ATE, as the sample mean of the CATE.
- Presence of heterogeneity (H0: τ constant), for the tests.

## Methods

| method | description |
|---|---|
| `causal_forest` | grf `causal_forest`, internal cross-fitting |
| `dr_random_forest` | DR-learner; outcome model fit per arm (T-learner RF, own-arm predictions OOB), RF propensity, RF second stage |
| `dr_superlearner` | DR-learner; SuperLearner per-arm outcome models and propensity, single leave-one-fold-out crossfit (V = 4, 5, 10 at n = 100, 250, ≥ 500), SuperLearner second stage |
| `dr_oracle` | DR-learner with the true outcome model and propensity 0.5 |
| `dr_semi_oracle` | DR-learner with the known propensity 0.5 and `dr_random_forest`'s per-arm outcome forests |

SuperLearner libraries, one per nuisance (`R/sl_library.R::sl_libraries`),
chosen without reference to the scenarios' HTE shapes. The outcome model is fit
in each arm separately, so its library is sized to about half the training
rows (~37 per arm at n = 100):

| | n = 100 | n = 250 | n ≥ 500 |
|---|---|---|---|
| propensity | mean, glm | same | same |
| outcome (per arm) | mean, lasso, ranger (min.node.size 25) | + glm, gam | + earth |
| CATE (stage 2) | mean, glm, lasso (lambda.min and lambda.1se), gam, ranger (min.node.size 25) | + lasso on pairwise interactions (both tunings), earth | same |

Each stage-2 lasso is included at both tunings and SuperLearner's CV weighs
them, rather than fixing the tuning in advance.

Learners that error or give non-finite predictions on a fold are dropped
before fitting (`pretest_superlearner`) and saved as `sl_dropped`.
Heterogeneity tests per method: BLP (GenericML) and independence tests on the
CATE and on the pseudo-outcome. Also run on the true CATE and nuisances
(`cts_true_cate_tests.RDS`).

## Performance measures

Per run, over the units in the sample (`R/metrics.R::cate_metrics`), then
averaged over runs:

- bias (mean of estimate − truth), ATE bias, relative ATE bias, relative CATE bias
- MSE, RMSE, MAE
- Pearson and Spearman correlation with the true CATE (set to 0 in scenario 1)
- sign accuracy
- HTE tests: p-values `BLP_p`, `indep_cate`, `indep_po` → rejection rates
  (type I error in scenario 1, power in 2–10)
