# ADEMP — binary outcome, sample size study

## Aims

Assess how well DR-learner and causal forest estimators recover the CATE, on
the risk-difference scale, in a two-arm RCT with a binary outcome, across CATE
structures and sample sizes, and the size and power of the heterogeneity tests.
Same design as `continuous/`.

## Data-generating mechanisms

- Covariates as in `continuous/`: `W ~ Bernoulli(0.5)`, `X1 ~ Bernoulli(0.4)`,
  `X2 ~ N(0, 1)`, `X3 ~ Bernoulli(0.7)`, `X4, X5 ~ N(0, 1)`, noise `X01`–`X05`
- `Y ~ Bernoulli(m0(x) + W·τ(x))`
- Control risk `m0(x) = 0.34 + 0.36·plogis(−1.72 − 0.5·X1 + X2)`; E[m0] = 0.40
- `τ(x) = bW + g(x)`; `g` is the continuous modifier with signs reversed,
  scaled by `RD_SCALE`, X4 and X5 through `tanh` (t4 = tanh(X4), t5 = tanh(X5))
- `bW` set so the true ATE (marginal RD) gives 80% power in a two-proportion
  test (RD −0.248, −0.164, −0.118, −0.085 at n = 100, 250, 500, 1000)
- Every treated risk stays within [0.01, 0.99]

| scenario | g(x) |
|---|---|
| 1 | 0 (no HTE) |
| 2 | 0.082·t4 |
| 3 | −0.102·X3 − 0.026·t4 + 0.026·t4·t5 |
| 4 | −0.204·cos(X4) |
| 5 | −0.274·X3 |
| 6 | −0.023·X3 + 0.075·t4 |
| 7 | −0.082·X3·t4 |
| 8 | −0.274·X3 − 0.069·t4 + 0.069·X3·t4 |
| 9 | 0.082·t4·t5 |
| 10 | −0.179·X3 − 0.060·exp(−\|X4\|) |

Full factorial: scenario (1–10) × n (100, 250, 500, 1000).
Repetitions: 100 per cell; 500 for scenarios 1–4.

## Estimands

- Unit-level CATE `τ(x_i) = P(Y=1 | x, W=1) − P(Y=1 | x, W=0)` for every unit in the sample.
- ATE (marginal risk difference), as the sample mean of the CATE.
- Presence of heterogeneity (H0: τ constant), for the tests.

## Methods

As `continuous/` (`family = binomial()` for SuperLearner):

| method | description |
|---|---|
| `causal_forest` | grf `causal_forest`, internal cross-fitting |
| `dr_random_forest` | DR-learner; S-learner RF nuisances, OOB predictions, RF second stage |
| `dr_superlearner` | DR-learner; SuperLearner, single leave-one-fold-out crossfit (V = 4, 5, 10 at n = 100, 250, ≥ 500) |
| `dr_oracle` | DR-learner with the true risk model and propensity 0.5 |
| `dr_semi_oracle` | DR-learner with the known propensity 0.5 only |

Heterogeneity tests per method: BLP and independence tests on the CATE and on
the pseudo-outcome. Also run on the true CATE (`bin_true_cate_tests.RDS`).

## Performance measures

As `continuous/`:

- bias, ATE bias, relative ATE bias, relative CATE bias
- MSE, RMSE, MAE
- Pearson and Spearman correlation with the true CATE
- sign accuracy
- HTE test rejection rates (`BLP_p`, `indep_cate`, `indep_po`)
