# ADEMP — model evaluation

## Aims

In a two-arm RCT, decide whether cheap proxy loss functions can choose among
candidate CATE learners as well as the true PEHE would. True PEHE cannot be
computed outside a simulation. The proxies are built from nuisance estimators
that are fit separately from the candidates.

The study compares proxies on three axes:

- the score family: influence-corrected PEHE, DR risk, and DR calibration;
- whether the proxy's propensity is estimated or fixed at the true 0.5;
- the data the scoring nuisance is fit on, relative to the data each
  candidate was trained on.

This is a model-*selection* study. The object of interest is how each proxy
ranks the 9 candidates within a replicate, not how accurate any one
candidate is.

## Data-generating mechanisms

The continuous outcome of the correlated sample-size studies,
`set = "continuous_corr_0.5"` in `R/dgm_scenarios.R` (see
`sample_size/correlated/ADEMP.md` for the copula and the calibration of bW):

- `Y = 0.4 − 0.5·X1 + X2 + W·τ(x) + ε`, `ε ~ N(0, 0.5²)`, `τ(x) = bW + g(x)`;
- `W ~ Bernoulli(0.5)`, independent of X, so the true propensity is 0.5;
- 10 covariates (X1–X5, X01–X05) from a Gaussian copula with exchangeable
  latent correlation ρ = 0.5 (`ME_RHO`); every learner sees all ten;
- bW is set so that the true ATE gives an unadjusted test 80% power at that n,
  planned as if the effect were homogeneous.

| scenario | description | g(x) |
|---|---|---|
| 1 | null | 0 |
| 4 | non-linear | cos(X4) |
| 6 | two variables | 0.3·X3 − X4 |
| 8 | single effects + interaction | 2·X3 + 0.5·X4 − 0.5·X3·X4 |

This is the same scenario subset as `crossfitting/`, which runs it on
independent covariates.

### Design

Full factorial, in three result trees:

| tree | script | scenarios | n | runs | array jobs |
|---|---|---|---|---|---|
| main | `me_analysis.R` | 1, 4, 6, 8 | 250, 500, 1000 | 100 | 1200 |
| strategies | `me_strategies.R` | the same | the same | the same | 1200 (reads the main tree) |
| split | `me_split.R` | 1, 4, 6, 8 | 500, 1000 | 100 | 800 |

- Each run index has its own L'Ecuyer-CMRG stream (`setup_rng_stream(run)`).
- One fold assignment (V = 10) is drawn per replicate. The candidates and the
  fold-based scoring arms share it.
- The strategies tree does not regenerate data or refit candidates. It reads
  the main tree's saved data, truth, folds and per-candidate τ̂, and adds two
  scoring arms. Every arm in it therefore scores the same candidate fits.
- The split tree leaves out n = 250, because 50 evaluation rows are too few to
  rank 9 candidates.

## Estimands

- The unit-level CATE τ(xᵢ) for the rows of the sample.
- For each candidate m, its true PEHE: `mean((τ̂_m(xᵢ) − τ(xᵢ))²)`. This is
  taken over all n rows, or over the 20% evaluation rows in the split tree.
- The estimand of the selection problem is the ordering of the 9 candidates
  by true PEHE, and in particular the candidate with the lowest true PEHE.

## Methods

### Candidates (what is being selected among)

Each candidate is a single-crossfit DR-learner (V = 10, `scf_scf`):

- outcome models are fit per arm (T-learner);
- propensities are trimmed to [0.05, 0.95];
- the stage-2 regression is fit on the same folds as the nuisances.

| id | family | outcome model and stage 2 | propensity |
|---|---|---|---|
| `rf1` | ranger | defaults (`mtry = floor(sqrt(p))`) | own forest |
| `rf2` | ranger | `mtry = p`, `max.depth = 5` | own forest |
| `rf3` | ranger | `mtry = ceiling(p/2)`, `max.depth = 3` | own forest |
| `net1` | glmnet | lasso (`alpha = 1`) | own glmnet |
| `net2` | glmnet | lasso on 10 main effects + 45 pairwise products | own glmnet |
| `net3` | glmnet | elastic net (`alpha = 0.5`) | own glmnet |
| `SL1` | SuperLearner | production `dr_superlearner` (`sl_libraries(n)`) | mean + glm |
| `SL2` | SuperLearner | glmnet, xgboost, cforest, earth, gam, mean | mean + glm |
| `SL3` | SuperLearner | svm, nnet, mean (deliberately weak) | mean + glm |

`p = 10`. In the split tree the candidates are refit on the 80% training
split, still crossfit within it and with the same hyperparameters, and predict
only on the 20%.

### Proxies (the selection methods)

A proxy is one combination of four factors. It selects the candidate with the
lowest score.

**1. Scoring-nuisance pipeline.** Each pipeline estimates `mu_DR` (one model
on X and W, read at W = 0 and W = 1) and `pi` (a model on X):

- **XGBoost** (`xgb`): a 36-point grid over `eta`, `max_depth`, `subsample`
  and `colsample_bytree`, tuned by 5-fold CV with up to 100 rounds and early
  stopping;
- **H2O AutoML** (`automl`): up to 20 models, with DeepLearning and XGBoost
  excluded.

From these, `calculate_pseudos()` forms the AIPW pseudo-outcome
`phi = mu1 − mu0 + (Y − mu_W)(W − pi) / (pi(1 − pi))`, where `mu_W` is
`mu_DR` read at the observed W. `phi05` is the same
quantity with `pi` fixed at 0.5. Neither is trimmed.

**2. Nuisance arm.** The arm sets which rows the scoring nuisance is fit on,
relative to the rows each candidate was trained on:

| arm | tree | nuisance fit on | predicted on | row-honest | decoupled from candidate |
|---|---|---|---|---|---|
| `whole` | main, strategies | all n rows | all n rows | no | no |
| `cv_shared` | strategies | the candidate's V − 1 training folds | its held-out fold | yes | no |
| `holdout` | strategies | the candidate's held-out fold only | that fold | no | yes |
| `split` | split | the 20% evaluation rows | the same 20% | no | yes, and the candidates never see those rows |

- At n = 250, `holdout` pools adjacent pairs of folds into 5 blocks of 50
  rows (`holdout_blocks()`). A pooled block is therefore only half decoupled
  from any one candidate fold model.

**3. Score family.** Each family is a loss (lower is better), computed on the
evaluation rows:

| family | score | definition |
|---|---|---|
| influence | `infl` | the influence-corrected PEHE estimate (`calc_infl_score`) |
| DR risk | `dr` | `mean((τ̂ − phi)²)` (`calc_dr_risk`) |
| calibration | `calqK` | `Σ_k \|G_k\| · \|mean_{G_k}(τ̂) − mean_{G_k}(phi)\|`, over K rank groups of τ̂, for K = 5 and 10 (`calc_cal_score`) |

- Calibration measures calibration, not discrimination. A candidate that
  predicts a constant scores close to zero on it.

**4. Propensity.** The proxy uses either the pipeline's estimated `pi` or the
true 0.5. The true value is used through `phi05` and `pi = 0.5`, and gives the
`*05` columns.

This makes 8 score types per (arm × pipeline): 17 score columns per
candidate-run in the main and split trees, and 49 in the strategies tree. Each
tree also carries `true_pehe` (`me_metrics.R`).

## Performance measures

Computed per (scenario, n, run, proxy) over the candidates scored in that run
(`me_results.qmd`). Each one is then averaged over runs and reported with its
Monte Carlo SE, `sd / sqrt(#non-missing)`.

- **Rank agreement:** the Spearman and Kendall correlations between the
  proxy's scores and true PEHE across the 9 candidates.
- **Top-1 selection accuracy:** the proportion of runs in which the proxy's
  pick is the candidate with the lowest true PEHE.
- **True rank of the pick:** where the selected candidate sits in the
  true-PEHE ordering, from 1 to 9.
- **Regret:** `true_pehe[pick] − min(true_pehe)`.
- **Relative regret:** regret divided by `min(true_pehe)`.

Regret has two reference points:

- **random pick:** `mean(true_pehe) − min(true_pehe)`, the expected regret of
  ignoring the scores;
- **best fixed choice:** the mean regret of always using the single candidate
  with the lowest mean regret in that (scenario, n). It is chosen with
  hindsight, so it is an optimistic bound.

Two descriptive summaries are also reported: how often each candidate is truly
best, and the spread `max − min` of true PEHE.

The report reads the strategies tree. It covers the `whole`, `cv_shared` and
`holdout` arms and leaves out the calibration family. It also reports
completeness: the runs missing per cell, and any ranking that scored fewer
than 9 candidates. The split tree is scored by `me_metrics.R split`.
