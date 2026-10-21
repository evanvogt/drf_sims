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
`sample_size/ADEMP.md` for the copula and the calibration of bW):

- `Y = 0.4 − 0.5·X1 + X2 + W·τ(x) + ε`, `ε ~ N(0, 0.5²)`, `τ(x) = bW + g(x)`;
- `W ~ Bernoulli(0.5)`, independent of X, so the true propensity is 0.5;
- 10 covariates (X1–X5, X01–X05) from a Gaussian copula with exchangeable
  latent correlation ρ = 0.5 (`ME_RHO`); every learner sees all ten;
- bW is set so that the true ATE gives an unadjusted test 80% power at that n,
  planned as if the effect were homogeneous.

| scenario | description | g(x) |
|---|---|---|
| 1 | null | 0 |
| 2 | simple, one continuous modifier | −X4 |
| 3 | complex: single effects + interaction | 2·X3 + 0.5·X4 − 0.5·X4·X5 |
| 4 | non-linear | cos(X4) |

These are the main four scenarios that the other studies report; the question
does not need every CATE structure. The grid lists them in the order 1, 4, 2,
3, which keeps the array indices of scenarios 1 and 4 unchanged
(`me_config.R`).

### Design

Full factorial, in three result trees:

| tree | script | scenarios | n | runs | array jobs |
|---|---|---|---|---|---|
| main | `me_analysis.R` | 1–4 | 250, 500, 1000 | 100 | 1200 |
| strategies | `me_strategies.R` | the same | the same | the same | 1200 (reads the main tree) |
| split | `me_split.R` | 1–4 | 500, 1000 | 100 | 800 |

- 100 runs, not `sample_size/`'s 500: each replicate fits 9 candidates plus
  two nuisance pipelines (a 36-point XGBoost CV grid and up to 20 AutoML
  models) across the arms and both nuisance targets.
- Each run index has its own L'Ecuyer-CMRG stream (`setup_rng_stream(run)`).
- One fold assignment (V = 10) is drawn per replicate. The candidates and the
  fold-based scoring arms share it.
- The strategies tree does not regenerate data or refit candidates. It reads
  the main tree's saved data, truth, folds and per-candidate τ̂ and fits all
  three scoring arms, `whole` included, so every arm in it scores the same
  candidate fits. It is the tree the report scores; the main tree's own
  `whole` arm is superseded.
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
| `SL3` | SuperLearner | svm, nnet, mean (a weak baseline) | mean + glm |

`p = 10`. In the split tree the candidates are refit on the 80% training
split, still crossfit within it and with the same hyperparameters, and predict
only on the 20%.

Why the candidates are set up this way:

- **One crossfitting scheme for all 9.** `rf1–3` could take whole-sample OOB
  predictions, but glmnet has no bagging and SuperLearner's internal CV
  selects weights rather than giving an honest prediction of a training row
  (`crossfitting/README.md`, "DR-learner, SuperLearner"). Mixing schemes would
  confound learner choice with crossfitting in the ranking. V = 10 at every n:
  stage 2 trains on 90% of the data, so there is no need to shrink V at
  n = 250.
- **Per-arm outcome models.** One model on `cbind(W, X)` read at W = 0 and 1
  would make the `net*` candidates' `Y1.hat − Y0.hat` a constant (W enters as
  a main effect only).
- **Propensity.** `rf*` and `net*` use their own learner. The SuperLearner
  candidates use the production library, mean + glm: in an RCT at 0.5 a
  flexible learner only fits noise that the weight `1/(e(1 − e))` amplifies.
- **`SL1` is the production estimator, not a copy.** Each fold calls
  `sl_split_fit()` and stage 2 is `stage_2_sl()` (`R/cate_models.R`), so the
  study answers directly whether a proxy would pick the estimator the other
  studies report. Its libraries are sized to the rows it is fit on, so
  `SL.earth` enters the outcome library only from 500 rows. Dropped learners
  are saved as `SL1$sl_dropped`.
- **`net2` is the interaction lasso** (the expansion of `SL.glmnet.int`, at
  `lambda.min`). It is the only linear-family candidate that can represent
  scenario 3's X4·X5, and it keeps the `net*` candidates apart: lasso, ridge
  and elastic net on 10 covariates give near-identical τ̂, and near-ties make
  top-1 selection a coin toss.
- **`SL3` is a weak baseline.** `SL.nnet`'s defaults (2 hidden units, no
  weight decay) are unstable on a noisy pseudo-outcome. It gives the set a
  candidate a good proxy should avoid, which is what regret measures.

### Proxies (the selection methods)

A proxy is one combination of four factors. It selects the candidate with the
lowest score.

**1. Scoring-nuisance pipeline.** Each pipeline estimates `mu_DR` per arm, as
a T-learner (`mu0` on the W = 0 rows and `mu1` on the W = 1 rows, each a model
on X, matching the candidates), and `pi` (a model on X). A per-arm fit with
fewer than `MIN_ARM_ROWS` = 10 training rows in either arm is refused.

- **XGBoost** (`xgb`): a 36-point grid over `eta`, `max_depth`, `subsample`
  and `colsample_bytree`, tuned by 5-fold CV with up to 100 rounds and early
  stopping;
- **H2O AutoML** (`automl`): up to 20 models, with DeepLearning and XGBoost
  excluded, and no `max_runtime_secs` cap, since a cap would change what is
  measured.

From these, `calculate_pseudos()` forms the AIPW pseudo-outcome
`phi = mu1 − mu0 + (Y − mu_W)(W − pi) / (pi(1 − pi))`, where `mu_W` is
the observed arm's model, `mu1` if W = 1 and `mu0` if W = 0, and `pi` is the
estimated propensity trimmed to [0.05, 0.95] (`trim_ps()`, as for the
candidates). `phi05` is the same quantity with `pi` fixed at 0.5, which needs
no trimming. The IF score uses the same trimmed `pi`.

**2. Nuisance arm.** The arm sets which rows the scoring nuisance is fit on,
relative to the rows each candidate was trained on. That is the only axis: the
candidate fits are identical across arms.

| arm | tree | nuisance fit on | predicted on | row-honest | decoupled from candidate |
|---|---|---|---|---|---|
| `whole` | main, strategies | all n rows | all n rows | no | no |
| `cv_shared` | strategies | the candidate's V − 1 training folds | its held-out fold | yes | no |
| `holdout` | strategies | the candidate's held-out fold only | that fold | no | yes |
| `split` | split | the 20% evaluation rows | the same 20% | no | yes, and the candidates never see those rows |

- **Why there is no `cv_indep` arm.** A nuisance fit over a second,
  independent fold draw buys nothing on row-level honesty (whether row i's own
  `(Y_i, W_i)` entered the model predicting row i), which holds under the
  shared draw too. It only cuts the overlap of the two training sets for row i
  from 100% to 90% (`crossfitting/README.md`: "Re-randomising the stage-2
  split cannot remove that dependence"), at the cost of a full nuisance
  pipeline per replicate.
- **What `holdout` gives up.** It fits and predicts on the same block, so each
  row's own Y is in the model predicting it, `mu_DR` interpolates it, and
  phi's correction `(Y − mu_DR)` collapses toward zero. That is the price of
  being the only fold-based arm decoupled from the candidate: one small block
  cannot be both row-honest and decoupled. `cv_shared` is the row-honest
  comparator.
- At n = 250 the 25-row folds are too small to fit on, so `holdout` pools
  adjacent pairs of the candidate's folds into 5 blocks of 50
  (`holdout_blocks()`). A pooled block is only half decoupled from any one
  candidate fold model, since the partner fold was in its training set.
  `cv_shared` stays at V = 10.
- `holdout`'s per-arm fits are thin: about 25 rows per treatment arm at
  n = 250 and 500, 50 at n = 1000. Expect degenerate AutoML fits there more
  often than in the other arms.
- **The 80:20 arm** is the only one that refits the candidates, and that is
  the point: every other arm scores crossfit fits that have seen all n rows,
  wherever the nuisance is fit. It is its own tree, not another column,
  because τ̂, the evaluation rows and the true PEHE all live on the 20%.

**3. Score family.** Each family is a loss (lower is better), computed on the
evaluation rows:

| family | score | definition |
|---|---|---|
| influence | `inflt` | the influence-corrected PEHE estimate (`calc_infl_score`), with the T-learner plug-in `mu1_DR − mu0_DR` as its surrogate; the reported IF proxy |
| influence | `infl` | the same, with `phi` as its surrogate |
| DR risk | `dr` | `mean((τ̂ − phi)²)` (`calc_dr_risk`) |
| calibration | `calqK` | `Σ_k \|G_k\| · \|mean_{G_k}(τ̂) − mean_{G_k}(phi)\|`, over K rank groups of τ̂, for K = 5 and 10 (`calc_cal_score`) |

The calibration score has three choices that are easy to "correct" into
something else:

- **Absolute, not signed.** The signed sum telescopes to
  `n · (mean(τ̂) − mean(phi))`, an ATE-bias measure that cannot see
  miscalibration.
- **Weights are group counts, not proportions,** so the score scales with n.
  Within a run that factor is common to the 9 candidates and moves no ranking,
  but raw scores are not comparable across n.
- **Groups come from ranks** (`ceiling(K · rank(τ̂) / n)`), not `quantile()`
  cut points, which collapse or give NA when τ̂ is constant or tied, as in
  scenario 1 and a fully shrunk `net*`.

Both K are carried: at n = 250, K = 5 puts 50 rows in a group and K = 10 puts
25, trading a stable group mean against finer resolution. Calibration measures
calibration, not discrimination: a candidate predicting a constant scores
close to zero whatever its PEHE, so read it as a complement to the PEHE
proxies.

**4. Propensity.** The proxy uses either the pipeline's estimated `pi` or the
true 0.5, through `phi05` and `pi = 0.5`, which gives the `*05` columns.
Treatment is a fair coin in every scenario, so the `*05` scores are oracle-π
scores, and `dr` against `dr05` isolates what it costs to estimate a
propensity that never needed estimating.

Trimming the estimated `pi` matters: untrimmed, AutoML's propensity is
sometimes exactly 0 or 1 where a model predicts rows it was fit on (`holdout`,
`whole`), which makes phi non-finite and every estimated-π score undefined
for the run (`docs/DEVLOG.md`, bug S). Trimmed `holdout` propensities still
pile up at the bounds, where the weight is about 21 against 4 at π = 0.5, so
those proxies are defined but noisy. `me_metrics.R` recomputes phi from each
arm's saved predictions and keeps only finite scores.

This makes 10 score types per (arm × pipeline): 20 score columns per
candidate-run in the main and split trees, and 60 in the strategies tree,
named `<score_type>_<arm>_<pipeline>`. `"05"` appears in a score type only
when the propensity is fixed. Each tree also carries `true_pehe` and
`mu1_mse`, the candidate's MSE against the true `p1`: with `pi = 0.5` the IF
score reduces in expectation to that plus a constant shared by the candidates,
so `mu1_mse` is what it selects on (NA in the split tree, whose truth keeps
only τ).

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
