# ADEMP — interim-analysis validation

Two studies on one design, differing only in the outcome:

| study | folder | prefix | outcome | scenario set |
|---|---|---|---|---|
| continuous | `continuous/` | `cts_val_*` | continuous | `continuous_corr_0.5` |
| binary | `binary/` | `bin_val_*` | binary, risk-difference scale | `binary_corr_0.5` |

Both use the shared code in `val_common.R`. The binary arm differs from the
continuous one in its DGM and in the DR SuperLearner's outcome family, and
nowhere else.

## Aims

A CATE estimator fit at an interim analysis of a two-arm RCT produces
findings: responder subgroups, a spread of estimated effects, a ranking of
covariates by importance, and a top-ranked covariate. The study asks whether
those findings replicate in the participants recruited after the interim.

It is a robustness check on CATE-based subgroup discovery, not an estimator
comparison. The object of interest is whether one estimator's findings on
chunk 1 hold up on chunk 2, not which estimator is most accurate.

Secondary aims:

- how replication changes with the interim's share of the trial
  (`interim_prop`), which trades chunk-1 size against chunk-2 size;
- whether the two importance measures (TE-VIM and surrogate TreeSHAP) rank
  covariates as stably, and nominate the same top covariate;
- which form of the top-covariate test can confirm a nominated modifier,
  including the non-monotone modifier of scenario 4, and which forms are
  fooled by a correlated proxy;
- how much classical and HC3 standard errors change the interaction tests;
- how much less power the binary arm has than the continuous one.

## Data-generating mechanisms

The correlated sample-size studies' DGM at ρ = 0.5 (`R/dgm_scenarios.R`,
CORRELATED COVARIATES; see `sample_size/ADEMP.md` for the copula
and the calibration of bW).

- `E[Y | x, W] = m0(x) + W·τ(x)`, `τ(x) = bW + g(x)`, `W ~ Bernoulli(0.5)`
  independent of X.
- **Continuous:** `Y = 0.4 − 0.5·X1 + X2 + W·τ(x) + ε`, `ε ~ N(0, 0.5²)`.
- **Binary:** `Y ~ Bernoulli(m0(x) + W·τ(x))`, with
  `m0(x) = 0.34 + 0.36·plogis(−1.72 − 0.5·X1 + X2)`. g is the continuous g
  with its sign reversed, scaled by `RD_SCALE`, and X4, X5 enter through tanh
  (t4, t5).
- 10 covariates. X1–X5 and X01–X03 come from a Gaussian copula with
  exchangeable latent correlation 0.5; X1 ~ Bernoulli(0.4) and
  X3 ~ Bernoulli(0.7) are thresholded, the rest are N(0, 1). X04/X05 are a
  3-level factor, independent of everything else. Every estimator sees all
  ten.
- bW is set so the true ATE gives an unadjusted test 80% power at n = 1000,
  the whole trial, planned as if the effect were homogeneous. Both chunks share
  that ATE.

| scenario | description | modifiers | continuous g(x) | binary g(x) |
|---|---|---|---|---|
| 1 | null | none | 0 | 0 |
| 2 | simple | X4 | −X4 | 0.082·t4 |
| 3 | complex | X3, X4, X5 | 2·X3 + 0.5·X4 − 0.5·X4·X5 | −0.102·X3 − 0.0255·t4 + 0.0255·t4·t5 |
| 4 | non-linear | X4 | cos(X4) | −0.204·cos(X4) |

Features of the scenarios that the interpretation depends on:

- **Scenario 1 is the null.** τ(x) = bW for everyone, so nothing chunk 1
  finds is real. Every significant chunk-2 test there is a type I error, and
  the replication rates are the tests' size.
- **Proxies.** Under ρ = 0.5, X1, X2 and X01–X03 are correlated with the
  modifiers, so E[τ | X] moves with them even though they modify nothing.
- **Scenario 3.** X3 is the largest single source of heterogeneity. X5 acts
  only through its product with X4, so its marginal linear interaction is weak.
- **Scenario 4.** cos is even and X4 is centred on 0, so τ has zero linear
  slope in X4 and equal means either side of its median.
- **Binary heterogeneity is small.** The SD of τ at n = 1000 is about 0.05,
  0.06 and 0.09 in scenarios 2, 3 and 4.

### Design

Full factorial within each arm:

| factor | levels |
|---|---|
| scenario | 1, 2, 3, 4 |
| ρ | 0.5 |
| n | 1000 |
| `interim_prop` | 0.25 to 0.75 by 0.05 (11 points) |
| runs | 100; in scenario 1, 1000 at `interim_prop` 0.25, 0.5 and 0.75 |

That is 4 × 11 × 100 = 4400 array jobs per arm, plus 3 × 900 = 2700 for the
null's extra runs: 7100 in all (`<prefix>_val_config.R`). Scenario 1 was
added on 2026-10-12, and its extension to 1000 runs on 2026-10-25. The
extension is appended to the grid, as rows 4401–7100, so the first 4400 array
indices are unchanged. It covers only the three landmark points the report's
tables use, because the null is the size reference for every replication rate.

- Each run draws **one** trial of 1000 (`split_trial()`). Chunk 1 is its first
  `round(1000 · interim_prop)` rows and chunk 2 the rest. Rows are iid, so
  the first n1 are a valid "enrolled by the interim" cohort.
- The seed is the run alone (`setup_rng_stream(run)`). A run's chunk 1 at
  0.25 is therefore nested in its chunk 1 at 0.30, so every curve over
  `interim_prop` is paired within run. The same run draws the same
  participants in every scenario; only τ changes.
- Chunks range from 250 to 750 rows.

## Estimands

There is no single parameter being estimated. The study's targets are the
truths that each chunk-1 finding is checked against:

- **The CATE** τ(xᵢ), and its variance `Var(τ(X))`. The reports compute
  `Var(τ)` from a draw of 2 × 10⁵ from the DGM. It is the same in both chunks
  and at every `interim_prop`.
- **The set of true modifiers** M: ∅ in scenario 1, {X4} in scenarios 2 and
  4, {X3, X4, X5} in scenario 3. A working estimator should rank a member of M
  most important in both chunks.
- **Chunk-2 interaction contrasts,** conditional on what chunk 1 found:
  - for a subgroup S learned on chunk 1, the difference in the average
    treatment effect between S and its complement;
  - for chunk 1's top covariate X_top, the `W:v` coefficient of each test's
    working model in the chunk-2 population, where v is X_top in one of five
    forms (see Methods).

  W is randomised, so each `W:v` coefficient is the slope of the linear
  projection of τ on v, and on all covariates for the adjusted form. These
  contrasts depend on the run, through chunk 1. Since 2026-10-12 the study
  computes each one per run from chunk 2's true τ (`true2`), and the same
  contrast from chunk 1's estimated τ (`ref1`, what chunk 1 claimed). All
  `true2` are 0 in scenario 1. A significant test is read as a replication.
  Whether that is power or a type I error depends on the scenario and on what
  chunk 1 nominated: the reports separate the two by whether X_top is a member
  of M, and by the sign of `true2`.

## Methods

### Estimators

Each is fit separately on chunk 1 and on chunk 2 (`fit_val_methods()`, from
`R/cate_models.R`):

| estimator | nuisances | CATE |
|---|---|---|
| `causal_forest` | grf's own internal cross-fitting | grf causal forest, OOB predictions |
| `dr_random_forest` | `nuisance_rf()`: per-arm regression forests and a propensity forest, whole-sample OOB, propensity trimmed to [0.05, 0.95] | regression forest of the AIPW pseudo-outcome, OOB predictions |
| `dr_superlearner` | `nuisance_sl()`: per-arm SuperLearners and a propensity SuperLearner (mean + glm), single crossfit | SuperLearner of the pseudo-outcome, crossfit on the same folds |

- The DR SuperLearner uses 5 folds for a chunk under 500 rows and 10 from 500
  (`chunk_folds()`), and `sl_libraries()` sized to the chunk.
- The binary arm fits its outcome SuperLearners with `binomial()`
  (`method.NNloglik`). The forests take the 0/1 outcome as it is.
- No oracle or semi-oracle arms are run.

### Variable importance

Two measures per estimator per chunk, both larger-is-more-important:

- **TE-VIM:** refit the stage 2 with each covariate dropped and take the
  increase in squared error against the pseudo-outcome (`te_vim_scores()`).
  The refit is OOB for the forests and crossfit on the same folds for the DR
  SuperLearner.
- **Surrogate TreeSHAP:** a tuned xgboost surrogate of the estimated CATEs
  (`cvboost_cate()`), then mean |SHAP| per covariate (`get_shap_vims()`).

### Chunk comparisons

`chunk_validations()` makes four comparisons per estimator:

| comparison | what is done |
|---|---|
| `subgroups` | cut chunk 1's estimated CATEs into its top 10%, middle and bottom 10%; fit a classification tree (rpart) of those groups on X; predict the groups into chunk 2; test `W × 1[top10]` and `W × 1[bottom10]` there |
| `variances` | `var(τ̂)` in chunk 1 (`vt1`) and chunk 2 (`vt2`) |
| `var_imps` | each measure's covariate ranks in each chunk, and the change |
| `top_var_tests` | take each measure's chunk-1 top covariate X_top and test its interaction with W in chunk 2, five ways (below) |

Forms of the top-covariate test, each `lm` in chunk 2:

| form | v in `Y ~ W * v` | shape learned on chunk 1? |
|---|---|---|
| `p_cts` | X_top itself | no |
| `p_cts_adj` | X_top, in `Y ~ W * (all covariates)` | no |
| `p_split` | X_top split at its chunk-2 median; a binary X_top is its own split | no |
| `p_dd` | high/low side of a depth-2 tree of chunk 1's τ̂ on X_top (`data_driven_split()`) | yes |
| `p_proj` | a 4-df smoothing spline of chunk 1's τ̂ on X_top, evaluated in chunk 2 (`projection_score()`) | yes |

- Only chunk 1 decides the shape of `p_dd` and `p_proj`, so both are honest
  tests in chunk 2.
- `p_cts_adj` conditions on the other covariates, so it should not flag a
  proxy of a modifier that `p_cts` does flag.
- In scenario 4, `p_cts`, `p_cts_adj` and `p_split` have no power against X4
  even when it is nominated. They are kept as the contrast: `p_dd` and `p_proj` are
  the forms that can see it.

**Standard errors.** Every interaction test is run twice
(`chunk_validations_se()`): classical, and HC3 (`sandwich::vcovHC`). HC3 is
primary. For a binary Y the model is a linear probability model, so the errors
are heteroskedastic by construction. For a continuous Y the working model
leaves out the prognostic covariates and the spread of τ within each level of
v, so the residual variance still differs across the `W × v` cells.

**Estimates and reference contrasts.** Each test saves its `W:v` estimate
(`est_<test>`) and standard error (`se_<test>`) next to its p-value, and the
two reference contrasts of the Estimands, `ref1_<test>` and `true2_<test>`
(`tau_slope()`, `tau_slope_adj()`). For `ref1`, the test's v is built in
chunk 1 with the same helper: the subgroup tree's own groups in chunk 1, chunk
1's median split, and the tree and spline evaluated at chunk 1's X_top. By
construction `ref1` is positive for the top-10% subgroup and the data-driven
split, and negative for the bottom 10%.

**Not compared.** A chunk-1 vs chunk-2 comparison of the HTE tests (BLP,
independence tests) is not implemented, although every estimator carries the
inputs.

## Performance measures

Computed per (arm, scenario, estimator, `interim_prop`) over all the runs the
cell has, HC3 unless stated (`val_results.qmd`, one render per outcome and
scenario). That is 100, or 1000 at the null's landmark points. A proportion's
Monte Carlo SE is `sqrt(p(1 − p) / runs)`; a mean's is `sd / sqrt(runs)`. At
100 runs that is at most 0.05 for a proportion, and about 0.022 at p = 0.05.
At 1000 runs it is about 0.007 at p = 0.05.

- **Completeness:** runs present per cell, against the number the grid
  expects.
- **Subgroup replication:** the proportion of runs with the chunk-2 interaction
  p < 0.05 and p < 0.10, for the top-10% and bottom-10% subgroups separately,
  with MCSE; and the empirical CDF of those p-values. Runs whose tree predicts
  no member of a group into chunk 2 are NA and leave the denominator.
- **CATE variance:** the distributions of `vt1` and `vt2` against the true
  `Var(τ)`; per-run `vt1` vs `vt2`; and the mean of `vt2 − vt1`, with MCSE.
- **Importance-rank stability,** per measure:
  - the mean Spearman correlation between the chunk-1 and chunk-2 ranks, with
    MCSE;
  - the mean rank of each covariate in each chunk;
  - the mean absolute rank change of each covariate;
  - the proportion of runs in which a true modifier is ranked top in chunk 1,
    in chunk 2, and in both.
- **Top-covariate replication,** per measure and test form: the proportion of
  runs with chunk-2 p < 0.05, with MCSE, against the nominal 0.05. It is read
  alongside:
  - the proportion of runs in which X_top is a true modifier;
  - the proportion in which both chunks nominate the same covariate, right or
    wrong;
  - for `p_dd`, the share of runs with no split, one cut, and two or more cuts,
    and the mean share of chunk 2 on the high side.
- **Direction:** per test and form, the proportion of runs with p < 0.05 and
  `est` of the same sign as `ref1` (replicated chunk 1's direction), and as
  `true2` (right in truth; not defined in scenario 1). These sit alongside the
  plain p < 0.05 proportion, so the gap is the share of "replications" that
  went the wrong way. Also the means of `ref1`, `est` and `true2`: `ref1` above
  `true2` is chunk 1 overstating its finding.
- **Classical vs HC3:** each test's rejection proportion under both kinds of
  standard error, and the difference.
- **Scenario 1:** every rejection proportion above is a type I error rate,
  read against the nominal 0.05.

A headline table collects these at three landmark interim points: 0.25, 0.50
and 0.75.
