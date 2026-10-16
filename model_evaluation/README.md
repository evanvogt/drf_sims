# Model evaluation — does cheap proxy scoring pick the right CATE learner?

Given a basket of candidate CATE-learner configurations, can cheap proxy loss
functions — built from *independent* nuisance estimators — rank them the way
the (normally unobservable) true PEHE would? This is a model-*selection*
study, not another "does estimator X recover the true CATE" study: the
object of interest is the 12-column-wide proxy-vs-truth comparison
`me_metrics.R` produces, not any one candidate's own accuracy.

Ported from a single-commit, never-successfully-run prototype that used the
external `benchtm` package for data generation. This version uses the same
shared DGM every other study in this repo does
(`R/dgm_scenarios.R`, `set = "continuous_corr_0.5"` since 2026-10-08 - see
"Data" below) and the same 7-role file shape
(`config`/`dgms`/`models`/`analysis`/`check`/`collect`/`metrics`) `continuous/`
uses, plus study-specific extras the same way `crossfitting/` has extras
beyond that floor.

## Design

Four scenarios varying the structure of the CATE, crossed with three sample
sizes, 100 runs each — **1200 array jobs** (30 before 2026-10-08).

| | |
|---|---|
| scenarios | 1–4, the main four (see `R/dgm_scenarios.R`, `DESC_10`); 1, 4, 6, 8 before 2026-10-11 |
| covariates | correlated, latent ρ = 0.5 (`me_dgms.R`'s `ME_RHO`) |
| n | 250, 500, 1000 |
| runs | 100 |
| folds | 10 (all n — see "Crossfitting strategy" below) |
| results | `../results/correlated/model_evaluation/scenario_<k>/<n>/res_sim_<run>.RDS` |

### Data

Since 2026-10-08 the study runs on `continuous_corr_0.5`: `continuous/`'s
scenario table, with the covariates drawn from the copula at latent
correlation 0.5 — the primary arm of `sample_size/`. τ(x) and m0(x)
are unchanged; only the covariates' joint distribution moves, and with it bW,
the ATE and SD(τ) wherever g multiplies two modifiers (scenario 3 here,
X4·X5), and the correlation between the prognostic m0 and τ. See
`sample_size/README.md` for the numbers. Only ρ = 0.5 is run, so ρ
is a constant rather than a grid column and the grid is unchanged.

**Why only 4 of the 10 scenarios, not all of them like `continuous/`.** This
study's research question — does a cheap proxy loss rank 9 candidate models
the way true PEHE would — doesn't need every CATE-structure re-litigated to
answer. The four it runs are the main four every other study reports: 1 null,
2 simple (−X4), 3 complex (2·X3 + 0.5·X4 − 0.5·X4·X5), 4 non-linear (cos X4).

Until 2026-10-11 it ran 1, 4, 6, 8 instead: `crossfitting/`'s subset, chosen
before the main four existed and carried through the 2026-09-26 renumbering
number for number. The grid lists the scenarios as `c(1, 4, 2, 3)`, not `1:4`,
so that scenarios 1 and 4 keep their array indices and 2 and 3 take 6 and 8's
(`me_config.R`'s header). That matters because `me_strategies.R` and
`me_split.R` seed on the array index: the runs already made for 1 and 4 stay
valid and reproducible in all three trees, and only 2 and 3 need running (see
"Status").

**100 runs.** Each replicate is far more expensive than a `continuous/` one -
9 single-crossfit candidate-model fits, *plus* two independent
nuisance-evaluation pipelines (XGBoost with a 36-combination CV grid search,
and H2O AutoML with up to 20 auto-tuned models), each across the nuisance arms
and 2 nuisance targets - so the study used to stop at 30. It runs 100 since
2026-10-08, with 10 cores per task and 10 tasks at once (see "Sizing the array
job" below); the real per-replicate cost has still never been measured.

## Crossfitting strategy

All 9 candidates use a single leave-one-fold-out crossfit ("`scf_scf`" in
`crossfitting/`'s naming): nuisances trained on `V-1` folds and predicted on
the held-out fold, with the stage-2 regression trained on the *same* folds —
not the double-crossfit-over-fold-pairs scheme (`C(V,2)` nuisance fits at
`V=10`) this study originally ported with. Brought in line with
`crossfitting/`'s comparison of alternatives, the same comparison that moved
`R/cate_models.R`'s production estimators off double crossfitting (forests to
whole-sample OOB, causal forest to `grf`'s internal cross-fitting,
SuperLearner to `scf_scf`) — see `crossfitting/README.md` and
`R/cate_models.R`'s header.

Single crossfitting applies to all 9 candidates uniformly, not a per-learner
OOB/crossfit split, even though `rf1-3` (`ranger`) could take whole-sample OOB
predictions natively. `net1-3` and `SL1-3` have no OOB analogue — `glmnet` has
no bagging, and `SuperLearner`'s internal CV selects ensemble weights rather
than producing an honest prediction of a training row
(`crossfitting/README.md`'s "DR-learner, SuperLearner" section). Since this
study's question is whether cheap proxy losses rank the 9 candidates the way
true PEHE would, they need one shared honesty regime — mixing OOB forests with
crossfit everything-else would confound learner choice with crossfitting
scheme in the ranking. All `n` now use `V=10`: single crossfitting trains each
fold's stage-2 model on `(V-1)/V` of the data (vs. double crossfitting's
`(V-2)/V`), so unlike `continuous/` there's no need to shrink `V` at `n=250`.

The candidates' crossfitting is fixed at `V=10` single crossfit everywhere and
is **not** what this study varies. The one exception is the 80:20 arm, where
refitting them on a training split is the point (see below).

`me_analysis.R` draws one fold assignment per replicate, for the 9
candidates. The fold-based scoring arms reuse it. A second, *independent* draw
for the scoring pipelines (the `cv_indep` arm) was removed before the
2026-10-08 rerun — see the next section.

## Candidate models

9 configurations, each fit as its own single-crossfit DR-learner
(`me_models.R`, `run_all_candidate_models()`):

| id | family | outcome model and stage 2 | propensity |
|---|---|---|---|
| `rf1` | random forest | ranger defaults (`mtry = floor(sqrt(p))`) | the candidate's own forest |
| `rf2` | random forest | `mtry = p, max.depth = 5` (every covariate each split) | the candidate's own forest |
| `rf3` | random forest | `mtry = ceiling(p/2), max.depth = 3` | the candidate's own forest |
| `net1` | elastic net | `alpha = 1` (lasso) | the candidate's own glmnet |
| `net2` | elastic net | `alpha = 1` over main effects **and all pairwise interactions** (interaction lasso) | the candidate's own glmnet |
| `net3` | elastic net | `alpha = 0.5` | the candidate's own glmnet |
| `SL1` | SuperLearner | **the production `dr_superlearner`**: `sl_libraries(n)$Y` per arm, `$tau` at stage 2, pretested (`R/sl_library.R`) | `sl_libraries(n)$W` = mean + glm |
| `SL2` | SuperLearner | `glmnet, xgboost, cforest, earth, gam, mean` | mean + glm, as SL1 |
| `SL3` | SuperLearner | `svm, nnet, mean` — deliberately weak | mean + glm, as SL1 |

`p = ncol(X)` = 10 in every scenario (X1–X5, X01–X05). `rf2`/`rf3`'s `mtry` is
scaled to `p` rather than fixed, unlike the ported prototype's original
`mtry = 30`/`mtry = 10` — those were calibrated to `benchtm`'s much wider
covariate set and error out of `ranger` ("mtry can not be larger than number
of variables in data") against this DGM's smaller one. Found by
`me_testing.R full`, not by inspection.

**Outcome models are per arm (T-learner) in all 9** (since 2026-10-08, before
the correlated rerun), as in `R/cate_models.R`'s DR-learners: one model per
arm, fit on `X` with that arm's training rows, predicting every held-out row.
Until then each candidate fit one model on `cbind(W, X)` and read it at
`W = 0` and `1`. For the `net*` candidates that made `Y1.hat − Y0.hat` a single
constant (W entered as a main effect only), so all the heterogeneity in the
pseudo-outcome came from its noisy residual term; for `rf1`, `W` was offered at
only ~30% of splits.

**Propensity.** `rf*` and `net*` fit it with the candidate's own learner, as
they always have. The SuperLearner candidates all use the production
propensity library, mean + glm: every study here is an RCT at P(W = 1) = 0.5,
and a flexible learner only fits noise that the DR weight `1/(e(1 − e))`
amplifies (`R/sl_library.R`'s header). Every candidate trims to
[0.05, 0.95].

**`SL1` is the production estimator, not a copy of it.** Each fold calls
`R/cate_models.R::sl_split_fit()` and stage 2 is `stage_2_sl()`, with
`sl_libraries(n)` on this study's `fold_indices` — the same `scf_scf` scheme
`nuisance_sl()` uses — so the study answers directly whether a proxy would
pick the estimator the other studies report. Its dropped learners are saved
as `SL1$sl_dropped`, as `cate_methods()` saves them. The libraries are sized
to the rows the candidate is fit on (`n`, or the 80% in the split arm), so
`SL.earth` is in the outcome library only from 500 rows.

**`net2` is the interaction lasso** (replaced ridge on 2026-10-08): `X`
expanded to its 10 main effects and 45 pairwise products (`net_design()`, the
same expansion as `R/sl_library.R`'s `SL.glmnet.int`), at `lambda.min`, in all
three nuisances and stage 2. It is the only linear-family candidate that can
represent scenario 3's X4·X5 effect, and it keeps the three `net*` candidates
apart: lasso, ridge and elastic net at `lambda.min` on 10 covariates are
expected to give near-identical τ, and near-ties make top-1 selection a coin
toss.

**`SL3` is deliberately weak.** `SL.nnet`'s defaults (2 hidden units, no
weight decay) are unstable on a noisy pseudo-outcome. It is kept so the set
has a candidate a good proxy should avoid, which is what regret measures.

`rf*`, `net*`, `SL2` and `SL3` are this study's own code path (e.g. the
elastic-net arms wrap a single learner via `create.Learner("SL.glmnet", ...)`
inside `SuperLearner()`, a code path `R/cate_models.R` doesn't have) — not
unified into `R/`, since only this study needs it.

## Nuisance-evaluation pipelines

A second, *independent* estimate of 2 nuisance targets (`mu_DR`, `pi`) —
used only to build proxy scores for the 9 candidates above, never to fit
them:

- **XGBoost** — a hand-tuned CV grid (`eta`, `max_depth`, `subsample`,
  `colsample_bytree`).
- **H2O AutoML** — up to 20 auto-tuned models per target
  (`exclude_algos = c("DeepLearning", "XGBoost")`).

**Outcome models are per arm (T-learner)** since 2026-10-14, in both
pipelines and every arm: `mu0_DR` is fit on the training set's `W = 0` rows
and `mu1_DR` on its `W = 1` rows, each with the pipeline's own tuning (the
same XGBoost grid, or its own AutoML run). Until then both pipelines fit one
outcome model on `(W, X)` and read it at `W = 0` and `1` — the S-learner the
candidates had dropped on 2026-10-08 (see "Candidate models"); the scoring
pipelines were missed by that change (`docs/DEVLOG.md`, R). The change is for
**consistency**, so the candidates and the nuisances that score them share
one outcome-model structure. It is not a fix for a demonstrated accuracy
problem: a local XGBoost check (n = 500, scenarios 2–4, two seeds each) found
`cor(mu1_DR − mu0_DR, tau)` about the same under both (0.73–0.96), and a
little lower per arm in scenario 4, where each arm's model has to learn the
full baseline from ~250 rows. AutoML and n = 1000 were not checked. Every
result tree written before 2026-10-14 is S-learner scored; see "Status". The propensity
model is unchanged. `me_nuisance.R` refuses a per-arm fit with fewer than
`MIN_ARM_ROWS` = 10 training rows in either treatment arm.

### The three arms

What varies is **what data the evaluation nuisance sees relative to the data
the candidate it is scoring was trained on**. That is the only axis; the
candidate fits are identical across all three arms, because `me_strategies.R`
reads them from the completed results rather than refitting them.

| arm | nuisance trained on | predicted on | row-honest | decoupled from candidate |
|---|---|---|---|---|
| `whole` | all `n` rows | all `n` rows | ✗ | ✗ |
| `cv_shared` | the candidate's own `V-1` training folds | the candidate's held-out fold | ✓ | ✗ (identical training set) |
| `holdout` | the candidate's held-out fold **only** | that same fold | ✗ | ✓ |

All three are fit by `me_strategies.R`. `me_analysis.R` still fits `whole`
for the main tree, but since 2026-10-14 the strategies pass recomputes it
rather than copying it, so the main tree's `nuisances` (and its
`me_metrics.RDS`) are superseded: the strategies tree is the scored one, and
the main tree supplies only the candidates, data, truth and `fold_info`.

**Why `cv_indep` was removed (2026-10-08).** A fourth arm fit the nuisance
leave-one-fold-out over a second, *independent* fold draw, justified as
stopping a candidate's `tau_hat` and the nuisance at the same row being fit on
the identical training set. Two things are wrong with that:

- Row-level honesty — whether row `i`'s own `(Y_i, W_i)` entered the model
  predicting at row `i` — holds under **both** a shared and an independent
  draw. Neither ever lets a nuisance see the row it predicts, so the second
  draw buys nothing on the axis it was justified by.
- What it actually changes is the *overlap* of the two training sets for row
  `i`, from 100% down to `(V-1)/V` = **90%**. It removes a tenth of the
  dependence.

`crossfitting/README.md` makes the general version of this point — *"Re-
randomising the stage-2 split cannot remove that dependence"*. The arm was
carried while it was already computed and free; the correlated rerun would
have paid for a full 10-fold nuisance pipeline per replicate to recompute it,
and the report already discarded it, so it was removed with its fold draw.

**What `holdout` gives up.** It is resubstitution — fit and predict on the
same block — so every row's own `Y` is in the model that predicts it, `mu_DR`
interpolates it, and `phi`'s AIPW correction `(Y - mu_DR)` collapses toward
zero. That is the arm's defining property, not a defect: it is the only
fold-based arm fully decoupled from the candidate, and row-level honesty and
decoupling cannot both be had from a single small block. `cv_shared` is the
row-honest comparator, and the pair makes the cost of each visible.

This restores the `infold` regime the ported prototype had and that was
removed. Its AutoML branch shipped with a `W`-recoding line
(`as.numeric(W) - 1`) flagged at the time as unverified; the flag was right,
and the line is deliberately **not** restored — `W` in that scope is the
caller's numeric 0/1 vector and was never converted to a factor, so the line
would send it to `-1/0` and break `calculate_pseudos()`. Only the fold filter
came back. See `run_automl_holdout()`.

**Block size at `n=250`.** `V=10` gives 25-row folds, and `mu_DR`/`pi` fit on
just the rows inside one — not estimable. So at
`n=250` only, `holdout_blocks()` pools *adjacent pairs of the candidate's
folds* into 5 blocks of 50. Pooling rather than drawing a fresh 5-fold split
is what keeps each block tied to the candidate's own partition without
refitting any candidate. The cost, stated plainly: a pooled block is only
*half*-decoupled from any single candidate fold model, because the partner
fold was in that model's training set. `cv_shared` stays at `V=10` at every
`n` — it trains on 225 rows at `n=250`, so nothing forces a change there.

**Per-arm fits are thin in `holdout`.** With per-arm outcome models each
block's `mu0_DR` and `mu1_DR` see only one treatment arm's rows: about 25 at
`n=250` (pooled blocks of 50) and `n=500` (unpooled folds of 50), 50 at
`n=1000`. The block rule was deliberately left as it is (2026-10-14) so the
arm keeps its definition; expect more degenerate AutoML fits, and more
`automl`/`holdout` NA, than the S-learner runs had.

### The 80:20 arm

A fifth position on the same axis, run as its own study (`me_split.R`,
`results/model_evaluation_split`) rather than a fifth column, because three
things move together: `tau_hat` exists only on the 20%, the evaluation rows
*are* that 20%, and `true_pehe` therefore has to be computed on those same
rows or the reference is not matched to the proxies. Reporting it as another
column would silently compare scores over different row sets against
differently-computed truths.

It is the **only** arm that refits the candidates, and that is the point:
every other arm scores `me_analysis.R`'s crossfit fits, whose training data
covers all `n` rows, so no matter where the nuisance is fit the candidate has
already seen those rows. Here the candidates see only the 80% (still
crossfit within it, still the same hyperparameters via
`candidate_hyperparams()`) and the nuisance sees only the 20%. `n=250` is
excluded — a 50-row evaluation set is too thin to rank 9 models.

### Scores

From these, `calculate_pseudos()` builds two AIPW pseudo-outcomes: `phi`,
using the pipeline's estimated `pi`, and `phi05`, with the propensity fixed
at 0.5. `me_metrics.R` scores every candidate's `tau`
against both pipelines under every arm, plus true PEHE (available because this
is simulated data).

**Three score families, each in two propensity regimes** — 10 score types per
(arm × pipeline):

| family | estimated `pi` | fixed `pi = 0.5` | what it is |
|---|---|---|---|
| influence, T-learner surrogate | `inflt` | `inflt05` | influence-corrected PEHE proxy (`calc_infl_score`) with the plug-in `mu1_DR − mu0_DR` as its surrogate. The IF proxy the thesis figures and tables report (`me_selection.R`) |
| influence, φ surrogate | `infl` | `infl05` | the same score with `phi` / `phi05` as its surrogate |
| DR risk | `dr` | `dr05` | MSE against the AIPW pseudo-outcome (`calc_dr_risk`) |
| calibration | `calq5`, `calq10` | `cal05q5`, `cal05q10` | DR calibration over K quantile groups (`calc_cal_score`) |

Column names stay `<score_type>_<arm>_<pipeline>`; the propensity regime and
the group count are folded into that first token (`cal05q10` = calibration,
fixed `pi`, K = 10) rather than added as a fourth field, because
`me_results.qmd` recovers the design axes by splitting the name. `"05"` never
appears in a score type unless the propensity is fixed — `calq10` carries a
`"10"`, never an `"05"` — so the regime is recoverable with `grepl("05", .)`.

#### Fixing the propensity at 0.5 is an oracle, not an approximation

`R/dgm_scenarios.R` assigns treatment with `W <- rbinom(n, 1, 0.5)` — a fair
coin, independent of `X`, in **every** scenario. So 0.5 *is* the true
propensity, the `*05` scores are oracle-π scores, and the `dr` / `dr05`
contrast isolates exactly one thing: what it costs to estimate a propensity
that never needed estimating. This is the RCT setting the study is about, so
the answer is not incidental.

**None of this required a rerun.** `phi05` has been computed by
`calculate_pseudos()` since this study's first commit, so every completed
replicate already carries it; `calc_infl_score()` already took `pi` as an
argument; and the calibration score needs only each candidate's saved per-row
`tau` and the arm's saved `phi`. Scores are a pure post-hoc function of
`<prefix>_all.RDS`, so adding a score family means re-running `me_metrics.R`
and nothing else — the same argument `me_strategies.R`'s header makes for its
own pass. This is separate from the propensity *inside the candidates*:
`me_models.R` builds each candidate's own pseudo-outcome from a trimmed
estimated `W.hat` (its own learner for `rf*`/`net*`, mean + glm for `SL*` —
see "Candidate models"), and changing that changes every stored `tau`. Every
arm here scores the identical candidate fits, which is the controlled
comparison the study rests on.

#### The calibration score

Split the evaluation rows into K quantile groups `G_1..G_K` of the
*candidate's own* `tau_hat`. Each group has the effect the candidate claims and
the effect the DR scores imply:

```
GATE_k^hat = (1/|G_k|) Σ_{i∈G_k} tau_hat(x_i)
GATE_k^DR  = (1/|G_k|) Σ_{i∈G_k} phi(x_i)
M^CAL-DR   = Σ_k |G_k| · | GATE_k^hat − GATE_k^DR |
```

Three implementation choices worth stating, because each is easy to "correct"
into something that measures a different thing:

- **Absolute, not signed.** The signed sum `Σ_k |G_k| (GATE_k^hat − GATE_k^DR)`
  telescopes to `n · (mean(tau_hat) − mean(phi))`: every group boundary cancels
  and what survives is an ATE-bias measure that cannot see miscalibration at
  all. Over- and under-estimation in different groups must not net out.
- **Weights are group counts, not proportions** — `sum(w)` is `n`, not 1, so
  the score scales with `n`. Within a run that factor is identical across the 9
  candidates, so it moves no ranking, correlation, pick or regret; but the raw
  magnitude is not comparable across sample sizes.
- **Groups come from ranks, not `quantile()` cut points.** Cut points collapse
  to fewer than K groups — or emit `NA` — the moment `tau_hat` is constant or
  heavily tied, which is exactly what scenario 1 and a fully-shrunk `net*`
  produce. `ceiling(K · rank(tau_hat) / n)` always gives K non-empty,
  near-equal groups, so the column never silently goes `NA` for the flattest
  candidates.

`CAL_QUANTILES` (`me_config.R`) carries **both** K = 5 and K = 10 rather than
picking one: at `n = 250` the first puts 50 rows in a group and the second 25,
trading a stable `GATE^DR` against finer resolution, and emitting both makes
the sensitivity to K a result instead of a hidden choice.

**What it cannot see.** Unlike `infl` and `dr` this is not an estimate of PEHE;
it measures calibration, not discrimination. A candidate predicting the
constant ATE everywhere has arbitrary quantile groups, so each group's
`GATE^DR` is about the overall mean and every discrepancy is near zero — a
near-perfect score at whatever PEHE the true heterogeneity implies. Expect the
shrunk-to-a-constant candidates to look good here, and read the family as a
complement to the PEHE proxies rather than a replacement.

#### Column counts

The column set is *derived* from whatever arms the nuisance list carries, never
enumerated, which is why one `me_per_model()` serves all three result trees:

| tree | arms | score columns |
|---|---|---|
| `model_evaluation` | `whole` | 2 + 10×1×2 = **22** |
| `model_evaluation_strategies` | `whole`, `cv_shared`, `holdout` | 2 + 10×3×2 = **62** |
| `model_evaluation_split` | `split` | 2 + 10×1×2 = **22** |

The 2 are `true_pehe` and `mu1_mse`, the candidate's MSE against the true
treated-arm mean `p1`. With `pi = 0.5` the IF score reduces, in expectation, to
`mean((tau_hat − mu1)²)` plus a constant shared by the candidates, whichever
surrogate it uses, so `mu1_mse` is what it selects on. It is `NA` in the split
tree, whose truth keeps only `tau`. `inflt`, `inflt05` and `mu1_mse` were added
on 2026-10-15; they need only `me_metrics.R` rerun, not the fits.

The split tree needs no scoring variant because `me_split.R` stores `data` and
`truth` already restricted to its 20% evaluation rows, so every vector
`me_per_model()` touches is the same length.

### Propensity trimming — [0.05, 0.95], every arm and pipeline (2026-10-14)

`calculate_pseudos()` trims the scoring propensity to [0.05, 0.95] with
`trim_ps()`, as the 9 candidates' propensities always were (`me_models.R`) and
as every other study's are. It applies to every arm and both pipelines alike,
which keeps them comparable; the earlier objection was only to trimming some
arms and not others. The true propensity is 0.5, so trimming can only move an
estimate that is already wrong toward the truth. What the estimated- vs
fixed-π contrast measures is now the cost of estimating a *trimmed*
propensity. `pi` keeps the raw prediction, `pi_trim` the trimmed one, and
`me_strategies.R`'s per-run `pi_diagnostics` (min/max/quantiles, max weight,
`n_extreme` = rows trimmed) are on the raw `pi`.

**Why (DEVLOG bug S).** Untrimmed, AutoML's propensity was **exactly 0 or 1**
on some rows of 294 of 1200 `holdout` runs, 58 `whole` and 6 `cv_shared`
(2026-10-14 cluster check; XGBoost never, and never an NA prediction). Most
likely from H2O's forests, whose vote-fraction probabilities reach 0 or 1
when a model predicts the rows it was fit on (`holdout`, `whole`). One such
row made phi NaN (π = 1 where W = 1) or ±Inf (π = 1 where W = 0) and every
estimated-π score undefined for the whole run, so some AutoML proxies'
summaries rested on as few as 65% of runs. Where phi was +Inf, the DR score was +Inf for all 9
candidates, and the argmin returned the first candidate as if it were a pick.
`me_selection.R` now keeps only finite scores, so that cannot recur silently.

The trimming is applied **when scoring**: `me_metrics.R`'s `arm_scores()`
recomputes phi from each arm's saved `mu0_DR`/`mu1_DR`/`pi`, so results fit
before 2026-10-14 need only collect and metrics re-run, no refit. It checks
the recomputed phi against the saved one on the rows the trimming did not
move, which catches Y/W that are not the nuisance's rows.
`me_strategies_verify.R` checks completeness the same way (no NA prediction,
finite trimmed phi) and tallies the runs that had a 0/1 propensity. Trimmed
`holdout` propensities still pile up at the bounds, where the weight
`1/(π(1−π))` is about 21 against 4 at π = 0.5, so those proxies are defined
in every run but noisy.

## Files

| file | role |
|---|---|
| `me_config.R` | the parameter grid, results path, and candidate-model list — **the** definition |
| `me_dgms.R` | names this study's slice of `R/dgm_scenarios.R` (`set = "continuous_corr_0.5"`) |
| `me_utils.R` | design-matrix prep and fold-splitting (DGM-agnostic, unchanged by the port) |
| `me_models.R` | the 9 candidate CATE-learner configurations and their fitting logic |
| `me_nuisance.R` | the two independent nuisance-evaluation pipelines — see below for why this exists outside the usual 7-file shape |
| `me_analysis.R` | array entry point; one row of the grid per index |
| `me_strategies.R` | second pass over completed runs — fits all three arms (`whole`, `cv_shared`, `holdout`) without refitting the candidates, writes to `model_evaluation_strategies` |
| `me_strategies_verify.R` | proves that pass carried the candidates, data, truth and `fold_info` through bit-identically and that all three arms are complete (no NA prediction, finite phi with π trimmed), and tallies the runs with a propensity of exactly 0 or 1 (see "Propensity trimming" above) |
| `me_split.R` | the 80:20 arm — the only script that refits the candidates |
| `me_check.R` | finds missing runs, writes `jobscripts/failed_ids.txt`, and updates `-J` and the resource request in the rerun jobscript. Takes a tree: `main` (default) / `strategies` / `split` |
| `me_collect.R` | gathers per-run files into `<prefix>_all.RDS`. Same tree argument |
| `me_metrics.R` | computes `<prefix>_metrics.RDS` (reuses `R/metrics.R::compute_metrics()`). Same tree argument |
| `me_results.qmd` | the results report — see below for why it derives its own quantities |
| `me_testing.R` | verification checks — run before submitting anything |

**Why `me_nuisance.R` exists outside the `config`/`dgms`/`models`/`analysis`/
`check`/`collect`/`metrics` shape**: no other study in this repo has a
second, independent nuisance-estimation pipeline used purely to *score*
candidate models rather than fit them — it doesn't map onto any of those 7
roles. `crossfitting/` is the precedent for a study needing files beyond that
floor (`cf_testing.R`, `cf_results.R`/`.qmd`); the 7-file shape is a floor,
not a ceiling.

**Why `me_results.qmd` derives its own quantities**: every other study's
report summarises a per-model metric straight out of its `*_metrics.RDS`
(mean bias, mean MSE, coverage). This one can't — the object of interest is
the *ranking* of the 9 candidates within a single (scenario, n, run), so the
report first has to reduce each run's 9x9 score matrix to per-run rank
agreement, top-1 selection accuracy and regret before anything can be
averaged. That derivation lives inline in the `.qmd`, as
the retired `sample_size/continuous/cts_results.qmd` and `sample_size/binary/bin_results.qmd` kept theirs; it is
not in `R/figures.R`, whose `summarise_metrics()` is built for the
bias/MSE/correlation columns this study doesn't have.

**Note on `me_metrics.R`**: `R/metrics.R::compute_metrics()` always does
`true_tau <- sim_res$truth$tau` — there's no equivalent of the ported
prototype's `calc_metrics(..., truth_avail = FALSE)` branch (for scoring
without ground truth). Since `me_analysis.R` only ever runs against
`dgm_scenarios.R`-generated data, which always carries truth, this is a real
but low-risk scope reduction from the original code.

## Running it

**Locally**, run only the fast check. `full` mode is the slow half and `SL2`
fails on the Windows machine for a documented package-version reason (see
"Known local-environment limitation"), so it is not a useful local signal —
run it on the cluster instead:

```bash
Rscript model_evaluation/me_testing.R               # structure + regression checks + arm plumbing
qsub    model_evaluation/jobscripts/me_testing.sh   # the full suite, with real XGB-CV/H2O AutoML
```

`me_testing.sh` also settles two things about the cluster environment that
nothing else does: whether `sim-env` really carries `h2o`/`xgboost`/`caret`
with a working Java runtime, and whether `SL2`'s local failure reproduces
there (it should not).

The main study:

```bash
qsub model_evaluation/jobscripts/me_1.sh            # the study itself - 1-1200
Rscript model_evaluation/me_check.R                 # writes failed_ids.txt if any are missing
qsub model_evaluation/jobscripts/me_collect.sh
qsub model_evaluation/jobscripts/me_metrics.sh
```

The nuisance-arm pass. **This does not require rerunning the study** — it
reads the completed results and adds arms to them (see `me_strategies.R`'s
header for why that is sound):

```bash
qsub    model_evaluation/jobscripts/me_strategies.sh      # 1-1200, reads the main tree
Rscript model_evaluation/me_check.R strategies            # progress / completion
qsub    model_evaluation/jobscripts/me_strategies_verify.sh
qsub -v TREE=strategies model_evaluation/jobscripts/me_collect.sh
qsub -v TREE=strategies model_evaluation/jobscripts/me_metrics.sh
```

The 80:20 arm, independent of the above:

```bash
qsub    model_evaluation/jobscripts/me_split.sh           # 1-800 (n = 500, 1000 only)
Rscript model_evaluation/me_check.R split
qsub -v TREE=split model_evaluation/jobscripts/me_collect.sh
qsub -v TREE=split model_evaluation/jobscripts/me_metrics.sh
```

```bash
quarto render model_evaluation/me_results.qmd   # the report - needs me_metrics.RDS
```

**Re-scoring an already-collected tree.** The score families in `me_metrics.R`
are a pure post-hoc function of `<prefix>_all.RDS` — adding one means
re-running only the last step, never the study:

```bash
qsub model_evaluation/jobscripts/me_metrics.sh   # rewrites me_metrics.RDS in place
```

`me_collect.R` output is untouched, so `me_all.RDS` does not need regenerating
either. Because this overwrites `me_metrics.RDS`, it is worth keeping a copy of
the old one and checking that the pre-existing columns come back bit-identical
and only new ones were appended — the same kind of inertness proof
`me_strategies_verify.R` gives for its pass.

**Checking progress of the derived trees.** `Rscript me_check.R strategies`
runs `check_failed(..., write = FALSE)` and reports how many runs are done.
Note the completion criterion is *not* zero missing: `study_strat`'s grid is
all 1200 rows (the array index has to keep meaning the same grid row), so any
runs excluded from the main study are reported missing forever. `me_check.R`
diffs the two trees and separates "missing because there is no source run"
from "missing because this pass failed" — only the latter needs resubmitting,
which is also why `failed_ids_strat.txt` is not written automatically.

Results land in `../results/correlated/model_evaluation{,_strategies,_split}/`
(outside the repo, as elsewhere). `me_results.qmd` reads `me_metrics.RDS` from
there, so it renders wherever the results are — not on a machine that only has
the repo. The rendered `.html` is gitignored, as every other study's report is.

## Sizing the array job

`me_1.sh`'s `#PBS -l` lines and trailing `Rscript` args (`workers`/`n_cores`)
are set by hand, never measured: the `syrup` profiling sweep meant to replace
them didn't work for this study (see `docs/DEVLOG.md`'s "Resource profiling
(removed)"). Since 2026-10-08 they are:

| jobscript | select | trailing args |
|---|---|---|
| `me_1.sh` | `ncpus=10:ompthreads=10:mem=24gb` | `10 10` |
| `me_rerun.sh` | `ncpus=11:ompthreads=11:mem=29gb` | `10 10` |
| `me_split.sh`, `me_split_rerun.sh` | `ncpus=10:ompthreads=10:mem=24gb` | `10 10` |
| `me_strategies.sh`, `me_strat_rerun.sh` | `ncpus=10:ompthreads=10:mem=16gb` | `10` |

Why 10: each candidate is crossfit over 10 folds in parallel, so 10 workers
run every fold loop in one round; 5 would take two rounds, and so would 8.
`ompthreads` must equal `ncpus`, or R does not see the extra cores.
`me_rerun.sh` is one core and 1.2x memory above `me_1.sh`, which is what
`check_failed()` writes there anyway. `me_strategies.sh` has no workers, and its fold loops run one after another,
so its 10 cores only feed XGBoost/H2O threads on small data - expect much less
speed-up there. Every walltime is still a placeholder. Keep in mind:

- **Two knobs, two sequential phases.** `workers` (the `future` multisession
  backend, controlling the 9 candidate models' single-crossfit fold-wise
  fitting) and `n_cores` (XGBoost's `nthread` / H2O's `nthreads`, controlling
  the nuisance-evaluation pipelines) parallelise two *sequential* phases of
  one replicate, not one combined computation. So `ncpus` needs to cover
  `max(workers, n_cores)`, not the sum, since the two are never both active at
  once.
- **H2O's JVM is a separate Java process.** Each task starts its own H2O JVM
  with a `mem = "10G"` heap, so `mem=` has to cover that on top of R - and in
  `me_1.sh`/`me_split.sh` on top of the 10 worker sessions too, which stay
  alive through the nuisance phase. Check the request against
  `qstat -fx <jobid> | grep resources_used` on the first real subjobs.

**The array's concurrency throttle (`-J 1-1200%N`) is a separate problem from
the per-task request.** Each concurrent task starts its own H2O JVM
cluster with a `mem = "10G"` heap, and too many at once makes the H2O calls
fail — nothing like `continuous/`'s `%190` or `crossfitting/`'s `%380` is safe
here; every array runs at `%10`. `me_check.R` passes `throttle = 10` to
`check_failed()`, so the `-J` it writes into `me_rerun.sh` stays at `%10` too
(its default is `%100`). With `%10` the main array takes about 120 times one
run's walltime (1200 / 10).

**Package availability.** `h2o`, `xgboost`, `caret`, `tidyverse` are
confirmed present in this machine's local ambient R library — `benchtm` is
confirmed absent, consistent with the old prototype never having run. That
does **not** confirm the *cluster's* `R/4.3.2-gfbf-2023a` module + `sim-env`
conda env has them: none of these four packages are dependencies of any
other study in this repo, and H2O additionally needs a working Java runtime
on the compute node. Verify against the actual cluster environment before
the first `qsub`.

## Known local-environment limitation

`SL2` (`SL.glmnet, SL.xgboost, SL.cforest, SL.earth, SL.gam, SL.mean`) fails
locally: the installed `xgboost` package has a redesigned API (`eta` renamed
to `learning_rate`, `data` renamed to `x`, an explicit `y` argument now
required) that `SuperLearner`'s bundled `SL.xgboost` wrapper predates.
`SuperLearner()`'s *fit* step handles that failure gracefully (gives
`SL.xgboost` weight 0, warns, continues) — but `predict.SuperLearner()`
unconditionally calls `predict()` on every library learner regardless of its
weight, and crashes on the `NULL` fit object that failure left behind:

```
Error in UseMethod("predict") :
  no applicable method for 'predict' applied to an object of class "NULL"
Calls: predict -> predict.SuperLearner -> do.call -> predict
```

Only the outcome and stage-2 fits use `SL2`'s library (its propensity is
mean + glm since 2026-10-08), but those still predict through
`predict.SuperLearner()`, so the failure stands. On 2026-10-08 a cluster
`me_testing.sh` run (on the candidates as they were before that day's changes)
was passing, so `SL.xgboost` is kept.

`SL2`'s library is unchanged from the ported prototype - this is a local
package-version mismatch, not something the DGM swap introduced, and exactly
the situation `.claude/CLAUDE.md` already documents ("local R is 4.5.3 ...
the cluster runs R 4.3.2-gfbf-2023a ... if package-version issues come up,
this is why"). Confirm whether this reproduces on the cluster's actual
`R/4.3.2-gfbf-2023a` module + `sim-env` conda env before relying on a clean
local `me_testing.R full` run as a sign this study is ready to submit.

## Status

### Scoring propensity trimmed: metrics owed, no refit (2026-10-14)

The scoring pseudo-outcome now uses the propensity trimmed to [0.05, 0.95]
(see "Propensity trimming"; `docs/DEVLOG.md`, bug S). `me_metrics.R`
recomputes phi from each run's saved predictions, so nothing is refitted:
after the strategies (and split) passes are done, re-run collect and metrics
for each tree, then the thesis tables and figures. Every `*_metrics.RDS`
written before this has undefined estimated-π scores in the runs where AutoML
predicted a propensity of exactly 0 or 1.

### Per-arm scoring nuisances: strategies and split trees owed in full (2026-10-14)

The scoring pipelines moved to per-arm outcome models (see "Nuisance-evaluation
pipelines"; `docs/DEVLOG.md`, bug R). Every nuisance on disk is S-learner
scored, so every score in every tree is stale. The main tree's candidates are
not: they were already per arm, and the strategies pass refits no candidate.
What reruns:

- **strategies**: all 1200 runs, not just scenarios 2/3 — `me_strategies.R`
  now recomputes `whole` as well as `cv_shared` and `holdout`;
- **split**: all 800 runs, candidates included — `me_split.R` fits its
  candidates and its nuisance in one pass;
- **main**: nothing beyond the scenario 2/3 runs owed below. Its `nuisances`
  are superseded by the strategies tree, and its `me_metrics.RDS` is not used.

From the repo root on the cluster, after `git pull`:

1. Finish the main tree first (step 1 of the section below) if its scenario
   2/3 runs are not done: both other trees read it.
2. `cd model_evaluation/jobscripts`, then `qsub me_testing.sh` (`full` mode,
   the only place the per-arm AutoML code runs before the array).
3. `qsub me_strategies.sh` — the full `1-1200%10`. The pass now makes about
   1.6× the model fits it did, in both pipelines (two outcome models per fit
   instead of one, plus the `whole` arm; the per-arm fits are on half the
   rows), so check the last strategies run's `resources_used`
   against the 1-hour walltime before submitting.
4. `qsub me_split.sh` — the full `1-800%10`. It does not wait for step 3.
5. `qsub me_strategies_verify.sh`. It fails only on an NA prediction or a
   phi that is non-finite even with π trimmed; propensities of exactly 0 or 1
   are tallied per pipeline and arm, not failed (see its header).
6. Check, collect and metrics for `strategies` and `split` as in step 4 below,
   then render `me_results.qmd` and `me_infl_investigation.qmd` (which reads
   `me_strat_metrics.RDS` directly, and needs `population.RDS` from
   `Rscript model_evaluation/me_infl_investigation.R`).

### Scenarios 2 and 3 owed (2026-10-11)

Superseded for the strategies and split trees by the section above; step 1
(main tree) still applies.

The grid moved from scenarios 1, 4, 6, 8 to the main four, 1–4 (see "Data").
Runs already made for 1 and 4 keep their array indices and stay valid in all
three trees, so only 2 and 3 need running. In each tree they are the indices
that used to be 6 and 8: 3, 4, 7, 8, 11, 12, …. That gives 600 main, 600
strategies and 400 split runs. Do these steps on the cluster, from the repo
root unless a step says otherwise. Every `qsub` is run from inside
`model_evaluation/jobscripts/`, because the jobscripts `cd "${PBS_O_WORKDIR}/.."`
and read their id file from `${PBS_O_WORKDIR}`.

If a tree has not been run at all on the current design, skip its step below
and submit its full array (`me_1.sh`, `me_strategies.sh`, `me_split.sh`). The
grid already covers 1–4.

0. **Set up.** `git pull`, then `Rscript make_log_dirs.R`, which creates
   `logs_split_rerun/` for the new `me_split_rerun.sh`. If
   `scenario_6/` or `scenario_8/` exist under any of
   `../results/correlated/model_evaluation{,_strategies,_split}/`, move them
   out of the tree. Nothing reads them now, because check, collect and metrics
   walk the grid's combinations, but they are not part of the study.

1. **Main tree.** `me_check.R` finds every scenario 2/3 row missing. It writes
   all 600 to `failed_ids.txt` and sets `me_rerun.sh`'s `-J` to `1-600%10`.
   Any genuinely failed 1/4 runs are included as well, and they need
   rerunning anyway.

   ```bash
   Rscript model_evaluation/me_check.R
   cd model_evaluation/jobscripts
   qsub me_rerun.sh
   ```

2. **Strategies tree.** Wait until step 1 has finished: `me_strategies.R`
   reads the main tree's runs. `me_check.R strategies` doesn't write an id
   file (see its header), so write one from the grid:

   ```bash
   Rscript -e 'source(here::here("model_evaluation/me_config.R")); cat(grid_indices(study_strat, scenario = c(2, 3)), file = study_strat$failed_file, sep = "\n")'
   cd model_evaluation/jobscripts
   qsub -J 1-600%10 me_strat_rerun.sh   # -J on the command line overrides the script's
   ```

3. **Split tree.** This also waits for step 1, because `me_split.R` reads the
   main tree's data, but it does not wait for step 2. `me_split_rerun.sh`'s
   `-J` is already `1-400%10`:

   ```bash
   Rscript -e 'source(here::here("model_evaluation/me_config.R")); cat(grid_indices(study_split, scenario = c(2, 3)), file = study_split$failed_file, sep = "\n")'
   cd model_evaluation/jobscripts
   qsub me_split_rerun.sh
   ```

4. **Check and collect every tree again.** Each tree's `<prefix>_all.RDS`
   has to be rebuilt so that it holds 1–4:

   ```bash
   Rscript model_evaluation/me_check.R
   Rscript model_evaluation/me_check.R strategies
   Rscript model_evaluation/me_check.R split
   cd model_evaluation/jobscripts
   qsub me_strategies_verify.sh
   qsub me_collect.sh
   qsub -v TREE=strategies me_collect.sh
   qsub -v TREE=split me_collect.sh
   # once each tree's collect has finished:
   qsub me_metrics.sh
   qsub -v TREE=strategies me_metrics.sh
   qsub -v TREE=split me_metrics.sh
   ```

   Then render `me_results.qmd`.

**Re-run owed — every tree.** The runs described below predate bug O, which
changed the continuous DGM this study generates from (`docs/DEVLOG.md`),
and use the pre-2026-09-26 scenario numbers (1/4/6/9, now 1/4/6/8). The main,
strategies and split trees are archived by `R/archive_old_results.R` (root
`README.md`, Status, step 0), and all three re-run from empty.

**Changed for the rerun (2026-10-08)**, none of it comparable with the
archived runs:

- the data: correlated covariates, ρ = 0.5 (see "Data"), written under
  `results/correlated/`;
- the candidates: per-arm outcome models in all 9, mean + glm propensity for
  `SL*`, `SL1` = the production `dr_superlearner`, `net2` = interaction lasso
  (see "Candidate models");
- the arms: `cv_indep` and its second fold draw removed, so `me_analysis.R`
  fits only `whole` and the strategies tree carries three arms.

Re-run `jobscripts/me_testing.sh` on the cluster before `me_1.sh`: the
2026-10-08 cluster run that was passing predates these candidate changes.

**Changed 2026-10-14**, not comparable with anything scored before it: the
scoring pipelines' outcome models are per arm, and `me_strategies.R`
recomputes `whole` instead of copying it (see "Status").

What follows is the state of the archived runs.

**The main study completed 358 of 360 runs.** Two runs failed repeatedly
and are permanently excluded. They have no `res_sim_*.RDS`, so every derived
pass skips them by design (exiting 0 with a message rather than erroring) and
they stay excluded consistently across all three trees — see `me_check.R` for
why that means the derived trees' completion criterion is "exactly those 2
missing", not zero.

**Of the 358 strategies-tree runs, 69 had `NA` phi in the automl `holdout`
arm**, from an untrimmed propensity of exactly 0 or 1 (see "Propensity
trimming" above; then put down to H2O returning NA predictions).

**The nuisance-arm reconfiguration does not require rerunning any of it.**
Everything the new arms consume — `data$Y/W/X`, `truth`, `fold_info`, and each
candidate's `tau` — is already saved per replicate. `me_strategies.R` reads
those, so the DGM is never re-run and the 9 candidates are never refit: all
three arms score the *identical* candidate fits, which is what makes the
comparison controlled rather than confounded with fit-to-fit variation. The
only new fitting is the two new nuisance arms themselves, and (separately) the
80:20 arm, which refits candidates because that is its entire purpose.

Cost, roughly: `cv_shared` re-runs the full 10-fold nuisance pipeline, about
what the old `cv` arm cost on its own; `holdout` fits on 25–100-row blocks, so
its cost is dominated by H2O's per-call JVM overhead across 10 blocks rather
than by model size. Budget on the order of the original job's nuisance half.
Don't trust `me_strategies.sh`'s placeholder walltime: check the first
subjobs' `resources_used`, since that per-call overhead is exactly what a
placeholder gets wrong.

**A second local/cluster version mismatch, found while adding the arms.**
`run_xgb_cv()` — fixed an xgboost 3.x API location bug; now reads whichever location the installed version uses, falling back to `which.min()` on the evaluation log. Same class of problem as the `SL2` limitation below, and the same caution applies: local xgboost and cluster xgboost are not the same package.

### Port history (bugs found and fixed)

Unlike every other study in this repo, the ported prototype never completed a run — it referenced an undefined variable and other bugs. Bugs fixed during the port:

- Nonexistent `here("src", ...)` paths
- Undefined variable reference
- A `future` plan leak
- A `create_rf_hyperparams()` typo (`NUL` for `NULL`)
- Missing `on.exit()` call parens
- A dead defensive check (`fix_automl()`)
- Duplicate `collate_predictions()` definition
- Missing `n_cores` argument in `run_all_xgb_nuisance()`
- `mtry` scaling needed for `rf2`/`rf3` after porting (calibrated to `benchtm`'s wider covariate set, crashed `ranger` against this DGM's smaller one)

The 9 candidates' crossfitting scheme *also* changed after the port, from
double crossfitting to single crossfitting — see "Crossfitting strategy"
above. The first 16 `res_sim_*.RDS` files this study produced predate that
change and are not comparable to anything produced after it; they need
deleting before the study is re-run.

H2O AutoML has no `max_runtime_secs` cap in the current code — a real
walltime-uncertainty risk (see "Sizing the array job" above), left alone
since capping it would change what's being measured. Worth revisiting if it
starts causing walltime failures.
