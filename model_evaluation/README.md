# Model evaluation — does cheap proxy scoring pick the right CATE learner?

Given a basket of candidate CATE-learner configurations, can cheap proxy loss
functions — built from *independent* nuisance estimators — rank them the way
the (normally unobservable) true PEHE would? This is a model-*selection*
study, not another "does estimator X recover the true CATE" study: the
object of interest is the proxy-vs-truth comparison `me_metrics.R` produces,
not any one candidate's own accuracy.

It uses the shared DGM (`R/dgm_scenarios.R`, `set = "continuous_corr_0.5"`)
and the usual `config`/`dgms`/`models`/`analysis`/`check`/`collect`/`metrics`
file shape, plus study-specific extras.

## Design

Four scenarios varying the structure of the CATE, crossed with three sample
sizes, 100 runs each — **1200 array jobs**.

| | |
|---|---|
| scenarios | 1–4, the main four (see `R/dgm_scenarios.R`, `DESC_10`) |
| covariates | correlated, latent ρ = 0.5 (`me_dgms.R`'s `ME_RHO`) |
| n | 250, 500, 1000 |
| runs | 100 |
| folds | 10 (all n — see "Crossfitting strategy" below) |
| results | `../results/correlated/model_evaluation/scenario_<k>/<n>/res_sim_<run>.RDS` |

### Data

`continuous_corr_0.5` is `sample_size/continuous/`'s scenario table with the
covariates drawn from the copula at latent correlation 0.5 — the primary arm of
`sample_size/`. See `sample_size/README.md` for bW, the ATE, SD(τ) and the
correlation between m0 and τ. Only ρ = 0.5 is run, so ρ is a constant rather
than a grid column.

**Why only 4 of the 10 scenarios.** The question — does a cheap proxy loss rank
9 candidate models the way true PEHE would — doesn't need every CATE structure.
The four are the main four every other study reports: 1 null, 2 simple (−X4),
3 complex (2·X3 + 0.5·X4 − 0.5·X4·X5), 4 non-linear (cos X4). The grid lists
them as `c(1, 4, 2, 3)`, not `1:4`, to keep array indices: see `me_config.R`'s header.

**100 runs.** Each replicate is far more expensive than a `sample_size/` one —
9 single-crossfit candidate fits, *plus* two independent nuisance-evaluation
pipelines (XGBoost with a 36-combination CV grid search, and H2O AutoML with up
to 20 auto-tuned models), each across the nuisance arms and 2 nuisance targets.

## Crossfitting strategy

All 9 candidates use a single leave-one-fold-out crossfit ("`scf_scf`" in
`crossfitting/`'s naming): nuisances trained on `V-1` folds and predicted on
the held-out fold, with the stage-2 regression trained on the *same* folds.
This follows `crossfitting/`'s comparison of alternatives — see
`crossfitting/README.md` and `R/cate_models.R`'s header.

Single crossfitting applies to all 9 candidates uniformly, not a per-learner
OOB/crossfit split, even though `rf1-3` (`ranger`) could take whole-sample OOB
predictions natively. `net1-3` and `SL1-3` have no OOB analogue — `glmnet` has
no bagging, and `SuperLearner`'s internal CV selects ensemble weights rather
than producing an honest prediction of a training row
(`crossfitting/README.md`'s "DR-learner, SuperLearner" section). Mixing OOB
forests with crossfit everything-else would confound learner choice with
crossfitting scheme in the ranking. All `n` use `V=10`: single crossfitting
trains each fold's stage-2 model on `(V-1)/V` of the data, so there's no need
to shrink `V` at `n=250`.

The candidates' crossfitting is fixed and is **not** what this study varies.
The one exception is the 80:20 arm, where refitting them on a training split is
the point (see below). `me_analysis.R` draws one fold assignment per
replicate; the fold-based scoring arms reuse it.

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
| `SL3` | SuperLearner | `svm, nnet, mean` — a weak baseline | mean + glm, as SL1 |

`p = ncol(X)` = 10 in every scenario (X1–X5, X01–X05); `rf2`/`rf3`'s `mtry` is
scaled to `p`.

**Outcome models are per arm (T-learner) in all 9**, as in `R/cate_models.R`'s
DR-learners: one model per arm, fit on `X` with that arm's training rows,
predicting every held-out row. A single model on `cbind(W, X)` read at
`W = 0` and `1` would make the `net*` candidates' `Y1.hat − Y0.hat` a constant
(W enters as a main effect only), and offers `W` at only ~30% of `rf1`'s splits.

**Propensity.** `rf*` and `net*` fit it with the candidate's own learner. The
SuperLearner candidates all use the production propensity library, mean + glm:
every study here is an RCT at P(W = 1) = 0.5, and a flexible learner only fits
noise that the DR weight `1/(e(1 − e))` amplifies (`R/sl_library.R`'s header).
Every candidate trims to [0.05, 0.95].

**`SL1` is the production estimator, not a copy of it.** Each fold calls
`R/cate_models.R::sl_split_fit()` and stage 2 is `stage_2_sl()`, with
`sl_libraries(n)` on this study's `fold_indices` — the same `scf_scf` scheme
`nuisance_sl()` uses — so the study answers directly whether a proxy would
pick the estimator the other studies report. Its dropped learners are saved
as `SL1$sl_dropped`. The libraries are sized to the rows the candidate is fit
on (`n`, or the 80% in the split arm), so `SL.earth` is in the outcome library
only from 500 rows.

**`net2` is the interaction lasso**: `X` expanded to its 10 main effects and 45
pairwise products (`net_design()`, the same expansion as `R/sl_library.R`'s
`SL.glmnet.int`), at `lambda.min`, in all three nuisances and stage 2. It is
the only linear-family candidate that can represent scenario 3's X4·X5 effect,
and it keeps the three `net*` candidates apart: lasso, ridge and elastic net at
`lambda.min` on 10 covariates give near-identical τ, and near-ties make top-1
selection a coin toss.

**`SL3` is a weak baseline.** `SL.nnet`'s defaults (2 hidden units, no
weight decay) are unstable on a noisy pseudo-outcome. It is kept so the set
has a candidate a good proxy should avoid, which is what regret measures.

`rf*`, `net*`, `SL2` and `SL3` are this study's own code path (e.g. the
elastic-net arms wrap a single learner via `create.Learner("SL.glmnet", ...)`
inside `SuperLearner()`), not in `R/`, since only this study needs it.

## Nuisance-evaluation pipelines

A second, *independent* estimate of 2 nuisance targets (`mu_DR`, `pi`) —
used only to build proxy scores for the 9 candidates above, never to fit
them:

- **XGBoost** — a hand-tuned CV grid (`eta`, `max_depth`, `subsample`,
  `colsample_bytree`).
- **H2O AutoML** — up to 20 auto-tuned models per target
  (`exclude_algos = c("DeepLearning", "XGBoost")`).

**Outcome models are per arm (T-learner)** in both pipelines and every arm:
`mu0_DR` is fit on the training set's `W = 0` rows and `mu1_DR` on its `W = 1`
rows, each with the pipeline's own tuning. This is for **consistency** with the
candidates, not a fix for a demonstrated accuracy problem (`docs/DEVLOG.md`,
bug R). `me_nuisance.R` refuses a per-arm fit with fewer than
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

All three are fit by `me_strategies.R`. `me_analysis.R` also fits `whole` for
the main tree, but the strategies pass recomputes it: the strategies tree is
the scored one, and the main tree supplies only the candidates, data, truth and
`fold_info`.

**Why there is no `cv_indep` arm.** An arm fitting the nuisance over a second,
*independent* fold draw buys nothing on row-level honesty — whether row `i`'s
own `(Y_i, W_i)` entered the model predicting at row `i` — which holds under
both a shared and an independent draw. It only cuts the overlap of the two
training sets for row `i` from 100% to `(V-1)/V` = 90% (`crossfitting/README.md`:
*"Re-randomising the stage-2 split cannot remove that dependence"*), at the
cost of a full 10-fold nuisance pipeline per replicate.

**What `holdout` gives up.** It is resubstitution — fit and predict on the
same block — so every row's own `Y` is in the model that predicts it, `mu_DR`
interpolates it, and `phi`'s AIPW correction `(Y - mu_DR)` collapses toward
zero. That is the arm's defining property, not a defect: it is the only
fold-based arm fully decoupled from the candidate, and row-level honesty and
decoupling cannot both be had from a single small block. `cv_shared` is the
row-honest comparator, and the pair makes the cost of each visible. In
`run_automl_holdout()`, `W` is the caller's numeric 0/1 vector: do not add an
`as.numeric(W) - 1` recode, which would send it to −1/0.

**Block size at `n=250`.** `V=10` gives 25-row folds, too few to fit
`mu_DR`/`pi` on. So at `n=250` only, `holdout_blocks()` pools *adjacent pairs
of the candidate's folds* into 5 blocks of 50, which keeps each block tied to
the candidate's own partition without refitting any candidate. The cost: a
pooled block is only *half*-decoupled from any single candidate fold model,
because the partner fold was in that model's training set. `cv_shared` stays
at `V=10` at every `n`.

**Per-arm fits are thin in `holdout`.** Each block's `mu0_DR` and `mu1_DR` see
only one treatment arm's rows: about 25 at `n=250` (pooled blocks of 50) and
`n=500` (unpooled folds of 50), 50 at `n=1000`. Expect degenerate AutoML fits,
and `automl`/`holdout` NA, more often than in the other arms.

### The 80:20 arm

A fourth position on the same axis, run as its own study (`me_split.R`,
`results/model_evaluation_split`) rather than another column, because three
things move together: `tau_hat` exists only on the 20%, the evaluation rows
*are* that 20%, and `true_pehe` therefore has to be computed on those same
rows or the reference is not matched to the proxies.

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
at 0.5. `me_metrics.R` scores every candidate's `tau` against both pipelines
under every arm, plus true PEHE (available because this is simulated data).

**Three score families, each in two propensity regimes** — 10 score types per
(arm × pipeline):

| family | estimated `pi` | fixed `pi = 0.5` | what it is |
|---|---|---|---|
| influence, T-learner surrogate | `inflt` | `inflt05` | influence-corrected PEHE proxy (`calc_infl_score`) with the plug-in `mu1_DR − mu0_DR` as its surrogate. The IF proxy the thesis figures and tables report (`me_selection.R`) |
| influence, φ surrogate | `infl` | `infl05` | the same score with `phi` / `phi05` as its surrogate |
| DR risk | `dr` | `dr05` | MSE against the AIPW pseudo-outcome (`calc_dr_risk`) |
| calibration | `calq5`, `calq10` | `cal05q5`, `cal05q10` | DR calibration over K quantile groups (`calc_cal_score`) |

Column names are `<score_type>_<arm>_<pipeline>`; the propensity regime and
the group count are folded into that first token (`cal05q10` = calibration,
fixed `pi`, K = 10) rather than added as a fourth field, because
`me_results.qmd` recovers the design axes by splitting the name. `"05"` never
appears in a score type unless the propensity is fixed — `calq10` carries a
`"10"`, never an `"05"` — so the regime is recoverable with `grepl("05", .)`.

Scores are a pure post-hoc function of `<prefix>_all.RDS`: adding a score
family means re-running `me_metrics.R` and nothing else. This is separate from
the propensity *inside the candidates*, which builds each candidate's own
pseudo-outcome from a trimmed estimated `W.hat`; changing that changes every
stored `tau`.

#### Fixing the propensity at 0.5 is an oracle, not an approximation

`R/dgm_scenarios.R` assigns treatment with `W <- rbinom(n, 1, 0.5)` — a fair
coin, independent of `X`, in **every** scenario. So 0.5 *is* the true
propensity, the `*05` scores are oracle-π scores, and the `dr` / `dr05`
contrast isolates exactly one thing: what it costs to estimate a propensity
that never needed estimating. This is the RCT setting the study is about, so
the answer is not incidental.

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
tree, whose truth keeps only `tau`.

The split tree needs no scoring variant because `me_split.R` stores `data` and
`truth` already restricted to its 20% evaluation rows, so every vector
`me_per_model()` touches is the same length.

### Propensity trimming — [0.05, 0.95], every arm and pipeline

`calculate_pseudos()` trims the scoring propensity to [0.05, 0.95] with
`trim_ps()`, as the 9 candidates' propensities are (`me_models.R`) and as
every other study's are. It applies to every arm and both pipelines alike,
which keeps them comparable. The true propensity is 0.5, so trimming can only
move an estimate that is already wrong toward the truth; the estimated- vs
fixed-π contrast measures the cost of estimating a *trimmed* propensity. `pi`
keeps the raw prediction, `pi_trim` the trimmed one, and `me_strategies.R`'s
per-run `pi_diagnostics` (min/max/quantiles, max weight, `n_extreme` = rows
trimmed) are on the raw `pi`.

**Why.** Untrimmed, AutoML's propensity is sometimes **exactly 0 or 1** when a
model predicts the rows it was fit on (`holdout`, `whole`), which makes phi
NaN or ±Inf and every estimated-π score undefined for the run (`docs/DEVLOG.md`,
bug S). `me_selection.R` keeps only finite scores, so that cannot recur
silently.

The trimming is applied **when scoring**: `me_metrics.R`'s `arm_scores()`
recomputes phi from each arm's saved `mu0_DR`/`mu1_DR`/`pi`, and checks it
against the saved phi on the rows the trimming did not move, which catches Y/W
that are not the nuisance's rows. `me_strategies_verify.R` checks completeness
the same way (no NA prediction, finite trimmed phi) and tallies the runs that
had a 0/1 propensity. Trimmed `holdout` propensities still pile up at the
bounds, where the weight `1/(π(1−π))` is about 21 against 4 at π = 0.5, so
those proxies are defined in every run but noisy.

## Files

| file | role |
|---|---|
| `me_config.R` | the parameter grid, results path, and candidate-model list — **the** definition |
| `me_dgms.R` | names this study's slice of `R/dgm_scenarios.R` (`set = "continuous_corr_0.5"`) |
| `me_utils.R` | design-matrix prep and fold-splitting (DGM-agnostic) |
| `me_models.R` | the 9 candidate CATE-learner configurations and their fitting logic |
| `me_nuisance.R` | the two independent nuisance-evaluation pipelines |
| `me_analysis.R` | array entry point; one row of the grid per index |
| `me_strategies.R` | second pass over completed runs — fits all three arms (`whole`, `cv_shared`, `holdout`) without refitting the candidates, writes to `model_evaluation_strategies` |
| `me_strategies_verify.R` | proves that pass carried the candidates, data, truth and `fold_info` through bit-identically and that all three arms are complete (no NA prediction, finite phi with π trimmed), and tallies the runs with a propensity of exactly 0 or 1 |
| `me_split.R` | the 80:20 arm — the only script that refits the candidates |
| `me_check.R` | finds missing runs, writes `jobscripts/failed_ids.txt`, and updates `-J` and the resource request in the rerun jobscript. Takes a tree: `main` (default) / `strategies` / `split` |
| `me_collect.R` | gathers per-run files into `<prefix>_all.RDS`. Same tree argument |
| `me_metrics.R` | computes `<prefix>_metrics.RDS` (reuses `R/metrics.R::compute_metrics()`). Same tree argument |
| `me_results.qmd` | the results report |
| `me_testing.R` | verification checks — run before submitting anything |

`me_results.qmd` derives its own quantities rather than summarising a
per-model metric: the object of interest is the *ranking* of the 9 candidates
within a single (scenario, n, run), so each run's score matrix is first reduced
to rank agreement, top-1 selection accuracy and regret. That derivation is
inline in the `.qmd`, not in `R/figures.R`, whose `summarise_metrics()` is
built for bias/MSE/correlation columns.

`R/metrics.R::compute_metrics()` always reads `sim_res$truth$tau`; there is no
scoring-without-truth path, which is fine because every run here is simulated.

## Running it

**Locally**, run only the fast check. `full` mode is the slow half and `SL2`
fails on the Windows machine for a package-version reason (see "Known
local-environment limitation"), so run it on the cluster instead:

```bash
Rscript model_evaluation/me_testing.R               # structure + regression checks + arm plumbing
qsub    model_evaluation/jobscripts/me_testing.sh   # the full suite, with real XGB-CV/H2O AutoML
```

`me_testing.sh` also checks that `sim-env` carries `h2o`/`xgboost`/`caret`
with a working Java runtime, and that `SL2`'s local failure does not reproduce
there.

Every `qsub` is run from inside `model_evaluation/jobscripts/` (the jobscripts
`cd "${PBS_O_WORKDIR}/.."` and read their id file from `${PBS_O_WORKDIR}`);
the commands below are shown from the repo root for brevity.

The main study:

```bash
qsub model_evaluation/jobscripts/me_1.sh            # the study itself - 1-1200
Rscript model_evaluation/me_check.R                 # writes failed_ids.txt if any are missing
qsub model_evaluation/jobscripts/me_collect.sh
qsub model_evaluation/jobscripts/me_metrics.sh
```

The nuisance-arm pass. It reads the completed main tree and adds arms to it
(see `me_strategies.R`'s header), so run it after the main study:

```bash
qsub    model_evaluation/jobscripts/me_strategies.sh      # 1-1200, reads the main tree
Rscript model_evaluation/me_check.R strategies            # progress / completion
qsub    model_evaluation/jobscripts/me_strategies_verify.sh
qsub -v TREE=strategies model_evaluation/jobscripts/me_collect.sh
qsub -v TREE=strategies model_evaluation/jobscripts/me_metrics.sh
```

The 80:20 arm, which reads the main tree's data but not the strategies tree:

```bash
qsub    model_evaluation/jobscripts/me_split.sh           # 1-800 (n = 500, 1000 only)
Rscript model_evaluation/me_check.R split
qsub -v TREE=split model_evaluation/jobscripts/me_collect.sh
qsub -v TREE=split model_evaluation/jobscripts/me_metrics.sh
```

```bash
quarto render model_evaluation/me_results.qmd   # the report - needs me_metrics.RDS
```

**Rerunning part of a derived tree.** `me_check.R strategies` / `split` runs
`check_failed(..., write = FALSE)`: `study_strat`'s grid is all 1200 rows (the
array index has to keep meaning the same grid row), so runs missing from the
main tree are reported missing forever, and `me_check.R` separates "no source
run" from "this pass failed". Write the id file yourself, e.g. for scenarios
2 and 3:

```bash
Rscript -e 'source(here::here("model_evaluation/me_config.R")); cat(grid_indices(study_strat, scenario = c(2, 3)), file = study_strat$failed_file, sep = "\n")'
cd model_evaluation/jobscripts
qsub -J 1-600%10 me_strat_rerun.sh   # -J on the command line overrides the script's
```

`study_split` and `me_split_rerun.sh` work the same way.

**Re-scoring an already-collected tree** means re-running only
`me_metrics.sh`, which rewrites `me_metrics.RDS` in place; `me_all.RDS` is
untouched. Keep a copy of the old file and check that the existing columns
come back bit-identical.

Results land in `../results/correlated/model_evaluation{,_strategies,_split}/`
(outside the repo, as elsewhere). `me_results.qmd` reads `me_metrics.RDS` from
there, so it renders only where the results are.

## Sizing the array job

The `#PBS -l` lines and trailing `Rscript` args (`workers`/`n_cores`) are set
by hand (see `docs/DEVLOG.md`'s "Resource profiling (removed)"):

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
`check_failed()` writes there anyway. `me_strategies.sh` has no workers, and
its fold loops run one after another, so its 10 cores only feed XGBoost/H2O
threads on small data. Keep in mind:

- **Two knobs, two sequential phases.** `workers` (the `future` multisession
  backend, for the 9 candidates' fold-wise fitting) and `n_cores` (XGBoost's
  `nthread` / H2O's `nthreads`, for the nuisance-evaluation pipelines)
  parallelise two *sequential* phases of one replicate. So `ncpus` needs to
  cover `max(workers, n_cores)`, not the sum.
- **H2O's JVM is a separate Java process.** Each task starts its own H2O JVM
  with a `mem = "10G"` heap, so `mem=` has to cover that on top of R - and in
  `me_1.sh`/`me_split.sh` on top of the 10 worker sessions too, which stay
  alive through the nuisance phase. Check the request against
  `qstat -fx <jobid> | grep resources_used`.
- **H2O AutoML has no `max_runtime_secs` cap**, since capping it would change
  what is measured; it is the main source of walltime uncertainty.

**The array's concurrency throttle (`-J 1-1200%N`) is a separate problem from
the per-task request.** Too many concurrent H2O JVMs make the H2O calls fail,
so every array runs at `%10`. `me_check.R` passes `throttle = 10` to
`check_failed()`, so the `-J` it writes into `me_rerun.sh` stays at `%10` too
(its default is `%100`). With `%10` the main array takes about 120 times one
run's walltime (1200 / 10).

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

On the cluster's xgboost 1.7.11.1 `me_testing.sh` passes, so `SL.xgboost` is
kept. Don't read a clean or failing local `me_testing.R full` run as a sign
about the cluster. `run_xgb_cv()` reads the best round from whichever location
the installed xgboost uses, for the same reason.
