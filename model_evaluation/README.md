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
(`R/dgm_scenarios.R`, `set = "continuous"`) and the same 7-role file shape
(`config`/`dgms`/`models`/`analysis`/`check`/`collect`/`metrics`) `continuous/`
uses, plus study-specific extras the same way `crossfitting/` has extras
beyond that floor.

## Design

Four scenarios varying the structure of the CATE, crossed with three sample
sizes, 30 runs each — **360 array jobs**.

| | |
|---|---|
| scenarios | 1, 4, 6, 9 (see `R/dgm_scenarios.R`, `DESC_10`) |
| n | 250, 500, 1000 |
| runs | 30 |
| folds | 5 at n=250, else 10 |
| results | `../results/model_evaluation/scenario_<k>/<n>/res_sim_<run>.RDS` |

**Why only 4 of the 10 scenarios, not all of them like `continuous/`.** This
study's research question — does a cheap proxy loss rank 9 candidate models
the way true PEHE would — doesn't need every CATE-structure re-litigated to
answer. `crossfitting/cf_analysis.R` made the same call for the same reason
(`scenario = c(1, 4, 6, 9)`); this study reuses its exact subset.

**Why 30 runs, not `continuous/`'s 100.** Each replicate here is far more
expensive than a `continuous/` replicate: 9 double-crossfit candidate-model
fits, *plus* two independent nuisance-evaluation pipelines (XGBoost with a
36-combination CV grid search, and H2O AutoML with up to 20 auto-tuned
models), each across 2 fold regimes and 4 nuisance targets. See "Sizing the
array job" below — the real per-replicate cost has never been measured, so
even 30 is provisional.

## Candidate models

9 configurations, each fit as its own double-crossfit DR-learner
(`me_models.R`, `run_all_candidate_models()`):

| id | family | hyperparameters |
|---|---|---|
| `rf1` | random forest | ranger defaults (`mtry = floor(sqrt(p))`) |
| `rf2` | random forest | `mtry = p, max.depth = 5` (every covariate each split) |
| `rf3` | random forest | `mtry = ceiling(p/2), max.depth = 3` |
| `net1` | elastic net | `alpha=1` (lasso) |
| `net2` | elastic net | `alpha=0` (ridge) |
| `net3` | elastic net | `alpha=0.5` |
| `SL1` | SuperLearner | `glmnet, ranger, earth, gam, mean` |
| `SL2` | SuperLearner | `glmnet, xgboost, cforest, earth, gam, mean` |
| `SL3` | SuperLearner | `svm, nnet, mean` |

`p = ncol(X)` (7-9 across this study's scenarios). `rf2`/`rf3`'s `mtry` is
scaled to `p` rather than fixed, unlike the ported prototype's original
`mtry = 30`/`mtry = 10` — those were calibrated to `benchtm`'s much wider
covariate set and error out of `ranger` ("mtry can not be larger than number
of variables in data") against this DGM's smaller one. Found by
`me_testing.R full`, not by inspection.

This is a *different* estimator family from `R/cate_models.R`'s DR-learner
(e.g. the elastic-net arms wrap a single learner via
`create.Learner("SL.glmnet", ...)` inside `SuperLearner()`, a code path
`R/cate_models.R` doesn't have) — not unified into `R/`, since only this
study needs it.

## Nuisance-evaluation pipelines

A second, *independent* estimate of 4 nuisance targets (`mu0_T`, `mu1_T`,
`mu_DR`, `pi`) — used only to build proxy scores for the 9 candidates above,
never to fit them:

- **XGBoost** — a hand-tuned CV grid (`eta`, `max_depth`, `subsample`,
  `colsample_bytree`).
- **H2O AutoML** — up to 20 auto-tuned models per target
  (`exclude_algos = c("DeepLearning", "XGBoost")`).

Each runs under two fold regimes: **`cv`** (leave-one-fold-out) and
**`whole`** (no split). A third regime, **`infold`** (resubstitution — fit
and predict on the same held-in fold), was in the ported prototype and has
been removed rather than left unused. `me_nuisance.R` keeps a comment
documenting exactly what the removed code did, including a documented risk
in the AutoML version's `infold` branch (a `W`-recoding line that was never
confirmed correct, since this pipeline has never completed a run) — a future
re-add should verify that rather than restoring it unchanged.

From these, `calculate_pseudos()` builds a T-learner CATE (`tau_T`), and an
AIPW pseudo-outcome (`phi`; `phi05` also computed with propensity fixed at
0.5, currently unused by any score). `me_metrics.R` scores every candidate's
`tau` against both pipelines under both fold regimes two ways — an
influence-corrected PEHE proxy (`calc_infl_score`) and a DR-risk proxy
(`calc_dr_risk`, MSE against the AIPW pseudo-outcome) — plus true PEHE
(available because this is simulated data). That's 2 pipelines × 2 fold
regimes × 2 score types + `true_pehe` = **9 score columns** per candidate per
run.

## Files

| file | role |
|---|---|
| `me_config.R` | the parameter grid, results path, and candidate-model list — **the** definition |
| `me_dgms.R` | names this study's slice of `R/dgm_scenarios.R` (`set = "continuous"`) |
| `me_utils.R` | design-matrix prep and fold-splitting (DGM-agnostic, unchanged by the port) |
| `me_models.R` | the 9 candidate CATE-learner configurations and their fitting logic |
| `me_nuisance.R` | the two independent nuisance-evaluation pipelines — see below for why this exists outside the usual 7-file shape |
| `me_analysis.R` | array entry point; one row of the grid per index |
| `me_check.R` | finds missing runs, writes `jobscripts/failed_ids.txt` |
| `me_collect.R` | gathers per-run files into `me_all.RDS` |
| `me_metrics.R` | computes `me_metrics.RDS` (reuses `R/metrics.R::compute_metrics()`) |
| `me_testing.R` | verification checks — run before submitting anything |
| `me_profile.R` | timing / memory / CPU sweep over `(workers, n_cores)`, instrumented with `syrup` |
| `me_profile_summary.R` | turns the sweep into PBS directives and writes them into `me_1.sh` |

**Why `me_nuisance.R` exists outside the `config`/`dgms`/`models`/`analysis`/
`check`/`collect`/`metrics` shape**: no other study in this repo has a
second, independent nuisance-estimation pipeline used purely to *score*
candidate models rather than fit them — it doesn't map onto any of those 7
roles. `crossfitting/` is the precedent for a study needing files beyond that
floor (`cf_testing.R`, `cf_profile.R`, `cf_profile_summary.R`,
`cf_diagnose_*.R`, `cf_results.R`/`.qmd`); the 7-file shape is a floor, not a
ceiling.

**Note on `me_metrics.R`**: `R/metrics.R::compute_metrics()` always does
`true_tau <- sim_res$truth$tau` — there's no equivalent of the ported
prototype's `calc_metrics(..., truth_avail = FALSE)` branch (for scoring
without ground truth). Since `me_analysis.R` only ever runs against
`dgm_scenarios.R`-generated data, which always carries truth, this is a real
but low-risk scope reduction from the original code.

## Running it

```bash
Rscript model_evaluation/me_testing.R              # quick: structure + regression checks
Rscript model_evaluation/me_testing.R full          # adds the real XGBoost-CV/H2O AutoML pipelines (slow)

Rscript model_evaluation/me_profile.R 1             # smoke-test the profiler locally
qsub model_evaluation/jobscripts/me_profile.sh      # 16 profiling jobs
Rscript model_evaluation/me_profile_summary.R       # writes measured directives into me_1.sh

qsub model_evaluation/jobscripts/me_1.sh            # the study itself - 1-360
Rscript model_evaluation/me_check.R                 # writes failed_ids.txt if any are missing
qsub model_evaluation/jobscripts/me_collect.sh
qsub model_evaluation/jobscripts/me_metrics.sh
```

Results land in `../results/model_evaluation/` (a sibling of the repo, as
elsewhere).

## Sizing the array job

`me_1.sh` ships with **placeholder** `#PBS -l` lines and trailing `Rscript`
args (`1 1`), the same way `crossfitting/jobscripts/cf_1.sh` did before it
was profiled. `me_profile.R` / `me_profile_summary.R` settle it the same way
`cf_profile.R` / `cf_profile_summary.R` do for `cf_1.sh` (`syrup`,
process-tree filtering, the peak-memory upper bound — see that folder's
README for the shared mechanics) — with two differences specific to this
study:

- **Two knobs, two sequential phases.** `workers` (the `future` multisession
  backend, controlling the 9 candidate models' double-crossfit fold-pair
  fitting) and `n_cores` (XGBoost's `nthread` / H2O's `nthreads`, controlling
  the nuisance-evaluation pipelines) parallelise two *sequential* phases of
  one replicate, not one combined computation — `me_profile.R` times them
  separately so the phase breakdown is visible, and every place that would
  otherwise add the two knobs together instead takes `max(workers, n_cores)`,
  since they're never both active at once.
- **H2O's JVM may not be a tracked child process.** `syrup` tracks R
  processes via `ps` parent/child relationships; whether H2O's own Java
  process shows up in that tree is unconfirmed. `me_profile.R` prints the
  full, unfiltered process list on every run so this is visible rather than
  silently under-counting H2O's memory footprint — treat the peak-memory
  recommendation with real caution (cross-check against
  `qstat -fx <jobid> | grep resources_used`) until that's settled.

**The array's concurrency throttle (`-J 1-360%N`) is a separate problem the
profiler doesn't solve.** Each concurrent task starts its own H2O JVM
cluster with a `mem = "10G"` heap — nothing like `continuous/`'s `%190` or
`crossfitting/`'s `%380` is safe here. `N` needs setting from the specific
HPC queue's real memory/fair-share limits, which needs consulting whoever
manages the allocation, not just measuring one replicate's own footprint.

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

`SL2`'s library is unchanged from the ported prototype - this is a local
package-version mismatch, not something the DGM swap introduced, and exactly
the situation `.claude/CLAUDE.md` already documents ("local R is 4.5.3 ...
the cluster runs R 4.3.2-gfbf-2023a ... if package-version issues come up,
this is why"). Confirm whether this reproduces on the cluster's actual
`R/4.3.2-gfbf-2023a` module + `sim-env` conda env before relying on a clean
local `me_testing.R full` run as a sign this study is ready to submit.

## Status

**First run, not a re-run.** Unlike every other study in this repo, there is
no prior successful run of this study to compare against or to regenerate
after a bug fix — the ported prototype never completed a run (it referenced
an undefined variable, among other bugs; see the git history around the
`me_*` files for what changed during the port). Bugs fixed during the port:
nonexistent `here("src", ...)` paths, the undefined-variable bug just
mentioned, a `future` plan that was never restored on exit, a
`create_rf_hyperparams()` typo (`NUL` for `NULL`), an `on.exit()` missing its
call parens (so the H2O cluster's shutdown safety net never actually fired),
a dead defensive check (`fix_automl()`), a duplicated `collate_predictions()`
now sourced from `R/utils.R` instead, and `run_all_xgb_nuisance()` being
called without `n_cores` (so XGBoost silently never used the intended thread
count).

One thing *did* change beyond pure plumbing, found by `me_testing.R full`
rather than by inspection: `rf2`/`rf3`'s fixed `mtry` values (30, 10) were
calibrated to `benchtm`'s wider covariate set and crash `ranger` against this
DGM's smaller one (7-9 columns) — see the candidate-models table above for
the scaled replacement. Every other candidate model's and both nuisance
pipelines' modeling choices are unchanged.

H2O AutoML has no `max_runtime_secs` cap in the current code — a real
walltime-uncertainty risk (see "Sizing the array job" above), left alone
since capping it would change what's being measured. Worth revisiting once
real profiling numbers are in hand.
