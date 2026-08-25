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
| folds | 10 (all n — see "Crossfitting strategy" below) |
| results | `../results/model_evaluation/scenario_<k>/<n>/res_sim_<run>.RDS` |

**Why only 4 of the 10 scenarios, not all of them like `continuous/`.** This
study's research question — does a cheap proxy loss rank 9 candidate models
the way true PEHE would — doesn't need every CATE-structure re-litigated to
answer. `crossfitting/cf_analysis.R` made the same call for the same reason
(`scenario = c(1, 4, 6, 9)`); this study reuses its exact subset.

**Why 30 runs, not `continuous/`'s 100.** Each replicate here is far more
expensive than a `continuous/` replicate: 9 single-crossfit candidate-model
fits, *plus* two independent nuisance-evaluation pipelines (XGBoost with a
36-combination CV grid search, and H2O AutoML with up to 20 auto-tuned
models), each across 2 fold regimes and 4 nuisance targets. See "Sizing the
array job" below — the real per-replicate cost has never been measured, so
even 30 is provisional.

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

`me_analysis.R` draws two *independent* fold assignments per replicate — one
for the 9 candidates, one for `me_nuisance.R`'s scoring pipelines. That second
draw was justified as preserving honesty. **It does not, and the justification
has been retired** — see the next section. The draw itself is left in place so
`me_analysis.R` keeps producing the shape the 358 completed runs already have;
what changed is that it is now one arm among four rather than the default.

## Candidate models

9 configurations, each fit as its own single-crossfit DR-learner
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

### The four arms

What varies is **what data the evaluation nuisance sees relative to the data
the candidate it is scoring was trained on**. That is the only axis; the
candidate fits are identical across all four arms, because `me_strategies.R`
reads them from the completed results rather than refitting them.

| arm | nuisance trained on | predicted on | row-honest | decoupled from candidate |
|---|---|---|---|---|
| `whole` | all `n` rows | all `n` rows | ✗ | ✗ |
| `cv_indep` | 90%, from a second *independent* fold draw | held-out 10% | ✓ | ✗ (90% overlap) |
| `cv_shared` | the candidate's own `V-1` training folds | the candidate's held-out fold | ✓ | ✗ (identical training set) |
| `holdout` | the candidate's held-out fold **only** | that same fold | ✗ | ✓ |

**Why `cv_indep` is retired as the default.** The second independent draw was
justified as stopping a candidate's `tau_hat` and the nuisance at the same row
being fit on the identical training set. Two things are wrong with that:

- Row-level honesty — whether row `i`'s own `(Y_i, W_i)` entered the model
  predicting at row `i` — holds under **both** a shared and an independent
  draw. Neither ever lets a nuisance see the row it predicts, so the second
  draw buys nothing on the axis it was justified by.
- What it actually changes is the *overlap* of the two training sets for row
  `i`, from 100% down to `(V-1)/V` = **90%**. It removes a tenth of the
  dependence.

`crossfitting/README.md` makes the general version of this point — *"Re-
randomising the stage-2 split cannot remove that dependence"* — which is why
`scf_scf_new` is an arm to be *tested* there rather than the default. The arm
is kept here (it is already computed, so it costs nothing) precisely so the
claim becomes empirical: if `cv_indep` and `cv_shared` rank the candidates
alike, that is the direct demonstration that the extra draw did nothing.

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

**Block size at `n=250`.** `V=10` gives 25-row folds, and `mu0_T`/`mu1_T` fit
on the ~12 control / ~12 treated rows inside one — not estimable. So at
`n=250` only, `holdout_blocks()` pools *adjacent pairs of the candidate's
folds* into 5 blocks of 50. Pooling rather than drawing a fresh 5-fold split
is what keeps each block tied to the candidate's own partition without
refitting any candidate. The cost, stated plainly: a pooled block is only
*half*-decoupled from any single candidate fold model, because the partner
fold was in that model's training set. `cv_shared` stays at `V=10` at every
`n` — it trains on 225 rows at `n=250`, so nothing forces a change there.

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

From these, `calculate_pseudos()` builds a T-learner CATE (`tau_T`), and an
AIPW pseudo-outcome (`phi`; `phi05` also computed with propensity fixed at
0.5, currently unused by any score). `me_metrics.R` scores every candidate's
`tau` against both pipelines under every arm two ways — an
influence-corrected PEHE proxy (`calc_infl_score`) and a DR-risk proxy
(`calc_dr_risk`, MSE against the AIPW pseudo-outcome) — plus true PEHE
(available because this is simulated data).

The column set is *derived* from whatever arms the nuisance list carries, never
enumerated, which is why one `me_per_model()` serves all three result trees:

| tree | arms | score columns |
|---|---|---|
| `model_evaluation` | `cv`, `whole` | 1 + 2×2×2 = **9** |
| `model_evaluation_strategies` | `whole`, `cv_indep`, `cv_shared`, `holdout` | 1 + 2×4×2 = **17** |
| `model_evaluation_split` | `split` | 1 + 2×1×2 = **5** |

The split tree needs no scoring variant because `me_split.R` stores `data` and
`truth` already restricted to its 20% evaluation rows, so every vector
`me_per_model()` touches is the same length.

### Propensity trimming — an open question, deliberately left open

`calculate_pseudos()` divides by `pi * (1 - pi)` with **no** trimming, while
the candidates trim via `trim_ps()` (`me_models.R`). The `whole` arm has always
carried that exposure and completed 358/360 runs. `holdout` and `split` fit
`pi` on 25–100 rows, so their predictions sit much closer to 0/1 and the
weight `1/(pi(1-pi))` can dominate `phi`. The formula is **not** changed —
trimming one arm and not the others would make them non-comparable — but
`me_strategies.R` records per-arm `pi` min/max/quantiles and the max weight in
each run's `pi_diagnostics`, so the decision can be made from measured numbers.

## Files

| file | role |
|---|---|
| `me_config.R` | the parameter grid, results path, and candidate-model list — **the** definition |
| `me_dgms.R` | names this study's slice of `R/dgm_scenarios.R` (`set = "continuous"`) |
| `me_utils.R` | design-matrix prep and fold-splitting (DGM-agnostic, unchanged by the port) |
| `me_models.R` | the 9 candidate CATE-learner configurations and their fitting logic |
| `me_nuisance.R` | the two independent nuisance-evaluation pipelines — see below for why this exists outside the usual 7-file shape |
| `me_analysis.R` | array entry point; one row of the grid per index |
| `me_strategies.R` | second pass over completed runs — adds the `cv_shared` and `holdout` arms, writes to `model_evaluation_strategies` |
| `me_strategies_verify.R` | proves that pass carried the candidates, data, truth and `whole`/`cv_indep` through bit-identically |
| `me_split.R` | the 80:20 arm — the only script that refits the candidates |
| `me_check.R` | finds missing runs, writes `jobscripts/failed_ids.txt`, and updates `-J` and the resource request in the rerun jobscript. Takes a tree: `main` (default) / `strategies` / `split` |
| `me_collect.R` | gathers per-run files into `<prefix>_all.RDS`. Same tree argument |
| `me_metrics.R` | computes `<prefix>_metrics.RDS` (reuses `R/metrics.R::compute_metrics()`). Same tree argument |
| `me_results.qmd` | the results report — see below for why it derives its own quantities |
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

**Why `me_results.qmd` derives its own quantities**: every other study's
report summarises a per-model metric straight out of its `*_metrics.RDS`
(mean bias, mean MSE, coverage). This one can't — the object of interest is
the *ranking* of the 9 candidates within a single (scenario, n, run), so the
report first has to reduce each run's 9x9 score matrix to per-run rank
agreement, top-1 selection accuracy and regret before anything can be
averaged. That derivation lives inline in the `.qmd`, as
`continuous/cts_results.qmd` and `binary/bin_results.qmd` keep theirs; it is
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

The main study (already complete — 358/360):

```bash
Rscript model_evaluation/me_profile.R 1             # smoke-test the profiler locally
qsub model_evaluation/jobscripts/me_profile.sh      # 16 profiling jobs
Rscript model_evaluation/me_profile_summary.R       # writes measured directives into me_1.sh

qsub model_evaluation/jobscripts/me_1.sh            # the study itself - 1-360
Rscript model_evaluation/me_check.R                 # writes failed_ids.txt if any are missing
qsub model_evaluation/jobscripts/me_collect.sh
qsub model_evaluation/jobscripts/me_metrics.sh
```

The nuisance-arm pass. **This does not require rerunning the study** — it
reads the completed results and adds arms to them (see `me_strategies.R`'s
header for why that is sound):

```bash
qsub    model_evaluation/jobscripts/me_strategies.sh      # 1-360, reads the main tree
Rscript model_evaluation/me_check.R strategies            # progress / completion
qsub    model_evaluation/jobscripts/me_strategies_verify.sh
qsub -v TREE=strategies model_evaluation/jobscripts/me_collect.sh
qsub -v TREE=strategies model_evaluation/jobscripts/me_metrics.sh
```

The 80:20 arm, independent of the above:

```bash
qsub    model_evaluation/jobscripts/me_split.sh           # 1-240 (n = 500, 1000 only)
Rscript model_evaluation/me_check.R split
qsub -v TREE=split model_evaluation/jobscripts/me_collect.sh
qsub -v TREE=split model_evaluation/jobscripts/me_metrics.sh
```

```bash
quarto render model_evaluation/me_results.qmd   # the report - needs me_metrics.RDS
```

**Checking progress of the derived trees.** `Rscript me_check.R strategies`
runs `check_failed(..., write = FALSE)` and reports how many runs are done.
Note the completion criterion is *not* zero missing: `study_strat`'s grid is
all 360 rows (the array index has to keep meaning the same grid row), so the
2 runs excluded from the main study are reported missing forever. `me_check.R`
diffs the two trees and separates "missing because there is no source run"
from "missing because this pass failed" — only the latter needs resubmitting,
which is also why `failed_ids_strat.txt` is not written automatically.

Results land in `../results/model_evaluation{,_strategies,_split}/` (siblings
of the repo, as elsewhere). `me_results.qmd` reads `me_metrics.RDS` from
there, so it renders wherever the results are — not on a machine that only has
the repo. The rendered `.html` is gitignored, as every other study's report is.

## Sizing the array job

`me_1.sh` ships with **placeholder** `#PBS -l` lines and trailing `Rscript`
args (`1 1`), the same way `crossfitting/jobscripts/cf_1.sh` did before it
was profiled. `me_profile.R` / `me_profile_summary.R` settle it the same way
`cf_profile.R` / `cf_profile_summary.R` do for `cf_1.sh` (`syrup`,
process-tree filtering, the peak-memory upper bound — see that folder's
README for the shared mechanics) — with two differences specific to this
study:

- **Two knobs, two sequential phases.** `workers` (the `future` multisession
  backend, controlling the 9 candidate models' single-crossfit fold-wise
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

**The main study has completed: 358 of 360 runs.** Two runs failed repeatedly
and are permanently excluded. They have no `res_sim_*.RDS`, so every derived
pass skips them by design (exiting 0 with a message rather than erroring) and
they stay excluded consistently across all three trees — see `me_check.R` for
why that means the derived trees' completion criterion is "exactly those 2
missing", not zero.

**The nuisance-arm reconfiguration does not require rerunning any of it.**
Everything the new arms consume — `data$Y/W/X`, `truth`, `fold_info`, and each
candidate's `tau` — is already saved per replicate. `me_strategies.R` reads
those, so the DGM is never re-run and the 9 candidates are never refit: all
four arms score the *identical* candidate fits, which is what makes the
comparison controlled rather than confounded with fit-to-fit variation. The
only new fitting is the two new nuisance arms themselves, and (separately) the
80:20 arm, which refits candidates because that is its entire purpose.

Cost, roughly: `cv_shared` re-runs the full 10-fold nuisance pipeline, about
what the old `cv` arm cost on its own; `holdout` fits on 25–100-row blocks, so
its cost is dominated by H2O's per-call JVM overhead across 10 blocks rather
than by model size. Budget on the order of the original job's nuisance half —
and re-profile rather than trusting `me_strategies.sh`'s placeholder walltime,
since that per-call overhead is exactly what a placeholder gets wrong.

**A second local/cluster version mismatch, found while adding the arms.**
`run_xgb_cv()` read `best_iteration` from the top level of `xgb.cv()`'s
result. xgboost 3.x moved it into an `early_stop` sub-list, so on local
xgboost (3.2.1.1) it is `integer(0)` for **every** fit at every data size —
`data.frame()` then errors with *"arguments imply differing number of rows:
0, 1"*. The cluster's `R/4.3.2` module still has an older xgboost, which is
why the 358 runs completed: they went through this exact line. `run_xgb_cv()`
now reads whichever location the installed version uses, falling back to
`which.min()` on the evaluation log.

That fallback is an equivalent, not a degradation — early stopping selects the
minimum of the very column being scanned, and the two agree exactly where both
are available (verified: both give iteration 31 on the same fit). So it
changes no number the cluster has already produced. Same class of problem as
the `SL2` limitation below, and the same caution applies: local xgboost and
cluster xgboost are not the same package.

**Port history.** Unlike every other study in this repo, there was no prior
successful run of this study to compare against or to regenerate
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

The 9 candidates' crossfitting scheme *also* changed after the port, from
double crossfitting to single crossfitting — see "Crossfitting strategy"
above. The first 16 `res_sim_*.RDS` files this study produced predate that
change and are not comparable to anything produced after it; they need
deleting before the study is re-run.

H2O AutoML has no `max_runtime_secs` cap in the current code — a real
walltime-uncertainty risk (see "Sizing the array job" above), left alone
since capping it would change what's being measured. Worth revisiting once
real profiling numbers are in hand.
