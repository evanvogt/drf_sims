# Crossfitting comparison

Does the double crossfitting procedure used throughout this study earn its cost?

Everywhere else in `drf_sims`, two-stage estimators fit nuisances over all `C(V,2)`
fold *pairs*; `collate_predictions` (`R/utils.R`) assembles an `n x V` pseudo-outcome
matrix whose column `k` is untouched by fold `k`, and the second stage trains on
`po_matrix[, k]` before predicting fold `k`. At `V = 10` that is 45 nuisance fits
instead of 10 — 4.5x the cost — and it has never been benchmarked.

**The hypothesis.** Under ordinary crossfitting, a stage-2 training row `i` carries a
pseudo-outcome from a nuisance model trained on every fold except `i`'s, which
*includes* the held-out test fold. So if double crossfitting matters, `scf_scf`
should differ from `dcf`. If it doesn't, the cheap procedure is enough and the
rest of the study can be sped up 4.5x.

**Every DR arm fits its outcome model separately in each treatment arm** (a
T-learner), as the production DR-learners in `R/cate_models.R` do, so the arms
differ only in how they split the sample. The per-arm fits are production's own
code, not copies: `t_learner_rf_split` / `nuisance_rf` for the forests and
`sl_split_fit` / `nuisance_sl` for SuperLearner. (Until 2026-09-27 most arms used
an S-learner outcome model on `cbind(W, X)`; that design, and the S-vs-T arms
`scf_oob_t` / `oob_oob_s` / `oob_oob_manual` that went with it, are gone.)

## Design

Fixed at **n = 500, V = 10, scenarios 1 / 4 / 6 / 8, 100 runs** (400 array jobs).
All variants within a replicate share one fold assignment, so differences are
attributable to the procedure rather than to the fold draw.

Every arm is scored twice: against the known truth on the **training** sample, and
on an independent **test** sample of 2000 drawn from the same DGP.

**There is no optimism to detect, and the train-to-test gap is not an overfitting
diagnostic here.** Optimism is what you see when a model is scored against the
labels it was fit to. Every arm here is scored against the *known true CATE*, while
the label the stage-2 model saw is a noisy pseudo-outcome.

What the two scoring sets are actually for: the **training** score is the estimand
the study cares about (CATEs for the units you have), and the **test** score
measures how the fitted surface generalises to fresh covariate draws, which is the
only like-for-like comparison of the fitted final models across procedures.

**Two test scores, because crossfit and whole-sample arms predict differently.**
A crossfit arm ends up with `V` fitted models and averages their test predictions —
that is the estimator you would actually deploy, so `mse` uses it — but the
averaging is a variance-reducing ensemble stacked on top of the honesty effect
being studied. `mse_test_single` scores each fold model on the test set separately
and averages the `V` scores, which is the like-for-like reading against a
whole-sample arm's single model. For whole-sample arms the two coincide, so any
gap between them is the ensembling effect alone (`cf_ensemble_effect.png`).

### DR-learner, random forest (4 arms)

| id | stage 1 (nuisances) | stage 2 (final model) |
|---|---|---|
| `dcf` | double CF over fold pairs | crossfit, same folds, column `k` (**status quo**) |
| `scf_scf` | single CF, leave-one-fold-out | crossfit, same folds |
| `scf_oob` | single CF | whole sample, **OOB** predictions |
| `oob_oob` | whole sample, **OOB** | whole sample, **OOB** (production `dr_random_forest`) |

The crossfit arms fit one forest per treatment arm on each split's training rows
and predict the held-out rows (`t_learner_rf_split`, via `rf_split_fit`, which
adds the propensity forest and the marginal `Y.hat.cf` forest the causal forest
arms use). `oob_oob` calls production's `nuisance_rf` itself. Per-arm forests make
the OOB arm honest without any workaround: a unit's own-arm prediction is
out-of-bag, and its other-arm prediction comes from a forest that never saw it.
So `scf_oob` vs `oob_oob` isolates crossfit vs OOB nuisances directly.

### DR-learner, SuperLearner (2 arms)

`dcf`, `scf_scf`. Each split's nuisances are production's `sl_split_fit`: one
SuperLearner per treatment arm plus the propensity, with the per-nuisance
libraries from `R/sl_library.R`. `scf_scf` calls production's `nuisance_sl`
itself, so it is `dr_superlearner`'s stage 1 by construction; `dcf` runs the same
per-split fit over the 45 fold pairs. At n = 500 each per-arm outcome fit sees
~200 rows under `dcf` and ~225 under `scf_scf`, which is what `sl_libraries(n)$Y`
is sized for. No OOB analogue exists for SuperLearner, so the OOB arms are
dropped.

### Causal forest (4 arms)

| id | `Y.hat` / `W.hat` | forest |
|---|---|---|
| `cf_dcf` | double-CF matrices, column `k` | fold-wise (**status quo**) |
| `cf_scf` | single-CF vectors | fold-wise |
| `cf_full_oob` | single-CF vectors | whole sample, OOB `tau` |
| `cf_default` | grf's own internal OOB | whole sample, OOB `tau` — plain `causal_forest(X, Y, W)` |

## Files

| file | role |
|---|---|
| `cf_models.R` | DGP wrapper, nuisance producers (over production's per-arm fits), stage-2 consumers, `run_all_crossfit_variants` |
| `cf_analysis.R` | array entry point, one replicate per index |
| `cf_testing.R` | verification checks — run before submitting anything |
| `cf_check.R` | finds missing runs, writes `jobscripts/failed_ids.txt`, and updates `-J` and the resource request in the rerun jobscript |
| `cf_metrics.R` | metric definitions (functions only, no side effects) |
| `cf_collect.R` | streams the per-run files through `cf_metrics.R` into `cf_metrics.RDS` |
| `cf_results.R` | figures |
| `confidence_intervals/cf_ci_analysis.R` | confidence-interval pilot, all 8 RF/CF arms — see below |
| `confidence_intervals/cf_ci_testing.R` | verification checks for the CI pilot (`full` adds the production-parity check) |
| `confidence_intervals/cf_ci_check.R` / `cf_ci_metrics.R` / `cf_ci_collect.R` | CI pilot's own check/metrics/collect, parallel to the files above |

## Half-sample bootstrap CI pilot

`cf_ci_analysis.R` adds confidence intervals to **all 8 non-SuperLearner arms**.
The 2 SuperLearner arms stay out of scope (not RF-based), which is why the pilot
calls `run_all_crossfit_variants(sl_lib = NULL)`.

The 8 split by stage-2 structure, and the bootstrap differs between them:

| arms | bootstrap | half sample | `tau_half` |
|---|---|---|---|
| `dcf`, `scf_scf`, `cf_dcf`, `cf_scf` | `rf_half_boot` / `cf_half_boot` | stratified by fold | refit per fold, predict the held-out fold |
| `scf_oob`, `oob_oob`, `cf_full_oob`, `cf_default` | `rf_oob_half_boot` / `cf_oob_half_boot` | unstratified `floor(n/2)` | one refit, OOB for in-half rows and `newdata` for the rest |

Nuisances are held fixed and sliced in every case — only the final-stage forest
is refit, which is `R/bootstrap_ci.R`'s existing design. That includes
`cf_default`, the one arm with no nuisance stage of its own: `cf_whole` now
returns the `Y.hat`/`W.hat` grf cross-fit internally so its bootstrap can hold
them fixed too, rather than letting grf re-cross-fit each half sample and put
nuisance variability into that arm's band and nobody else's.

### `tau_half` for an OOB arm, and why there are two of them

A whole-sample refit has no held-out fold to predict, so there are two defensible
ways to score the `n` units. `oob_bands()` produces **both**, off one set of
forest fits — a paired contrast, since masking costs nothing:

- **`half_boot`** (`hb_lb`/`hb_ub`) — in-half rows take the half forest's own OOB
  predictions, out-of-half rows take `newdata` predictions. Every unit gets a
  root in every draw, so `S_star` is a supremum over all `n` units. For in-half
  rows the functional matches the point estimate exactly: both are OOB
  predictions, at `n` vs `n/2`. Neither branch is contaminated by the row's own
  outcome. The one asymmetry is tree count — an OOB prediction averages the
  `1 - sample.fraction` share of trees (~1000 of 2000 at grf's defaults) while a
  `newdata` prediction averages all of them, which is second-order Monte Carlo
  noise beside the statistical variance at `n = 500`.
- **`half_boot_out`** (`hb_out_lb`/`hb_out_ub`) — only rows the half forest never
  saw, in-half cells masked to `NA`. Uniform in tree count, but each unit gets
  ~`B/2` roots and `S_star` becomes a supremum over ~250 rather than 500 units.
  The sup of `|N(0,1)|` over `m` units grows like `sqrt(2 log m)` — `3.32` vs
  `3.53` — so the prediction is a band ~6% narrower, a systematic downward bias
  in the critical value. Whether that survives the correlation between roots at
  this `n` is a question for the pilot, which is why both are computed rather
  than one being argued for on paper.

`simultaneous_band()` gained an `na.rm` argument for the masked variant. It
defaults to `FALSE`, so every pre-existing caller is bit-for-bit unaffected —
`cf_ci_testing.R` check 5 asserts that directly.

### grf's own variance, for free

The OOB arms carry a **third** interval, `grf_normal`. For a whole-sample forest
grf returns OOB variance estimates (bootstrap of little bags) alongside the
predictions at no extra compute, and `R/metrics.R`'s `normal_interval()` — the
same function `sample_size/confidence_intervals/binary/` uses for its `causal_forest_inbuilt`
method — turns those into a CI. `stage2_whole_rf`/`cf_whole` now return
`var_oob`, and `arm()` carries it; the 4 crossfit arms have none, and downstream
code keys off exactly that.

This is natural for an OOB arm and awkward for a crossfit one, whose `tau` is
stitched from `V` fold models predicting quantities grf's variance theory does
not cover — the exact mirror image of the bootstrap. So the OOB arms are the
first place in this study where both methods apply to the *same* arm, which
turns "does the bootstrap extend to OOB arms?" into the sharper "does an
expensive bootstrap buy anything over a free closed-form interval?"

**`grf_normal` is a pointwise interval; both bootstrap bands are simultaneous.**
Its near-zero `simultaneous_coverage` is the method working as designed, not a
defect — compare it on `marginal_coverage`.

Neither method reflects first-stage nuisance uncertainty: the bootstrap holds
nuisances fixed and grf's variance treats the pseudo-outcomes as known outcomes.
Whether that is second-order is the Neyman-orthogonality claim this study probes,
so it is at least apples-to-apples across the three.

### Scale and reproducibility

It is a **pilot, not the production run**: 3 scenarios (`1, 4, 8`) × 50 runs
= 150 replicates, `CI_boot = 200`, `CI_sf` fixed at 0.5 (grf's default
`sample.fraction`, no sweep).

Because the pilot now makes the same orchestrator call as `cf_analysis.R` minus
a SuperLearner block that sits strictly *after* every RF/CF arm in the RNG
stream, its point estimates are **bit-identical** to the production study's for
the same `(scenario, run)` — `cf_ci_testing.R` check 2 asserts that arm by arm
(under `full`, since it needs SuperLearner). This replaces the earlier
`run_crossfit_structured_arms()` orchestrator, which trimmed the out-of-scope
nuisance fits and consequently produced a *different, equally valid* draw of
`dcf`/`cf_dcf` that could only be checked by correlation. That function is gone;
its trimming bought nothing once the OOB arms came into scope, since they need
the very nuisances it was skipping.

`cf_half_boot` accepts single-crossfit vector nuisances (`cf_scf`) alongside its
original double-crossfit matrix nuisances (`cf_dcf`) — a shape-detection change,
backward compatible with its existing caller in `R/cate_models.R`.

Results land in `../results/crossfitting_ci/`, a wholly separate tree from
`../results/crossfitting/` — the production study's 400 replicates are
never read or touched. Per-run files drop the bootstrap `draws` matrices
before saving (only the bounds and `var_oob` are needed downstream), following
this folder's existing small-file convention.

Coverage from this method is already known (from `confidence_intervals/`) to
run below nominal — that's the pilot's actual research question, so
`cf_ci_testing.R` checks band well-formedness (finite, brackets `tau` for
most units, non-degenerate width) rather than gating on ~95% coverage.

`cf_ci_metrics.R` emits one row per **(arm, `ci_method`)**, so a replicate
produces `4 + 4 × 3 = 16` rows. That multi-row shape is the convention
`R/metrics.R`'s `compute_metrics` already documents for the CI studies.

```bash
Rscript crossfitting/confidence_intervals/cf_ci_testing.R              # structure + regression checks
Rscript crossfitting/confidence_intervals/cf_ci_testing.R full         # adds the "identical to production" check (needs SuperLearner)
Rscript crossfitting/confidence_intervals/cf_ci_analysis.R 1 10 2 1    # local smoke test: index 1, CI_boot=10

qsub crossfitting/confidence_intervals/jobscripts/cf_ci_1.sh           # the pilot itself (150 jobs)
Rscript crossfitting/confidence_intervals/cf_ci_check.R                # 150/150?
qsub crossfitting/confidence_intervals/jobscripts/cf_ci_collect.sh
```

### Sizing the CI pilot's array job

`cf_ci_1.sh`'s `#PBS -l` lines and trailing `CI_boot`/`workers`/`grf_threads`
args are **placeholders** set by hand (currently 2 cores, 5gb, 1h, `200 2 1`).
The `syrup` sweep meant to measure them didn't work for this study (see the root
README's "Resource profiling (removed)"). Size them by hand, and check the first
real subjobs with `qstat -fx <jobid> | grep resources_used`. What drives the
cost:

- **Bootstrap refits, linear in `CI_boot`.** Each bootstrap draw (`future_map()`
  in `R/bootstrap_ci.R`) is an independent refit: `V` forests on a half-sample for
  a crossfit arm, one for an OOB arm. So elapsed time scales roughly linearly in
  `CI_boot`, and a timed run at a small `CI_boot` (the smoke test above uses 10)
  extrapolates to the real 200.
- **A fixed cost.** The 4 OOB arms add little to the bootstrap total (`B x 1`
  refit each against a crossfit arm's `B x V`), and `half_boot_out` adds
  nothing since it reuses the same refits. But every replicate first fits all
  the RF/CF point estimates, including the 45 fold-pair double-crossfit
  nuisances, so extrapolate with an intercept, not just a per-draw rate.
- **Memory grows with `CI_boot` too.** `future_map()` accumulates all `CI_boot`
  result vectors before the draws matrix is assembled, so peak memory measured at
  a small `CI_boot` underestimates the real `CI_boot = 200`.

Change the trailing `workers` and the `ncpus` in `#PBS -l select` together, so
the two can't drift apart.

Nothing is forked: `R/utils.R` supplies `setup_rng_stream` and
`collate_predictions`, `sample_size/continuous/cts_dgms.R` supplies the DGP, and
`R/cate_models.R` supplies the per-arm outcome models (`t_learner_rf_split`,
`nuisance_rf`, `sl_split_fit`, `nuisance_sl`) and, via `R/sl_library.R`,
`pretest_superlearner`, the per-nuisance SuperLearner libraries and
`sl_fit_predict`. `cf_testing.R` section 1 checks that `oob_oob` reproduces
production's `dr_random_forest` bit-for-bit (production's `stage2_whole_rf` is
masked by this folder's, so the check sources `R/cate_models.R` into an
environment of its own), and section 1b that each arm's outcome forest never
sees the other arm's outcomes.

This folder was the model for the repo-wide `R/` refactor: it was already
sourcing shared code rather than copying it, at a time when the same CATE
estimators existed in seven files. The reference implementations it compares
against moved from `sample_size/continuous/cts_models.R` into `R/cate_models.R`, which is
now the only copy - `cts_models.R` is a thirteen-line profile shim.

## Status

**Re-run owed** - the arms were rebuilt around per-arm outcome models
(2026-09-27, above), so no existing result reflects the current design. Before
that, bug O changed the continuous DGM this study generates from
(`sample_size/continuous/README.md`), and the 2026-09-26 renumbering made its scenarios
1 / 4 / 6 / 8 (were 1 / 4 / 6 / 9; the pilot's 1 / 6 / 9 are now 1 / 4 / 8).
The comparison's conclusions were drawn on the old DGM. Archive both old trees
(`../results/crossfitting/`, `../results/crossfitting_ci/`) first with
`R/archive_old_results.R` (root `README.md`, Status, step 0) - `cf_results.R`
and `cf_ci_results.R` label the new numbers, so they would mislabel the old
`cf_metrics.RDS` / `cf_ci_metrics.RDS`.

## Running it

```bash
Rscript crossfitting/cf_testing.R              # structure + regression checks (fast)
Rscript crossfitting/cf_testing.R full         # adds the SuperLearner family

qsub crossfitting/jobscripts/cf_1.sh        # the study itself
Rscript crossfitting/cf_check.R             # 400/400?
qsub crossfitting/jobscripts/cf_collect.sh
Rscript crossfitting/cf_results.R
```

Results land in `../results/crossfitting/` (a sibling of the repo, as elsewhere).

## Sizing the array job

`cf_1.sh`'s `#PBS -l` lines and trailing `workers`/`grf_threads` args are set by
hand (currently 1 core, 2gb, 30 minutes, `1 1`). They were meant to be measured by
a `syrup` profiling sweep, but that didn't work for this study (see the root
README's "Resource profiling (removed)"). When changing them:

- keep `workers x grf_threads` within `ncpus`. Otherwise the grf threads
  oversubscribe the allocated cores
- change the trailing args and the `#PBS -l` line together, so they can't drift
  apart
- keep `ompthreads=ncpus` (see "Deviations" below)
- check the memory request against `qstat -fx <jobid> | grep resources_used` on
  the first real subjobs

## Deviations from the rest of the study, on purpose

- **Propensities are trimmed to `[0.05, 0.95]` in every arm**, including the RF ones
  (no longer a deviation: `R/cate_models.R` now trims its RF propensities too, where
  it used to trim only for SuperLearner). With `W ~ Bernoulli(0.5)` this is a no-op
  for the double-crossfit nuisances — `cf_testing.R` asserts it — but it stops any
  arm from producing exploding pseudo-outcomes and losing on a technicality.
- **`bias` is `estimate - truth`** (bug G, fixed repo-wide).
- **`stage2_crossfit_sl` uses the pretested library in both branches** (bug F, fixed repo-wide).
- **The per-run files carry no `data` and no nuisance matrices** — only `tau`,
  `tau_test` and timings per arm, plus the truth vectors. Replicates are
  reproducible from `run` via `setup_rng_stream`. This is why `cf_collect.sh` asks
  for 8gb where the other studies' collect jobs need 15gb.
- **`cf_metrics.R` holds definitions, not a script.** `cf_collect.R` applies them
  while streaming, so the large nested intermediate the other studies write to disk
  is never materialised and there is no separate `cf_metrics.sh` to submit.
- **`num.threads` is passed explicitly to grf** and `OMP_NUM_THREADS` is set before
  the `multisession` workers spawn. Elsewhere in the repo neither is set, which
  means `ompthreads=2` alongside `workers <- 2` lets each worker claim 2 threads
  against 2 allocated cores.
- **The select lines now request `ompthreads=ncpus`, like every other study's
  jobscripts.** PBS Pro sets `NCPUS` from `ompthreads`, and
  `parallelly::availableCores()` reads `NCPUS` — so `ompthreads` below `ncpus`
  understates the cores available to `plan(multisession, workers = ...)` and can
  make it fail outright once `workers` exceeds `ompthreads`. Per-worker thread
  control is unaffected: it's still done in R, via the `Sys.setenv(OMP_NUM_THREADS
  = grf_threads)` call before `plan()` and via grf's `num.threads`, both of which
  run before or independently of whatever PBS put in the environment.
  `cf_1.sh` requests `ompthreads=ncpus` accordingly.
