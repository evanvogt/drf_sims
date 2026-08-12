# Crossfitting comparison

Does the double crossfitting procedure used throughout this study earn its cost?

Everywhere else in `drf_sims`, two-stage estimators fit nuisances over all `C(V,2)`
fold *pairs*; `collate_predictions` (`R/utils.R`) assembles an `n x V` pseudo-outcome
matrix whose column `k` is untouched by fold `k`, and the second stage trains on
`po_matrix[, k]` before predicting fold `k`. At `V = 10` that is 45 nuisance fits
instead of 10 — 4.5x the cost — and it has never been benchmarked.

**The hypothesis.** Under ordinary crossfitting, a stage-2 training row `i` carries a
pseudo-outcome from a nuisance model trained on every fold except `i`'s, which
*includes* the held-out test fold. Re-randomising the stage-2 split cannot remove
that dependence. So if double crossfitting matters, `scf_scf` and `scf_scf_new`
should perform alike and both should differ from `dcf`. If they don't, the cheap
procedure is enough and the rest of the study can be sped up 4.5x.

## Design

Fixed at **n = 500, V = 10, scenarios 1 / 4 / 6 / 9, 500 runs** (2000 array jobs).
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

### DR-learner, random forest (8 arms)

| id | stage 1 (nuisances) | stage 2 (final model) |
|---|---|---|
| `dcf` | double CF over fold pairs | crossfit, same folds, column `k` (**status quo**) |
| `scf_scf` | single CF, leave-one-fold-out | crossfit, same folds |
| `scf_scf_new` | single CF | crossfit, **fresh independent** split |
| `scf_oob` | single CF | whole sample, **OOB** predictions |
| `oob_oob` | whole sample, **OOB**, **T-learner** | whole sample, **OOB** |
| `oob_oob_s` | whole sample, **OOB**, S-learner (`X.orig` workaround) | whole sample, **OOB** |
| `oob_oob_manual` | whole sample, **OOB**, S-learner (manual tree-loop) | whole sample, **OOB** |
| `scf_oob_t` | single CF, **T-learner** | whole sample, OOB (**control**) |

**Why `oob_oob` uses a T-learner, and what `oob_oob_s`/`oob_oob_manual` add.** The
nuisance model elsewhere is an S-learner — one forest on `cbind(W, X)`, predicted at
`W = 0` and `W = 1` (`cts_models.R:70`). grf's public API only returns OOB
predictions at each unit's *observed* covariate row, so a plain S-learner has no OOB
counterfactual; `oob_oob` sidesteps that with separate arm forests instead (a
treated unit's `Y1.hat` is OOB and its `Y0.hat` comes from a forest that never saw
it, so both arms are honest) — but that confounds "OOB vs crossfit" with "T vs S
learner", which is what `scf_oob_t` is for.
`oob_oob_s` and `oob_oob_manual` remove that confound directly, by getting a real
S-learner OOB counterfactual: `oob_oob_s` uses a maintainer-endorsed but unsupported
workaround ([grf-labs/grf#307](https://github.com/grf-labs/grf/issues/307)) —
point the fitted forest's `X.orig` at the perturbed covariates, clear its cached
`predictions`/`debiased.error`, and call `predict(forest)` so it re-reads `X.orig`
and returns OOB-for-row-`i` predictions at the perturbed point. `oob_oob_manual`
gets the same answer through grf's *documented* `get_tree()`/`get_leaf_node()` API
instead — loop over every tree, and for each of row `i`'s out-of-bag trees, average
the training `Y` in the leaf the counterfactual point falls into — and exists mainly
as a check on the shortcut (see `cf_testing.R`'s section 4). Both are kept alongside
`oob_oob` rather than replacing it.

### DR-learner, SuperLearner (3 arms)

`dcf`, `scf_scf`, `scf_scf_new`. No OOB analogue exists for SuperLearner, so the OOB
arms and the T-learner control are dropped.

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
| `cf_models.R` | DGP wrapper, nuisance producers, stage-2 consumers, `run_all_crossfit_variants` |
| `cf_analysis.R` | array entry point, one replicate per index |
| `cf_testing.R` | verification checks — run before submitting anything |
| `cf_profile.R` | timing / memory / CPU sweep over `(workers, grf_threads)`, instrumented with `syrup` |
| `cf_profile_summary.R` | turns the sweep into PBS directives and writes them into `cf_1.sh` |
| `cf_check.R` | finds missing runs, writes `jobscripts/failed_ids.txt`, and updates `-J` and the resource request in the rerun jobscript |
| `cf_metrics.R` | metric definitions (functions only, no side effects) |
| `cf_collect.R` | streams the per-run files through `cf_metrics.R` into `cf_metrics.RDS` |
| `cf_results.R` | figures |
| `confidence_intervals/cf_ci_analysis.R` | confidence-interval pilot, all 12 RF/CF arms — see below |
| `confidence_intervals/cf_ci_testing.R` | verification checks for the CI pilot (`full` adds the production-parity check) |
| `confidence_intervals/cf_ci_check.R` / `cf_ci_metrics.R` / `cf_ci_collect.R` | CI pilot's own check/metrics/collect, parallel to the files above |
| `confidence_intervals/cf_ci_profile.R` | timing / memory / CPU sweep over `(workers, grf_threads, CI_boot)` for the CI pilot, instrumented with `syrup` |
| `confidence_intervals/cf_ci_profile_summary.R` | turns the sweep into PBS directives (extrapolating to the pilot's real `CI_boot`) and writes them into `cf_ci_1.sh` |

## Half-sample bootstrap CI pilot

`cf_ci_analysis.R` adds confidence intervals to **all 12 non-SuperLearner arms**.
The 3 SuperLearner arms stay out of scope (not RF-based), which is why the pilot
calls `run_all_crossfit_variants(sl_lib = NULL)`.

The 12 split by stage-2 structure, and the bootstrap differs between them:

| arms | bootstrap | half sample | `tau_half` |
|---|---|---|---|
| `dcf`, `scf_scf`, `scf_scf_new`, `cf_dcf`, `cf_scf` | `rf_half_boot` / `cf_half_boot` | stratified by fold | refit per fold, predict the held-out fold |
| `scf_oob`, `scf_oob_t`, `oob_oob`, `oob_oob_s`, `oob_oob_manual`, `cf_full_oob`, `cf_default` | `rf_oob_half_boot` / `cf_oob_half_boot` | unstratified `floor(n/2)` | one refit, OOB for in-half rows and `newdata` for the rest |

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
same function `confidence_intervals/binary/` uses for its `causal_forest_inbuilt`
method — turns those into a CI. `stage2_whole_rf`/`cf_whole` now return
`var_oob`, and `arm()` carries it; the 5 crossfit arms have none, and downstream
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

It is a **pilot, not the production run**: 3 scenarios (`1, 6, 9`) × 50 runs
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
`../results/crossfitting/` — the production study's 2000 replicates are
never read or touched. Per-run files drop the bootstrap `draws` matrices
before saving (only the bounds and `var_oob` are needed downstream), following
this folder's existing small-file convention.

Coverage from this method is already known (from `confidence_intervals/`) to
run below nominal — that's the pilot's actual research question, so
`cf_ci_testing.R` checks band well-formedness (finite, brackets `tau` for
most units, non-degenerate width) rather than gating on ~95% coverage.

`cf_ci_metrics.R` emits one row per **(arm, `ci_method`)**, so a replicate
produces `5 + 7 × 3 = 26` rows. That multi-row shape is the convention
`R/metrics.R`'s `compute_metrics` already documents for the CI studies.

```bash
Rscript crossfitting/confidence_intervals/cf_ci_testing.R              # structure + regression checks
Rscript crossfitting/confidence_intervals/cf_ci_testing.R full         # adds the "identical to production" check (needs SuperLearner)
Rscript crossfitting/confidence_intervals/cf_ci_analysis.R 1 10 2 1    # local smoke test: index 1, CI_boot=10

Rscript crossfitting/confidence_intervals/cf_ci_profile.R 1            # smoke-test the CI profiler locally
qsub crossfitting/confidence_intervals/jobscripts/cf_ci_profile.sh     # 48 profiling jobs
Rscript crossfitting/confidence_intervals/cf_ci_profile_summary.R      # writes measured directives into cf_ci_1.sh

qsub crossfitting/confidence_intervals/jobscripts/cf_ci_1.sh           # the pilot itself (150 jobs)
Rscript crossfitting/confidence_intervals/cf_ci_check.R                # 150/150?
qsub crossfitting/confidence_intervals/jobscripts/cf_ci_collect.sh
```

### Sizing the CI pilot's array job

`cf_ci_1.sh` shipped with **placeholder** `#PBS -l` lines (a hand-sizing step the
pilot's original scope deferred rather than skipped). `cf_ci_profile.R` /
`cf_ci_profile_summary.R` settle it the same way `cf_profile.R` / `cf_profile_summary.R`
do for `cf_1.sh` — see "Sizing the array job" below for the shared mechanics
(`syrup`, process-tree filtering, the peak-memory upper bound) — with one addition
specific to the bootstrap: the sweep profiles `CI_boot` at **20 and 60**, not at the
pilot's real 200. Each bootstrap draw (`future_map()` in `R/bootstrap_ci.R`) is an
independent refit — `V` forests on a half-sample for a crossfit arm, one for an OOB
arm — so elapsed time scales ~linearly in `CI_boot`. Profiling directly at
`CI_boot = 200` across the sweep would cost as much as the production array it
exists to size. `cf_ci_profile_summary.R` fits `elapsed ~ CI_boot` per
`(workers, grf_threads)`, pooled over scenario/run so the R² is a real diagnostic,
and warns if any configuration's fit falls below R² = 0.9 rather than trusting a bad
extrapolation silently. It also reports a per-arm bootstrap cost breakdown over all
12 arms, extrapolated to `CI_boot = 200` — the "why is it slow" answer, alongside
the resource sizing.

Expect the fixed cost to dominate more than it used to. The 7 OOB arms add only
~15% to the bootstrap total (`B x 1` refit each against a crossfit arm's `B x V`),
and `half_boot_out` adds nothing at all since it reuses the same refits — but the
pilot now pays for the 4 nuisance objects the old trimmed orchestrator skipped, and
`nuisance_oob_rf_manual` is a pure-R double tree loop over `num.trees x 2`
counterfactual passes. That is likely the largest single non-bootstrap cost in a
replicate, and it inflates the intercept of the `elapsed ~ CI_boot` fit, so the
R² < 0.9 warning is worth reading rather than skimming.

> **Known-wrong, deferred: the memory figure.** `cf_ci_profile_summary.R`
> extrapolates only *elapsed time* in `CI_boot`; for memory it applies a flat
> `mem_factor = 1.5` to the peak RSS observed at `CI_boot = 20/60`. The
> justification for that (still written in `cf_ci_profile.R`'s header) is that peak
> memory is governed by how many draws run concurrently across `workers` rather
> than by `CI_boot` itself. **That reasoning does not hold** — `future_map()`
> accumulates all `CI_boot` result vectors before the draws matrix is assembled, so
> memory does grow with `CI_boot`, and the `mem=` written into `cf_ci_1.sh` is an
> underestimate for the real `CI_boot = 200`. Not fixed in the change that added
> the OOB arms. Until it is, cross-check the request against
> `qstat -fx <jobid> | grep resources_used` on the first real subjobs rather than
> trusting the written figure.

The `Rscript` line it writes carries `CI_boot`, `workers` and `grf_threads`
together, fixing the mismatch the placeholder version had (`workers=2` on the
Rscript line vs. `ncpus=1` in `#PBS -l select` — the same kind of drift
`cf_profile_summary.R`'s header comment warns about for `cf_1.sh`).

Nothing is forked: `R/utils.R` supplies `setup_rng_stream` and
`collate_predictions`, `continuous/cts_dgms.R` supplies the DGP, and
`R/cate_models.R` supplies `pretest_superlearner` plus the reference
implementation the regression check in `cf_testing.R` compares against — though
that last one no longer holds, which is what breaks section 1 of `cf_testing.R`
(see "Known issue" below).

This folder was the model for the repo-wide `R/` refactor: it was already
sourcing shared code rather than copying it, at a time when the same CATE
estimators existed in seven files. The reference implementations it compares
against moved from `continuous/cts_models.R` into `R/cate_models.R`, which is
now the only copy - `cts_models.R` is a thirteen-line profile shim.

## Known issue: `cf_testing.R` section 1 aborts

**`cf_testing.R` currently fails on its first check** and takes the whole script
down with it, so none of sections 2–5 run:

```
=== 1. regression check: dcf against cts_models.R ===
Error in validate_num_threads(num.threads) :
  'list' object cannot be coerced to type 'double'
Calls: nuisance_rf -> regression_forest -> validate_num_threads
```

This is a **stale test, not a broken estimator**. Section 1 was left behind by
the repo-wide crossfitting simplification (`b1bcb16`), which changed the shared
reference implementations out from under it:

| `cf_testing.R` line | calls | problem |
|---|---|---|
| 64 | `nuisance_rf(X, Y, W, fold_indices, fold_pairs)` | signature is now `(X, Y, W, ipw, num.threads)`, so `fold_indices` binds to `ipw` and the `fold_pairs` **list** binds to `num.threads` — hence the error |
| 64 | `nz_old$po_matrix`, `nz_old$W.hat_matrix` | `nuisance_rf` returns vectors now; there are no `*_matrix` fields |
| 92 | `stage_2_rf(...)` | no longer exists in `R/cate_models.R` — replaced by `stage2_whole_rf` |

The check's premise is also gone: it asserted that this folder's
`nuisance_double_rf` reproduces the *production* implementation, but production
is no longer double crossfitting — `dcf` is now only an arm of this study, with
no shared counterpart to regress against. Fixing it means either pointing
section 1 at a pinned copy of the old double-crossfit code, or deleting it and
relying on sections 2–5.

Sections 2–5 (structure checks, the OOB S-learner equivalence check against the
manual `get_tree()` reimplementation, and the SuperLearner family under `full`)
are unaffected in principle but currently unreachable, so **this folder has no
working verification** until section 1 is fixed. Discovered while migrating
`competing_risk/` onto the shared strategy; confirmed to reproduce on pristine
`HEAD`, so it predates that work.

## Running it

```bash
Rscript crossfitting/cf_testing.R              # structure + regression checks (fast) - SEE ABOVE, currently aborts
Rscript crossfitting/cf_testing.R full         # adds the SuperLearner family

Rscript crossfitting/cf_profile.R 1         # smoke-test the profiler locally
qsub crossfitting/jobscripts/cf_profile.sh  # 36 profiling jobs
Rscript crossfitting/cf_profile_summary.R   # writes measured directives into cf_1.sh

qsub crossfitting/jobscripts/cf_1.sh        # the study itself
Rscript crossfitting/cf_check.R             # 2000/2000?
qsub crossfitting/jobscripts/cf_collect.sh
Rscript crossfitting/cf_results.R
```

Results land in `../results/crossfitting/` (a sibling of the repo, as elsewhere).

## Sizing the array job

`cf_1.sh` ships with **placeholder** `#PBS -l` lines. `cf_profile.R` measures what
one replicate actually costs across a `(workers, grf_threads)` sweep, and
`cf_profile_summary.R` turns that into directives and writes them in.

Profiling is instrumented with [`syrup`](https://simonpcouch.github.io/syrup/),
which runs a separate R session that samples `ps` every second and reports every R
process — so it sees the `multisession` workers, which the parent's own `gc()`
cannot. Nothing depends on the scheduler, so `Rscript crossfitting/cf_profile.R 1`
runs on a laptop as a smoke test of the harness (but not for the final numbers: a
machine with fewer cores than `workers x grf_threads` self-oversubscribes, and the
summary flags any cell where that happened).

Two things to know about the memory figure:

- `syrup` reports every R session `ps` can see, which on a shared node can include
  another of your own jobs. `cf_profile.R` filters to `pid == Sys.getpid() |
  ppid == Sys.getpid()` — `multisession` workers are PSOCK children — and warns if
  the tracked process count is not `workers + 1`.
- Peak memory is the summed RSS across the tree at the worst snapshot. That
  overcounts shared library pages, so it is an upper bound — the safe direction for
  a `mem=` request. Cross-check it once against `qstat -fx <jobid> | grep
  resources_used` on the first real subjobs; PBS's cgroup figure counts shared
  pages once and should sit below the request.

`syrup`'s `pct_cpu` also settles the oversubscription question outright rather than
inferring it from wall-clock differences: the summary prints observed CPU% against
`workers x grf_threads x 100` per configuration. This matters because
`continuous/jobscripts/cts_1.sh` asks for `ompthreads=2` alongside
`cts_analysis.R`'s `workers <- 2`, and `num.threads` is set nowhere in the repo — so
each worker may be claiming two threads against two allocated cores.

`syrup` must be installed in the `sim-env` conda environment; `cf_profile.R` fails
in the first second with a clear message if it is not.

## Deviations from the rest of the study, on purpose

- **Propensities are trimmed to `[0.05, 0.95]` in every arm**, including the RF ones.
  `cts_models.R` only trims for SuperLearner. With `W ~ Bernoulli(0.5)` this is a
  no-op for the crossfit arms — `cf_testing.R` asserts it — but it stops the in-sample
  nuisances from producing exploding pseudo-outcomes and losing on a technicality.
- **`bias` is `estimate - truth`**, the usual convention. This used to be a
  deviation: `cts_metrics.R` computed `bias` as `mean(true - est)` while its
  `ate_bias` used the opposite sign. That was bug G, and the whole repo now uses
  `est - true` (`R/metrics.R`), so these numbers ARE comparable with the other
  studies once their metrics are regenerated.
- **`stage2_crossfit_sl` uses the pretested library in both branches.** This was
  bug F, and it was worse than it looked: every caller of the shared `stage_2_sl`
  passed a matrix, so the pretested library had never been used in stage 2 in any
  study and the correct vector branch was dead code. Fixed repo-wide
  (`PRETEST_STAGE2` in `R/cate_models.R`); this folder was right all along.
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
  against 2 allocated cores. The profiling sweep measures whether that matters.
- **The select lines now request `ompthreads=ncpus`, like every other study's
  jobscripts.** PBS Pro sets `NCPUS` from `ompthreads`, and
  `parallelly::availableCores()` reads `NCPUS` — so `ompthreads` below `ncpus`
  understates the cores available to `plan(multisession, workers = ...)` and can
  make it fail outright once `workers` exceeds `ompthreads`. Per-worker thread
  control is unaffected: it's still done in R, via the `Sys.setenv(OMP_NUM_THREADS
  = grf_threads)` call before `plan()` and via grf's `num.threads`, both of which
  run before or independently of whatever PBS put in the environment.
  `cf_profile_summary.R` writes `ompthreads=ncpus` into `cf_1.sh` accordingly.
