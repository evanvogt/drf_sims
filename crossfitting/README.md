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
| `cf_check.R` | finds missing runs, writes `jobscripts/failed_ids.txt` |
| `cf_metrics.R` | metric definitions (functions only, no side effects) |
| `cf_collect.R` | streams the per-run files through `cf_metrics.R` into `cf_metrics.RDS` |
| `cf_results.R` | figures |
| `confidence_intervals/cf_ci_analysis.R` | half-sample bootstrap CI pilot — see below |
| `confidence_intervals/cf_ci_testing.R` | verification checks for the CI pilot |
| `confidence_intervals/cf_ci_check.R` / `cf_ci_metrics.R` / `cf_ci_collect.R` | CI pilot's own check/metrics/collect, parallel to the files above |
| `confidence_intervals/cf_ci_profile.R` | timing / memory / CPU sweep over `(workers, grf_threads, CI_boot)` for the CI pilot, instrumented with `syrup` |
| `confidence_intervals/cf_ci_profile_summary.R` | turns the sweep into PBS directives (extrapolating to the pilot's real `CI_boot`) and writes them into `cf_ci_1.sh` |

## Half-sample bootstrap CI pilot

`cf_ci_analysis.R` adds half-sample bootstrap confidence intervals
(`R/bootstrap_ci.R`'s `cf_half_boot`/`rf_half_boot`, the same machinery
`confidence_intervals/` uses) to the 5 arms whose stage 2 is a genuine
per-fold crossfit: `dcf`, `scf_scf`, `scf_scf_new` (family `dr_rf`) and
`cf_dcf`, `cf_scf` (family `causal_forest`). Out of scope: the 6 whole-sample
/ OOB arms (no per-fold structure for the bootstrap to refit against) and the
3 SuperLearner arms (not RF-based).

It is a **pilot, not the production run**: 3 scenarios (`1, 6, 9`) × 50 runs
= 150 replicates, `CI_boot = 200`, `CI_sf` fixed at 0.5 (grf's default
`sample.fraction`, no sweep). Nuisances are computed once per replicate via a
new trimmed orchestrator, `run_crossfit_structured_arms()` in `cf_models.R`
(mechanically extracted from `run_all_crossfit_variants`). Its `nz_double`/
`nz_single` nuisance objects are bit-identical to the production study's
under the same `setup_rng_stream(run)` seed, since nothing precedes them in
the RNG stream in either orchestrator — but its `dcf`/`cf_dcf` stage-2
estimates are **not** bit-identical to the saved production arms: skipping
the 4 out-of-scope nuisance fits that `run_all_crossfit_variants` performs
between `nz_single` and its own stage-2 calls shifts every later
`future_map()` forest fit to a different position in the RNG stream. That's
by design (replaying the skipped fits just to stay in lockstep would defeat
the point of trimming them); `dcf`/`cf_dcf` here are a different, equally
valid draw of the same estimator for the same `(scenario, run)`, not a
reproduction of the saved production values — see `cf_ci_testing.R` check 1,
which verifies nuisance-level identity and stage-2 estimator-level agreement
(high correlation, not exact equality) instead. Only the final-stage forest
is refit per half-sample-per-fold, matching `R/bootstrap_ci.R`'s existing
design. `cf_half_boot` was generalized to accept single-crossfit vector
nuisances (`cf_scf`) alongside its original double-crossfit matrix nuisances
(`cf_dcf`) — a shape-detection change, backward compatible with its existing
caller in `R/cate_models.R`.

Results land in `../results/crossfitting_ci/`, a wholly separate tree from
`../results/crossfitting/` — the production study's 2000 replicates are
never read or touched. Per-run files drop the bootstrap `draws` matrices
before saving (only `hb_lb`/`hb_ub` are needed downstream), following this
folder's existing small-file convention.

Coverage from this method is already known (from `confidence_intervals/`) to
run below nominal — that's the pilot's actual research question, so
`cf_ci_testing.R` checks band well-formedness (finite, brackets `tau` for
most units, non-degenerate width) rather than gating on ~95% coverage.

```bash
Rscript crossfitting/confidence_intervals/cf_ci_testing.R              # structure + regression checks
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
pilot's real 200. Each bootstrap draw (`future_map()` in `rf_half_boot`/
`cf_half_boot`, `R/bootstrap_ci.R`) is an independent refit of `V` forests on a
half-sample, and the `n x CI_boot` draws matrix is only assembled after
`future_map()` returns — so elapsed time scales ~linearly in `CI_boot` while peak
memory (governed by how many draws run concurrently across `workers`, not by
`CI_boot` itself) stays roughly flat. Profiling directly at `CI_boot = 200` across
the sweep would cost as much as the production array it exists to size.
`cf_ci_profile_summary.R` fits `elapsed ~ CI_boot` per `(workers, grf_threads)`,
pooled over scenario/run so the R² is a real diagnostic, and warns if any
configuration's fit falls below R² = 0.9 rather than trusting a bad extrapolation
silently. It also reports a per-arm bootstrap cost breakdown (which of `dcf`,
`scf_scf`, `scf_scf_new`, `cf_dcf`, `cf_scf` dominates, extrapolated to
`CI_boot = 200`) — the "why is it slow" answer, alongside the resource sizing.

The `Rscript` line it writes carries `CI_boot`, `workers` and `grf_threads`
together, fixing the mismatch the placeholder version had (`workers=2` on the
Rscript line vs. `ncpus=1` in `#PBS -l select` — the same kind of drift
`cf_profile_summary.R`'s header comment warns about for `cf_1.sh`).

Nothing is forked: `R/utils.R` supplies `setup_rng_stream` and
`collate_predictions`, `continuous/cts_dgms.R` supplies the DGP, and
`R/cate_models.R` supplies `pretest_superlearner` plus the reference
implementation the regression check in `cf_testing.R` compares against.

This folder was the model for the repo-wide `R/` refactor: it was already
sourcing shared code rather than copying it, at a time when the same CATE
estimators existed in seven files. The reference implementations it compares
against moved from `continuous/cts_models.R` into `R/cate_models.R`, which is
now the only copy - `cts_models.R` is a thirteen-line profile shim.

## Running it

```bash
Rscript crossfitting/cf_testing.R              # structure + regression checks (fast)
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
