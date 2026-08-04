# Crossfitting comparison

Does the double crossfitting procedure used throughout this study earn its cost?

Everywhere else in `drf_sims`, two-stage estimators fit nuisances over all `C(V,2)`
fold *pairs*; `collate_predictions` (`utils.R:16`) assembles an `n x V` pseudo-outcome
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
on an independent **test** sample of 2000 drawn from the same DGP. The test set is
the only like-for-like comparison of the fitted final models, and the train-to-test
gap is the overfitting measure for the whole-sample arms.

### DR-learner, random forest (8 arms)

| id | stage 1 (nuisances) | stage 2 (final model) |
|---|---|---|
| `dcf` | double CF over fold pairs | crossfit, same folds, column `k` (**status quo**) |
| `scf_scf` | single CF, leave-one-fold-out | crossfit, same folds |
| `scf_scf_new` | single CF | crossfit, **fresh independent** split |
| `scf_full` | single CF | whole sample, in-sample predictions |
| `scf_oob` | single CF | whole sample, **OOB** predictions |
| `oob_oob` | whole sample, **OOB** | whole sample, **OOB** |
| `naive` | whole sample, in-sample | whole sample, in-sample |
| `scf_oob_t` | single CF, **T-learner** | whole sample, OOB (**control**) |

`scf_full` / `scf_oob` are two views of the same fitted forest, as are
`cf_full_oob` / `cf_naive`, so those contrasts isolate the prediction type rather
than also picking up forest-to-forest randomness.

**Why `oob_oob` uses a T-learner.** The nuisance model elsewhere is an S-learner —
one forest on `cbind(W, X)`, predicted at `W = 0` and `W = 1` (`cts_models.R:70`).
grf only returns OOB predictions at each unit's *observed* covariate row, so an
S-learner has no OOB counterfactual and `oob_oob` cannot be built from one. It uses
separate arm forests instead: a treated unit's `Y1.hat` is OOB and its `Y0.hat`
comes from a forest that never saw it, so both arms are honest. That confounds
"OOB vs crossfit" with "T vs S learner", which is what `scf_oob_t` is for — it
differs from `scf_oob` only in learner structure and from `oob_oob` only in
splitting, making both effects identifiable.

### DR-learner, SuperLearner (5 arms)

`dcf`, `scf_scf`, `scf_scf_new`, `scf_full`, `naive`. No OOB analogue exists for
SuperLearner, so the OOB arms and the T-learner control are dropped.

### Causal forest (5 arms)

| id | `Y.hat` / `W.hat` | forest |
|---|---|---|
| `cf_dcf` | double-CF matrices, column `k` | fold-wise (**status quo**) |
| `cf_scf` | single-CF vectors | fold-wise |
| `cf_full_oob` | single-CF vectors | whole sample, OOB `tau` |
| `cf_default` | grf's own internal OOB | whole sample, OOB `tau` — plain `causal_forest(X, Y, W)` |
| `cf_naive` | single-CF vectors | whole sample, in-sample `tau` |

## Files

| file | role |
|---|---|
| `cf_models.R` | DGP wrapper, nuisance producers, stage-2 consumers, `run_all_crossfit_variants` |
| `cf_analysis.R` | array entry point, one replicate per index |
| `cf_test.R` | verification checks — run before submitting anything |
| `cf_profile.R` | timing / memory sweep over `(workers, grf_threads)` |
| `cf_profile_summary.R` | turns the sweep into PBS directives and writes them into `cf_1.sh` |
| `cf_check.R` | finds missing runs, writes `jobscripts/failed_ids.txt` |
| `cf_metrics.R` | metric definitions (functions only, no side effects) |
| `cf_collect.R` | streams the per-run files through `cf_metrics.R` into `cf_metrics.RDS` |
| `cf_results.R` | figures |

Nothing is forked: `utils.R` supplies `setup_rng_stream` and `collate_predictions`,
`continuous/cts_dgms.R` supplies the DGP, and `continuous/cts_models.R` supplies
`pretest_superlearner` plus the reference implementation the regression check in
`cf_test.R` compares against.

## Running it

```bash
Rscript crossfitting/cf_test.R              # structure + regression checks (fast)
Rscript crossfitting/cf_test.R full         # adds the SuperLearner family

qsub crossfitting/jobscripts/cf_profile.sh  # 36 profiling jobs
Rscript crossfitting/cf_profile_summary.R   # writes measured directives into cf_1.sh

qsub crossfitting/jobscripts/cf_1.sh        # the study itself
Rscript crossfitting/cf_check.R             # 2000/2000?
qsub crossfitting/jobscripts/cf_collect.sh
Rscript crossfitting/cf_results.R
```

Results land in `../results/crossfitting/` (a sibling of the repo, as elsewhere).

## Deviations from the rest of the study, on purpose

- **Propensities are trimmed to `[0.05, 0.95]` in every arm**, including the RF ones.
  `cts_models.R` only trims for SuperLearner. With `W ~ Bernoulli(0.5)` this is a
  no-op for the crossfit arms — `cf_test.R` asserts it — but it stops the in-sample
  nuisances from producing exploding pseudo-outcomes and losing on a technicality.
- **`bias` is `estimate - truth`**, the usual convention. `cts_metrics.R:37` uses
  `mean(true - est)` for `bias` while its `ate_bias` at line 43 uses the opposite
  sign, so these numbers are not sign-comparable with the main continuous study.
- **`stage2_crossfit_sl` uses the pretested library in both branches.**
  `cts_models.R:387` computes `po_lib` and then passes the untested `sl_lib` in the
  matrix branch; line 384 uses `po_lib` correctly in the vector branch.
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
