# Competing risks

CATE estimation when event 1 can be pre-empted by a competing event 2. This
is the study the rest of the repo builds towards, and the only one whose
estimators are genuinely different rather than a configuration of the shared
ones: 19 arms across IPCW, causal survival forest, pseudo-value and RSF
families, on restricted-mean time scales.

Scenarios 1–7 × censoring on/off × ρ ∈ {0.5 primary, 0 sensitivity}, n = 500,
500 runs — **14,000 array jobs**: rows 1–7000 are ρ = 0 (`surv_1.sh`),
7001–14000 ρ = 0.5 (`surv_2.sh`). The full design — DGM, estimands, every arm,
the pseudo-value comparison and the performance measures — is in
[`ADEMP.md`](ADEMP.md).

On the cluster, from `competing_risk/jobscripts/` (run `Rscript make_log_dirs.R`
from the repo root first, for `logs_1/` and `logs_2/`):

```bash
qsub surv_1.sh               # 1-7000, rho = 0 (sensitivity)
qsub surv_2.sh               # 7001-14000, rho = 0.5 (primary)
Rscript ../../study.R check competing_risk/surv_config.R      # writes failed_ids.txt, points surv_rerun.sh at it
qsub surv_rerun.sh           # only if the check found failures
qsub surv_collect.sh
qsub surv_metrics.sh
```

Results go to `../results/competing_risk/rho_<ρ>/scenario_<k>/<n>/censor_<TRUE|FALSE>/`.

Local smoke test, from `competing_risk/`: `Rscript surv_analysis.R 7003`
(scenario 3, censoring on, run 1, ρ = 0.5, the primary), then
`Rscript surv_analysis.R 3` (the same run at ρ = 0).

## Without the queue

Any subset of the grid can be produced inside one interactive RStudio session
instead of the array: request an RStudio session with 8 cores and 64gb, set
`ids` in `surv_run.R` (e.g. `grid_indices(study, rho = 0.5, scenario = 7)`),
then

```r
source(here::here("competing_risk", "surv_run.R"))
```

With `ids <- NULL` it targets every row without a results file, which for a
fresh grid is all 14,000, far more than a session gets through. Same
`surv_analysis.R`, same results, 4 rows at a time (2 cores each, matching
`surv_1.sh`'s `ncpus=2`). See `surv_run.R`'s header for why 4 and not 8.

## Implementation notes

- **Threads.** The RSF fold workers set `options(rf.cores = 1)`. The OOB RSF
  fits run in the main process and follow `OMP_NUM_THREADS`, which
  `surv_analysis.R` sets to `grf_threads` (default 1) before starting its
  workers. The grf forests take `num.threads = grf_threads` directly.
- **The RSF arms are not bitwise reproducible**: they differ run to run on
  unchanged code, so leave them out of exact equivalence checks.
- `surv_analysis.R` sources `R/cate_models.R` before `surv_models.R`, so the
  study's own definitions win where they still exist. That is what supplies
  `t_learner_rf`, `stage2_whole_rf`, `stage_2_sl`, `dr_pseudo`, and (via
  `R/sl_library.R`) `pretest_superlearner`, `sl_fit_predict` and `sl_libraries`.
- `surv_metrics.R` does not use `compute_metrics()`: its results nest as
  framework × target rather than one entry per model. It does use
  `cate_metrics()`, so the bias convention matches the rest of the repo, and
  `framework_truth_map` to score each arm against its own truth.

## Known issues

**Two unguarded failure routes.** Both were diagnosed on an earlier DGM, where
they failed 19 runs deterministically; under the current parameters neither
should recur, but neither code path is guarded, so check `failed_ids.txt` after
a rerun. `surv_analysis.R`'s error handler prints the call stack to stdout, and
the jobscripts use `#PBS -j oe`, so keep the `.o` logs.

- **Whole-sample pseudo-values NA.** When nobody is observed past the horizon,
  `pseudoyl()` returns NA for the max-time individual (the `pseudo:::ci.omit`
  bug below) and `pseudo_all()` has no fallback, so the run dies at
  `pseudo_cf_whole_oob` with `The vector of observations (W, Y, Z or D)
  contains at least one NA.` The estimand is not identified past the last
  observed time, so any fix needs to count what it drops or imputes.
- **SL.mean fallback dropped.** On a constant treated-arm RMTL2 cell,
  `pretest_superlearner()` falls back to `SL.mean`, and the live
  `SuperLearner()` then drops that too: `All algorithms dropped from library`,
  in `sl_t_split`.

**SuperLearner DR-learner, split pseudo-obs (disabled).** `results$sl_dr_split`
is commented out in `all_cate_surv_models()`; `nuisance_pseudo_sl_split()` and
`stage_2_sl_vec()` are still defined. The cause is in the `pseudo` package:
`pseudoyl()`'s `pseudo:::ci.omit()` computes leave-one-out cumulative incidence
over a shared risk-set matrix, and for the individual with the **maximum
observed time** the leave-one-out risk set can reach zero at the tail, giving
`0/0 = NaN` in both the KM term (`km`) and the cause-specific ratio (`cumi`).
The package forward-fills `km` but not `cumi`, so `cumsum(cumi * km)` is NA for
that individual. `compute_split_pseudoyl()`/`compute_split_pseudomean()` hit
this whenever a validation observation's `Y` exceeds the max `Y` in that fold's
KM set — reliably, at the production n = 500 and 10 folds. `pseudo_crossfit()`
guards the same NA with a `pseudo_whole` fallback; the split construction has
none. Fixing it needs a design decision: fall back (to the training-fold
pseudo-value or `pseudo_whole`) or drop the observation from that fold's
stage-2 fit. The T-learner split arm has the same latent risk and is guarded
(`keep <- !is.na(y)`).

To work on either:

```bash
Rscript R/regression_check.R baseline competing_risk
```

The diagnostic scripts (`surv_failed_diagnose.R`, `surv_dr_split_na_diagnose.R`)
are at git tag `diagnostics-2026-10`; restore one with
`git show diagnostics-2026-10:competing_risk/surv_failed_diagnose.R > surv_failed_diagnose.R`
and pass it `--out`, since its stages 2 and 3 overwrite their CSVs in place.

## Files

`surv_config.R` (grid), `surv_dgm.R`, `surv_models.R`, `surv_analysis.R`,
`surv_run.R` (the no-queue RStudio-session alternative, see "Without the
queue" above), `surv_metrics.R`; check and collect are the root `study.R`.
`surv_dgm_check.R` produces ADEMP's "What the DGM implies" numbers (event mix,
population truths, X3 separation) for each ρ.
`single_event/` is a separate, smaller study: event 1 with event 2 removed,
reusing this study's arms (see its README).
`snr/` is a small signal-to-noise check: scenario 7 with the X3 effect scaled
(`generate_surv_data(hte_mult = )`) or n raised, the pseudo-value CF against the
production DR-learners and their semi-oracle twins (see its README).
