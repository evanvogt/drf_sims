# Model evaluation — does cheap proxy scoring pick the right CATE learner?

Given 9 candidate CATE learners, can cheap proxy losses built from separately
fitted nuisance estimators rank them the way the (normally unobservable) true
PEHE would? A model-*selection* study: the object of interest is how each
proxy ranks the candidates, not any one candidate's accuracy.

Scenarios 1–4 of `continuous_corr_0.5` × n 250, 500, 1000 × 100 runs
(**1200 array jobs**), plus a strategies pass and an 80:20 split arm. The full
design — candidates, nuisance arms, score families and performance measures —
is in [`ADEMP.md`](ADEMP.md).

## Files

| file | role |
|---|---|
| `me_config.R` | the parameter grid, results paths, candidate list and arm names — **the** definition |
| `me_dgms.R` | names this study's slice of `R/dgm_scenarios.R` (`set = "continuous_corr_0.5"`) |
| `me_utils.R` | design-matrix prep and fold-splitting (DGM-agnostic) |
| `me_models.R` | the 9 candidate CATE-learner configurations and their fitting logic |
| `me_nuisance.R` | the two independent nuisance-evaluation pipelines |
| `me_analysis.R` | array entry point; one row of the grid per index |
| `me_strategies.R` | second pass over completed runs — fits all three arms (`whole`, `cv_shared`, `holdout`) without refitting the candidates, writes to `model_evaluation_strategies` |
| `me_strategies_verify.R` | proves that pass carried the candidates, data, truth and `fold_info` through bit-identically and that all three arms are complete (no NA prediction, finite phi with π trimmed), and tallies the runs with a propensity of exactly 0 or 1 |
| `me_split.R` | the 80:20 arm — the only script that refits the candidates |
| `me_cs_learners.R` | the case-study learners (DR-RF, CF, DR-SL) scored by the DR score: reads the strategies tree and the sample-size study's matching runs, fits nothing, writes `me_cs_metrics.RDS` (`ADEMP.md`, "The case-study learners") |
| `me_check.R` | finds missing runs, writes `jobscripts/failed_ids.txt`, and updates `-J` and the resource request in the rerun jobscript. Takes a tree: `main` (default) / `strategies` / `split` |
| `me_collect.R` | gathers per-run files into `<prefix>_all.RDS`. Same tree argument |
| `me_metrics.R` | computes `<prefix>_metrics.RDS` (reuses `R/metrics.R::compute_metrics()`). Same tree argument |
| `me_selection.R` | per-run selection, ranks, regret and references, shared by the report and the thesis figures/tables |
| `me_strategies_run.R` | runs the strategies pass inside one interactive session instead of the array |
| `me_results.qmd` | the results report |
| `me_infl_investigation.R`, `.qmd` | the IF score's population check with oracle nuisances |
| `me_testing.R` | verification checks — run before submitting anything |

Scores are a pure post-hoc function of `<prefix>_all.RDS`: adding a score
family means re-running `me_metrics.R` and nothing else.

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

The nuisance-arm pass. It reads the completed main tree, so run it after the
main study:

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

The case-study learners, after the strategies tree's metrics exist. It is one
job, not an array: nothing is fit, it only reads each strategies run and the
sample-size study's matching continuous run (ρ = 0.5), and it exits non-zero
if a dataset differs or SL1's scores do not reproduce the strategies tree's:

```bash
qsub model_evaluation/jobscripts/me_cs_learners.sh                 # writes me_cs_metrics.RDS
Rscript results_processing/thesis_tables/me_cs_tables.R            # me_cs_learners.tex, me_cs_ranks.tex
Rscript results_processing/thesis_tables/me_cs_tables.R "" "" "" extra   # + me_cs_learners_extra.tex
```

`me_cs_tables.R` always prints the estimated-π version and the other two
pairwise accuracies, and writes them as `me_cs_learners_extra.tex` only when
its fourth argument is `extra` (an empty argument keeps that path's default).

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
- **H2O AutoML has no `max_runtime_secs` cap**; it is the main source of
  walltime uncertainty.

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
