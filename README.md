# drf_sims

Simulation studies evaluating CATE (Conditional Average Treatment Effect)
estimation across outcome types and data settings. Part of WP2 of a PhD at
Imperial College London on causal machine learning methods in survival and
competing-risks settings.

The question throughout: **when the treatment effect varies between people, which
estimators recover that variation, and under what conditions do they stop
working?** Each folder answers one version of it — as sample size shrinks, as the
outcome becomes binary or time-to-event, as covariates go missing, and when an
interval rather than a point estimate is needed.

## Layout

```
R/                    shared library - every study sources from here
├── utils.R           RNG streams, crossfitting fold bookkeeping
├── dgm_scenarios.R   scenario tables + the data generator
├── missingness.R     amputation and missing-data handling
├── cate_models.R     the estimators + cate_methods()
├── bootstrap_ci.R    half-sample bootstrap, MI pooling, sf calibration
├── metrics.R         metric definitions + the metrics pipeline
├── pipeline.R        study configs, collect, check
├── figures.R         display labels, palette, plot helpers
├── bin_verify_hte.R  checks the binary risk-difference scenarios, re-derives RD_SCALE
└── regression_check.R  old-vs-new equivalence harness

sample_size/
├── correlated/           continuous and binary outcomes, sample size sweep,
│                         correlated covariates at rho 0 / 0.5
│                         (continuous, binary, confidence_intervals/)
├── continuous/           RETIRED 2026-10-06 - README only (method notes)
├── binary/               RETIRED 2026-10-06 - README only (method notes)
└── confidence_intervals/ RETIRED 2026-10-06 - READMEs only (method notes)
competing_risk/       competing risks - the target setting
missing/              missing covariates (continuous, binary, CI example)
crossfitting/         compared double crossfitting against cheaper alternatives
model_evaluation/    does cheap proxy scoring pick the right CATE learner?
validation/           do CATE subgroups/variance/importance found at an interim
                      analysis replicate on the rest of the trial? (continuous)
results_processing/   thesis figures
scratch/              unmaintained exploratory code
```

Each study folder has the same shape:

| file | |
|---|---|
| `<s>_config.R` | the parameter grid and results path — **the** definition |
| `<s>_dgms.R` | names this study's scenario set |
| `<s>_models.R` | this study's configuration of the shared estimators |
| `<s>_analysis.R` | array entry point; one grid row per index |
| `<s>_check.R` | finds missing runs |
| `<s>_collect.R` | gathers per-run files |
| `<s>_metrics.R` | computes metrics |
| `jobscripts/` | PBS submission scripts |

Results are written **outside the repo**, to `../results/<study>/...`.

Every folder has its own README with that study's design, gotchas and status.

## Why `R/` exists

The repo grew by copy-paste: each new study started as a duplicate of
`sample_size/continuous/`. By the time that stopped, `sample_size/continuous/cts_models.R` and
`sample_size/binary/bin_models.R` differed in **two** places out of 438 lines, the same DGM
existed in four files, and the same collect/check boilerplate in eight.

Consolidating removed about 3,000 lines. The more useful outcome is that the
copies can no longer drift — which is where most of the bugs below came from.
`crossfitting/` was the model: it already sourced shared code rather than
forking it.

See `R/README.md` for the axes of `cate_methods()`, the orchestration
profiles, and the grid contract.

## Scenarios

The binary and continuous studies share one set of ten treatment-effect
scenarios (`R/dgm_scenarios.R`). The four the chapters report come first:

| scenario | treatment effect | before 2026-09-26 |
|---|---|---|
| 1 | none (null) | 1 |
| 2 | simple — continuous X4 | 3 |
| 3 | complex — X3 + X4 + X4·X5 | 8 |
| 4 | non-linear — cos(X4) | 9 |
| 5 | binary X3 | 2 |
| 6 | X3 + X4 | 4 |
| 7 | X3·X4 | 5 |
| 8 | X3 + X4 + X3·X4 | 6 |
| 9 | X4·X5 | 7 |
| 10 | exponential in X4 | 10 |

The last column is the old number. The renumbering changed no data: each
scenario generates exactly what it did under its old number. Everything saved
before it uses the old numbers: `../collected_metrics/`, results directories,
figures and the `results_processing/` notebooks.
The missing-data studies use scenarios 1–6 under the same numbers. The subset
studies run: `crossfitting/` and `model_evaluation/` 1, 4, 6, 8;
`crossfitting/confidence_intervals/` 1, 4, 8; `validation/continuous/` 2.

Both outcomes are one model, `E[Y | x, W] = m0(x) + W·τ(x)`. For a binary
outcome the treatment effect is on the **risk-difference** scale, the scale of
the estimand, with the control risk bounded in [0.34, 0.70] (a 40% event rate).
The binary modifiers are the continuous ones with X4 and X5 through `tanh`,
signs reversed, and scaled per scenario to the largest effect the risk bounds
allow. Until 2026-09-26 the binary effect was on the logit scale, which made X1
and X2 effect modifiers of the risk difference too. See `sample_size/binary/README.md`.

## Methods

| method | |
|---|---|
| Causal Forest | `grf::causal_forest` with cross-fitted nuisances |
| DR-RF | doubly-robust learner, random forest second stage |
| DR-SL | doubly-robust learner, SuperLearner throughout |
| DR-Oracle | true outcome model, known propensity |
| DR-Semi-Oracle | known propensity, estimated outcome model |

Competing risks adds IPW-transformed RMST, cause-specific and subdistribution
causal survival forests, and pseudo-value approaches — see
`competing_risk/README.md`.

`crossfitting/` compared double crossfitting (fitting nuisances over all
`C(V,2)` fold pairs, 45 fits at `V = 10` rather than 10) against cheaper
alternatives for DR-RF, DR-SL and causal forest. Based on that comparison,
`R/cate_models.R` now uses, per method:

- **DR-RF, DR-Oracle, DR-Semi-Oracle** — whole-sample, out-of-bag: outcome
  forests fit separately in each arm (a T-learner, `t_learner_rf`; each unit's
  own-arm prediction OOB), no sample splitting, and an OOB stage-2 regression
  forest (`nuisance_rf` / `stage2_whole_rf`). Until 2026-09-27 the outcome
  forest was an S-learner on `cbind(W, X)` (crossfitting's `oob_oob_s`).
- **Causal Forest** — `grf`'s own internal cross-fitting (a plain
  `causal_forest(X, Y, W)` with no externally-supplied nuisances).
- **DR-SL** — a single leave-one-fold-out crossfit, with the *same* fold
  assignment shared by the nuisance stage and the stage-2 regression
  (`nuisance_sl` / `stage_2_sl`), rather than double-crossfit nuisances
  feeding a separately-split stage 2. The outcome model is one SuperLearner
  per arm, and each nuisance has its own library (`R/sl_library.R`).

See `crossfitting/README.md` for the full arm comparison behind this choice.

## Metrics

**Point estimation** — bias, ATE bias, MSE, RMSE, MAE, correlation, Spearman
correlation, sign accuracy. `bias` is `estimate - truth` throughout.

**HTE detection** — BLP test p-value (`GenericML`), independence test p-value
(`coin`).

**Intervals** — marginal coverage, simultaneous coverage, mean width.

## Running a study

`qsub` submits to Imperial's HPC cluster (PBS) — those lines only work
there. `Rscript` lines run anywhere, including as a local smoke test.

```bash
cd sample_size/correlated/continuous/jobscripts   # jobscripts cd to ${PBS_O_WORKDIR}/..
qsub cts_corr_1.sh            # the array job — cluster only
Rscript ../cts_corr_check.R   # any missing runs? — runs locally too
qsub cts_corr_collect.sh      # cluster only
qsub cts_corr_metrics.sh      # cluster only
```

The array index is a **row number** of `study$grid`. Never filter or reorder the
grid — that renumbers every job. To run a subset:

```r
idx <- grid_indices(study, method = "complete_data")
```

## Bug ledger

Found during the de-duplication. Each is written up in the relevant folder README.

| | what | where | fixed |
|---|---|---|---|
| A | ran on the **continuous coefficient table** on a logit scale | `confidence_intervals/binary` | yes — re-run |
| B | collect looked for `AUX` where everything else said `MNAR` | `missing/continuous` | yes |
| C | `rel_efficiency` was `NA` everywhere — the reference arm was never collected | `missing/*` | yes |
| D | grid filtered after `expand.grid`, so the array index meant two different things | `missing/*` | yes — structurally |
| E | `.gitignore` `*test*.*` hid the verification harness from git | repo | yes |
| F | `stage_2_sl` discarded the pretested SuperLearner library; the correct branch was dead code | shared | yes — re-run |
| G | `bias` and `ate_bias` had **opposite signs** | all metrics | yes |
| H | propensity trimming only on the SuperLearner path | shared | yes — resolved when `nuisance_rf` moved to whole-sample OOB and picked up `trim_ps` too |
| I | competing risks never validates its SuperLearner library | `competing_risk` | yes — resolved when it adopted the shared crossfitting strategy and picked up `pretest_superlearner`; the split-pseudo T-learner branch is the one remaining unvalidated call site |
| J | stale comments and filenames | various | yes |
| K | `pretest_superlearner()` could return an empty SL library, crashing `nuisance_sl`/`stage_2_sl` | shared | yes |
| L | `run_blp_whole()` had no `tryCatch`, crashed on a constant/degenerate CATE | shared | yes |
| M | `dr_oracle` handed log-odds as outcome predictions — `6b06db3` dropped the `plogis` from the oracle formulas but not `oracle_link = "identity"` | `missing/binary` | yes — no re-run: every finished result predates it |
| N | binary MNAR-Y truth evaluated at U = 0 rather than averaged over U — equal on the identity scale, not the logit one | `missing/binary` | yes — repaired at metrics time from the saved truth; no re-run |
| O | continuous `bW` calibration used `sd = s_err + s2` (ignoring b1, b2 and the heterogeneity variance, and adding SDs) and set `bW` rather than the ATE, so realised power ran 3–100% and scenarios 3, 5, 8 had a positive ATE; b0, b1, b2 also varied by scenario. Now one baseline (0.4, −0.5, 1), and each trial planned for 80% power under homogeneity with the true ATE equal to the planned effect — realised power 61–80% as heterogeneity grows | `R/dgm_scenarios.R` — every continuous study | yes — re-run |
| P | binary `bW` calibration planned at 75% power at the risk plogis(b0), ignoring X1 and X2 (0.401 against a population risk of 0.453), and set `bW` rather than the ATE (a marginal risk difference), so `bW` was the same in every scenario and realised power ran 5–99%; the modifiers in scenarios 2, 5, 6 and 7 also had the opposite sign to the continuous ones, and `b5` was a column no scenario used. Now each trial is planned for 80% power under homogeneity with the true RD equal to the planned effect in every scenario (realised power 79–81%), the signs match `sample_size/continuous/`, and `b5` is gone. *The risk-difference DGM (2026-09-26) has since reversed every binary modifier's sign relative to `sample_size/continuous/`, on purpose — see `sample_size/binary/README.md`* | `R/dgm_scenarios.R` — every binary study | yes — re-run |
| Q | `pretest_superlearner()` dropped a learner on **any warning** (`tryCatch(warning = )` aborts the fit at the first one), so benign glm/gam warnings — and SL.mean's "All algorithms have zero weight" for a propensity near 0.5 — cut folds' libraries to one or two learners; it also kept a learner if even one prediction was non-NA. Now warnings are recorded and muffled, a learner is dropped only on an error or non-finite predictions, and the drops are saved as `sl_dropped`. Fixed alongside the move to per-nuisance libraries (`R/sl_library.R`) | shared, `crossfitting`, `competing_risk` | yes — re-run (SuperLearner arms only) |

Three more surfaced along the way:

- `missing/binary` was a **half-converted fork**: continuous coefficients, a
  continuous power calculation, and truth on the **log-odds** scale while every
  estimator targets a risk difference
- `binary`'s grid was declared three ways and they disagreed. Submitting
  indices 1–1600 against the analysis script's ten-scenario grid ran runs 1–40
  of all ten scenarios, so the results on disk have 40 replicates per cell,
  not 100 (see `sample_size/binary/bin_config.R` at tag `independent-ss-final`). The study re-runs in full anyway
- `combine_mi()` in `missing/ci_example` read `alpha` as a **free variable** from
  the global environment

## Status

**Retired 2026-10-06:** the independent-covariate sample-size studies
(`continuous`, `binary`, and `confidence_intervals/{continuous, binary,
optimal_sf}`) are not used in the thesis. The correlated studies replace them,
with their ρ = 0 arm as the independent-covariate comparator. Their code is at
git tag `independent-ss-final`, and their READMEs stay for the method notes.
Their results are archived with
`R/archive_old_results.R --label retired_2026-10 --trees continuous binary confidence_intervals`.
The rows below that name them are historical.

**Step 0 — archive the old results.** Nothing under `../results` is what the
current code produces (bug O and the risk-difference DGM, below), and every
`scenario_<k>/` directory there uses the pre-2026-09-26 numbers. Re-running into
those trees would overwrite some results and leave others to be collected under
the wrong scenario, so pack them away first — on the cluster, from the repo
root:

```bash
Rscript R/archive_old_results.R              # dry run: what goes where
cd R/jobscripts; qsub archive_old_results.sh # one .tar per study, then deletes the tree
```

Each study's tree becomes `results/_archive/pre_2026-09-26/<path>.tar`
(uncompressed: the result files are gzipped already), and is deleted only once
its `.tar` lists every file; `tar -xf <that .tar>` from `results/` restores it.
Every study below then re-runs from empty. `competing_risk/` is left alone.

| study | |
|---|---|
| `continuous`, `binary`, `missing/continuous`, `missing/binary`, `missing/ci_example`, `confidence_intervals/continuous`, `confidence_intervals/binary`, `confidence_intervals/optimal_sf` (cts), `confidence_intervals/optimal_sf` (bin), `validation/continuous` | **re-run — crossfitting strategy changed** (see Methods above), on top of any bug-fix re-run already listed below |
| `crossfitting` | its own comparison arms are unchanged by the crossfitting change — but it re-runs for bug O, below |
| `continuous`, `missing/continuous` | also re-run for bug F (`dr_superlearner` only) |
| `binary` | also re-run for bug F (`dr_superlearner` only) |
| `missing/binary` | also re-run — the DGM was wrong three ways |
| `confidence_intervals/binary`, `confidence_intervals/optimal_sf` (bin) | also re-run — the DGM was wrong |
| every continuous study: `continuous`, `missing/continuous`, `missing/ci_example`, `confidence_intervals/continuous`, `confidence_intervals/optimal_sf` (cts), `validation/continuous`, `model_evaluation`, `crossfitting`, `crossfitting/confidence_intervals` | also re-run for bug O — the shared baseline and `bW` changed, which moves every dataset and the level of the true CATE. For `crossfitting`, `crossfitting/confidence_intervals` and `model_evaluation` this is new: any finished results they hold are stale. See `continuous/README.md`'s "Outcome model and `bW` calibration" |
| every binary study: `binary`, `missing/binary`, `confidence_intervals/binary`, `confidence_intervals/optimal_sf` (bin) | also re-run for bug P and the **risk-difference DGM**. Bug P changed the `bW` calibration and the modifier signs in scenarios 2, 5, 6 and 7. The risk-difference DGM then put the treatment effect on the risk-difference scale, with a bounded control risk and rescaled, sign-reversed modifiers. Together they move every dataset's outcome and the true CATE. Submit on the risk-difference code; any bug-P re-run submitted on the logit-scale DGM is superseded too. For `missing/binary` this is new: its finished rows 1–9,900 are superseded, so all 12,600 rows re-run. See `binary/README.md`'s "Outcome model and `bW` calibration" |
| `competing_risk` | **first run under the new strategy** — it has now adopted the crossfitting change (it was the last production study still double-crossfitting) and runs clean end-to-end. Its pseudo-value and SuperLearner frameworks each ship in several arms, because it also crosses a second factor — whole-sample vs crossfit pseudo-values. See its README |
| `model_evaluation` | independent estimator/nuisance code (see its README), whose 9 candidates moved onto the shared single-crossfit strategy. Its 358/360 runs on that strategy predate bug O, so they are archived with the rest and the study re-runs — main, strategies and split trees |

Roughly 32,000 array jobs. Bug G costs no cluster time: it is computed from the
saved `*_all.RDS` files, so only the metrics scripts and the figures rerun.

### Tracking the rerun

Each study still resubmits its own failures the usual way:
`Rscript <study>/<prefix>_check.R` writes `jobscripts/failed_ids.txt`, then
`qsub jobscripts/<prefix>_rerun.sh` resubmits them. The check script also
rewrites that rerun jobscript - `-J` to match the number of failures, and the
resource request to sit above `<prefix>_1.sh`'s (one more core, 1.2x the
memory, 2x the walltime, capped at 72 hours and never lowered below what the
script already asks for). There is no `-J` to set by hand any more.

Every jobscript writes its PBS output to a `jobscripts/logs*/` directory, and
`*logs*/` is gitignored, so none of them exist after a fresh clone and PBS
rejects the job. Create them all on the cluster before the first submit:

```bash
Rscript make_log_dirs.R            # --dry-run to just list them
```

For a bird's-eye view
across every study at once - how many jobs are expected, found and missing,
and why each study is being rerun - run:

```bash
Rscript check_all.R
```

on the HPC login node (results only exist there - too large to sync back).
It's read-only (never writes `failed_ids.txt`, never calls `qsub`), so it's
safe to run as often as useful during the campaign. It writes
`check_all_studies.csv`/`.md` next to itself; those two files are committed,
so `git push` on the HPC + `git pull` locally is how progress gets checked
from off the cluster. The registry of studies it scans lives in
`R/study_registry.R` - add a row there for any new study.

## Resource profiling (removed)

Profiling with the [`syrup`](https://simonpcouch.github.io/syrup/) package is
**not effective for this study**, and the profiling scripts have been deleted.
Each study used to have a `<prefix>_profile.R` sweep (timing / memory / CPU over
`workers` × `grf_threads`, instrumented with `syrup`), a `jobscripts/<prefix>_profile.sh`
array job, and a `<prefix>_profile_summary.R` that wrote the measured `#PBS -l`
directives into `<prefix>_1.sh`. None of it produced usable numbers:

- on the cluster, `syrup()` starts its sampler with `callr::r_session$new()` and a
  hardcoded 3 s timeout that it gives no way to raise. The profiling subjobs died
  with `Could not start R session, timed out` — most likely because starting an R
  session off the networked `/rds` filesystem takes longer than that
- on cells with `workers > 1` the `future::multisession` workers never showed up
  in the process tree `syrup` samples, so the CPU / memory figures missed them
- locally `syrup` measures nothing useful on Windows

So the `#PBS -l` lines and trailing `Rscript` arguments in the `*_1.sh`
jobscripts are hand-set, not measured. Any that still say *placeholder* stay
that way — adjust them by hand from how runs actually behave (the `_rerun.sh`
resource bump above does part of this automatically). The study READMEs, code
comments and jobscript headers have been updated to match.

The deleted files (the 26 `*_profile*.R`/`.sh` scripts, plus
`crossfitting/cf_diagnose_sampler.R`, `cf_diagnose_multisession.R` and
`jobscripts/cf_diagnose.sh`, which diagnosed the failures above) are in git
history:

```bash
git log --diff-filter=D --name-only -- '*_profile*' '*cf_diagnose*'   # find the deleting commit
git checkout <commit>^ -- crossfitting/cf_profile.R                  # restore a file from its parent
```

## Verifying a change

```bash
Rscript R/regression_check.R baseline   # before touching anything
Rscript R/regression_check.R verify     # after - must be 8/8
Rscript crossfitting/cf_testing.R       # independent check of the estimators
```

An intentional behaviour change (like the crossfitting strategy switch above)
is the one case where `verify` is *expected* to fail — re-baseline afterward
once the diffs are confirmed to be exactly the fields the change touched, so
the baseline is an equivalence reference for the next refactor again.

The harness fingerprints generated datasets, estimates **and** saved nuisance
structure, and runs each study in its own subprocess. It proves the code
reproduces the current behaviour on this machine; it is not a claim about the
cluster's numbers (R 4.5.3 locally vs 4.3.2 there).

## Dependencies

renv isn't pinned yet — there's an `renv/` folder but no `.Rprofile` or
committed `renv.lock`, so `renv::restore()` currently has nothing to restore
from. Until a lockfile exists, install directly:

```r
install.packages(c("grf", "SuperLearner", "GenericML", "pseudo", "coin",
  "furrr", "future", "ranger", "glmnet", "gam", "earth", "mice",
  "missForest", "VIM", "dplyr", "here"))
```

**TODO:** once the package set stabilizes, run `renv::snapshot()` and commit
`renv.lock`, then switch this section back to `renv::restore()`.
