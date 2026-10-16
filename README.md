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

## Where each result in the thesis comes from

Keyed by topic. Every study folder has a README with its design, and an
`ADEMP.md` where the design is written out in full (aims, data-generating
mechanisms, estimands, methods, performance measures).

| question | study (code) | what it produces |
|---|---|---|
| How do bias, accuracy and HTE tests behave as sample size shrinks, for continuous and binary outcomes with correlated covariates? | `sample_size/{continuous,binary}/` | figures `results_processing/thesis_figures/sample_size.R`; tables `thesis_tables/ss_tables.R`, `ss_test_tables.R` |
| Do CATE intervals reach their nominal coverage? | `sample_size/confidence_intervals/{continuous,binary}/`, with `optimal_sf/` choosing the half-sample fraction | figures `thesis_figures/ss_ci.R`; tables `thesis_tables/ss_ci_tables.R` |
| What do missing covariates do, and which handling method is best? | `missing/{continuous,binary}/` | figures `thesis_figures/missing.R`; tables `thesis_tables/miss_tables.R` |
| How should an interval be built after multiple imputation? | `missing/ci_example/` | figures `thesis_figures/miss_cts_ci.R` |
| Do the estimators work with competing risks (the target setting)? | `competing_risk/` (and `competing_risk/single_event/`, the same without the competing event) | figures `thesis_figures/surv.R`; tables `thesis_tables/surv_tables.R` |
| Why does the code use the crossfitting strategy it does? | `crossfitting/` (and `crossfitting/confidence_intervals/`) | report `crossfitting/cf_results.qmd` |
| Can a cheap proxy score pick the right CATE learner? | `model_evaluation/` | report `model_evaluation/me_results.qmd` |
| Do subgroups, CATE variance and variable importance found at an interim analysis replicate on the rest of the trial? | `validation/{continuous,binary}/` | report `validation/val_results.qmd` (Quarto `outcome` and `scenario` params) |

The shared code is in `R/`; the layout section below says what is where.

## Reproducing the results

There are three levels, from cheapest to dearest.

**1. Check the code (a laptop, about 15 seconds).**

```bash
Rscript tests/run_tests.R
```

It checks that the data-generating mechanisms deliver the effects they plan
for, that every study's jobscript ranges match its parameter grid, that the
bootstrap bands have their nominal coverage, and that the metrics have the
conventions the text states. It needs no results.

**2. Rebuild the figures and tables from the metrics files (a laptop).** Each
study's metrics step writes a small `*_metrics.RDS` into its results folder,
and the scripts in `results_processing/` read only those, for example
`../results/competing_risk/surv_metrics.RDS` and
`../results/correlated/continuous/cts_corr_metrics.RDS`. With those in place:

```bash
cd results_processing/thesis_figures && Rscript surv.R
cd ../thesis_tables && Rscript surv_tables.R
```

`results_processing/README.md` lists which script reads which study. The
metrics files are not in this repository (see "Where the results live").

**3. Re-run a study (the HPC cluster).** The simulations are about 128,000
independent array jobs in total, each a few minutes to a few hours. They were
run on Imperial's PBS cluster and are not meant to be re-run on a laptop. One
run of any study can be tried locally, as a smoke test, by giving the analysis
script a row number of the study's grid:

```bash
cd sample_size/continuous
Rscript cts_corr_analysis.R 1        # row 1 of the grid; writes to ../results, beside the repository
```

See "Running a study" below for the full cluster sequence (submit, check,
collect, metrics).

### Where the results live

Results are large, so they are written **outside the repository**, next to it:

```
<parent>/drf_sims/       this repository
<parent>/results/        per-run files, <study>/<parameters>/res_sim_<run>.RDS
                         and each study's collected *_all.RDS and *_metrics.RDS
```

Every path comes from the study's `<prefix>_config.R` (`res_path`), built from
`here()`, so the repository can sit anywhere as long as `results/` is its
sibling. <!-- add the location of the deposited *_metrics.RDS files here -->

### Software

The results were produced on Imperial's cluster with R 4.3.2 from the conda
environment `sim-env`, whose exact contents are in `environment.yml`; the
versions of the packages the code loads are in `sessionInfo_cluster.txt`. See
"Dependencies" at the end.

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

sample_size/          continuous and binary outcomes, sample size sweep,
                      correlated covariates at rho 0 / 0.5
                      (continuous, binary, confidence_intervals/, snr/)
competing_risk/       competing risks - the target setting
missing/              missing covariates (continuous, binary, CI example)
crossfitting/         compared double crossfitting against cheaper alternatives
model_evaluation/    does cheap proxy scoring pick the right CATE learner?
validation/           do CATE subgroups/variance/importance found at an interim
                      analysis replicate on the rest of the trial? (continuous)
results_processing/   thesis figures and tables
tests/                quick invariant checks (Rscript tests/run_tests.R)
docs/                 DEVLOG.md - bug ledger and re-run history
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

Every folder has its own README with that study's design and how to run it.
`R/README.md` covers the axes of `cate_methods()`, the orchestration profiles
and the grid contract.

## Scenarios

The binary and continuous studies share one set of ten treatment-effect
scenarios (`R/dgm_scenarios.R`). The four the chapters report come first:

| scenario | treatment effect |
|---|---|
| 1 | none (null) |
| 2 | simple — continuous X4 |
| 3 | complex — X3 + X4 + X4·X5 |
| 4 | non-linear — cos(X4) |
| 5 | binary X3 |
| 6 | X3 + X4 |
| 7 | X3·X4 |
| 8 | X3 + X4 + X3·X4 |
| 9 | X4·X5 |
| 10 | exponential in X4 |

Results saved before 2026-09-26 use older scenario numbers; the map is in
`docs/DEVLOG.md`. The missing-data studies use scenarios 1–6. The subset
studies run: `crossfitting/` 1, 4, 6, 8; `model_evaluation/` 1–4;
`crossfitting/confidence_intervals/` 1, 4, 8; `validation/continuous/` and
`validation/binary/` 2–4.

Both outcomes are one model, `E[Y | x, W] = m0(x) + W·τ(x)`. For a binary
outcome the treatment effect is on the **risk-difference** scale, the scale of
the estimand, with the control risk bounded in [0.34, 0.70] (a 40% event rate).
The binary modifiers are the continuous ones with X4 and X5 through `tanh`,
signs reversed, and scaled per scenario to the largest effect the risk bounds
allow. See `sample_size/README.md`.

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
`R/cate_models.R` uses, per method:

- **DR-RF, DR-Oracle, DR-Semi-Oracle** — whole-sample, out-of-bag: outcome
  forests fit separately in each arm (a T-learner, `t_learner_rf`; each unit's
  own-arm prediction OOB), no sample splitting, and an OOB stage-2 regression
  forest (`nuisance_rf` / `stage2_whole_rf`).
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
cd sample_size/continuous/jobscripts   # jobscripts cd to ${PBS_O_WORKDIR}/..
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

A run that fails is found by `Rscript <study>/<prefix>_check.R`, which writes
`jobscripts/failed_ids.txt` and points `jobscripts/<prefix>_rerun.sh` at it.
Every jobscript writes its PBS output to a `jobscripts/logs*/` folder, which is
not in git, so create them once on the cluster before the first submit
(`Rscript make_log_dirs.R`). `Rscript check_all.R` on the login node counts
expected, found and missing runs for every study in `R/study_registry.R`.

`docs/DEVLOG.md` holds the history: the bug ledger (A–S), the re-run campaign,
the scenario renumbering and each study's past status notes. Nothing there is
needed to read or run the code.

## Verifying a change

```bash
Rscript tests/run_tests.R               # quick invariants, no results needed
Rscript R/regression_check.R baseline   # before touching anything
Rscript R/regression_check.R verify     # after - every study identical
Rscript crossfitting/cf_testing.R       # independent check of the estimators
```

An intentional behaviour change (like the crossfitting strategy switch above)
is the one case where `verify` is *expected* to fail — re-baseline afterward
once the diffs are confirmed to be exactly the fields the change touched, so
the baseline is an equivalence reference for the next refactor again.

The harness fingerprints generated datasets, estimates **and** saved nuisance
structure, and runs each study in its own subprocess. It proves the code
reproduces the current behaviour on this machine; it is not a claim about the
cluster's numbers (R 4.5.3 and newer packages locally, R 4.3.2 there - see Dependencies).

## Dependencies

Every result was produced inside one conda environment on the cluster,
`sim-env`, recreated with

```bash
conda env create -f environment.yml
```

`environment.yml` is `conda env export --no-builds` of it, so it pins R and
every package, and `sessionInfo_cluster.txt` is the `sessionInfo()` of a session
with all the study packages loaded. The packages the studies depend on, at the
versions that ran:

| | version | | version |
|---|---|---|---|
| R | 4.3.2 | `randomForestSRC` | 3.2.2 |
| `grf` | 2.4.0 | `pseudo` | 1.4.3 |
| `SuperLearner` | 2.0-29 | `mice` | 3.17.0 |
| `GenericML` | 0.2.2 | `missForest` | 1.5 |
| `coin` | 1.4-3 | `xgboost` | 1.7.11.1 |
| `ranger`, `glmnet`, `gam`, `earth` | 0.17.0, 4.1-8, 1.22-5, 5.3.4 | `SHAPforxgboost` | 0.2.0 |
| `future`, `furrr` | 1.67.0, 0.3.1 | `h2o` | 3.44.0.3 |

**The local machine is not the same.** Smoke tests, `tests/run_tests.R` and
`R/regression_check.R` were run locally on R 4.5.3 with newer packages (for
example `grf` 2.5.0, `mice` 3.19.0, `randomForestSRC` 3.9.0, `xgboost` 3.2.1.1).
Anything checked only locally, such as the `grf` weight-scale behaviour noted in
`R/missingness.R`, may differ on the cluster. Use `environment.yml` to reproduce
the cluster's numbers.
