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
└── regression_check.R  old-vs-new equivalence harness

continuous/           continuous outcome, sample size sweep
binary/               binary outcome, sample size sweep
competing_risk/       competing risks - the target setting
missing/              missing covariates (continuous, binary, CI example)
confidence_intervals/ interval estimation (continuous, binary, optimal_sf)
crossfitting/         does double crossfitting earn its 4.5x cost?
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
`continuous/`. By the time that stopped, `continuous/cts_models.R` and
`binary/bin_models.R` differed in **two** places out of 438 lines, the same DGM
existed in four files, and the same collect/check boilerplate in eight.

Consolidating removed about 3,000 lines. The more useful outcome is that the
copies can no longer drift — which is where most of the bugs below came from.
`crossfitting/` was the model: it already sourced shared code rather than
forking it.

See `R/README.md` for the four axes of `cate_methods()`, the orchestration
profiles, and the grid contract.

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

Nuisances use **double crossfitting**: models are fit over all `C(V,2)` fold
pairs, so the pseudo-outcome column used to predict fold `k` was never touched by
fold `k`. At `V = 10` that is 45 fits rather than 10. Whether it earns that cost
is exactly what `crossfitting/` measures.

## Metrics

**Point estimation** — bias, ATE bias, MSE, RMSE, MAE, correlation, Spearman
correlation, sign accuracy. `bias` is `estimate - truth` throughout.

**HTE detection** — BLP test p-value (`GenericML`), independence test p-value
(`coin`).

**Intervals** — marginal coverage, simultaneous coverage, mean width.

## Running a study

```bash
qsub continuous/jobscripts/cts_1.sh      # the array job
Rscript continuous/cts_check.R           # any missing runs?
qsub continuous/jobscripts/cts_collect.sh
qsub continuous/jobscripts/cts_metrics.sh
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
| H | propensity trimming only on the SuperLearner path | shared | documented |
| I | competing risks never validates its SuperLearner library | `competing_risk` | open |
| J | stale comments and filenames | various | yes |

Three more surfaced along the way:

- `missing/binary` was a **half-converted fork**: continuous coefficients, a
  continuous power calculation, and truth on the **log-odds** scale while every
  estimator targets a risk difference
- `binary`'s grid was declared three ways and they disagreed. The runs
  themselves used the intended four-scenario design, so the results are sound;
  the stale `c(1:10)` in the analysis script was simply out of step
- `combine_mi()` in `missing/ci_example` read `alpha` as a **free variable** from
  the global environment

## Status

| study | |
|---|---|
| `confidence_intervals/continuous`, `optimal_sf/cts` | unchanged — must stay bit-identical |
| `missing/ci_example` | unchanged |
| `crossfitting` | unchanged |
| `continuous`, `missing/continuous` | re-run for bug F (`dr_superlearner` only) |
| `binary` | re-run for bug F (`dr_superlearner` only) |
| `missing/binary` | re-run — the DGM was wrong three ways |
| `confidence_intervals/binary`, `optimal_sf/bin` | re-run — the DGM was wrong |
| `competing_risk` | **currently fails to run** — see its README |

Roughly 32,000 array jobs. Bug G costs no cluster time: it is computed from the
saved `*_all.RDS` files, so only the metrics scripts and the figures rerun.

## Verifying a change

```bash
Rscript R/regression_check.R baseline   # before touching anything
Rscript R/regression_check.R verify     # after - must be 8/8
Rscript crossfitting/cf_testing.R       # independent check of the estimators
```

The harness fingerprints generated datasets, estimates **and** saved nuisance
structure, and runs each study in its own subprocess. It proves the code
reproduces the current behaviour on this machine; it is not a claim about the
cluster's numbers (R 4.5.3 locally vs 4.3.2 there).

## Dependencies

Managed with `renv`:

```r
renv::restore()
```

Key packages: `grf`, `SuperLearner`, `GenericML`, `pseudo`, `coin`, `furrr`,
`future`, `ranger`, `glmnet`, `gam`, `mice`, `missForest`, `VIM`, `dplyr`,
`here`.
