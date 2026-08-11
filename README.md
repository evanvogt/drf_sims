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

`crossfitting/` compared double crossfitting (fitting nuisances over all
`C(V,2)` fold pairs, 45 fits at `V = 10` rather than 10) against cheaper
alternatives for DR-RF, DR-SL and causal forest. Based on that comparison,
`R/cate_models.R` now uses, per method:

- **DR-RF, DR-Oracle, DR-Semi-Oracle** — whole-sample, out-of-bag: an
  S-learner nuisance forest with no sample splitting, and an OOB stage-2
  regression forest (`nuisance_rf` / `stage2_whole_rf`).
- **Causal Forest** — `grf`'s own internal cross-fitting (a plain
  `causal_forest(X, Y, W)` with no externally-supplied nuisances).
- **DR-SL** — a single leave-one-fold-out crossfit, with the *same* fold
  assignment shared by the nuisance stage and the stage-2 regression
  (`nuisance_sl` / `stage_2_sl`), rather than double-crossfit nuisances
  feeding a separately-split stage 2.

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
qsub continuous/jobscripts/cts_1.sh         # the array job — cluster only
Rscript continuous/cts_check.R              # any missing runs? — runs locally too
qsub continuous/jobscripts/cts_collect.sh   # cluster only
qsub continuous/jobscripts/cts_metrics.sh   # cluster only
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
| `continuous`, `binary`, `missing/continuous`, `missing/binary`, `missing/ci_example`, `confidence_intervals/continuous`, `confidence_intervals/binary`, `confidence_intervals/optimal_sf` (cts), `confidence_intervals/optimal_sf` (bin), `validation/continuous` | **re-run — crossfitting strategy changed** (see Methods above), on top of any bug-fix re-run already listed below |
| `crossfitting` | its own comparison arms are unchanged; only the production consumers of `R/cate_models.R` above moved |
| `continuous`, `missing/continuous` | also re-run for bug F (`dr_superlearner` only) |
| `binary` | also re-run for bug F (`dr_superlearner` only) |
| `missing/binary` | also re-run — the DGM was wrong three ways |
| `confidence_intervals/binary`, `confidence_intervals/optimal_sf` (bin) | also re-run — the DGM was wrong |
| `competing_risk` | **currently fails to run** — see its README; untouched by the crossfitting strategy change (independent estimator code) |
| `model_evaluation` | **first run, not a re-run** — independent estimator/nuisance code (see its README); unaffected by the crossfitting strategy change |

Roughly 32,000 array jobs. Bug G costs no cluster time: it is computed from the
saved `*_all.RDS` files, so only the metrics scripts and the figures rerun.

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
  "furrr", "future", "ranger", "glmnet", "gam", "mice", "missForest", "VIM",
  "dplyr", "here"))
```

**TODO:** once the package set stabilizes, run `renv::snapshot()` and commit
`renv.lock`, then switch this section back to `renv::restore()`.
