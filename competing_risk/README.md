# Competing risks

CATE estimation when event 1 can be pre-empted by a competing event 2. This
is the study the rest of the repo builds towards, and the only one whose
estimators are genuinely different rather than a configuration of the shared
ones.

| | |
|---|---|
| scenarios | 1–7 |
| n | 500 |
| censoring | TRUE, FALSE |
| covariate correlation ρ | 0.5 (primary), 0 (sensitivity) |
| runs | 500 |
| array | **14,000 jobs**: 1–7000 ρ = 0, sensitivity (`surv_1.sh`), 7001–14000 ρ = 0.5, primary (`surv_2.sh`) |
| horizon | 28 |
| results | `../results/competing_risk/rho_<ρ>/scenario_<k>/<n>/censor_<TRUE\|FALSE>/` |

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

## Design

Event 1 (E1) and event 2 (E2) are generated from Weibull
distributions via joint hazards (Beyersmann 2009). Scenarios vary where the
treatment acts:

| # | |
|---|---|
| 1 | ATE on E1 only |
| 2 | ATE on E2 only |
| 3 | HTE on E1, no ATE on E2 |
| 4 | HTE on E1, ATE on E2 |
| 5 | HTE on E2, no ATE on E1 |
| 6 | HTE on E2, ATE on E1 |
| 7 | HTE on both |

Scenarios 2, 5 and 6 are the interesting ones: the treatment moves event 2, so
a naive analysis of event 1 can show an apparent effect that is really a change
in who survives long enough to have one.

The event 1 baseline scale is 25 and the effects are log-HR-matched across
events (HR 1.42 / 2.86 on event 1, 0.70 / 0.35 on event 2, at X3 = 0 / 1), so
the scenarios separate on the RMTL scale as well as the hazard scale (ADEMP.md,
"Why these values").

**Covariates.** X1–X5 and the noise X01–X03 come from the copula `sample_size/`
uses (`correlated_covariates()` in `R/dgm_scenarios.R`), at ρ ∈ {0, 0.5}, so ρ
is a grid factor and a `rho_<ρ>/` level in the results path. Estimators see the
same 10 covariates as `sample_size/` (X1–X5, X01–X05); X4 and X5 enter no
hazard here, so they are noise.

- τ(x) is identical at both ρ, and the scenario labels separate on the RMTL
  scale at ρ = 0.5 too. `Rscript surv_dgm_check.R` gives the numbers (ADEMP
  "What the DGM implies").
- The two ρ share each run's seed, so compare them with paired per-run
  differences. `surv_results.qmd`, "Sensitivity analysis", does this. Every
  other section of the report is ρ = 0.5.

`get_ipw()` keeps a unit censored or with the competing event after the
horizon (Y ≥ 28), as in the grf tutorial; dropping them biases `ipw` badly.

### Scenario 5 flips the sign of the RMTL1 CATE — expect it

In scenario 5 almost every arm correlates *negatively* with `tau_RMTL1` while
recovering `tau_RMTL2` correctly. This is a property of the scenario, not a
bug:

- **The cause-specific truth is flat.** `tau_RMST1_cs` has zero variance in
  scenario 5 — the treatment has no cause-specific effect on event 1 at all.
  Every bit of `tau_RMTL1` heterogeneity is the competing-risk artefact: the
  treatment changes who is still at risk for event 1, nothing else.
- **The two truths are near-mirror images.** `cor(tau_RMTL1, tau_RMTL2)` is
  −0.92 in the truth at ρ = 0.5 (−0.91 at ρ = 0), because the RMTL1
  heterogeneity is mechanically induced by the RMTL2 heterogeneity. So an
  estimator that tracks the real (event-2) signal correctly lands
  anti-correlated with the *induced* RMTL1 estimand.
- **It is consistent across arms**, whole-sample OOB, single-crossfit and
  split-pseudo alike, which a crossfitting artefact would not be.

This is the phenomenon scenarios 2/5/6 exist to expose, so it belongs in the
write-up rather than being chased as a defect.

**Downstream consequence:** because `tau_RMST1_cs` is constant here, `cor()`
returns `NA` for the `ipw` and `csf_cs` Event 1 rows in scenario 5, and
`surv_metrics.R` emits `the standard deviation is zero` warnings. The `corr`,
`spearman` and `c_stat` columns are `NA` for exactly those cells. That is the
statistically correct answer — there is no heterogeneity to rank — but it means
those cells must be dropped rather than treated as zeros when summarising.

## Estimators

Five frameworks, each targeting restricted mean survival or lost time:

| framework | |
|---|---|
| `ipw` | causal forest on IPW-transformed RMST outcomes |
| `csf_cs` | causal survival forest, competing events treated as censoring (cause-specific) |
| `csf_sh` | causal survival forest, competing events kept in the risk set (subdistribution / Fine-Gray) |
| `pseudo_cf` | causal forest on pseudo-values |
| `pseudo_dr` | DR-learner on jackknife pseudo-values |

Plus SuperLearner T-learner and DR-learner variants. The DR-learner's split
pseudo-observation construction (Cwiling et al. 2025) is **disabled** — see
"Known issues" below.

The two pseudo-value frameworks and the SuperLearner ones each ship in
**several arms**, because this study also asks how the pseudo-values should be
built — see "Crossfitting" below.

**Semi-oracle DR-learners.** `pseudo_dr_whole_oob_semi`, `sl_dr_whole_semi` and
`rsf_dr_oob_semi` repeat the three production DR-learners with the true
propensity 0.5 (`TRUE_PS`) in the pseudo-outcome in place of the estimated one.
They reuse the parent arm's fitted outcome models (the nuisance functions
return `po_semi` next to `po`), so each pair differs in the propensity alone.
The only extra fit is stage 2. They run after the RSF arms, so the other arms
keep their RNG stream. In `surv_results.qmd` they are variant "Semi-oracle
(e = 0.5)" within their parent's family.

**The truth column depends on the framework**, which is why
`surv_metrics.R` carries `framework_truth_map`: `ipw` and `csf_cs` remove
competing events so they target the cause-specific (net) RMST, while `csf_sh`
keeps them in the risk set and targets the subdistribution RMST. Scoring one
against the other's truth would be a category error.

## Crossfitting

The production strategy matches `R/cate_models.R` — see that file's header
comment and `crossfitting/README.md` for the benchmark the choice rests on:

| arm family | strategy | |
|---|---|---|
| `ipw`, `csf_cs`, `csf_sh`, `pseudo_cf_whole_oob` | grf's own internal crossfitting, **`cf_default`** | plain `causal_forest`/`causal_survival_forest`, OOB `tau` |
| `pseudo_dr_whole_oob` | whole-sample OOB, T-learner (one forest per arm), **`oob_oob`** | `nuisance_pseudo_rf_oob` + `stage2_whole_rf` |
| the SuperLearner arms | single leave-one-fold-out, **`scf_scf`**; DR outcome model per arm | `nuisance_pseudo_sl` + `stage_2_sl`, sharing one `fold_indices` |
| `rsf_dr_oob` | whole-sample OOB, T-learner (one competing-risks RSF per arm) | `nuisance_rsf_oob` + `stage2_whole_rf` |
| `rsf_dr_scf` | single leave-one-fold-out, the OOB arm's crossfit twin | `nuisance_rsf_scf` + `stage_2_rf_scf` |

Every DR-learner outcome model (pseudo-value on X) is fit separately in each
arm - a T-learner, as in `R/cate_models.R`. Every DR arm clamps its propensity
with `R/utils.R::trim_ps`, and the SuperLearner arms pretest their library;
libraries are per nuisance, from `R/sl_library.R`.

### The pseudo-value comparison

Because the outcome here is a jackknife pseudo-value rather than an observed
`Y`, "which rows the nuisance model saw" and "which rows the outcome was
computed from" are two independent choices. The study varies both, crossed
as far as they can be:

- **pseudo-values** — `whole` (one Aalen–Johansen / Kaplan–Meier fit on all `n`,
  `pseudo_all`) vs `cvps` (leave-one-fold-out, `pseudo_crossfit`)
- **fitting** — `oob` (whole sample, grf-internal) vs `scf` (single crossfit)

`cvps + oob` is not a cell: column `k` of the crossfit matrix is `NA` on fold
`k`'s own rows, so it only exists inside a fold loop. That is why the
random-forest families get three arms rather than four, with `whole_scf` as the
control that separates the pseudo-value factor from the fitting factor.

| framework | pseudo-values | fitting | |
|---|---|---|---|
| `pseudo_cf_whole_oob`, `pseudo_dr_whole_oob` | whole | whole-sample / grf-internal | **production parity** |
| `pseudo_cf_whole_scf`, `pseudo_dr_whole_scf` | whole | single CF | control |
| `pseudo_cf_cvps_scf`, `pseudo_dr_cvps_scf` | crossfit | single CF | comparison |
| `sl_t_whole`, `sl_dr_whole` | whole | single CF | |
| `sl_t_cvps`, `sl_dr_cvps` | crossfit | single CF | |
| `rsf_dr_oob` | whole (correction term only) | whole-sample OOB | |
| `rsf_dr_scf` | whole (correction term only) | single CF | |

SuperLearner has no OOB analogue — `crossfitting/README.md` drops the OOB arms
from its SuperLearner comparison for the same reason — so the SL families stay
`scf` throughout and vary the pseudo-values alone.

### What the comparison actually measures

**For `pseudo_cf` it is complete**: the pseudo-value *is* the outcome, so
`whole_scf` vs `cvps_scf` is a clean read on the pseudo-value question.

**For the DR arms it is partial.** The DR correction term needs a pseudo-value
for the *held-out* rows, and `pseudo_crossfit` has none there by construction —
a jackknife pseudo-value for observation `j` requires `j` to be in the sample
being decomposed. Obtaining one honestly means appending the test row to the
training KM set, which is exactly the split-pseudo construction that is
disabled (see "Known issues"). **So every DR arm keeps whole-sample
pseudo-values in its correction term**, and the factor being varied is only
*the pseudo-values used to train the nuisance regressions*. Do not report it as
more than that.

A second qualification: `pseudo_crossfit` substitutes the whole-sample value
wherever `pseudoyl`/`pseudomean` returns `NA`, which leaks the held-out fold
into the `cvps` arms. `results$pseudos$cf_cv$n_na_fallback` counts how often
that fired, per estimand, so the comparison can be qualified rather than
assumed clean. It is seed-dependent — check it rather than assuming 0.

### The RSF DR-learner (`rsf_dr_oob`, `rsf_dr_scf`)

The outcome model is a **competing-risks random survival forest**
(randomForestSRC) fit to the observed `(Y, D)`, rather than a regression on
pseudo-values. There is one forest per arm (T-learner), with the default
composite Gray split rule (`splitrule = "logrankCR"`), so a single fit gives
both CIFs and all three estimands (RMSTc = 28 − RMTL1 − RMTL2). Ŵ and the
stage-2 forest are the grf ones the `pseudo_dr` arms use, so the outcome model
is the only thing that differs from `pseudo_dr_whole_*`.

- **Pseudo-values enter only the correction term**, as whole-sample values. The
  whole/cvps factor therefore has nothing to act on, and the two arms differ in
  fitting alone. `rsf_dr_oob` takes each unit's own-arm prediction from
  `cif.oob`/`predicted.oob`, with the same honesty argument as `t_learner_rf`.
  `rsf_dr_scf` is its crossfit twin, kept to check that OOB is adequate for
  rfsrc, since `crossfitting/` only showed it for grf.
- **Valid because censoring here is independent of X**, so E[θ | X, W] is the
  conditional RMTL that the forest estimates.
- **RMTL to the horizon.** For competing risks, rfsrc's `predicted` /
  `predicted.oob` is "expected life years lost due to cause j": the step CIF
  integrated over `time.interest`, but only up to the **last event time**.
  (Checked against `src/survival.c::getMortality` in 3.9.0.) τ can't be set,
  because `ntime` snaps to observed event times. So `rsf_fit` censors at the
  horizon and uses `ntime = 0`, and `rsf_rmtl` adds the tail
  (28 − t_T)·CIF_j(t_T). This matches a manual integration of `cif.oob` to
  1e-14. The package's own `get.rmst()` is not used: it is unexported, works
  for the `"surv"` family only, and takes the in-bag survival whenever the OOB
  one exists.
- **Single-cause guard.** If a training subset has only one event type, rfsrc
  silently becomes a plain survival forest, whose `predicted` is not years
  lost. `rsf_rmtl` then integrates 1 − S for that cause, sets the other to 0,
  and warns.
- **Threads.** The fold workers set `options(rf.cores = 1)`. The OOB fits run in
  the main process and follow `OMP_NUM_THREADS`, which `surv_analysis.R` sets
  to `grf_threads` (default 1) before starting its workers. The grf forests
  take `num.threads = grf_threads` directly.
- **Not bitwise reproducible**: the RSF arms differ run to run on unchanged
  code, so leave them out of exact equivalence checks.

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

`surv_analysis.R` sources `R/cate_models.R` before `surv_models.R`, so the
study's own definitions win where they still exist. That is what supplies
`t_learner_rf`, `stage2_whole_rf`, `stage_2_sl`, `dr_pseudo`, and (via
`R/sl_library.R`) `pretest_superlearner`, `sl_fit_predict` and `sl_libraries`.

`surv_metrics.R` does not use `compute_metrics()`: its results nest as
framework × target rather than one entry per model. It does use `cate_metrics()`,
so the bias convention matches the rest of the repo.
