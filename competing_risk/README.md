# Competing risks

CATE estimation when event 1 can be pre-empted by a competing event 2. This
is the study the rest of the repo builds towards, and the only one whose estimators are genuinely different rather than a configuration of the
shared ones.

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
Rscript ../surv_check.R      # writes failed_ids.txt, points surv_rerun.sh at it
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
`surv_1.sh`'s `ncpus=2`). See `surv_run.R`'s header for why 4 and not 8, and
`validation/continuous/cts_val_run.R` for the pattern this follows.

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

**Parameters retuned 2026-10-01.** The event 1 baseline scale went from 15 to 25,
and the effects became log-HR-matched across events (HR 1.42 / 2.86 on
event 1, 0.70 / 0.35 on event 2, at X3 = 0 / 1). Under the old values the scenarios
held on the hazard scale but ran together on the RMTL scale. ADEMP.md "Why
these values" has the numbers. **Every result, failure count and smoke-test
figure below that predates 2026-10-01 comes from the old DGM**, and the full
array has to be rerun. `R/archive_old_results.R` skips this study, so archive
the old tree by hand before resubmitting. From `results/`:
`tar -cf _archive/pre_2026-10-01/competing_risk.tar competing_risk`, then
remove `competing_risk/` so `surv_check.R` and `surv_collect.R` don't mix the
two DGMs.

**Correlated covariates added 2026-10-01.** X1, X2, X3 and the noise X01–X03 (X4, X5 too since 2026-10-10)
now come from the copula `sample_size/correlated/` uses
(`correlated_covariates()` in `R/dgm_scenarios.R`), at ρ ∈ {0, 0.5}. That
makes ρ a grid factor and a `rho_<ρ>/` level in the results path. Key points:

- ρ = 0.5 is the primary analysis and ρ = 0 (independent) the sensitivity
  analysis, since 2026-10-06 (ρ = 0 was primary until then). ρ = 0 has the
  same distribution as before 2026-10-01 but a new draw order, so it is not
  paired with any older run.
- τ(x) is identical at both ρ, and the scenario labels still separate on the
  RMTL scale at ρ = 0.5. `Rscript surv_dgm_check.R` gives the numbers (ADEMP
  "What the DGM implies").
- The two ρ share each run's seed, so compare them with paired per-run
  differences. `surv_results.qmd`, "Sensitivity analysis", does this.
- Every other section of the report is ρ = 0.5.

**X4 / X5 kept and IPW fixed, 2026-10-10.** Two changes, both needing a full rerun:

- The dataset now carries the copula's X4 and X5, which used to be drawn and
  then dropped. Estimators see the same 10 covariates as `sample_size/` (X1–X5, X01–X05).
  X4 and X5 are noise here: they enter no hazard. No extra random draws, so each
  run's W, Y, D and truth are unchanged. Only the estimators' inputs change.
- `get_ipw()` used to drop every unit censored or with the competing event, even
  after the horizon. Now a unit with Y ≥ 28 is kept, as in the grf tutorial.
  This changes `ipw` (badly biased before, event 2 worst) and `csf_sh` when
  `censoring = TRUE`.

Archive the earlier results the same way as above (`_archive/pre_2026-10-10/`).

### Scenario 5 flips the sign of the RMTL1 CATE — expect it

On the smoke-tested scenario-5 replicate (old DGM, n = 500, censoring on), **almost every
arm correlates *negatively* with `tau_RMTL1`**, between −0.05 and −0.29, while
recovering `tau_RMTL2` correctly and positively (+0.10 to +0.42). This is a
property of the scenario, not a bug, and three things pin that down:

- **The cause-specific truth is flat.** `tau_RMST1_cs` has zero variance in
  scenario 5 — the treatment has no cause-specific effect on event 1 at all.
  Every bit of `tau_RMTL1` heterogeneity is the competing-risk artefact: the
  treatment changes who is still at risk for event 1, nothing else.
- **The two truths are near-mirror images.** `cor(tau_RMTL1, tau_RMTL2) = −0.93`
  in the truth itself (−0.91 under the current parameters), because the RMTL1 heterogeneity is mechanically induced by
  the RMTL2 heterogeneity. So an estimator that tracks the real (event-2) signal
  correctly lands anti-correlated with the *induced* RMTL1 estimand.
- **It is consistent across arms**, including `sl_t_split`, whose internals were
  not touched by the crossfitting change. A crossfitting artefact would not show
  up uniformly across whole-sample OOB, single-crossfit and split-pseudo arms
  alike.

This is the phenomenon scenarios 2/5/6 exist to expose, so it belongs in the
write-up rather than being chased as a defect. **Two caveats before quoting the
numbers**: they come from a single replicate, and the correlations are weak
enough that per-replicate noise at n = 500 is a real part of them — confirm the
pattern across the full 100 runs before making anything of it.

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
pseudo-observation construction (Cwiling et al. 2025) is currently **disabled**
— see "Known issues" below.

The two pseudo-value frameworks and the SuperLearner ones each ship in
**several arms**, because this study also asks how the pseudo-values should be
built — see "Crossfitting" below.

**The truth column depends on the framework**, which is why
`surv_metrics.R` carries `framework_truth_map`: `ipw` and `csf_cs` remove
competing events so they target the cause-specific (net) RMST, while `csf_sh`
keeps them in the risk set and targets the subdistribution RMST. Scoring one
against the other's truth would be a category error.

## Crossfitting

This study used to double-crossfit every nuisance over all `C(V,2)` fold pairs.
It no longer does. The production strategy now matches `R/cate_models.R` — see
that file's header comment and `crossfitting/README.md` for the benchmark the
choice rests on:

| arm family | strategy | |
|---|---|---|
| `ipw`, `csf_cs`, `csf_sh`, `pseudo_cf_whole_oob` | grf's own internal crossfitting, **`cf_default`** | plain `causal_forest`/`causal_survival_forest`, OOB `tau` |
| `pseudo_dr_whole_oob` | whole-sample OOB, T-learner (one forest per arm), **`oob_oob`** | `nuisance_pseudo_rf_oob` + `stage2_whole_rf` |
| the SuperLearner arms | single leave-one-fold-out, **`scf_scf`**; DR outcome model per arm | `nuisance_pseudo_sl` + `stage_2_sl`, sharing one `fold_indices` |
| `rsf_dr_oob` | whole-sample OOB, T-learner (one competing-risks RSF per arm) | `nuisance_rsf_oob` + `stage2_whole_rf` |
| `rsf_dr_scf` | single leave-one-fold-out, the OOB arm's crossfit twin | `nuisance_rsf_scf` + `stage_2_rf_scf` |

Every DR-learner outcome model (pseudo-value on X) is fit separately in each
arm - a T-learner, as in `R/cate_models.R`. Until 2026-09-27 they were
S-learners on `cbind(W, X)`, and `pseudo_dr_whole_oob` was the `oob_oob_s`
arm. The disabled split-pseudo DR arm (`nuisance_pseudo_sl_split`) followed on
2026-10-01, along with a trimmed propensity. SuperLearner libraries are per
nuisance, from `R/sl_library.R`.

### The pseudo-value comparison

Because the outcome here is a jackknife pseudo-value rather than an observed
`Y`, "which rows the nuisance model saw" and "which rows the outcome was
computed from" are two independent choices. The study now varies both, crossed
as far as they can be:

- **pseudo-values** — `whole` (one Aalen–Johansen / Kaplan–Meier fit on all `n`,
  `pseudo_all`) vs `cvps` (leave-one-fold-out, `pseudo_crossfit`)
- **fitting** — `oob` (whole sample, grf-internal) vs `scf` (single crossfit)

`cvps + oob` is not a cell: column `k` of the crossfit matrix is `NA` on fold
`k`'s own rows, so it only exists inside a fold loop. That is why the
random-forest families get three arms rather than four, with `whole_scf` as the
control that separates the pseudo-value factor from the fitting factor — the
role `scf_oob_t` plays in `crossfitting/`.

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
currently broken (see "Known issues"). **So every DR arm keeps whole-sample
pseudo-values in its correction term**, as the code always did, and the factor
being varied is only *the pseudo-values used to train the nuisance
regressions*. Do not report it as more than that.

A second qualification: `pseudo_crossfit` substitutes the whole-sample value
wherever `pseudoyl`/`pseudomean` returns `NA`, which leaks the held-out fold
into the `cvps` arms. `results$pseudos$cf_cv$n_na_fallback` counts how often
that fired, per estimand, so the comparison can be qualified rather than
assumed clean. It was `0` on the first smoke-tested replicate (scenario 1,
n = 500, V = 10, censoring on), but it is seed-dependent — check it rather
than trusting that.

### The RSF DR-learner (`rsf_dr_oob`, `rsf_dr_scf`)

Added 2026-10-01. The outcome model is a **competing-risks random survival
forest** (randomForestSRC) fit to the observed `(Y, D)`, rather than a
regression on pseudo-values. There is one forest per arm (T-learner), with the
default composite Gray split rule (`splitrule = "logrankCR"`), so a single fit
gives both CIFs and all three estimands (RMSTc = 28 − RMTL1 − RMTL2). Ŵ and the
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
- **Needs `randomForestSRC` in `sim-env`.**

### Two behaviour changes that came with this

- Propensities are now clamped by `R/utils.R::trim_ps` in every DR arm. This
  study previously did not clamp at all, so the DR denominator
  `W.hat * (1 - W.hat)` could explode.
- The SuperLearner arms now pretest their library (see "Known issues" below).

## Known issues

### 19 runs still fail — two routes the round-1 fix left open

> **Old DGM.** Everything in this section was measured under the
> pre-2026-10-01 parameters. Route C needed a run where nobody was observed
> past 28. Under the current parameters about 30 controls and at least about 8
> treated units are past 28 in every cell, so route C should not recur. Route
> A' needed a near-constant treated-arm RMTL2 cell, which is now much less
> likely (at least about 13% treated E2s, against 3% before). Neither code
> path is guarded, so both remain possible in principle. Confirm against the
> new array's `failed_ids.txt`.

**Open.** After the split-T-learner fix below, a full rerun left 398 missing
runs, and a rerun of those at 4h/20gb left **19**. All 19 were in the original
225. They are deterministic in the array index — the same 19 came back with 4x
the walltime and 10x the memory — so this is code, not resources.

They are **not** the bug that fix addressed. Both routes below were reproduced
end to end locally with `Rscript surv_analysis.R <index>`, whose error handler
prints the call stack; `surv_failed_diagnose.R` stages 2 and 3 separate them and
score them against 100 successful controls. Artefacts in
`<results>/competing_risk/diagnostics/round2/`.

| | count | dies in | error |
|---|---|---|---|
| whole-sample pseudo-value NA | 15 | `pseudo_cf_whole_oob` | `The vector of observations (W, Y, Z or D) contains at least one NA.` |
| pretest's own fallback dropped | 4 | `pseudo_sl_t_split` | `All algorithms dropped from library` |

Both discriminate perfectly: **0 of 19 unexplained, 0 of 100 controls flagged.**

#### Route C — the whole-sample pseudo-values are NA, and nothing guards them

15 indices, all `censoring = TRUE` (runs 22, 68, 86 x scenarios 1, 3, 4, 6, 7).

The fix below guarded `pseudo_sl_t_split`'s *training-fold* pseudo-values.
`pseudo_all()` itself was never guarded — it is the only pseudo-value producer in
`surv_models.R` without an `is.na()` fallback — and `pseudo_crossfit()` fills its
own NAs *from* `pseudo_all()`'s vector, so when that vector is NA the fallback
substitutes an NA and the `cvps` arms inherit it too.

The driver is **`at_risk_past_horizon == 0`**: when nobody is observed past the
horizon of 28, `pseudoyl()`'s leave-one-out risk set empties at or below `tmax`,
so the `pseudo:::ci.omit` bug already written up below reaches the estimand
instead of an unused tail. On index 295, `max(Y) = 27.83` and `pseudoyl()`
returns one NA on each of cause 1 and cause 2. It is a clean split: **15 of 19
failures have nobody past the horizon, against 0 of 100 controls.**

`pseudomean()` is unaffected, and that is the internal control — stage 3 on index
295 errors in **20 arms, every one of them RMTL1 or RMTL2, and none of the RMSTc
or `csf_*`/`ipw` arms**. The failure tracks the estimand whose pseudo-values are
NA, exactly. `sl_t_split` is among the survivors: its `keep <- !is.na(y)`
backstop does its job here.

The run dies at `pseudo_cf_whole_oob`, the *first* pseudo-value arm — so for
these 15 the SuperLearner arms blamed in round 1 are never even reached.

#### Route A' — `pretest_superlearner`'s SL.mean fallback is not safe either

4 indices, all `censoring = FALSE` (84, 462, 923, 924).

Round 1's degenerate-library *condition*, but not its failure, and the difference
is the point: these four are more extreme than anything round 1 saw —
`min_n_unique` is **1**, a literally constant treated-arm RMTL2 cell, against 3
on index 67. On a constant `y`, `pretest_superlearner()` drops all four
candidates and takes its own "every candidate failed/warned" branch
(`R/cate_models.R:415-425`), returning `SL.mean`. The live `SuperLearner()` call
then drops `SL.mean` too and errors before there is any fit at all:

```
Removed libraries due to NA/error: SL.glm SL.glmnet SL.ranger SL.gam
Error in SuperLearner(): All algorithms dropped from library
```

Reproduced in isolation on a constant vector, so it needs none of this study's
data. **`onlySL = TRUE` is irrelevant here** — the failure is at fit time, not
predict time — so the round-1 fix could not have caught it, and it is live in
every SuperLearner call site, not just the ones that still omit `onlySL`. Stage 2
measures that directly: `onlysl_would_fix = 0` of 4.

Stage 3 on index 923 errors in **exactly one arm**, `sl_t_split`. The other
SuperLearner arms survive because they train on `V-1` folds; the split arm's
3-way split gives the most degenerate cell.

#### Fixing it needs two decisions, not one

- **Route C.** There is no `pseudo_whole`-style fallback available *for* the
  whole sample, so this is the same design decision the disabled split-DR arm
  faces: drop the affected observation, or impute it. Note the estimand is
  genuinely not identified past the last observed time when
  `at_risk_past_horizon == 0`, so dropping silently would be worse than it looks
  — whatever is chosen needs an `n_na_fallback`-style counter like the others.
- **Route A'.** `pretest_superlearner` should not return a library it has not
  shown to survive the *caller's* `cvControl`, and a constant-outcome cell has no
  meaningful CATE contribution anyway. Both point at handling the degenerate cell
  explicitly rather than at another library-level patch.

#### The 398 were probably walltime, and that matters for the next submission

The intermediate batch is a separate story: 318 of its 398 were *new* indices,
spread across all seven scenarios **including 2 and 5**, which neither route
above can touch (`bW_1 = 0`). The likely cause is the same commit trimming
`surv_1.sh` from 2h to 1h while `pretest_superlearner` made every SuperLearner
arm markedly more expensive. The 4h rerun cleared 379 of them, which fits. This
is inference from the failure pattern, not proof — the logs that would settle it
were deleted. **Do not re-submit the full array at 1h.**

#### The logs were deleted, and that cost most of this

`surv_analysis.R`'s error handler writes the message and call stack to **stdout**
precisely so a lost `.e` file does not lose the cause. Both the `.e` and the `.o`
files for this rerun were deleted before anyone read them, so stage 4 had no
input and the whole diagnosis had to be rebuilt from local reproduction under
R 4.5.3 rather than the cluster's 4.3.2. Keep `logs_rerun/` next time.

Since then the jobscripts use `#PBS -j oe`, so there is no separate `.e` file:
stderr (R's own error line and `all_cate_surv_models()`'s `message()`s) is
merged into the `.o` file alongside the handler's call stack.

### SuperLearner DR-learner, split pseudo-obs (disabled)

`all_cate_surv_models()` used to abort with:

```
SuperLearner(): missing data is currently not supported. Check Y, X, and newX
```

This is now root-caused, and the offending branch (`results$sl_dr_split`,
`surv_models.R` ~lines 307–339) is **commented out** until it's fixed
properly — `nuisance_pseudo_sl_split()`/`stage_2_sl_vec()` themselves are
still defined, just not wired into `all_cate_surv_models()`.

Root cause: it's a structural bug in the `pseudo` package, not a data
problem. `pseudoyl()`'s internal `pseudo:::ci.omit()` computes leave-one-out
cause-specific cumulative incidence for every individual at once via a
shared risk-set matrix. For whichever individual holds the **maximum
observed time** in the sample being decomposed, the leave-one-out risk set
can hit exactly zero at the tail, producing `0/0 = NaN` in two places: the
shared KM survival term (`km`) and the cause-specific numerator ratio
(`cumi`). The package patches only `km` (forward-fills its trailing NaNs);
it never patches the matching NaNs in `cumi`. Since the cause-specific
incidence is `cumsum(cumi * km)`, that one unpatched `NaN` poisons every
later entry, and `pseudoyl()` returns `NA` for exactly the max-time
individual.

`compute_split_pseudoyl()`/`compute_split_pseudomean()` hit this every time a
validation observation's `Y` exceeds the max `Y` already in that fold's KM set
— confirmed across every scenario/censoring/seed combination tested, so **this
is not a small-sample artefact**: it reproduces reliably at n = 500 with 10
folds, the production settings. Unlike `pseudo_crossfit()`, which already
guards against the same underlying `pseudoyl()`/`pseudomean()` NA behaviour
with an `is.na()` → `pseudo_whole` fallback, the split-pseudo-obs construction
has no such guard, so the NA reaches `stage_2_sl_vec()`'s `SuperLearner()` call
and aborts the run.

The crossfitting change above did **not** touch this. The split-pseudo arms use
a 3-way train/KM/validation split that is the estimator's definition rather
than a crossfitting choice, so they were deliberately left alone.

Diagnosed and reproduced end-to-end in `surv_dr_split_na_diagnose.R`
(confirmed on R 4.5.3 / `pseudo` 1.4.3 / `SuperLearner` 2.0.29; the cluster
runs R 4.3.2 with different package versions, so behaviour there is
unverified). Fixing it needs a design decision, since there's no
`pseudo_whole`-style equivalent available for the split construction today:
either extend the fallback pattern to the split pseudo-obs (e.g. falling
back to the standard training-fold pseudo-value, or to `pseudo_whole`), or
drop the affected observation from that fold's stage-2 fit.

### T-learner split pseudo-obs — same latent risk as the SuperLearner branch; confirmed and now guarded. Fixed.

### `pretest_superlearner` is now used in this study too — resolved by the crossfitting change.

To work on it:

```bash
Rscript R/regression_check.R baseline competing_risk
```

## Files

`surv_config.R` (grid), `surv_dgm.R`, `surv_models.R`, `surv_analysis.R`,
`surv_run.R` (the no-queue RStudio-session alternative, see "Without the
queue" above), `surv_check.R`, `surv_collect.R`, `surv_metrics.R`.
`surv_dgm_check.R` produces ADEMP's "What the DGM implies" numbers (event mix,
population truths, X3 separation) for each ρ.
`scratch_dgm_params_check.R` is exploratory.
`single_event/` is a separate, smaller study: event 1 with event 2 removed,
reusing this study's arms (see its README).
`surv_dr_split_na_diagnose.R` reproduces and traces the split-DR-learner
NA bug documented in "Known issues" above.
`surv_failed_diagnose.R` is the diagnostic that found both the 225-run failure
and the 19-run one: it decodes `jobscripts/failed_ids.txt` against the grid
(stage 1), probes every failed index against a control sample for each known
failure route (stage 2), reruns each arm of `all_cate_surv_models()` separately
and keeps going past a failure so one run reports all of them (stage 3), and
triages PBS `.o` files (stage 4). Point it at the next batch of missing runs
rather than starting from scratch:

```bash
Rscript surv_failed_diagnose.R --out=<results>/diagnostics/round3
Rscript surv_failed_diagnose.R --stage=3 --ids=295,923 --out=...
```

**Pass `--out`.** Stage 2 and stage 3 overwrite `failed_probe.csv` and
`failed_arms.csv` in place, so a fresh run silently destroys the previous round's
evidence — which is how round 1's 325-row probe was lost (its `stage2.log` still
carries the per-index lines, the CSV does not). Round 2's artefacts are under
`diagnostics/round2/`.

Its stage-2 mode labels are versioned against the code, not against history: the
round-1 labels would have relabelled route C as "the bug we already fixed" and
hidden it. If a future fix closes a route, retire its label rather than leaving
it able to claim a failure it can no longer cause.

`surv_analysis.R` sources `R/cate_models.R` before `surv_models.R`, so the
study's own definitions win where they still exist. That is what supplies
`t_learner_rf`, `stage2_whole_rf`, `stage_2_sl`, `dr_pseudo`, and (via
`R/sl_library.R`) `pretest_superlearner`, `sl_fit_predict` and `sl_libraries`.

`surv_metrics.R` does not use `compute_metrics()`: its results nest as
framework × target rather than one entry per model. It does use `cate_metrics()`,
so the bias convention matches the rest of the repo.
