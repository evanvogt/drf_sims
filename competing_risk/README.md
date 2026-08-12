# Competing risks

CATE estimation when the event of interest can be pre-empted by a competing
event. This is the study the rest of the repo builds towards, and the only one
whose estimators are genuinely different rather than a configuration of the
shared ones.

| | |
|---|---|
| scenarios | 1–7 |
| n | 500 |
| censoring | TRUE, FALSE |
| runs | 100 |
| array | **1,400 jobs** |
| horizon | 28 |
| results | `../results/competing_risk/scenario_<k>/<n>/censor_<TRUE\|FALSE>/` |

## Design

Event of interest (EOI) and competing event (CE) are generated from Weibull
distributions via joint hazards (Beyersmann 2009). Scenarios vary where the
treatment acts:

| # | |
|---|---|
| 1 | ATE on EOI only |
| 2 | ATE on CE only |
| 3 | HTE on EOI, no ATE on CE |
| 4 | HTE on EOI, ATE on CE |
| 5 | HTE on CE, no ATE on EOI |
| 6 | HTE on CE, ATE on EOI |
| 7 | HTE on both |

Scenarios 2, 5 and 6 are the interesting ones: the treatment moves the competing
event, so a naive analysis of the EOI can show an apparent effect that is really
a change in who survives long enough to have one.

### Scenario 5 flips the sign of the RMTL1 CATE — expect it

On the smoke-tested scenario-5 replicate (n = 500, censoring on), **almost every
arm correlates *negatively* with `tau_RMTL1`**, between −0.05 and −0.29, while
recovering `tau_RMTL2` correctly and positively (+0.10 to +0.42). This is a
property of the scenario, not a bug, and three things pin that down:

- **The cause-specific truth is flat.** `tau_RMST1_cs` has zero variance in
  scenario 5 — the treatment has no cause-specific effect on the EOI at all.
  Every bit of `tau_RMTL1` heterogeneity is the competing-risk artefact: the
  treatment changes who is still at risk for event 1, nothing else.
- **The two truths are near-mirror images.** `cor(tau_RMTL1, tau_RMTL2) = −0.93`
  in the truth itself, because the RMTL1 heterogeneity is mechanically induced by
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
| `pseudo_dr_whole_oob` | whole-sample OOB, S-learner, **`oob_oob_s`** | `nuisance_pseudo_rf_oob` + `stage2_whole_rf` |
| the SuperLearner arms | single leave-one-fold-out, **`scf_scf`** | `nuisance_pseudo_sl` + `stage_2_sl`, sharing one `fold_indices` |

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

### Two behaviour changes that came with this

- Propensities are now clamped by `R/utils.R::trim_ps` in every DR arm. This
  study previously did not clamp at all, so the DR denominator
  `W.hat * (1 - W.hat)` could explode.
- The SuperLearner arms now pretest their library (see "Known issues" below).

## Known issues

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

### T-learner split pseudo-obs — same latent risk, unverified

`pseudo_sl_t_split()` (feeding `results$sl_t_split`) computes its
training-fold pseudo-values with the identical unguarded
`pseudoyl()`/`pseudomean()` call as the DR-learner branch above (same
3-way fold split). `surv_dr_split_na_diagnose.R`'s breadth-check already
found scenario/seed combinations where that exact call produces NAs, so
this branch is likely exposed to the same crash under different conditions
— it just hasn't been directly observed to abort yet, and is left running.
If `all_cate_surv_models()` still aborts with the same
`SuperLearner(): missing data` error after the split-DR-learner fix above,
this is the next place to look.

### ~~No `pretest_superlearner`~~ — fixed

**Resolved** by the crossfitting change. `surv_analysis.R` now sources
`R/cate_models.R`, and every `SuperLearner()` call in the non-split arms takes
its library from `pretest_superlearner()` (`gaussian()` for the pseudo-value
fits, `binomial()` for the propensity). The stage-2 regression gets it too, via
`R/cate_models.R::stage_2_sl`. On the first smoke-tested replicate this
immediately dropped `SL.gam`, which had been failing silently.

The `pseudo_sl_t_split()` branch is the exception: it is out of scope for the
crossfitting change and still passes `sl_library` unvalidated.

Consequences:

- there is still **no regression baseline** for this study, so it stays excluded
  from the default `R/regression_check.R` sweep (`deferred = TRUE`). The full
  pipeline now does run clean end-to-end on the smoke-tested indices, so
  capturing one is finally possible — the deferred status is a decision to
  make, not a blocker
- wiring `surv_models.R` into `R/` is partly done: the local `stage_2_rf` and
  `stage_2_sl` are gone, replaced by `R/cate_models.R`'s `stage2_whole_rf` and
  `stage_2_sl`. What remains local is genuinely study-specific (the pseudo-value
  nuisances, the fold-wise `stage_2_rf_scf`, and the split-pseudo arms)

To work on it:

```bash
Rscript R/regression_check.R baseline competing_risk
```

## Files

`surv_config.R` (grid), `surv_dgm.R`, `surv_models.R`, `surv_analysis.R`,
`surv_check.R`, `surv_collect.R`, `surv_metrics.R`.
`scratch_dgm_params_check.R` is exploratory.
`surv_dr_split_na_diagnose.R` reproduces and traces the split-DR-learner
NA bug documented in "Known issues" above.

`surv_analysis.R` sources `R/cate_models.R` before `surv_models.R`, so the
study's own definitions win where they still exist. That is what supplies
`oob_predict_counterfactual`, `stage2_whole_rf`, `stage_2_sl`,
`pretest_superlearner` and `dr_pseudo`.

`surv_metrics.R` does not use `compute_metrics()`: its results nest as
framework × target rather than one entry per model. It does use `cate_metrics()`,
so the bias convention matches the rest of the repo.
