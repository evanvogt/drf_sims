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

## Estimators

Five frameworks, each targeting restricted mean survival or lost time:

| framework | |
|---|---|
| `ipw` | causal forest on IPW-transformed RMST outcomes |
| `csf_cs` | causal survival forest, competing events treated as censoring (cause-specific) |
| `csf_sh` | causal survival forest, competing events kept in the risk set (subdistribution / Fine-Gray) |
| `pseudo_cf` | causal forest on pseudo-values |
| `pseudo_dr` | DR-learner on jackknife pseudo-values |

Plus SuperLearner T-learner variants (standard and split pseudo-obs,
Cwiling et al. 2025) and a DR-learner variant on standard pseudo-obs. The
DR-learner's split pseudo-observation construction is currently **disabled**
— see "Known issues" below.

**The truth column depends on the framework**, which is why
`surv_metrics.R` carries `framework_truth_map`: `ipw` and `csf_cs` remove
competing events so they target the cause-specific (net) RMST, while `csf_sh`
keeps them in the risk set and targets the subdistribution RMST. Scoring one
against the other's truth would be a category error.

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

`compute_split_pseudoyl()`/`compute_split_pseudomean()` (`surv_models.R`
~lines 929–965) hit this every time a validation observation's `Y` exceeds
the max `Y` already in that fold's KM set — confirmed across every
scenario/censoring/seed combination tested, so **this is not a small-sample
artefact**: it reproduces reliably at n = 500 with 10 folds, the production
settings. Unlike `pseudo_crossfit()`/`pseudo_double_crossfit()`, which
already guard against the same underlying `pseudoyl()`/`pseudomean()` NA
behaviour with an `is.na()` → `pseudo_whole` fallback, the split-pseudo-obs
construction has no such guard, so the NA reaches `stage_2_sl_vec()`'s
`SuperLearner()` call and aborts the run.

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

### No `pretest_superlearner`

This folder never calls `pretest_superlearner`. `surv_models.R` passes
`SL.library = sl_library` at eleven remaining call sites with no
validation, so one bad algorithm can still take the entire run down instead
of being dropped. The `continuous`/`binary` studies pretest and (before bug
F was fixed) threw the result away; this study does not pretest at all.
Orthogonal to the two issues above and still unresolved.

Consequences:

- there is **no regression baseline** for this study, so it is excluded from the
  default `R/regression_check.R` sweep (`deferred = TRUE`) — disabling the
  split-DR-learner branch reduces crash risk but hasn't been verified to make
  the full pipeline run clean end-to-end (see the T-learner risk above), so
  the deferred status is not lifted yet
- wiring `surv_models.R` into `R/` is on hold until it runs — the obvious first
  step is dropping its local `stage_2_rf` (lines ~890–917), which is
  `R/cate_models.R`'s matrix branch verbatim

To work on it:

```bash
Rscript R/regression_check.R baseline competing_risk
```

## Files

`surv_config.R` (grid), `surv_dgm.R`, `surv_models.R` (1,395 lines — the largest
in the repo), `surv_analysis.R`, `surv_check.R`, `surv_collect.R`,
`surv_metrics.R`. `scratch_dgm_params_check.R` is exploratory.
`surv_dr_split_na_diagnose.R` reproduces and traces the split-DR-learner
NA bug documented in "Known issues" above.

`surv_metrics.R` does not use `compute_metrics()`: its results nest as
framework × target rather than one entry per model. It does use `cate_metrics()`,
so the bias convention matches the rest of the repo.
