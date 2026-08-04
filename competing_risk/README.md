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

Plus SuperLearner T-learner and DR-learner variants, including the split
pseudo-observation construction of Cwiling et al. (2025).

**The truth column depends on the framework**, which is why
`surv_metrics.R` carries `framework_truth_map`: `ipw` and `csf_cs` remove
competing events so they target the cause-specific (net) RMST, while `csf_sh`
keeps them in the risk set and targets the subdistribution RMST. Scoring one
against the other's truth would be a category error.

## ⚠ Currently broken

`all_cate_surv_models()` aborts with:

```
SuperLearner(): missing data is currently not supported. Check Y, X, and newX
```

An NA is reaching a stage-2 SuperLearner fit — most likely a leave-fold-out
pseudo-value matrix, which carries NAs by construction. **This is not a
small-sample artefact**: it reproduces at n = 500 with both 5 and 10 folds, the
production settings.

It is made worse by a second issue: **this folder never calls
`pretest_superlearner`**. `surv_models.R` passes `SL.library = sl_library` at
twelve call sites with no validation, so one bad algorithm takes the entire run
down instead of being dropped. The `continuous`/`binary` studies pretest and
(before bug F was fixed) threw the result away; this study does not pretest at
all.

Caveat: reproduced locally on R 4.5.3, while the cluster runs 4.3.2 with
different package versions, so cluster behaviour may differ.

Consequences:

- there is **no regression baseline** for this study, so it is excluded from the
  default `R/regression_check.R` sweep (`deferred = TRUE`)
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

`surv_metrics.R` does not use `compute_metrics()`: its results nest as
framework × target rather than one entry per model. It does use `cate_metrics()`,
so the bias convention matches the rest of the repo.
