# Continuous outcome — sample size study

How well do doubly-robust and forest-based CATE estimators recover heterogeneous
treatment effects as sample size grows, for a continuous outcome?

This is the reference study: the binary, missing-data and confidence-interval
studies are all variations on it.

## Design

Ten scenarios varying the structure of the CATE, crossed with four sample sizes,
100 runs each — **4,000 array jobs** — plus runs 101–500 for scenarios 1–4,
another 6,400 appended as grid rows 4001–10400.

| | |
|---|---|
| scenarios | 1–10 (see `R/dgm_scenarios.R`, `DESC_10`) |
| n | 100, 250, 500, 1000 |
| runs | 100; **500 for scenarios 1–4** |
| folds | 4 at n=100, 5 at n=250, else 10 |
| results | `../results/continuous/scenario_<k>/<n>/res_sim_<run>.RDS` |

Folds are reduced at small n because the double-crossfitting procedure fits
nuisances over all `C(V,2)` fold pairs — 45 fits at V=10 — and the training
sets become too small otherwise.

The SuperLearner libraries (one each for the propensity, outcome and CATE
models, `R/sl_library.R::sl_libraries`) also shrink at n=100. Until the
library change the single shared library dropped `SL.earth` and `SL.ranger`
there, which is why `dr_superlearner` is sometimes filtered out of the older
n=100 figures. The CATE library now keeps ranger at n=100, so it still has an
interaction-capable learner (see `../ADEMP.md`).

### Scenarios

| # | CATE structure |
|---|---|
| 1 | no HTE (ATE only) |
| 2 | simple, continuous variable (X4) |
| 3 | single effects + a different interaction |
| 4 | cosine |
| 5 | simple, binary variable (X3) |
| 6 | two variables, additive |
| 7 | continuous × binary interaction |
| 8 | single effects + interaction |
| 9 | continuous × continuous interaction |
| 10 | exponential |

Scenarios 1–4 are the ones the chapter reports: null, simple, complex and
non-linear. They were 1, 3, 8 and 9 before 2026-09-26; the root README has the
full old → new table.

Each dataset also carries five deliberately unrelated covariates (`X01`–`X05`),
so the estimators have to find the signal rather than being handed it.

### Outcome model and `bW` calibration

Every scenario shares one outcome model and differs only in its treatment
effect:

`Y = b0 + b1·X1 + b2·X2 + W·(bW + g(x)) + ε`, with `X1 ~ Bernoulli(0.4)`,
`X2 ~ N(0, 1)` and `ε ~ N(0, 0.5²)`.

| | |
|---|---|
| baseline | `b0 = 0.4`, `b1 = −0.5`, `b2 = 1`, in every scenario |
| `g(x)` | the scenario's heterogeneity term: its `te_expr` in `R/dgm_scenarios.R`, with `bW = 0` |
| `bW` | set so the true **ATE**, `bW + E[g]`, equals the effect the trial was planned to detect |

Each simulated RCT is planned the way trials usually are: to detect an ATE,
assuming the effect is the same for everyone. The planned effect δ gives 80%
power (`TARGET_POWER`) in an unadjusted two-sample t-test with n/2 per arm, using
the outcome SD with no heterogeneity: sqrt(b1²·0.4·0.6 + b2² + 0.5²) = 1.145.
`bW` is then set to −δ − E[g], so the true ATE is −δ. `te_moments()` computes
E[g] by quadrature, with no random draws, so the draw order is untouched.

| n | 100 | 250 | 500 | 1000 |
|---|---|---|---|---|
| true ATE, every scenario | −0.65 | −0.41 | −0.29 | −0.20 |

The plan gets the average effect right but knows nothing of the heterogeneity
around it, which makes the treated arm noisier by Var(g). So only scenario 1
realises the planned 80%; the others realise less, falling as heterogeneity
grows:

| scenario | 1, 10 | 4, 9 | 7 | 2, 5, 6, 8 | 3 |
|---|---|---|---|---|---|
| realised power | 0.78–0.81 | 0.75–0.77 | 0.69–0.71 | 0.65–0.69 | 0.61–0.63 |

**Before bug O** (root README), the baseline varied by scenario: b0 ran from 0.2
to 1, b1 was −0.05, and b2 was 1 in scenario 4 and 2 elsewhere. The calibration
used `sd = s_err + s2 = 1.5`, which adds SDs and ignores b1 and b2, and it set
`bW` rather than the ATE to the planned effect. The true ATE drifted by E[g]
(to a *positive* value in scenarios 3, 5 and 8), and power ran from 3% to 100%.
Each scenario's heterogeneity around its mean, g(x) − E[g], is unchanged; the
level of the true CATE moved, so `sign_acc` and the relative metrics are not
comparable with earlier results.

`Rscript R/calibration_report.R` prints `bW`, the true ATE, and the planned and
realised power for every continuous scenario at every n a study uses, including
the missing-data and validation studies, then the same for the binary scenarios
(see `binary/README.md`). It takes seconds and needs no
simulation. After any change to the scenario tables or `calibrate_bW()`, run it
on the cluster too: its `bW` tables should match a local run exactly, since
R 4.3.2 there could round a borderline value differently from 4.5.3 here.

## Estimators

`causal_forest`, `dr_random_forest`, `dr_oracle`, `dr_semi_oracle`,
`dr_superlearner` — all defined in `R/cate_models.R`, all using double
crossfitting for the nuisances. The oracle uses the true outcome model and a
known propensity of 0.5; the semi-oracle knows only the propensity.

`t_random_forest` and `t_superlearner` are T-learners, μ̂₁(x) − μ̂₀(x), built
from `dr_random_forest`'s and `dr_superlearner`'s own per-arm outcome models
(`nuisances_rf`, own-arm predictions OOB; `nuisances_sl`, out-of-fold). They
are not fitted: `Y1.hat` isn't saved, but the saved pseudo-outcome is
`(Y1.hat − Y0.hat) + (Y − Y.hat)(W − W.hat)/(W.hat(1 − W.hat))`, so
`R/metrics.R::add_t_learners()` subtracts the residual term at metrics time
and gets the contrast back exactly. It runs their BLP and `indep_cate` tests
there too (`BLP_p_os` via `blp_inputs()` as usual). They have no variance
estimate, and the n = 100 SuperLearner caveat above applies to
`t_superlearner` as to `dr_superlearner`. Only this study and `binary/` ask
for them (`compute_metrics(augment = add_t_learners)`); the CI and
missing-data studies don't.

## Files

| file | role |
|---|---|
| `cts_config.R` | the parameter grid and results path — **the** definition |
| `cts_dgms.R` | names this study's slice of `R/dgm_scenarios.R` |
| `cts_models.R` | `family = gaussian()`, `profile = "base"` |
| `cts_analysis.R` | array entry point; one row of the grid per index |
| `cts_check.R` | finds missing runs, writes `jobscripts/failed_ids.txt`, and updates `-J` and the resource request in the rerun jobscript |
| `cts_collect.R` | gathers per-run files into `cts_all.RDS` |
| `cts_metrics.R` | computes `cts_metrics.RDS`, plus `cts_true_cate_tests.RDS` — see below |
| `results_cts.R` | summaries |

### True-CATE HTE test evaluation

`cts_true_cate_tests.RDS` reruns the BLP and independence tests
(`run_true_cate_tests()`, `R/cate_models.R`) against the *true* CATE and true
nuisances (`truth$tau`, `truth$p0`, `W.hat = 0.5`) instead of an estimator's
fitted ones — one `BLP_p`/`BLP_p_os`/`indep_cate` row per (scenario, n, run), with no
per-model dimension, since nothing here is estimated. This isolates the
tests' own size/power from any estimator's error: scenario 1 is the null
(no heterogeneity), scenarios 2-10 the alternative.

Scenario 1's true CATE is constant, so neither test is defined there and
`BLP_p`, `BLP_p_os` and `indep_cate` are all `NA` for every scenario-1 run. The
true-CATE rows therefore have no null scenario; the estimators' scenario-1 rows
give the tests' size. Two details, both handled by `is_constant()`
(`R/cate_models.R`):

- `truth$tau` is computed as `p1 − p0 = (p0 + bW) − p0`, which leaves
  floating-point noise (SD ~1e-17) rather than an exact constant. Before the
  guard, `GenericML::BLP()` fitted beta.2 ≈ ±1e15 on that noise and returned a
  meaningless p-value; it was only `NA` when the noise happened to vanish.
- `coin::independence_test()` does not fail on a constant response: it warns
  ("zero diagonal elements") and returns p ≈ 0. Before the guard, the
  scenario-1 true-CATE `indep_cate` rejected in every run, and so did any
  estimator run whose CATE estimate collapsed to a constant.

### HTE tests

Every model's fit carries three post-estimation tests (`R/cate_models.R`,
"post-estimation heterogeneity tests"), read into the metrics by
`hte_test_metrics()` (`R/metrics.R`):

| Column | Test |
|---|---|
| `BLP_p` | `GenericML::BLP()`: `Y` on the baseline prediction `Y0.hat`, `(W − W.hat)` and `(W − W.hat)(τ̂ − mean τ̂)`, weighted by 1/(W.hat(1 − W.hat)). p-value for the interaction coefficient beta.2, **two-sided, homoskedastic SEs** (GenericML's default `vcovHC(type = "const")`), as computed when the models were fitted. |
| `BLP_p_os` | The same regression, **one-sided (H1: beta.2 > 0) with HC3 SEs**, recomputed at metrics time from each run's saved nuisances (`blp_inputs()`), so it needed no simulation re-run. |
| `indep_cate` | `coin::independence_test(τ̂ ~ X, teststat = "quadratic")`, asymptotic reference distribution: whether τ̂ is linearly associated with any covariate. |
| `indep_po` | The same test on the DR pseudo-outcome — the global heterogeneity test of [WATCH](https://github.com/Novartis/WATCH) (`src/04_explore_TEH.R`). |

Why two BLP columns. beta.2 = 1 means τ̂ is well calibrated and beta.2 = 0 that
it carries no heterogeneity signal; a τ̂ that ranks units in reverse gives
beta.2 < 0, which is not evidence of heterogeneity in the sense the test is
defined for (Chernozhukov et al.; `grf::test_calibration()` is one-sided for
the same reason) but gets the same two-sided p as a correct τ̂. And the
residual variance is not constant — for a binary outcome it is p(1 − p) by
construction — so the homoskedastic SEs are not the right ones. `BLP_p` is
kept unchanged so earlier figures stay reproducible.

Properties of the independence tests, evaluated as they are rather than changed:

- The quadratic statistic is built from linear correlations between the
  response and each covariate, so it has essentially no power against an
  effect with no linear component, however strong: scenario 4 (`cos(X4)`),
  scenario 9 (`X4·X5`; `tanh(X4)·tanh(X5)` for binary) and the `exp(−|X4|)`
  part of scenario 10 are uncorrelated with every covariate. Conversely, for a
  true CATE that is linear in X the statistic saturates (χ² = n − 1), so the
  true-CATE `indep_cate` shows whether a scenario's CATE has a linear
  component, not a power ceiling.
- WATCH applies the test to the pseudo-outcome, which has independent noise
  per unit. `indep_cate` applies it to τ̂, a fitted function of X, whose
  values are not exchangeable across units even with no heterogeneity, so its
  scenario-1 rejection rate is not guaranteed to be 5%.
- `causal_forest` and `dr_random_forest` share `nuisances_rf$po`, so their
  `indep_po` would be the same number; it is reported once, under
  `dr_random_forest`, and `NA` for `causal_forest`. Likewise `NA` for
  `t_random_forest` and `t_superlearner`, whose pseudo-outcomes are
  `dr_random_forest`'s and `dr_superlearner`'s.

Power and mean p-values leave out runs whose p-value is `NA` (a constant τ̂,
or a test that wasn't run); `n_na_BLP_p`, `n_na_BLP_p_os` and
`n_na_indep_cate` in the results summaries count them.

## Running it

```bash
qsub sample_size/continuous/jobscripts/cts_1.sh     # 1-4000
qsub sample_size/continuous/jobscripts/cts_extra.sh # 4001-10400: runs 101-500, scenarios 1-4
Rscript sample_size/continuous/cts_check.R          # writes failed_ids.txt if any are missing
qsub sample_size/continuous/jobscripts/cts_collect.sh
qsub sample_size/continuous/jobscripts/cts_metrics.sh
```

## Sizing the array job

`cts_1.sh` shipped with **unmeasured** `#PBS -l` lines: `ncpus=2:ompthreads=2` was
requested to match `cts_analysis.R`'s `workers <- 2`, but `num.threads` was never
set anywhere in `R/cate_models.R` (the shared estimator engine this study, and
six others, source), so grf always defaulted to using every visible core per
`multisession` worker regardless of what PBS allocated — the 2x2 pairing was a
guess, not a measurement. `R/cate_models.R`'s `cate_methods()` now accepts an
explicit `num.threads` argument (default `NULL`, so every other study using it
is unaffected), and `cts_analysis.R` forwards it from an optional trailing CLI
arg, the same pattern `crossfitting/cf_analysis.R` uses.

`cts_1.sh` now asks for one core and 2gb for an hour, and runs
`cts_analysis.R <i> 1 1` (one worker, one grf thread). Those figures are set by
hand. The `syrup` profiling sweep that was meant to measure them didn't work
for this study (see the root README's "Resource profiling (removed)"). When
changing them, size for the `n = 1000` cell, not a middle value: the array
spans `n ∈ {100, 250, 500, 1000}` under one `#PBS -l` line. Keep
`workers × grf_threads` within `ncpus`, and change the trailing args and the
`#PBS -l` line together.

## Status

**Archive the old results first** - `R/archive_old_results.R` (root `README.md`, Status, step 0). They predate the current DGM and use the pre-2026-09-26 scenario numbers, so running into that tree would mix old and new results.

**Re-runs required** — for the crossfitting strategy change to
`R/cate_models.R` (see root README Methods/Status), which moves all five
estimator arms, and separately for bug F, now fixed permanently, which changes
`dr_superlearner`: the second-stage SuperLearner library was pretested and
the result discarded, so failing algorithms were never dropped. Only the
`dr_superlearner` arm moves for bug F specifically — the harness can confirm
the other four are unchanged there.

**Also re-run for bug O** — the shared baseline and the `bW` calibration
changed (see "Outcome model and `bW` calibration" above), which moves every
dataset and the level of the true CATE in all ten scenarios.

**Also re-run for the per-arm outcome models and the SuperLearner libraries**
(2026-09-27). The DR-learners now fit the outcome model separately in each
treatment arm, which moves `dr_random_forest`, `dr_semi_oracle` and
`dr_superlearner` (`R/cate_models.R`). `dr_superlearner` also moves for its
per-nuisance libraries and the fixed pretest (bug Q; `R/sl_library.R`). A
smoke run at n = 1000 (index 34) on one core, as `jobscripts/cts_1.sh` runs it,
took 6.2 minutes (sharing the machine with other runs) against the 1h walltime. Results now carry `sl_dropped`: the
learners the pretest dropped, and any SuperLearner fit that failed and fell
back to the mean.

Nothing else in this study was affected by the *bug ledger*. `bias` also
changes sign when the metrics are regenerated (bug G), but that needs no
cluster time.
