# Sample-size studies

How well do doubly-robust and forest-based CATE estimators recover a
heterogeneous treatment effect as the trial gets smaller, for continuous and
binary outcomes? And do their confidence bands still cover?

Every study here uses one data-generating mechanism: ten treatment-effect
scenarios, with the covariates drawn from a Gaussian copula at latent
correlation ρ ∈ {0, 0.5}. ρ = 0.5 is the primary setting, where prognosis and
effect modification move together. ρ = 0 is the independent-covariate
comparison, paired run for run with ρ = 0.5. The full design is in
[`ADEMP.md`](ADEMP.md).

| study | folder | prefix | question | grid |
|---|---|---|---|---|
| continuous | `continuous/` | `cts_corr` | CATE estimation and HTE tests | scenarios 1–10 × n 100, 250, 500, 1000 × ρ; 500 runs for scenarios 1–4, 100 for 5–10 (20,800 rows) |
| binary | `binary/` | `bin_corr` | the same, binary outcome (risk-difference scale) | the same |
| CI, continuous | `confidence_intervals/continuous/` | `cts_corr_ci` | coverage of half-sample bootstrap bands across `CI_sf` | scenarios 1–4 × n 500, 1000 × ρ × 10 `CI_sf`; 100 runs (16,000 rows) |
| CI, binary | `confidence_intervals/binary/` | `bin_corr_ci` | the same, binary outcome | the same |
| optimal_sf | `confidence_intervals/optimal_sf/` | `cts_corr_ci_sf` | a data-driven choice of `CI_sf`, continuous outcome only | scenarios 1–4 × n 500, 1000 × ρ; 100 runs (1,600 rows) |
| SNR check | `snr/` | `bin_snr` | does the DR-learner catch the causal forest with more signal? | see its README |

Scenarios 1–4 (null, simple, complex, non-linear) are the main text's; 5–10 are
the appendix's. Results are written beside the repository, under
`../results/correlated/<outcome>/`, `../results/correlated/confidence_intervals/<outcome>/`
and `.../continuous/sf_calibration/` for optimal_sf.

## Data-generating mechanism

The implementation is `R/dgm_scenarios.R`. The scenario sets are
`continuous_corr_<ρ>` and `binary_corr_<ρ>` (`corr_set()`).

Both outcomes are one model, `E[Y | x, W] = m0(x) + W·τ(x)`, with
`τ(x) = bW + g(x)` and `W ~ Bernoulli(0.5)`. X1 and X2 are purely prognostic.
X3, X4 and X5 are the potential effect modifiers. All ten covariates (X1–X5,
X01–X05) are drawn in every scenario and every estimator sees all ten.

### Covariates

- A latent `Z ~ N(0, R)` over X1–X5 and X01–X03, with R exchangeable at ρ
  (`correlated_covariates()`).
- X1 and X3 are thresholded at prevalences 0.4 and 0.7. X2, X4, X5 and
  X01–X03 are the latent columns, N(0, 1).
- X01–X03 are noise for the outcome and the CATE, but correlated with X1–X5,
  so they are partial proxies for the modifiers.
- X04 and X05, indicators of a 3-level factor, are drawn independently.
- At ρ = 0.5 the observed correlations are about 0.5 between continuous
  covariates and 0.3–0.4 for pairs involving X1 or X3.

### Continuous outcome

`Y = 0.4 − 0.5·X1 + X2 + W·τ(x) + ε`, `ε ~ N(0, 0.5²)`.

### Binary outcome: the risk-difference scale

`Y ~ Bernoulli(m0(x) + W·τ(x))`, with control risk
`m0(x) = 0.34 + 0.36·plogis(−1.72 − 0.5·X1 + X2)`. The event is harmful and
treatment lowers its risk. Because τ adds to the risk, the true CATE is exactly
`bW + g(x)`: X1 and X2 are purely prognostic and scenario 1 is an exact null.

**Why the risk-difference scale.** The estimand is a risk difference. With
prognostic covariates an effect can be constant on at most one scale. Had the
effect been put on the logit scale, the link would make X1 and X2 modify the
risk difference in every scenario, scenario 1 would not be a null, and the
heterogeneity would change with n as bW is recalibrated.

**Why the control risk is bounded.** An additive effect needs room:
`m0(x) + τ(x)` must stay inside [0, 1] for every x. So m0 is a logistic scaled
into [0.34, 0.70] (E[m0] = 0.40, SD 0.048). The tight end is n = 100, where the
planned RD is −0.248 and every control risk must clear that plus the extra
benefit the HTE gives some patients. That leaves little room for m0 to spread,
so X1 and X2 are only weakly prognostic.

**The size of the HTE.** g is the continuous g with X4 and X5 entering through
`tanh`, its signs reversed, and scaled by `RD_SCALE_CORR[k]`:

- `RD_SCALE[k]` is the largest scale, floored to 3 dp, that keeps every treated
  risk inside [0.01, 0.99] (`RD_EPS`) at n = 100 to 1000. The floor at n = 100
  binds in every scenario but 4.
- Scenario 4 is held at 0.204, below the 0.209 the bounds would allow.
- `RD_SCALE_CORR` is `RD_SCALE` except scenario 9, which takes 0.139: under
  `RD_SCALE`, its lowest treated risk at ρ = 0.5, n = 100 would be −0.004. The
  same scale is used at both ρ, so the arms stay paired.
- At ρ = 0.5 two floors fall inside [0, 1] but below `RD_EPS`, and are accepted
  to keep the HTE at `RD_SCALE`'s size: scenario 3's (0.007) and 8's (0.003). Scales of 0.049 and 0.127 would restore them, at the cost of a
  smaller HTE.
- The signs are reversed so that each asymmetric scenario's larger swing goes
  towards *less* benefit, where the ceiling leaves room; only the floor is
  tight. With the continuous signs SD(τ) would be much smaller (at ρ = 0):

  | scenario | continuous signs | reversed |
  |---|---|---|
  | 3 | 0.034 | 0.050 |
  | 4 | 0.023 | 0.091 |
  | 5 | 0.053 | 0.126 |
  | 8 | 0.040 | 0.128 |
  | 10 | 0.044 | 0.083 |

  The price is that the opposite subgroup benefits more. It is always the
  X3 = 0 minority (30%) that benefits least; in scenarios 5, 8 and 10 it is
  harmed at n = 1000.

`Rscript R/bin_verify_hte.R` re-derives every scale and checks the floors.
Change `p0_lo`, `p0_hi`, b0–b2, `RD_EPS` or the studies' n, and the scales must
be recomputed.

### Calibration of bW

Each trial is planned the way trials usually are: to detect an ATE with 80%
power (`TARGET_POWER`), assuming the effect is the same for everyone
(`calibrate_bW()`). bW is then set so the true ATE equals the planned effect.

- **Continuous:** a two-sample t-test with n/2 per arm, at the outcome SD with
  no heterogeneity, `sqrt(b1²·0.4·0.6 + b2² + 2·b1·b2·Cov(X1, X2) + 0.5²)`.
  The plan ignores Var(g), so realised power falls below 0.80 as heterogeneity
  grows (0.61–0.63 in scenario 3 at ρ = 0).
- **Binary:** a two-proportion test at E[m0]. A Bernoulli outcome's variance is
  fixed by its marginal risk, so every scenario realises 0.80.
- `bW = δ − E[g]`, with E[g] and E[m0] taken by quadrature over the correlated
  latent (`te_moments()`, `baseline_grid()`), so no random draws are used.

At ρ = 0.5 bW and the ATE move through E[g] wherever g multiplies two modifiers
(scenarios 3, 7, 8, 9) and, for the planned SD, through Cov(X1, X2):

| | ρ = 0 | ρ = 0.5 |
|---|---|---|
| continuous true ATE, n = 100 / 1000 | −0.65 / −0.20 | −0.60 / −0.19 |
| binary true RD, n = 100 / 1000 | −0.248 / −0.085 | −0.247 / −0.085 |
| continuous SD(τ), scenario 3 | 1.16 | 1.36 |
| cor(m0, τ), continuous, scenarios 2 / 3 / 4 | 0 / 0 / 0 | −0.43 / 0.39 / 0.01 |
| cor(m0, τ), binary, scenarios 2 / 3 / 4 | 0 / 0 / 0 | 0.39 / −0.34 / 0.06 |

`Rscript R/calibration_report.R` prints bW, the true ATE, the power and (binary)
the treated-risk bounds for every scenario at every n. `Rscript
corr_truth_summary.R` prints SD(τ) and cor(m0, τ). Neither needs simulation
output. After any change to the scenario tables, run the report on the cluster
too: R 4.3.2 there could round a borderline bW differently.

## Pairing

Each run is seeded by its run index alone (`setup_rng_stream(run)`). The draw
order is W, the n × 8 latent block, the continuous noise, then X04/X05. So run
*r* at ρ = 0 and ρ = 0.5 shares W, the raw latent normals (multiplied by a
different Cholesky factor), the noise and X04/X05. Binary Y is coupled only
through `rbinom`'s shared uniforms; about 91–99% of outcomes agree. In the CI
studies run *r* is also the same dataset at every `CI_sf`.

How to compare:

- **ρ = 0.5 against ρ = 0:** paired per-run differences, MCSE = sd(diff)/√runs.
- **Estimators within one ρ:** against `dr_oracle` on the same run.
- **Across scenarios:** compute MCSEs per scenario, never pooled. Scenarios
  share seeds, so one fluctuation shows up in every scenario at once (per-run
  bias is about 0.9 correlated across scenarios). Sample sizes can be treated
  as independent: per-run errors are uncorrelated across n (|r| ≤ 0.02).
  These were checked on the independent-covariate design with
  `bin_bias_mc_checks.R`, at git tag `independent-ss-final`.

## Estimators

`causal_forest`, `dr_random_forest`, `dr_superlearner`, `dr_oracle` and
`dr_semi_oracle`, all from `cate_methods()` in `R/cate_models.R`. The oracle
knows the outcome mean and the propensity 0.5; the semi-oracle knows only the
propensity. `dr_superlearner` crossfits with 4, 5 and 10 folds at n = 100, 250
and ≥ 500, and its SuperLearner libraries are smaller at n = 100
(`sl_libraries()`, `R/sl_library.R`). The CI and optimal_sf studies run no
SuperLearner arm.

`t_random_forest` and `t_superlearner` are T-learners, μ̂₁(x) − μ̂₀(x), built
from the DR-learners' own per-arm outcome models. They are not fitted:
`add_t_learners()` (`R/metrics.R`) subtracts the residual term from the saved
pseudo-outcome at metrics time and gets the contrast back exactly. They have
no variance estimate.

## HTE tests

The continuous and binary studies run three tests on every model's fit
(`R/cate_models.R`), read into the metrics by `hte_test_metrics()`:

| column | test |
|---|---|
| `BLP_p` | `GenericML::BLP()`: Y on the baseline prediction, (W − ê) and (W − ê)(τ̂ − mean τ̂), weighted by 1/(ê(1 − ê)). p-value for β2, two-sided, homoskedastic SEs, as saved at estimation time |
| `BLP_p_os` | the same regression, one-sided (H1: β2 > 0) with HC3 SEs, recomputed at metrics time from the saved nuisances (`blp_inputs()`) |
| `indep_cate` | `coin::independence_test(τ̂ ~ X, teststat = "quadratic")`: is τ̂ linearly associated with any covariate? |
| `indep_po` | the same test on the DR pseudo-outcome (the global test of [WATCH](https://github.com/Novartis/WATCH)) |

**Why two BLP columns.** β2 = 1 means τ̂ is well calibrated and β2 = 0 that it
carries no signal. A τ̂ that ranks units in reverse gives β2 < 0, which is not
evidence of heterogeneity in the sense the test is defined for, but gets the
same two-sided p. And the residual variance is not constant (for a binary
outcome it is p(1 − p) by construction), so homoskedastic SEs are not the right
ones.

**The independence tests' blind spots.**

- The quadratic statistic is built from linear correlations, so it has almost
  no power against an effect with no linear component: scenario 4 (cos X4), 9
  (X4·X5) and the exp(−|X4|) part of 10. For a CATE linear in X it saturates
  (χ² = n − 1), so the true-CATE `indep_cate` shows whether a scenario has a
  linear component, not a power ceiling.
- `indep_cate` applies the test to τ̂, a fitted function of X whose values are
  not exchangeable across units even with no heterogeneity, so its scenario-1
  rejection rate is not guaranteed to be 5%.
- `causal_forest` shares `dr_random_forest`'s pseudo-outcome, so its
  `indep_po` is reported once, under `dr_random_forest`. Likewise the
  T-learners'.

**True-CATE tests.** `*_true_cate_tests.RDS` reruns the BLP and `indep_cate`
on the true CATE and true nuisances (`run_true_cate_tests()`), one row per
(ρ, scenario, n, run). This isolates the tests' own size and power from
estimation error. Scenario 1's true CATE is constant, so both tests are NA
there; the estimators' scenario-1 rows give the size. `is_constant()` catches
a constant up to floating-point noise: without it `GenericML::BLP()` fits β2
to that noise, and `coin` returns p ≈ 0 on a constant response instead of
failing.

Power and mean p-values leave out NA p-values (a constant τ̂, or a test that
wasn't run); `n_na_*` columns count them.

## Confidence intervals

**The method** (`R/bootstrap_ci.R`). For each of `CI_boot = 200` draws, refit
the second stage on a half sample with `sample.fraction = CI_sf` and form the
root `tau_full − tau_half`. Standardise each root by its bootstrap SD, take the
maximum absolute root over units within each draw, and use the 1 − α quantile
of that maximum as one critical value. That single value makes the band
simultaneous, and the absolute value makes it two-sided, so α = 0.05 gives a
95% band. `cf_half_boot` does this for the causal forest, `rf_half_boot` for
the DR-learners. The causal forest's own variance gives a pointwise normal
interval as well, scored as `causal_forest_inbuilt`.

`CI_sf` is swept from 0.05 to 0.5 because the right value is not known: too
small and the half-sample forests are noisy, too large and the halves stop
being independent enough. The CI studies use `profile = "ci"`: no SuperLearner
arm and no HTE tests, because the bootstrap dominates the runtime.

**Metrics** (`interval_metrics()`): marginal coverage (share of units their
own interval covers), simultaneous coverage (0/1: did the band cover every
unit) and mean length. Simultaneous coverage is the one the method controls.

### The query grid

The bands are also built over a fixed covariate grid (`build_query_grid()`),
scored as `<model>_grid` rows. It crosses whichever of X3 (0, 1), X4 and X5
(−2 to 2) the scenario uses, with every other covariate at 0. Because the grid
points are the same in every run, the grid rows also carry bias-eliminated
coverage: coverage of the across-run mean estimate rather than the true τ, so
bias drops out and only the width is judged.

At ρ = 0.5 a point such as (X4 = 2, X5 = −2) is about 4 SD from the centre of
the correlated latent, against 2 SD at ρ = 0, so the forests extrapolate there.
A fall in grid coverage at ρ = 0.5 can be extrapolation rather than
miscalibration. **The per-unit band is the primary ρ comparison.**

### optimal_sf

`find_optimal_sf()` picks `CI_sf` per dataset for `dr_random_forest`, using the
run's own estimate as a plug-in truth. For each candidate `CI_sf`, 50 times:
resample the pseudo-outcome residuals around τ̂, re-estimate the CATE, build a
band (`CI_boot = 100`) and record its coverage of τ̂. The candidate whose mean
coverage is closest to 1 − α wins, and the final band is built at it
(`CI_boot = 200`). Here α = 0.10, so the bands are 90% bands.

The calibration measures the bootstrap's self-consistency, not its true
coverage: it cannot see τ̂'s bias, and at ρ = 0.5, where m0 tracks τ, that bias
is where growth is expected. `cts_corr_ci_sf_metrics.R` records each run's pick
and the calibration's own coverage at it (`plugin_coverage`) beside the
true-τ coverage, so the gap is visible. There is no binary optimal_sf study.

## Files

Each study folder holds these files (the CI and optimal_sf studies use the
continuous and binary studies' `*_dgms.R`, and optimal_sf the CI studies'
`*_models.R`):

| file | |
|---|---|
| `<prefix>_config.R` | the parameter grid and results path, declared once. The array index is a row of `study$grid`: never filter or reorder it |
| `<prefix>_dgms.R` | the generator and oracle wrappers, with a `rho` argument |
| `<prefix>_models.R` | this study's call to `cate_methods()` |
| `<prefix>_analysis.R` | the array entry point: one grid row per index |
| `<prefix>_check.R` | lists missing runs in `jobscripts/failed_ids.txt` and points `<prefix>_rerun.sh` at it |
| `<prefix>_collect.R`, `<prefix>_metrics.R` | write `<prefix>_all.RDS` and `<prefix>_metrics.RDS` (and, for continuous and binary, `<prefix>_true_cate_tests.RDS`), all with a `rho` column |

In `sample_size/` itself:

| file | |
|---|---|
| `corr_results.qmd` | both outcomes: levels by ρ, paired differences, tests, ATE bias against `dr_oracle` |
| `corr_outcome_results.qmd` | one outcome (Quarto `outcome` param): every metric by ρ, paired differences, the true-CATE tests, NA tables, a headline table |
| `confidence_intervals/corr_ci_results.qmd` | the CI studies, one outcome per render (`outcome` param); the optimal_sf section, `corr_ci_results_sf.qmd`, is included for the continuous outcome only |
| `corr_truth_summary.R` | SD(τ) and cor(m0, τ) for every outcome, ρ and scenario; needs no simulation output |
| `corr_sl_failures.R` | summarises the SuperLearner arm's dropped learners and failed fits |
| `corr_blp_beta2.R` | pulls the saved BLP β2 estimates back out of the collected runs |

## Running a study

On the cluster, run `Rscript make_log_dirs.R` from the repo root once. Then
submit from the study's `jobscripts/` folder, because the scripts
`cd "${PBS_O_WORKDIR}/.."`:

```bash
cd sample_size/continuous/jobscripts
qsub cts_corr_1.sh            # rows 1-10000 (rho = 0 is 1-8000)
qsub cts_corr_2.sh            # rows 10001-16000, all rho = 0.5
qsub cts_corr_extra.sh        # rows 16001-20800, scenarios 5-10
Rscript ../cts_corr_check.R   # run after all three finish
qsub cts_corr_rerun.sh        # only if the check found failures
qsub cts_corr_collect.sh
qsub cts_corr_metrics.sh
```

Binary is the same with `bin_corr`. The CI studies have `*_1.sh` (1–10000) and
`*_2.sh` (10001–16000); optimal_sf has one, `cts_corr_ci_sf_1.sh` (1–1600), and
its todo list is `failed_cts_ids.txt`. The jobscripts split at the PBS array
cap, not at the ρ boundary.

Local smoke tests, from the study folder:

- `Rscript cts_corr_analysis.R 3` (scenario 3, n = 100, run 1, ρ = 0) and
  `8003` (the same run at ρ = 0.5);
- `Rscript bin_corr_analysis.R 18405` (binary scenario 9, n = 100, run 1,
  ρ = 0.5: the cell with its own scale);
- `Rscript cts_corr_ci_analysis.R 3` (scenario 3, n = 500, `CI_sf` = 0.05,
  run 1, ρ = 0).

A full optimal_sf calibration takes hours; shrink `n_sim`, `CI_boot` and
`sf_grid` in a scratch copy to check the wiring.

**Resources** are set by hand. The continuous and binary jobs take one core,
2 GB and 30 minutes; the CI jobs two cores, 3 GB and 30 minutes; optimal_sf two
cores, 7 GB and 8 hours. Size for the n = 1000 cell, since one `#PBS -l` line
covers every n. Keep `workers × grf_threads` (the trailing arguments on the
`Rscript` line) within `ncpus`, and change the two together.
