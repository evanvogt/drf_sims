# Binary outcome — sample size study

The continuous study with a binary outcome. Same estimators, same crossfitting,
same sample-size sweep. The outcome is `rbinom(n, 1, m0(x) + W·τ(x))`: the
treatment effect adds to the risk, so it is on the **risk-difference** scale,
and so is the estimand.

## Design

| | |
|---|---|
| scenarios | 1–10 |
| n | 100, 250, 500, 1000 |
| runs | 100; **500 for scenarios 1–4** (runs 101–500 are grid rows 4001–10400, `jobscripts/bin_extra.sh`) |
| array | **4,000 jobs** (`bin_1.sh`), plus **6,400** (`bin_extra.sh`) |
| results | `../results/binary/scenario_<k>/<n>/res_sim_<run>.RDS` |

Scenarios 1–4 are the ones the chapter reports: null, simple, complex
and non-linear (`SS_SCENARIO_LABELS` in `R/figures.R`).

### Outcome model and `bW` calibration

Every scenario shares one control-risk model and differs only in its treatment
effect, which adds to the risk:

`P(Y = 1 | x, W) = m0(x) + W·τ(x)`, with
`m0(x) = 0.34 + 0.36·plogis(b0 + b1·X1 + b2·X2)`, `τ(x) = bW + g(x)`,
`X1 ~ Bernoulli(0.4)` and `X2 ~ N(0, 1)`.

The event is harmful, and treatment lowers its risk.

| | |
|---|---|
| control risk `m0` | `b0 = −1.72`, `b1 = −0.5`, `b2 = 1`, in every scenario. b1 and b2 are `continuous/`'s; b0 sets the control event rate E[m0] to 0.400. m0 lies in [0.34, 0.70]: 5th–95th percentile 0.35–0.50, SD 0.048 |
| `g(x)` | the continuous scenario's heterogeneity with its signs reversed, scaled by `RD_SCALE`, and with X4 and X5 entering through `tanh`: the scenario's `te_expr` in `R/dgm_scenarios.R`, with `bW = 0` |
| `bW` | set so the true **ATE**, a marginal risk difference, equals the effect the trial was planned to detect |

Because τ adds to the risk, **the true CATE is exactly `bW + g(x)`**. X1 and X2
are purely prognostic, scenario 1 is an exact null, and g is the same at every
n; only the ATE moves with n. This is also the continuous model,
`E[Y | x, W] = m0(x) + W·τ(x)`, with a different `m0` and noise.

**Why the control risk is bounded.** A risk-difference effect can only add if
every control risk leaves room for it: `m0(x) + τ(x)` must stay inside [0, 1]
for every x. A logistic control risk reaches 0, so m0 is a logistic scaled into
[0.34, 0.70]. The floor is set by n = 100, where the planned RD is −0.248, so
every control risk must clear that plus the extra benefit the HTE gives some
patients. At a 40% control event rate that leaves little room for m0 to spread,
so X1 and X2 are only weakly prognostic: SD 0.048, against 0.13 under the old
logit-scale design.

**The HTE sizes.**

- `g = −RD_SCALE[k] × (the continuous coefficients)`.
- Each `RD_SCALE[k]` is the largest scale, floored to 3 dp, that keeps every
  treated risk inside [0.01, 0.99] (`RD_EPS`):
  - at n = 100 to 1000;
  - for scenarios 2–6, also under `missing/binary`'s MNAR-Y, which adds
    `bU·tanh(U)` with bU = 0.08 at n = 500.
- The binding bound:
  - the floor at n = 100 binds in every scenario except 4;
  - in scenario 4 the MNAR-Y ceiling binds.
- `bin_verify_hte.R` re-derives `RD_SCALE`. If `p0_lo`, `p0_hi`, b0–b2, bU,
  `RD_EPS` or the studies' n change, `RD_SCALE` has to be recomputed.

In the table below:

- "τ 5–95%" is the middle 90% of the true CATE.
- "harmed" is the share of patients with τ > 0 at n = 1000.
- "floor" and "ceiling" are the lowest and highest treated risk over the
  covariate support, with MNAR-Y included for scenarios 2–6.
- t4 = tanh(X4) and t5 = tanh(X5).

| # | g(x) | SD τ | τ 5–95%, n = 100 | τ 5–95%, n = 1000 | harmed | floor | ceiling |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 0 | −0.248 | −0.085 | 0% | 0.092 | 0.615 |
| 2 | 0.082·t4 | 0.051 | −0.324, −0.172 | −0.161, −0.009 | 0% | 0.010 | 0.744 |
| 3 | −0.102·X3 − 0.026·t4 + 0.026·t4·t5 | 0.050 | −0.308, −0.158 | −0.145, 0.005 | 6.5% | 0.011 | 0.784 |
| 4 | −0.204·cos(X4) | 0.091 | −0.328, −0.047 | −0.165, 0.116 | 16.8% | 0.012 | 0.989 |
| 5 | −0.274·X3 | 0.126 | −0.330, −0.056 | −0.167, 0.107 | 30% | 0.010 | 0.853 |
| 6 | −0.023·X3 + 0.075·t4 | 0.048 | −0.322, −0.176 | −0.159, −0.013 | 1.6% | 0.011 | 0.752 |
| 7 | −0.082·X3·t4 | 0.043 | −0.322, −0.174 | −0.159, −0.011 | 0% | 0.010 | 0.697 |
| 8 | −0.274·X3 − 0.069·t4 + 0.069·X3·t4 | 0.128 | −0.330, −0.005 | −0.167, 0.158 | 30% | 0.010 | 0.875 |
| 9 | 0.082·t4·t5 | 0.032 | −0.304, −0.192 | −0.141, −0.029 | 0% | 0.010 | 0.697 |
| 10 | −0.179·X3 − 0.060·exp(−\|X4\|) | 0.083 | −0.325, −0.107 | −0.162, 0.056 | 30% | 0.010 | 0.771 |

X3 subgroup means at n = 1000:

| scenario | X3 = 0 | X3 = 1 |
|---|---|---|
| 3 | −0.013 | −0.115 |
| 5 | +0.107 | −0.167 |
| 6 | −0.069 | −0.091 |
| 7 | −0.085 | −0.085 |
| 8 | +0.107 | −0.167 |
| 10 | +0.040 | −0.139 |

It is always the X3 = 0 minority (30%) that benefits least. In scenarios 5, 8
and 10 that group is harmed at n = 1000; at n = 100 every group still benefits.
In scenario 4 the harmed patients are those with large |X4|.

**The signs are reversed relative to `continuous/`.** Reversing them points
each asymmetric scenario's larger swing towards *less* benefit. The ceiling
leaves room in that direction; only the floor is tight. Under the continuous
signs, SD(τ) would be:

| scenario | continuous signs | reversed |
|---|---|---|
| 3 | 0.034 | 0.050 |
| 4 | 0.023 | 0.091 |
| 5 | 0.053 | 0.126 |
| 8 | 0.040 | 0.128 |
| 10 | 0.044 | 0.083 |

The price is that the opposite subgroup benefits more than in `continuous/`.
This undoes bug P's sign harmonisation (below). The relative sizes of the
scenarios come from the bounds, not from `continuous/`.

**Calibration.** The trial is planned as in `continuous/`: to detect an ATE,
assuming the effect is the same for everyone.

1. The control event rate is p̄0 = E[m0] = 0.400.
2. The planned treated rate p̄1 gives 80% power (`TARGET_POWER`) in a
   two-proportion test with n/2 per arm.
3. `bW = (p̄1 − p̄0) − E[g]`, the same form as the continuous
   `bW = −δ − E[g]`.

g enters only through E[g]. The plan assumes homogeneity, but the data aren't
homogeneous, so the `bW` that delivers the planned RD depends on g. Both
expectations are taken by quadrature (`baseline_grid()`, `te_grid()`), with no
random draws, so the draw order is untouched. `bW` is rounded to 3 dp; at 2 dp
the RD would move by up to 0.005.

| n | 100 | 250 | 500 | 1000 |
|---|---|---|---|---|
| true ATE (RD), every scenario | −0.248 | −0.164 | −0.118 | −0.085 |

Unlike `continuous/`, every scenario realises the planned 80% (0.796–0.803,
from rounding `bW`). A binary outcome's variance in each arm is p̄(1 − p̄),
fixed by the marginal risks, so the heterogeneity has no variance to add.

`Rscript R/calibration_report.R` prints `bW`, the true RD, the power and the
treated-risk floor and ceiling for every binary scenario at every n, and for
`missing/binary/`. `Rscript bin_verify_hte.R` (from `sample_size/binary/`) checks the
design end to end. See `continuous/README.md` for when to run them.

**Before bug P** (root README), `bW` was calibrated for 75% power at the risk
plogis(b0) = 0.401 — the risk at X1 = X2 = 0, not the population's 0.453 — and
set to the planned log-odds ratio rather than the ATE. So `bW` was the same in
every scenario, the true RD drifted with g, and power ran from 5% (scenario 4 at
n = 1000, RD −0.010) to 99% (scenarios 5 and 6). The modifiers in scenarios 2,
5, 6 and 7 also had the opposite sign to the continuous ones; bug P gave them
the continuous signs, and the risk-difference design has now reversed all of
them again, for the reason above.

**Before the risk-difference DGM** (2026-09-26), the effect was on the logit
scale: `P(Y = 1) = plogis(b0 + b1·X1 + b2·X2 + W·(bW + g))`, with b0 = −0.4,
b1 = b2 = 0.5, log-odds modifiers, and scenario 10's effect `b4·exp(X4)`. See
"Why the effect is on the risk-difference scale" below. Every binary dataset
generated before then is superseded.

## The grid was declared three ways

Before the configs existed, this study's grid appeared in three places and they
did not agree (scenario numbers in this section are the pre-2026-09-26 ones;
1, 3, 8 and 9 are now 1–4):

| | |
|---|---|
| `bin_analysis.R` | `scenario = c(1:10)` |
| `bin_check.R` / `bin_collect.R` | `scenario = c(1, 3, 8, 9)` |
| `jobscripts/bin_1.sh` | `#PBS -J 1-1600` |

4 scenarios × 4 sample sizes × 100 runs is exactly 1,600, which made
`c(1, 3, 8, 9)` look like the design. But `expand.grid` varies the first column
fastest, so submitting indices 1–1600 against the analysis script's 4,000-row
`c(1:10)` grid ran runs 1–40 of **all ten** scenarios. `bin_collect.R` then
looked for scenarios 1, 3, 8 and 9 and found 40 runs in each, so the results on
disk have 40 replicates per cell, not 100 (see `bin_config.R`'s header).

`bin_config.R` now declares the grid once and every script reads it from there,
so the drift cannot recur. The grid is all ten scenarios at 100 runs, with 1–4
taken to 500. The study re-runs in full anyway.

## Files

Same layout as `continuous/`. `bin_models.R` sets `family = binomial()`. The
oracle formula from `bin_dgms.R` returns the risk itself, so `dr_oracle` applies
no link — as in every study now (see `R/README.md`).

`bin_metrics.R` also writes `bin_true_cate_tests.RDS` — the BLP and
independence tests run on the true CATE and true nuisances instead of an
estimator's, one row per (scenario, n, run). See
`continuous/README.md`'s "True-CATE HTE test evaluation" for what it means
and why scenario 1's `BLP_p` is `NA`. That now holds here too: scenario 1's
true CATE is exactly constant. Under the logit-scale design it was not, so the
test rejected there.

## Running it

```bash
qsub sample_size/binary/jobscripts/bin_1.sh     # 1-4000
qsub sample_size/binary/jobscripts/bin_extra.sh # 4001-10400: runs 101-500, scenarios 1-4
Rscript sample_size/binary/bin_check.R          # writes failed_ids.txt if any are missing
qsub sample_size/binary/jobscripts/bin_collect.sh
qsub sample_size/binary/jobscripts/bin_metrics.sh
```

## Status

**Archive the old results first** - `R/archive_old_results.R` (root `README.md`, Status, step 0). They predate the current DGM and use the pre-2026-09-26 scenario numbers, so running into that tree would mix old and new results.

**Re-runs required** — for the crossfitting strategy change to
`R/cate_models.R` (see root README Methods/Status), which moves all five
estimator arms, and separately for bug F, which changes `dr_superlearner` as
in `continuous/`. Only the `dr_superlearner` arm moves for bug F specifically;
the harness can confirm the other four are unchanged there.

**Also re-run for bug P and the risk-difference DGM.** Bug P changed the `bW`
calibration. The risk-difference DGM then changed the outcome model, every
modifier, and scenario 10's covariates. Together they move every dataset's
outcome and the true CATE in all ten scenarios. Submit on the risk-difference
code. Any bug-P re-run submitted on the logit-scale DGM is superseded too, so
archive it with the older results. Bug A is the *confidence-interval* binary
study, not this one. `bias` also flips sign when the metrics are regenerated
(bug G), but that needs no cluster time.

**Also re-run for the per-arm outcome models and the SuperLearner libraries**
(2026-09-27), as in `continuous/`. At n = 100 a treated arm's training rows can
hold ~2 events (scenario 4, where the treated risk is lowest), so per-arm
SuperLearner fits there often drop learners, and a whole fit can fail; it then
falls back to the arm's event rate and is recorded in `sl_dropped`, rather
than aborting the run. A smoke run at n = 1000 (index 34) on one core, as
`jobscripts/bin_1.sh` runs it, took 2.8 minutes against the 1h walltime.

## Why the effect is on the risk-difference scale

Until 2026-09-26 the scenarios put the treatment effect on the **logit** scale,
but the estimand (`truth$tau`) is a **risk difference**. The logit-scale
descriptions held exactly, but on the RD scale the prognostic covariates X1 and
X2 modified the effect too, through the link, in every scenario. The
logit-scale version of `bin_verify_hte.R` (in git history) split `var(tau)` into
three parts:

- the part explained by the described modifiers (X3/X4/X5);
- the part explained by X1/X2;
- a remainder of at most 6%, the modifier × X1/X2 interaction the link creates.

After bug P:

| scenario | n | SD of true RD CATE | from described modifiers | from X1, X2 | ATE (RD) | marginal power |
|---|---|---|---|---|---|---|
| 1 Null | 100 | 0.048 | — | 100% | −0.26 | 0.80 |
| 1 Null | 1000 | 0.011 | — | 100% | −0.09 | 0.81 |
| 2 Simple | 100 | 0.057 | 27% | 71% | −0.26 | 0.80 |
| 2 Simple | 1000 | 0.045 | 93% | 6% | −0.09 | 0.80 |
| 3 Complex | 100 | 0.120 | 75% | 20% | −0.26 | 0.80 |
| 3 Complex | 1000 | 0.138 | 95% | 2% | −0.09 | 0.79 |
| 4 Non-linear | 100 | 0.058 | 27% | 70% | −0.26 | 0.80 |
| 4 Non-linear | 1000 | 0.048 | 92% | 6% | −0.09 | 0.81 |

That had two consequences:

- **Scenario 1 was not an RD null**, so its true-CATE test rejections were not
  type I error, and `BLP_p` was not `NA` there as it is for continuous
  scenario 1.
- **The HTE structure changed with n.** `bW` is recalibrated per n and was
  larger at small n (−1.31 at n = 100, −0.39 at n = 1000 in scenario 1), so
  link-induced HTE dominated small samples. Scenarios 10 and 5 were the
  extreme: only 20% and 22% of `var(tau)` came from the described modifiers at
  n = 100.

With prognostic covariates, an effect can be constant on at most one of the two
scales. The design above makes it the risk-difference scale, which is the scale
of the estimand. On the logit scale the same scenarios are now heterogeneous
through X1 and X2, which no estimator here targets.
