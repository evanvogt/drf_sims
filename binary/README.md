# Binary outcome — sample size study

The continuous study, on a logit scale. Same estimators, same crossfitting, same
sample-size sweep; the outcome is `rbinom(n, 1, plogis(lp))` and the estimand is
a **risk difference**.

## Design

| | |
|---|---|
| scenarios | 1–10 |
| n | 100, 250, 500, 1000 |
| runs | 100; **500 for scenarios 1, 3, 8, 9** (runs 101–500 are grid rows 4001–10400, `jobscripts/bin_extra.sh`) |
| array | **4,000 jobs** (`bin_1.sh`), plus **6,400** (`bin_extra.sh`) |
| results | `../results/binary/scenario_<k>/<n>/res_sim_<run>.RDS` |

Scenarios 1, 3, 8 and 9 are the ones the chapter reports: null, simple, complex
and non-linear (`SS_SCENARIO_LABELS` in `R/figures.R`).

### Outcome model and `bW` calibration

Every scenario shares one outcome model and differs only in its treatment
effect, which is on the logit scale:

`P(Y = 1) = plogis(b0 + b1·X1 + b2·X2 + W·(bW + g(x)))`, with
`X1 ~ Bernoulli(0.4)` and `X2 ~ N(0, 1)`.

| | |
|---|---|
| baseline | `b0 = −0.4`, `b1 = 0.5`, `b2 = 0.5`, in every scenario |
| `g(x)` | the scenario's heterogeneity term, in log-odds: its `te_expr` in `R/dgm_scenarios.R`, with `bW = 0` |
| `bW` | set so the true **ATE**, a marginal risk difference, equals the effect the trial was planned to detect |

The trial is planned as in `continuous/`: to detect an ATE, assuming the effect
is the same for everyone. The control-arm event rate is
p̄0 = E[plogis(b0 + b1·X1 + b2·X2)] = 0.453, and the planned treated-arm rate p̄1
gives 80% power (`TARGET_POWER`) in a two-proportion test with n/2 per arm.
`bW` then solves E[plogis(b0 + b1·X1 + b2·X2 + bW + g)] = p̄1, the logit-link
counterpart of the continuous `bW = −δ − E[g]`. g enters only that second step:
the plan assumes homogeneity, but the data aren't homogeneous, so the `bW` that
delivers the planned RD depends on g. Both steps use quadrature
(`baseline_grid()`, `te_grid()`, `marginal_risk()`), with no random draws, so
the draw order is untouched.

| n | 100 | 250 | 500 | 1000 |
|---|---|---|---|---|
| true ATE (RD), every scenario | −0.26 | −0.17 | −0.12 | −0.09 |

Unlike `continuous/`, every scenario realises the planned 80% (0.79–0.81, from
rounding `bW`). A binary outcome's variance in each arm is p̄(1 − p̄), fixed by
the marginal risks, so the heterogeneity has no variance to add.

The effect modifiers point the same way as in `continuous/`, so the same
subgroups benefit more. Their magnitudes are smaller because they are log-odds:
the continuous `b3 = 2` would be an odds ratio of 7.4. The other differences
from `continuous/` are the baseline, for the same reason, and scenario 10's
treatment effect, `b4·exp(X4)` rather than `b3·X3 + b4·exp(−|X4|)`.

**Before bug P** (root README), `bW` was calibrated for 75% power at the risk
plogis(b0) = 0.401 — the risk at X1 = X2 = 0, not the population's 0.453 — and
set to the planned log-odds ratio rather than the ATE. So `bW` was the same in
every scenario, the true RD drifted with g, and power ran from 5% (scenario 9 at
n = 1000, RD −0.010) to 99% (scenarios 2 and 4). The modifiers in scenarios 2–5
also had the opposite sign to the continuous ones (`b3 = −0.4` in 2 and 4,
`b4 = 0.2` in 3 and 0.3 in 4, `b34 = −0.5` in 5); only the signs changed. The
level of the true CATE moved in every scenario, so the metrics are not
comparable with earlier results.

`Rscript R/calibration_report.R` prints `bW`, the true RD and the power for every
binary scenario at every n, and for `missing/binary/`. See `continuous/README.md`
for when to run it.

## The grid was declared three ways

Before the configs existed, this study's grid appeared in three places and they
did not agree:

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
so the drift cannot recur. The grid is all ten scenarios at 100 runs, with 1, 3,
8 and 9 taken to 500. The study re-runs in full anyway.

## Files

Same layout as `continuous/`. `bin_models.R` sets `family = binomial()` and
`oracle_link = "logit"` — the oracle formula from `bin_dgms.R` is a linear
predictor, so the model code applies `plogis`. (`missing/binary/` uses the
opposite convention; see `R/README.md`.)

`bin_metrics.R` also writes `bin_true_cate_tests.RDS` — the BLP and
independence tests run on the true CATE and true nuisances instead of an
estimator's, one row per (scenario, n, run). See
`continuous/README.md`'s "True-CATE HTE test evaluation" for what it means
and why scenario 1's `BLP_p` is `NA`.

## Running it

```bash
qsub binary/jobscripts/bin_1.sh     # 1-4000
qsub binary/jobscripts/bin_extra.sh # 4001-10400: runs 101-500, scenarios 1, 3, 8, 9
Rscript binary/bin_check.R
```

**TODO:** `bin_metrics.sh` doesn't exist yet in `jobscripts/` — create it,
mirroring `continuous/jobscripts/cts_metrics.sh`, before the metrics step can
run on the cluster. (`bin_collect.sh` is there.)

## Status

**Re-runs required** — for the crossfitting strategy change to
`R/cate_models.R` (see root README Methods/Status), which moves all five
estimator arms, and separately for bug F, which changes `dr_superlearner` as
in `continuous/`. Only the `dr_superlearner` arm moves for bug F specifically;
the harness can confirm the other four are unchanged there.

**Also re-run for bug P** — the `bW` calibration and the modifier signs in
scenarios 2–5 changed (see "Outcome model and `bW` calibration" above), which
moves every dataset's outcome and the level of the true CATE in all ten
scenarios. Bug A is the *confidence-interval* binary study, not this one.
`bias` also flips sign when the metrics are regenerated (bug G), but that needs
no cluster time.

## Open idea: HTE on the risk-difference scale

*Parked 2026-09-26 for later exploration. Nothing implemented. Numbers updated
for bug P.*

### The problem

The scenarios put the treatment effect on the **logit** scale, but the estimand
(`truth$tau`) is a **risk difference**. The logit-scale descriptions hold
exactly — the oracle formula and the generator's truth agree to 0 — but on the
RD scale the prognostic covariates X1 and X2 modify the effect too, through the
link, in every scenario. `bin_verify_hte.R` splits `var(tau)` into the part
explained by the described modifiers (X3/X4/X5) and the part explained by X1/X2:

| scenario | n | SD of true RD CATE | from described modifiers | from X1, X2 | ATE (RD) | marginal power |
|---|---|---|---|---|---|---|
| 1 Null | 100 | 0.048 | — | 100% | −0.26 | 0.80 |
| 1 Null | 1000 | 0.011 | — | 100% | −0.09 | 0.81 |
| 3 Simple | 100 | 0.057 | 27% | 71% | −0.26 | 0.80 |
| 3 Simple | 1000 | 0.045 | 93% | 6% | −0.09 | 0.80 |
| 8 Complex | 100 | 0.120 | 75% | 20% | −0.26 | 0.80 |
| 8 Complex | 1000 | 0.138 | 95% | 2% | −0.09 | 0.79 |
| 9 Non-linear | 100 | 0.058 | 27% | 70% | −0.26 | 0.80 |
| 9 Non-linear | 1000 | 0.048 | 92% | 6% | −0.09 | 0.81 |

The remainder (≤ 6% everywhere) is the modifier × X1/X2 interaction the link
creates. Two consequences:

- **Scenario 1 is not an RD null**, so its true-CATE test rejections in
  `bin_true_cate_tests.RDS` are not type I error, and `BLP_p` is not `NA` there
  as it is for continuous scenario 1.
- **The HTE structure changes with n.** `bW` is recalibrated per n and is
  larger at small n (−1.31 at n = 100, −0.39 at n = 1000 in scenario 1), so
  link-induced HTE dominates small samples. Scenarios 10 and 2 are the extreme:
  20% and 22% from the described modifiers at n = 100.

Before bug P there was a third: `bW` was calibrated on the log-odds scale, so
the modifiers moved the ATE and the power (5% in scenario 9 at n = 1000). The
calibration now targets the marginal RD, so the ATE and power are the same in
every scenario (see "Outcome model and `bW` calibration").

The shape of the described HTE survives in direction: E[tau | modifiers] is a
monotone transform of the logit-scale effect, though not a linear one.

### A proposed design for a true RD null

With prognostic covariates, an effect can be constant on at most one of the two
scales. For a constant RD, bound the control risk and add δ:

```r
p0 <- L + (U - L) * plogis(b0 + b1 * X1 + b2 * X2)   # control risk in [L, U]
p1 <- p0 + delta                                      # tau = delta for every x
Y  <- rbinom(n, 1, ifelse(W == 1, p1, p0))

# calibrate on the RD scale; p0_bar by quadrature over X1, X2 (GH_NODES, no RNG)
delta <- power.prop.test(n = n / 2, p2 = p0_bar, power = TARGET_POWER)$p1 - p0_bar
```

This needs L ≥ |δ|. With the current `b0`/`b1`/`b2`, **L = 0.3, U = 0.95**
(p̄0 = 0.595):

| n | 100 | 250 | 500 | 1000 |
|---|---|---|---|---|
| δ | −0.276 | −0.176 | −0.125 | −0.088 |

Treated risk is then ≥ 0.02 for every x — a bound, not an empirical check.
L = 0.25 is not enough: its bound is L + δ = −0.03 at n = 100. Simulated power,
2,000 datasets per n with an uncorrected two-proportion test, was 0.79–0.81.

Costs: control-risk SD drops from 0.13 to 0.08 (scaling `b1`, `b2` by about 1.6
restores it — recheck L + δ ≥ 0 afterwards), and the control-arm model is a
scaled logistic rather than a logistic regression. On the logit scale the same
scenario is heterogeneous: log-OR(x) = qlogis(p0 + δ) − qlogis(p0).

### Fitting it into the pipeline

- **Oracle:** `"qlogis(L + (U-L)*plogis(b0 + b1*X$X1 + b2*X$X2) + W*bW)"` works
  with the existing `oracle_link = "logit"`, since plogis(qlogis(p)) = p —
  `run_dr_oracle` needs no change.
- **Generator:** `generate_scenario_data()` and `truth_at()` need a branch,
  e.g. driven by an `effect_scale = "logit"/"rd"` column on the scenario table.
- **Draw order:** Y is still one `rbinom(n, 1, ·)` in the same position, so the
  existing scenarios' regression fingerprints don't move.
- **RD-scale HTE:** τ(x) = δ + g(modifiers) works the same way if g is bounded
  and L ≥ max |τ(x)|. Binary X3 and cos(X4) are fine as they are; a term linear
  in a normal X4 is unbounded and needs bounding (e.g. `tanh(X4)`).

### Still to decide

Whether to add this as a new scenario, or keep the current scenarios and
describe them as logit-scale HTE, reporting the variance-share table above as a
property of the design.

Reproduce the table from `binary/` with `Rscript bin_verify_hte.R` (writes
nothing).
