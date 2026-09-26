# Binary outcome — sample size study

The continuous study, on a logit scale. Same estimators, same crossfitting, same
sample-size sweep; the outcome is `rbinom(n, 1, plogis(lp))` and the estimand is
a **risk difference**.

## Design

| | |
|---|---|
| scenarios | **1, 3, 8, 9** — a subset, not all ten |
| n | 100, 250, 500, 1000 |
| runs | 100; **500 for scenarios 1, 3, 8, 9** (runs 101–500 are grid rows 4001–10400, `jobscripts/bin_extra.sh`) |
| array | **1,600 jobs** |
| results | `../results/binary/scenario_<k>/<n>/res_sim_<run>.RDS` |

The four scenarios are the ones the chapter reports: null, simple, complex and
non-linear (`SS_SCENARIO_LABELS` in `R/figures.R`).

The coefficient table differs from the continuous study — `b0 = -0.4`,
`b1 = 0.5`, `b2 = 0.5` rather than the continuous values — because the same
numbers on a logit scale would saturate `plogis`. Scenario 10's treatment effect
also differs (`exp(X4)` rather than `exp(-abs(X4))`), and `bW` is calibrated with
`power.prop.test` rather than `power.t.test`. The power targets differ too: here
`bW` is calibrated to 75%, while since bug O the continuous studies calibrate the
ATE to 80%. Those are the only differences from `continuous/`.

## The grid was declared three ways

Before the configs existed, this study's grid appeared in three places and they
did not agree:

| | |
|---|---|
| `bin_analysis.R` | `scenario = c(1:10)` |
| `bin_check.R` / `bin_collect.R` | `scenario = c(1, 3, 8, 9)` |
| `jobscripts/bin_1.sh` | `#PBS -J 1-1600` |

4 scenarios × 4 sample sizes × 100 runs is exactly 1,600, so `c(1, 3, 8, 9)` is
the design and the analysis script's `c(1:10)` was the stale line.

**The results on disk are the intended ones** — four scenarios, 100 runs each.
The stale `c(1:10)` did not corrupt them; it was simply out of step with the
grid the runs were actually launched from.

`bin_config.R` now declares `c(1, 3, 8, 9)` once and every script reads it from
there, so the three-way drift cannot recur.

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
qsub binary/jobscripts/bin_1.sh     # 1-1600
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

The DGM is unaffected by the bug ledger — bug A is the *confidence-interval*
binary study, not this one. `bias` also flips sign when the metrics are
regenerated (bug G), but that needs no cluster time.

## Open idea: HTE on the risk-difference scale

*Parked 2026-09-26 for later exploration. Nothing implemented.*

### The problem

The scenarios put the treatment effect on the **logit** scale, but the estimand
(`truth$tau`) is a **risk difference**. The logit-scale descriptions hold
exactly — the oracle formula and the generator's truth agree to 0 — but on the
RD scale the prognostic covariates X1 and X2 modify the effect too, through the
link, in every scenario. `bin_verify_hte.R` splits `var(tau)` into the part
explained by the described modifiers (X3/X4/X5) and the part explained by X1/X2:

| scenario | n | SD of true RD CATE | from described modifiers | from X1, X2 | ATE (RD) | marginal power |
|---|---|---|---|---|---|---|
| 1 Null | 100 | 0.044 | — | 100% | −0.24 | 0.74 |
| 1 Null | 1000 | 0.009 | — | 100% | −0.08 | 0.72 |
| 3 Simple | 100 | 0.055 | 34% | 63% | −0.24 | 0.73 |
| 3 Simple | 1000 | 0.045 | 94% | 5% | −0.08 | 0.71 |
| 8 Complex | 100 | 0.126 | 85% | 10% | −0.20 | 0.57 |
| 8 Complex | 1000 | 0.140 | 97% | 1% | −0.04 | 0.25 |
| 9 Non-linear | 100 | 0.050 | 59% | 37% | −0.19 | 0.51 |
| 9 Non-linear | 1000 | 0.050 | 99% | 0% | −0.01 | 0.05 |

The remainder (≤ 6% everywhere) is the modifier × X1/X2 interaction the link
creates. Three consequences:

- **Scenario 1 is not an RD null**, so its true-CATE test rejections in
  `bin_true_cate_tests.RDS` are not type I error, and `BLP_p` is not `NA` there
  as it is for continuous scenario 1.
- **The HTE structure changes with n.** `bW` is recalibrated per n and is
  larger at small n (−1.21 at n = 100, −0.35 at n = 1000), so link-induced HTE
  dominates small samples. Scenarios 2 and 10 are the extreme: 18% from the
  described modifiers at n = 100.
- **75% power only holds when the modifier terms average zero.** In scenario 9,
  E[0.5·cos(X4)] ≈ +0.30 almost cancels `bW` at n = 1000. Scenarios 6 and 8
  fall to 31% and 25%; scenario 2 rises to 99%.

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
delta <- power.prop.test(n = n / 2, p2 = p0_bar, power = 0.75)$p1 - p0_bar
```

This needs L ≥ |δ|. With the current `b0`/`b1`/`b2`, **L = 0.3, U = 0.95**:

| n | 100 | 250 | 500 | 1000 |
|---|---|---|---|---|
| δ | −0.26 | −0.17 | −0.12 | −0.08 |

Treated risk is then ≥ 0.04 for every x — a bound, not an empirical check.
L = 0.25 looks fine in 200,000 draws but its bound is L + δ = −0.01 at n = 100.
Simulated two-proportion-test power was 0.69–0.74 (the n = 100 shortfall is
`power.prop.test`'s normal approximation and random arm sizes).

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
