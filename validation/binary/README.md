# Interim-analysis validation — binary outcome

The continuous study (`../continuous/`) with a binary outcome: the same design,
estimators, importance measures and four chunk comparisons, all from the shared
`../val_common.R`. This README covers what differs. For everything else —
the one-trial split, the DR SuperLearner and its crossfit TE-VIMs, the two
importance measures, `p_cts_adj`, the no-queue runner and job sizing — see
`../continuous/README.md`.

## Design

| | |
|---|---|
| outcome | binary, treatment effect on the **risk-difference** scale |
| covariates | correlated, rho = 0.5 — the `binary_corr_0.5` set (`R/dgm_scenarios.R`), the same set `sample_size/binary/` runs as its primary arm |
| scenarios | 1: the null, `tau = bW`; 2: `tau = bW + b4 * tanh(X4)`, `b4 = 0.082`; 3: `bW + b3 * X3 + b4 * tanh(X4) + b45 * tanh(X4) * tanh(X5)`, `b3 = -0.102`, `b4 = -0.0255`, `b45 = 0.0255`; 4: `bW + b4 * cos(X4)`, `b4 = -0.204` (`R/dgm_scenarios.R`, `TE_RD`) |
| n | 1000, **one** trial split into its first `n * interim_prop` participants and the rest |
| interim_prop | 0.25 to 0.75 in steps of 0.05 — 11 interim points |
| runs | 100 per scenario — **4400 array jobs**, 1 h walltime (as continuous; check with `bin_val_testing.R full`) |
| folds | DR SuperLearner only: 5 for a chunk below 500 rows, else 10 |
| results | `../results/validation/binary/rho_0.5/scenario_<k>/1000/<interim_prop>/res_sim_<run>.RDS` |

## What differs from the continuous arm

**The DGM.** `P(Y = 1 | x, W) = m0(x) + W * tau(x)`: the binary sets put the
effect on the risk-difference scale (`R/dgm_scenarios.R`, BINARY OUTCOMES ARE
ON THE RISK-DIFFERENCE SCALE). X1 and X2 are prognostic only, through the
control risk `m0`. X4 and X5 enter through `tanh` rather than as they are (cos
is bounded already), each scenario's coefficients are the continuous ones
scaled by `RD_SCALE[k]` with the sign reversed, so every treated risk stays
inside [0.01, 0.99]. The modifiers are the continuous scenarios': none in 1, X4 in 2 and 4,
X3, X4 and X5 in 3.

**The heterogeneity is small.** At n = 1000 the CATE's SD is ≈ 0.051 in
scenario 2, 0.056 in scenario 3 and 0.094 in scenario 4 (in-sample, run 1).
`../val_results.qmd` computes each scenario's `Var(tau)` reference line from a
large draw of the DGM rather than hard-coding it.
With 250–750 binary outcomes per chunk, every interaction test here has much
less power than its continuous counterpart. Low replication rates partly
reflect that power, not only an interim finding that fails to hold.

**The DR SuperLearner's outcome model** is `binomial()` (`method.NNloglik`) —
`bin_val_models.R`'s whole difference from `cts_val_models.R`. Its propensity
model and stage 2 are as in the continuous arm. The forests take the 0/1
outcome as it is, so their outcome regressions estimate the risk. The CATE is
a risk difference throughout, so the TE-VIMs and the TreeSHAP surrogate need
nothing outcome-specific.

**HC3 interaction tests.** Every interaction test in `chunk_validations()` (the
subgroup tests and the five top-covariate forms, `p_cts`, `p_cts_adj`,
`p_split`, `p_dd`, `p_proj`) is `lm(Y ~ W * v)`. With a binary `Y` that is a
linear probability model. It is linear in the risk, the scale the effect is
on, so `W:v` is a risk-difference interaction. But its errors are
heteroskedastic by construction (`Var = p(1 - p)`), so the classical standard
errors are wrong. `bin_val_analysis.R` calls `chunk_validations_se()`, which
saves every test with classical and with `sandwich::vcovHC(type = "HC3")`
standard errors (`coef_pval()`, `../val_common.R`). HC3 is the primary kind;
the classical one is kept so the report can show how much the choice matters.
The continuous arm does the same, with HC3 primary for a different reason (see
`../continuous/README.md`). `sandwich` is already in `sim-env` (GenericML
depends on it), and `bin_val_testing.R` check 1 confirms it.

**The new top-covariate tests, the `p_split` fix and scenario 4's symmetry**
are shared code and are the same as in the continuous arm — see
`../continuous/README.md`, What is compared across chunks. In scenario 4
`b4 < 0`, so the top-10% responders are in X4's tails and the bottom 10% near
0, the reverse of the continuous arm; the tests do not care which.

## Files

| file | role |
|---|---|
| `bin_val_config.R` | the parameter grid and results path — **the** definition |
| `bin_val_dgms.R` | names this study's slice of `R/dgm_scenarios.R` (`binary_corr_<rho>`) |
| `bin_val_models.R` | `run_all_cate_methods()`: `../val_common.R`'s `fit_val_methods()` with `family = binomial()` |
| `bin_val_analysis.R` | array entry point; splits the trial, fits both chunks, runs `chunk_validations_se()` |
| `bin_val_testing.R` | pre-submission verification — dependencies, grid, the HC3 tests, the binary DGM in each scenario, the split, one chunk's fit, both standard-error kinds; `full` adds one replicate end to end |
| `bin_val_check.R` | finds missing runs, writes `jobscripts/failed_ids.txt`, and updates the rerun jobscript |
| `bin_val_collect.R` | gathers per-run files into `bin_val_all.RDS` |
| `bin_val_metrics.R` | flattens `validations` into tidy `bin_val_metrics.RDS` |
| `../val_results.qmd` | the written-up report, shared with `continuous/`, one render per scenario (Quarto `outcome` and `scenario` params) |

There is no `bin_val_run.R` (the continuous arm's no-queue runner) and no
`results_bin_val.R` (its standalone plot script). The report carries every
plot, and the runner could be copied from `cts_val_run.R` with its three paths
changed if the queue is ever slower than the work.

## Running it

```bash
Rscript validation/binary/bin_val_testing.R full   # do this first
qsub validation/binary/jobscripts/bin_val_1.sh     # 1-4400
Rscript validation/binary/bin_val_check.R          # writes failed_ids.txt if any are missing
qsub validation/binary/jobscripts/bin_val_collect.sh
qsub validation/binary/jobscripts/bin_val_metrics.sh
cd validation/binary
quarto render ../val_results.qmd -P outcome:binary -P scenario:1 --output bin_val_results_sc1.html   # and 2, 3, 4
```

If old results become incomparable, archive the tree first — see
`../continuous/README.md`, Running it; the one `archive_old_results.R` command
there covers both arms.

`bin_val_1.sh` runs `bin_val_analysis.R <i> 1 1` on `ncpus=1`, as the
continuous arm does. Its 1 h walltime is copied from `cts_val_1.sh`, not
measured for a binary row: the binomial SuperLearner fits may run slower.
Read the replicate time `bin_val_testing.R full` reports against it before
submitting.
