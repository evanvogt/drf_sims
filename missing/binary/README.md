# Missing covariates — binary outcome

The `missing/continuous` design with a binary outcome, on the sample-size studies'
risk-difference DGM (`sample_size/ADEMP.md`, "The binary risk-difference scale").
See `missing/ADEMP.md` for the mechanisms and handling methods.

| | |
|---|---|
| array | **12,600 jobs**: `bin_miss_1.sh` (1–9900, scenarios 1, 3, 4, 5) and `bin_miss_extra.sh` (9901–12600, scenario 2) |
| results | `../results/missing/binary/scenario_<k>/<n>/<type>/<prop>/<mechanism>/<method>/` |
| figures | `../miss_results.R` / `.qmd` (`outcome = binary`) — every metric, to `results/all_figures/`; the diagnostic counterpart to the chapter script `results_processing/thesis_figures/miss_bin.R` |

## The DGM

The `binary_missing` table is built from `binary`'s scenarios 1–6 in
`R/dgm_scenarios.R`, so scenario k here is scenario k there by construction.

- **The effect is on the risk-difference scale.** `P(Y = 1) = m0(x) + W·τ(x)`,
  with the control risk m0 bounded in [0.34, 0.70] and a 40% control event
  rate. The modifiers are the continuous ones sign-reversed and scaled by
  `RD_SCALE_MISS` (the main study's `RD_SCALE` except scenario 4, 0.179 against
  0.204; see "This study's own HTE scale" below), with X4 and X5 through
  `tanh`. That makes X1 and X2 purely prognostic, and only **weakly**
  prognostic: m0's SD is 0.048.
- **`bW`** is set so that the true ATE is the same marginal RD in all six
  scenarios. At n = 500 it is −0.118, −0.118, −0.052, −0.010, 0.073 and
  −0.103 for scenarios 1–6; the true RD is −0.118 and the power 0.80.
  Correlated covariates move scenario 3's `bW` (E[tanh(X4)·tanh(X5)] ≠ 0) and
  leave the 40% control event rate (0.398) as it is.
- **MNAR-tau.** U enters the treatment effect as `bU·tanh(U)`, with bU = 0.12
  and U ~ N(0, 2²). That term is bounded, so every treated risk stays inside
  [0.01, 0.99]; `RD_SCALE_MISS` allows for it, and it binds scenario 4's. It
  has mean zero, so it leaves the RD and the power as they are, and the U-free
  truth *is* the average over U. It is calibrated without U, as in
  `missing/continuous`, so every mechanism shares one `bW` and one truth per
  scenario.
- **MNAR-Y0.** The same `bU·tanh(U)` added to the control risk instead, so in
  both arms: an SD of 0.095, about twice m0's own (0.048). The control risk
  stays in [0.22, 0.82] and the treated risk has exactly MNAR-tau's bound, so
  `RD_SCALE_MISS` needs nothing new.

`dr_oracle`'s outcome model returns the true outcome mean, so it is an exact
oracle under every mechanism.

### This study's own HTE scale

The main binary studies derive `RD_SCALE` at n = 100, where the planned RD
(−0.248) leaves almost no room above the risk floor. This study runs n = 500
only, where the RD is −0.118, so it has its own scale, `RD_SCALE_MISS`, and
spends the extra floor room on U:

- `bU = 0.12` (`BU_MISS`) is the largest value, to 2 dp, that leaves the
  floor-bound scenarios 2, 3, 5 and 6 exactly the HTE `sample_size/binary`
  gives them.
- Scenario 4 is bound by the ceiling instead, so its scale drops from 0.204
  to 0.179 (SD of τ 0.080 rather than 0.091).
- `sU = 2` pushes `tanh(U)` towards ±1. The amputation standardises U, so
  this strengthens the selection without touching the bounds.

`R/bin_verify_hte.R` re-derives both (its check 6) on the correlated set.

The mechanisms' strength, from `missing/miss_dgm_checks.R` (scenario 2,
n = 20,000):

| | binary | `missing/continuous` |
|---|---|---|
| U's share of Var(Y) | 4% | 45% |
| AUC of incompleteness from Y, MNAR-Y0 | 0.54 | 0.66 |
| E[U-term \| complete] | −0.023 | −0.25 |
| CC ATE bias under MNAR-tau, % of the ATE | ~20% | ~90% |

Still far weaker than the continuous mechanisms: a Bernoulli outcome's own
variance dominates, and the bounds cap bU. Read the binary MNAR results as a
mild version of the continuous ones, and binary MNAR-Y0 as nearly MCAR.

`Rscript R/calibration_report.R` prints `bW`, the RD and the power, with the
MNAR treated-risk floor and ceiling (shared by MNAR-Y0 and MNAR-tau).

## True-CATE HTE test evaluation

`bin_miss_metrics.R` also writes `bin_miss_true_cate_tests.RDS` — the BLP and
independence tests run on the true CATE and true nuisances (`truth$tau`,
`truth$p0`, `W.hat = 0.5`) instead of an estimator's fitted ones, one
`BLP_p`/`BLP_p_os`/`indep_cate` row per (scenario, n, type, prop, mechanism, method,
run). See `sample_size/ADEMP.md`'s "Heterogeneity tests" for what it
means, why every scenario-1 true-CATE test is `NA`, and the test columns
themselves (`BLP_p_os` is the one-sided, HC3 BLP).
`method == "multiple_imputation"` rows are `NA`/`NA` too: `data` is a list of
50 imputed data.frames there, with no single covariate matrix to test against
(see "Open decision" in `missing/ADEMP.md`).

## Running it

```bash
qsub missing/binary/jobscripts/bin_miss_1.sh        # 1-9900
qsub missing/binary/jobscripts/bin_miss_extra.sh    # 9901-12600, scenario 2
Rscript study.R check missing/binary/bin_miss_config.R
qsub missing/binary/jobscripts/bin_miss_collect.sh
qsub missing/binary/jobscripts/bin_miss_metrics.sh
```

Every model's HTE tests are written by the simulation itself
(`PROFILES$missing`), and the MI arm saves its per-imputation tests
(`mi_tests`, see `missing/ADEMP.md`).

Then, for the figures:

```bash
Rscript missing/miss_results.R binary              # every metric, to results/all_figures/
quarto render missing/miss_results.qmd -P outcome:binary --output bin_miss_results.html
```
