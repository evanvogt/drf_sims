# Missing covariates — binary outcome

The `missing/continuous` design with a binary outcome, on `binary/`'s
risk-difference DGM (`sample_size/binary/README.md`, "Outcome model and `bW` calibration").
See `missing/README.md` for the mechanisms, handling methods and the shared bug
fixes.

| | |
|---|---|
| array | **12,600 jobs**: `bin_miss_1.sh` (1–9900, scenarios 1, 3, 4, 5) and `bin_miss_extra.sh` (9901–12600, scenario 2) |
| results | `../results/missing/binary/scenario_<k>/<n>/<type>/<prop>/<mechanism>/<method>/` |
| figures | `bin_miss_results.R` / `.qmd` — every metric, to `results/all_figures/`; the diagnostic counterpart to the chapter script `results_processing/thesis_figures/miss_bin.R` |

## ⚠ This file was a half-converted fork

`bin_miss_dgms.R` was copied from the continuous version and only partly converted to a binary outcome, carrying three related defects (continuous coefficient table, wrong power-test calibration, un-plogis'd truth) — all fixed together; the code now always uses the corrected values. Its oracle link convention broke later, when the fix flags were removed — see bug M below.

### The coefficients

The `binary_missing` table is built from `binary`'s scenarios 1–6 in
`R/dgm_scenarios.R`, so scenario k here is scenario k there by construction
(scenario numbers in this README are the post-2026-09-26 ones;
`missing/README.md` has the old numbering). Until the risk-difference DGM this
table carried its own copy of the coefficients, and for scenarios 1 and 3–6 the
correspondence had been inferred from the scenario descriptions.

### Bug P and the risk-difference DGM

The binary calibration changed with bug P, and the whole binary DGM with the
move to the risk-difference scale; `sample_size/binary/README.md`'s "Outcome model and `bW`
calibration" has the design. Here it means:

- **The effect is on the risk-difference scale.** `P(Y = 1) = m0(x) + W·τ(x)`,
  with the control risk m0 bounded in [0.34, 0.70] and a 40% control event
  rate. The modifiers are the continuous ones sign-reversed and scaled by
  `RD_SCALE`, with X4 and X5 through `tanh`. That makes X1 and X2 purely
  prognostic, and only **weakly** prognostic: m0's SD is 0.048. So amputating
  X1 and X2 (they are two of the `both` covariates) costs an estimator less
  than it did under the logit design.
- **`bW`** is set so that the true ATE is the same marginal RD in all six
  scenarios. At n = 500 it is −0.118, −0.118, −0.047, 0.005, 0.073 and −0.103
  for scenarios 1–6; the true RD is −0.118 and the power 0.80.
- **MNAR-tau** (was MNAR-Y). U enters the treatment effect as `bU·tanh(U)`,
  with bU = 0.08. That term is bounded, so every treated risk stays inside
  [0.01, 0.99]. `RD_SCALE` allows for it, and it binds scenario 4's. It has
  mean zero, so it leaves the RD and the power as they are. Under the logit
  design, averaging over U pulled the treated risk towards 0.5 and shrank the
  RD to about −0.10, and the power to 0.61–0.63. It is still calibrated
  without U, as in `missing/continuous`, so every mechanism shares one `bW`
  and one truth per scenario.
- **MNAR-Y0** (since 2026-09-28). The same `bU·tanh(U)` added to the control
  risk instead, so in both arms: an SD of 0.050, about m0's own (0.048). The
  control risk stays in [0.26, 0.78] and the treated risk has exactly
  MNAR-tau's bound, so `RD_SCALE` needs nothing new.
- **Correlated covariates** (since 2026-09-28, `missing/ADEMP.md`). They move
  scenario 3's `bW` to −0.052 (E[tanh(X4)·tanh(X5)] ≠ 0) and leave the others
  and the 40% control event rate (0.398) as they were.
  `sample_size/binary/bin_verify_hte.R` re-derives the MNAR caps on the
  correlated set; `RD_SCALE` is unchanged.

`Rscript R/calibration_report.R` prints this table, with the MNAR
treated-risk floor and ceiling (shared by MNAR-Y0 and MNAR-tau).

## Bug M — `dr_oracle` on the log-odds scale (fixed; no finished result affected)

`dr_oracle` is the DR-learner handed the *true* outcome model. This study used
to keep its own copy of the oracle formulas, wrapped in `plogis(...)` (the legacy
`binary_missing` table), and passed `oracle_link = "identity"` to match.
`6b06db3` ("remove historical bug flags", 2026-09-03 13:32) deleted that table and
rebuilt `binary_missing_fixed` on `continuous_missing`, whose formulas are plain
linear predictors, but left the link at `"identity"`. From then on the oracle's
outcome predictions for a 0/1 outcome were log-odds. The propensity (0.5) is
known, so the pseudo-outcome stayed unbiased whatever outcome model it was given
— the arm was not biased, it just stopped being an oracle: pseudo-outcome
variance 4–6× the true oracle's, and on single runs a CATE that could correlate
*negatively* with the truth (−0.23 and −0.15 on grid rows 9901 and 9918).

**Fix:** `bin_miss_models.R` passes `oracle_link = "logit"`, as `binary/` does.
The oracle's outcome model now reproduces the true `p0`/`p1` exactly, and
scenario 2's `dr_oracle` is identical to `binary/` scenario 2's from the same seed.
*(Since the risk-difference DGM every oracle formula returns the outcome mean
and the `oracle_link` argument is gone, so there is no link left to mismatch.)*

**Nothing re-runs for it.** `check_all_studies.md` had this study at
9,900/9,900, HTE back-fill complete, at 11:50 on 2026-09-03 — before `6b06db3`
existed. And the fix restores the old behaviour exactly: grid rows 5, 15, 64 and
91 (null + complete cases, mean imputation, MNAR-Y + IPW, complete data) run
through `0df4a9b`, the last commit before the bug, and through the fixed code
give byte-identical `truth`, `tau` for all five arms, and `dr_oracle$po`. So rows
1–9900 stay valid and directly comparable with scenario 2 run now. The only way
that could be wrong is a result file written from code at or after `6b06db3`,
which the file times would show — on the cluster, from the repo root (those
results predate the renumbering, so scenario 2's directory is still
`scenario_6`):

```bash
find ../results/missing/binary -path '*scenario_6*' -prune -o \
  -name 'res_sim_*.RDS' -newermt '2026-09-03 13:32' -print | wc -l    # expect 0
```

## Bug N — the MNAR-Y truth was taken at U = 0 (repaired at metrics time; no re-run)

**Moot since the risk-difference DGM.** U now enters the treated risk as
`bU·tanh(U)`, which has mean zero, so the U-free truth *is* the average over U,
as on the continuous scale. `mnar_y_truth()` and `repair_mnar_y_truth()` are
deleted, `bin_miss_metrics.R` no longer repairs anything, and truths no longer
carry `tau_u0`. The account below is of the logit-scale design.

Under MNAR-Y the unobserved `U` enters the treated arm's linear predictor,
`lp = base + W·(te + bU·U)`. The truth removed it by evaluating at U = 0,
`p1 = plogis(base + te)`, on the reasoning that E[U] = 0. That holds on the
identity scale (`missing/continuous` is fine) but not on the logit scale. What
the data identify, and what every estimator targets since U is unobserved, is
the risk averaged over U, `E_U[plogis(base + te + bU·U)]`, which sits closer to
0.5. At n = 500:

| scenario | mean τ at U = 0 | mean τ averaged over U | mean \|gap\| | max \|gap\| |
|---|---|---|---|---|
| 2 | −0.110 | −0.087 | 0.025 | 0.039 |
| 3 | −0.066 | −0.052 | 0.024 | 0.039 |
| 4 | −0.044 | −0.030 | 0.021 | 0.039 |
| 5 | −0.165 | −0.137 | 0.029 | 0.039 |

The data agree. Pooling 4,000 generated datasets (2M rows), the treated-minus-
control difference in means matches the averaged truth (scenario 5: −0.1383 vs
−0.1383) and is 42 standard errors from the U = 0 one (−0.1667); scenario 2
likewise (32 SE). So every binary MNAR-Y result so far — every arm,
`complete_data` included — carries about +0.02 to +0.03 of bias that belongs to
the truth, not the estimator, and the other CATE metrics and the true-CATE tests
were scored against the same wrong target.

**Fix:** `generate_scenario_data()` builds the binary MNAR-Y truth through
`mnar_y_truth()` (`R/dgm_scenarios.R`), which averages over U by Gauss–Hermite
quadrature (80 nodes, within 5e-12 of `integrate()`). It draws no random
numbers, so the generated data and every fit are unchanged. The old U = 0 value
is kept as `truth$tau_u0`.

**No re-run.** The saved truth still pins down the right one: `p0` is unaffected
(U never enters the control arm) and `qlogis(p1)` is the U-free linear predictor,
so `repair_mnar_y_truth()` rebuilds it exactly. `bin_miss_metrics.R` applies it
to the collected results before computing anything. It skips truths that already
carry `tau_u0`, so a collection mixing runs from both sides of the fix is fine.
Re-run `bin_miss_metrics.sh` and the figures; no simulation jobs. Checked:
repaired pre-fix truths equal the fixed generator's (differences ~1e-16), and
with and without the repair every MAR and MNAR metric row is byte-identical.

**Not changed at the time:** `dr_oracle`'s outcome model under MNAR-Y was
still the U = 0 one. With the propensity known its pseudo-outcome was unbiased
for the averaged CATE regardless, only slightly noisier than a true oracle's.
Under the risk-difference DGM the U-free outcome model *is* the one averaged
over U, so `dr_oracle` is an exact oracle under MNAR-Y too.

## Every model carries the HTE tests

**Applies to `missing/continuous` too** — the cause was one field in
`profile = "missing"`, which both studies share. Written up once, here.

`dr_random_forest` used to record no BLP or independence test, so `BLP_p`,
`indep_cate` and `indep_po` were `NA` for that one model. The cause was a single
field in `PROFILES` (`R/cate_models.R`):

```r
missing = list(cf_variance = TRUE, tests = TRUE, dr_rf_tests = FALSE)
```

Under `profile = "missing"` the arm was built inline as
`list(tau = stage2_whole_rf(...)$tau)` and never reached the test block, while
the base studies called `run_dr_random_forest()` and got both tests. Nothing
marked this as deliberate — it was copy-paste drift from the file this study was
forked from.

**Decision: all models should carry the HTE tests in this study where possible.**
`PROFILES$missing` now sets `dr_rf_tests = TRUE`, so anything run from here on
produces them natively.

### History: the back-fill of the archived results

The results made before the flag changed were not re-run for it. Both tests are
deterministic and every input survived in the saved files, so a one-off patch
recomputed the three fields in place, and a verification script showed it
reproduced a re-run exactly. Those results are now archived: the re-run
writes the tests natively and has no patch step. The patch, its audit and its
jobscripts were removed after commit `e7b1d59`; check that commit out to read
them. Two things carry over for anyone comparing against the archived tree:

- **Patched files have no `dr_random_forest$variance`.** It could not be
  recovered without refitting. Nothing in this study reads it.
- **The archived `multiple_imputation` runs have no HTE tests for any model.**
  They kept no nuisances, so the patch refused them. The re-run saves each
  imputation's tests; how to pool them is still open. See "Open decision:
  pooling the `multiple_imputation` arm's HTE tests" in `missing/README.md`.

## True-CATE HTE test evaluation

`bin_miss_metrics.R` also writes `bin_miss_true_cate_tests.RDS` — the BLP and
independence tests run on the true CATE and true nuisances (`truth$tau`,
`truth$p0`, `W.hat = 0.5`) instead of an estimator's fitted ones, one
`BLP_p`/`BLP_p_os`/`indep_cate` row per (scenario, n, type, prop, mechanism, method,
run). See `sample_size/continuous/README.md`'s "True-CATE HTE test evaluation" for what
it means and why every scenario-1 true-CATE test is `NA`, and its "HTE tests"
for the test columns themselves (`BLP_p_os` is the one-sided, HC3 BLP). `method == "multiple_imputation"`
rows are `NA`/`NA` too, for the same reason as above: `data` is a list of 50
imputed data.frames there, with no single covariate matrix to test against.

## Running it

```bash
qsub missing/binary/jobscripts/bin_miss_1.sh        # 1-9900
qsub missing/binary/jobscripts/bin_miss_extra.sh    # 9901-12600, scenario 2
Rscript missing/binary/bin_miss_check.R
qsub missing/binary/jobscripts/bin_miss_collect.sh
qsub missing/binary/jobscripts/bin_miss_metrics.sh
```

There is no patch step: every model's HTE tests are written by the simulation
itself, and the MI arm saves its per-imputation tests (`mi_tests`, see
`missing/README.md`).

Then, for the figures:

```bash
Rscript missing/binary/bin_miss_results.R           # every metric, to results/all_figures/
quarto render missing/binary/bin_miss_results.qmd   # the same, as a browsable report
```

## Status

**Full re-run owed for bug P and the risk-difference DGM.** Bug P changed the
`bW` calibration, and the risk-difference DGM then changed the outcome model and
every modifier (see "Bug P and the risk-difference DGM" above). Together they
move every dataset's outcome and the truth in all six scenarios. That supersedes
the complete rows 1–9900 described below: all 12,600 rows (`bin_miss_1.sh`,
`bin_miss_extra.sh`) re-run on the risk-difference code, then collect and
metrics. Any rows already re-run on the logit-scale DGM are superseded too. New
runs carry the `dr_random_forest` HTE tests (`PROFILES$missing`), the MI arm's
per-imputation tests, and a truth that needs no bug N repair, so nothing needs
back-filling. There is no patch step. The old tree, including the back-fill's
`bin_miss_hte_patch/` manifests, was archived before this re-run with
`R/archive_old_results.R` (root `README.md`, Status, step 0). It uses the
pre-2026-09-26 numbers, so the new scenario 2 would otherwise have landed on
the old scenario 2's paths.

**Scenario 2 owed** — `bin_miss_extra.sh` (rows 9901–12600), then
collect/metrics. Rows 1–9900 are unchanged. *(Folded into the bug P
re-run above.)*

**Metrics owed (bug N)** — *superseded by the full re-run above, which needs no
repair.* Was: re-run `bin_miss_metrics.sh` so the MNAR-Y rows are scored against
the averaged truth.

**Rows 1–9900 were re-run** for the three DGM fixes, bug F and the crossfitting
change, and are complete — 9,900/9,900 with the HTE back-fill done, per
`check_all_studies.md` (2026-09-03). Bug M does not touch them; see above.

## Known issue found while profiling (fixed — see below)

A profiling run found `pretest_superlearner()` could return an empty SuperLearner library and crash downstream calls (bug K); fixed — it now falls back to `"SL.mean"` when every candidate algorithm fails a fold.
