# Missing covariates — binary outcome

The `missing/continuous` design on a logit scale. See `missing/README.md` for the
mechanisms, handling methods and the shared bug fixes.

| | |
|---|---|
| array | **12,600 jobs**: `bin_miss_1.sh` (1–9900, scenarios 1, 2, 4, 5) and `bin_miss_extra.sh` (9901–12600, scenario 6) |
| results | `../results/missing/binary/scenario_<k>/<n>/<type>/<prop>/<mechanism>/<method>/` |
| figures | `bin_miss_results.R` / `.qmd` — every metric, to `results/all_figures/`; the diagnostic counterpart to the chapter script `results_processing/thesis_figures/miss_bin.R` |

## ⚠ This file was a half-converted fork

`bin_miss_dgms.R` was copied from the continuous version and only partly converted to a binary outcome, carrying three related defects (continuous coefficient table, wrong power-test calibration, un-plogis'd truth) — all fixed together; the code now always uses the corrected values. Its oracle link convention broke later, when the fix flags were removed — see bug M below.

### The corrected coefficients

`b0`, `b1`, `b2` come straight from the binary table. `b3`/`b4`/`b45` are
taken from the binary scenario each reduced scenario corresponds to (1→1, 2→2,
3→4, 4→8, 5→9). **That mapping is an inference from the scenario descriptions,
not something the original code recorded** — worth a sanity check before
committing cluster time. Scenario 6 (→3) was added later specifically as binary
scenario 3, so its coefficients are copied, not inferred.

### Bug P — the `bW` calibration and the modifier signs

The binary calibration and coefficients changed with bug P; `binary/README.md`'s
"Outcome model and `bW` calibration" has the design. Here it means:

- **Signs.** Scenarios 2, 3 and 6 carry the sign flips of main-study scenarios
  2, 4 and 3: `b3 = 0.4` (was −0.4) in 2 and 3, `b4 = −0.3` (was 0.3) in 3, and
  `b4 = −0.2` (was 0.2) in 6. `b5`, which no treatment effect used, is gone.
- **`bW`** now differs by scenario, so that the true ATE is the same marginal
  RD in all six: at n = 500, −0.55, −0.84, −0.85, −0.75, −0.86 and −0.56 for
  scenarios 1–6 (−0.50 in every scenario before). The true RD is −0.122 and
  the power 0.80.
- **MNAR-Y** is calibrated without U, as in `missing/continuous`, so every
  mechanism shares one `bW` and one truth per scenario. Averaging the treated
  risk over U pulls it towards 0.5, so under MNAR-Y the RD shrinks to about
  −0.10 and the power to 0.61–0.63.

`Rscript R/calibration_report.R` prints this table.

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
scenario 6's `dr_oracle` is identical to `binary/` scenario 3's from the same seed.

**Nothing re-runs for it.** `check_all_studies.md` had this study at
9,900/9,900, HTE back-fill complete, at 11:50 on 2026-09-03 — before `6b06db3`
existed. And the fix restores the old behaviour exactly: grid rows 5, 15, 64 and
91 (null + complete cases, mean imputation, MNAR-Y + IPW, complete data) run
through `0df4a9b`, the last commit before the bug, and through the fixed code
give byte-identical `truth`, `tau` for all five arms, and `dr_oracle$po`. So rows
1–9900 stay valid and directly comparable with scenario 6 run now. The only way
that could be wrong is a result file written from code at or after `6b06db3`,
which the file times would show — on the cluster, from the repo root:

```bash
find ../results/missing/binary -path '*scenario_6*' -prune -o \
  -name 'res_sim_*.RDS' -newermt '2026-09-03 13:32' -print | wc -l    # expect 0
```

## Bug N — the MNAR-Y truth was taken at U = 0 (repaired at metrics time; no re-run)

Under MNAR-Y the unobserved `U` enters the treated arm's linear predictor,
`lp = base + W·(te + bU·U)`. The truth removed it by evaluating at U = 0,
`p1 = plogis(base + te)`, on the reasoning that E[U] = 0. That holds on the
identity scale (`missing/continuous` is fine) but not on the logit scale. What
the data identify, and what every estimator targets since U is unobserved, is
the risk averaged over U, `E_U[plogis(base + te + bU·U)]`, which sits closer to
0.5. At n = 500:

| scenario | mean τ at U = 0 | mean τ averaged over U | mean \|gap\| | max \|gap\| |
|---|---|---|---|---|
| 2 | −0.165 | −0.137 | 0.029 | 0.039 |
| 4 | −0.066 | −0.052 | 0.024 | 0.039 |
| 5 | −0.044 | −0.030 | 0.021 | 0.039 |
| 6 | −0.110 | −0.087 | 0.025 | 0.039 |

The data agree. Pooling 4,000 generated datasets (2M rows), the treated-minus-
control difference in means matches the averaged truth (scenario 2: −0.1383 vs
−0.1383) and is 42 standard errors from the U = 0 one (−0.1667); scenario 6
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

**Not changed:** `dr_oracle`'s outcome model under MNAR-Y is still the U = 0 one.
With the propensity known its pseudo-outcome is unbiased for the averaged CATE
regardless; it is only slightly noisier than a true oracle's. Making it exact
would change `dr_oracle` and mean redoing these runs.

## Patched: every model now carries the HTE tests

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

### The finished results did not need re-running

The obvious reading is that 19,800 completed runs are now wrong and owe ~10,000
CPU-hours of re-runs. They are not, and they do not. Both tests are
deterministic — `GenericML::BLP` is an OLS with a sandwich vcov, and
`coin::independence_test` under `teststat = "quadratic"` uses the asymptotic
null — and every input survives in the saved file (`nuisances_rf$W.hat`,
`$Y0.hat`, `$po`, `data`, and the arm's own `tau`). So the three fields were
recomputed exactly, in place, by `R/patch_hte_tests.R`:

```bash
Rscript missing/binary/bin_miss_patch.R dry     # report only, writes nothing
qsub    missing/binary/jobscripts/bin_miss_patch.sh
```

The array index is a row of `combos(study)` — one parameter combination and its
100 files, 99 of them — not a row of `study$grid`. The job is idempotent
(already-patched files are skipped) and each write is `saveRDS` to `.tmp`
followed by a rename, so it can be re-run and cannot leave a truncated result.
It writes a manifest to `../results/missing/binary/bin_miss_hte_patch/`, which
is what `check_all.R` reads for its `patch_status` column.

`missing/patch_hte_verify.R` is the evidence. It runs one grid row twice from
the same seed, with `dr_rf_tests` off and on, and checks that (1) **every** arm's
`tau` is byte-identical — so the added tests draw no random numbers and cannot
perturb `dr_oracle` / `dr_semi_oracle` / `dr_superlearner`, which are fitted
afterwards; (2) the patch reproduces the re-run's three fields exactly; and
(3) `dr_random_forest$independence_po` equals `causal_forest$independence_po`,
which catches a mis-rebuilt `X`. Run it before trusting the patch:

```bash
cd missing
Rscript patch_hte_verify.R binary 2      # complete_cases, all five arms
```

### Two things the patch does not do

**`dr_random_forest$variance` is gone from patched files.**
`run_dr_random_forest()` also returns `variance` from the stage-2 forest, which
the old inline branch discarded and which cannot be recovered without refitting.
Patched files therefore lack it while newly run ones have it. Nothing in this
study reads it — `bin_miss_metrics.R` calls only `cate_metrics()` and
`hte_test_metrics()`, and there are no interval metrics here — so the asymmetry
is recorded rather than hidden.

**The `multiple_imputation` arm still has no HTE tests, for any model.** That is
a separate and larger gap; see the multiple-imputation note in
`missing/README.md`. The patch detects those runs and refuses them, which is why
`check_all.R`'s `patchable_jobs` is 11,200 rather than 12,600.

## True-CATE HTE test evaluation

`bin_miss_metrics.R` also writes `bin_miss_true_cate_tests.RDS` — the BLP and
independence tests run on the true CATE and true nuisances (`truth$tau`,
`truth$p0`, `W.hat = 0.5`) instead of an estimator's fitted ones, one
`BLP_p`/`indep_cate` row per (scenario, n, type, prop, mechanism, method,
run). See `continuous/README.md`'s "True-CATE HTE test evaluation" for what
it means and why scenario 1's `BLP_p` is `NA`. `method == "multiple_imputation"`
rows are `NA`/`NA` too, for the same reason as above: `data` is a list of 50
imputed data.frames there, with no single covariate matrix to test against.

## Running it

```bash
qsub missing/binary/jobscripts/bin_miss_1.sh        # 1-9900
qsub missing/binary/jobscripts/bin_miss_extra.sh    # 9901-12600, scenario 6
Rscript missing/binary/bin_miss_check.R
qsub missing/binary/jobscripts/bin_miss_patch.sh    # 1-99, the HTE back-fill
Rscript missing/binary/bin_miss_patch_check.R       # did the back-fill land?
qsub missing/binary/jobscripts/bin_miss_collect.sh
qsub missing/binary/jobscripts/bin_miss_metrics.sh
```

The patch step goes **before** collect: collect reads the per-run files into
`bin_miss_all.RDS`, so patching afterwards would leave the collected copy
carrying the old, testless `dr_random_forest`. It is a one-off — once these
results are patched and `PROFILES$missing` is set, future runs need only the
usual four steps.

### Checking the back-fill landed

`bin_miss_patch_check.R` is to the repair what `bin_miss_check.R` is to the
simulation, and it exists because the first submission of `bin_miss_patch.sh`
lost ten of its 99 array elements without leaving a trace. `check_all.R` showed
the study at **patchable 8,800 / patched 7,800** — it counts manifest *rows*, and
combos 21, 22, 28, 29, 30, 32, 40, 44, 45 and 48 had written no manifest at all,
which from there is indistinguishable from never having been submitted. Neither
end of the job could say more: the HPC was returning no `.e` files and that run's
`.o` files were gone.

So the audit reconstructs the diagnosis from the result files, which do survive.
`patch_status_of()` says whether a file was patched; mtimes say *when*, so the
first and last patched file bracket how long the element ran before it stopped;
and an orphan `res_sim_<n>.RDS.tmp` names the file it was writing when it died.
Comparing the failed elements' runtimes against the ones that finished, and
against the job's own walltime, is what turns that into a cause.

```bash
Rscript missing/binary/bin_miss_patch_check.R     # writes failed_patch_ids.txt
qsub    missing/binary/jobscripts/bin_miss_patch_rerun.sh
Rscript missing/binary/bin_miss_patch_check.R     # confirm, then check_all.R
```

The re-run is safe over combinations that are already correct — the patch is
idempotent, so an element that only lost its manifest simply rewrites it. The
audit writes `bin_miss_patch_check.{csv,md}`, which are committed, so `git push`
from the HPC is how the answer leaves the cluster.

Three things changed so this cannot recur silently. `patch_hte_tests()` now
flushes each combination's manifest every `MANIFEST_FLUSH_EVERY` files instead
of only at the end, and prints one line per file to stdout, so a killed element
leaves both a short manifest and a log saying where it stopped. `bin_miss_patch.sh`
gained `-j oe`, merging stderr into the `.o` file that the HPC does return, and a
`%20` throttle on its array — it was the only job in this study without one, and
99 elements each doing 100 `readRDS` + 100 `saveRDS` against the shared
filesystem at once is the leading explanation for the original kills.

Then, for the figures:

```bash
Rscript missing/binary/bin_miss_results.R           # every metric, to results/all_figures/
quarto render missing/binary/bin_miss_results.qmd   # the same, as a browsable report
```

## Status

**Full re-run owed for bug P** — the `bW` calibration and the signs in
scenarios 2, 3 and 6 changed (see "Bug P" above), which moves every dataset's
outcome and the truth in all six scenarios. That supersedes the complete rows
1–9900 described below: all 12,600 rows (`bin_miss_1.sh`, `bin_miss_extra.sh`)
re-run, then collect and metrics. New runs carry the `dr_random_forest` HTE
tests (`PROFILES$missing`) and the averaged MNAR-Y truth, so nothing needs
back-filling and the bug N repair leaves them alone. `check_all.R`'s
`patch_status` still counts manifest rows, though, so, as for scenario 6 in
`missing/README.md`'s Status, a bookkeeping pass over every combination
(`qsub -J 1-126%20 jobscripts/bin_miss_patch.sh`) is what makes it read
complete. Archive the old `bin_miss_hte_patch/` manifests with the old results
first.

**Scenario 6 owed** — `bin_miss_extra.sh` (rows 9901–12600), then the
bookkeeping patch pass over combinations 100–126 (see `missing/README.md`
Status) and collect/metrics. Rows 1–9900 are unchanged. *(Folded into the bug P
re-run above.)*

**Metrics owed (bug N)** — re-run `bin_miss_metrics.sh` so the MNAR-Y rows are
scored against the averaged truth. No simulation jobs; doing it once, after
scenario 6's collect, covers both.

**Rows 1–9900 were re-run** for the three DGM fixes, bug F and the crossfitting
change, and are complete — 9,900/9,900 with the HTE back-fill done, per
`check_all_studies.md` (2026-09-03). Bug M does not touch them; see above.

## Known issue found while profiling (fixed — see below)

A profiling run found `pretest_superlearner()` could return an empty SuperLearner library and crash downstream calls (bug K); fixed — it now falls back to `"SL.mean"` when every candidate algorithm fails a fold.
