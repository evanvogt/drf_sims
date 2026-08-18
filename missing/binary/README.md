# Missing covariates — binary outcome

The `missing/continuous` design on a logit scale. See `missing/README.md` for the
mechanisms, handling methods and the shared bug fixes.

| | |
|---|---|
| array | **9,900 jobs** |
| results | `../results/missing/binary/scenario_<k>/<n>/<type>/<prop>/<mechanism>/<method>/` |
| figures | `bin_miss_results.R` / `.qmd` — every metric, to `results/all_figures/`; the diagnostic counterpart to the chapter script `results_processing/thesis_figures/miss_bin.R` |

## ⚠ This file was a half-converted fork

`bin_miss_dgms.R` was copied from the continuous version and only partly
converted to a binary outcome. **Three things came across unchanged**, all fixed
together behind `LEGACY_BIN_MISS` in `R/dgm_scenarios.R`:

**1. The continuous coefficient table.** `b0 = c(0.4, 0.2, 0.4, 1, 0.4)`,
`b1 = -0.05`, `b2 = c(2, 2, 2, 2, 1)` — the continuous values — fed through
`plogis`. With `b2 = 2` on a logit scale the outcome model saturates. The binary
study uses `b0 = -0.4`, `b1 = 0.5`, `b2 = 0.5`.

**2. The wrong power calculation.** `bW` was calibrated with `power.t.test` on
`s_err + s2`, the two-sample t-test for a *continuous* outcome, rather than
`power.prop.test`. The table even carries an `s_err` column despite the outcome
being binary and no error term ever being drawn.

**3. Truth on the wrong scale — the most consequential.** `truth$p0` and `p1`
were computed as `b0 + b1*X1 + b2*X2` with **no `plogis`**, making `truth$tau` a
difference in **log-odds**. But the outcome is `rbinom(n, 1, plogis(lp))` and
every estimator targets a **risk difference**, and `bin_miss_metrics.R` compared
the two directly. Because `plogis` compresses, the two differ by a lot, not by
rounding — bias and MSE for this study were computed against a mis-scaled target.

The oracle was *not* affected: `get_binary_oracle_info()` here wraps its formula
in `plogis(...)`, which is why `bin_miss_models.R` passes
`oracle_link = "identity"` — the opposite convention to `binary/`. That is
correct, not a fourth bug.

### The corrected coefficients

`b0`, `b1`, `b2` come straight from the binary table. `b3`/`b4`/`b5`/`b45` are
taken from the binary scenario each reduced scenario corresponds to (1→1, 2→2,
3→4, 4→8, 5→9). **That mapping is an inference from the scenario descriptions,
not something the original code recorded** — worth a sanity check before
committing cluster time.

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
`check_all.R`'s `patchable_jobs` is 8,800 rather than 9,900.

## Running it

```bash
qsub missing/binary/jobscripts/bin_miss_1.sh        # 1-9900
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

**Everything under `../results/missing/binary/` is superseded.** All three fixes
change the generated data or the target, so the whole study re-runs. Bug F
affects `dr_superlearner` here as elsewhere, but that is moot given the DGM
changes.

## Known issue found while profiling (fixed — see below)

While smoke-testing `bin_miss_profile.R` (see that file's header) the DR
SuperLearner arm crashed for `scenario = 1, mechanism = MAR, run = 1` at both
`method = complete_cases` (n = 331 after complete-case removal) and
`method = mean_imputation` (n = 500, no rows dropped) - same seed, same crash,
two different sample sizes, so it is not purely a small-`n` artifact:

```
Removed libraries due to NA/error:
[1] "SL.glm"    "SL.glmnet" "SL.earth"  "SL.gam"    "SL.mean"   "SL.ranger"
Error in (function (.x, .f, ..., .progress = FALSE)  : ℹ In index: 7.
Caused by error in `data.frame()`:
! arguments imply differing number of rows: 0, 1
```

**Reproduces on the unmodified `bin_miss_analysis.R`** (confirmed by running
`Rscript bin_miss_analysis.R 1`, i.e. `study$grid`'s row 1) - it predates and is
unrelated to the profiling patches added alongside `bin_miss_profile.R`.

**Root cause**, read from `R/cate_models.R`: `pretest_superlearner()` drops
any candidate algorithm that errors, warns, or returns all-`NA` on a 2-fold
inner CV, and returns whatever survives - which could be `character(0)` if
every candidate failed on a given fold (exactly what the "Removed libraries"
list above shows: all six did). Both callers fed that result straight into
`SuperLearner(..., SL.library = <possibly empty>)` with no guard against the
empty case:
- `nuisance_sl()`, for `Y.hat`/`W.hat` per fold
- `stage_2_sl()`, for the pseudo-outcome regression - this is the one that
  crashed above (`PRETEST_STAGE2 = TRUE` means it always uses the pretested
  library now, per bug F in `R/README.md`)

This is a different failure mode from bug F, which was about *which* library
variant gets used once pretesting leaves at least one survivor. This is what
happens when pretesting leaves zero — **fixed as bug K** in `R/README.md`'s
ledger: `pretest_superlearner()` now falls back to `"SL.mean"` (asserted
directly rather than re-run through the pretest loop) whenever nothing else
survives, so it can never return `character(0)` again.

**Confirmed trigger:** `scenario = 1` (no `X3`/`X4`/`X5`, the covariate-poorest
scenario in `binary_missing`), `mechanism = MAR`, `run = 1`, `n_folds = 10`.
Also independently confirmed in `continuous/` — 23 of the 24 array IDs that
failed on the cluster there hit this exact bug, for several different
scenarios and sample sizes (not `missing/binary`-specific, or
small-`n`-specific).
