# Missing covariates — binary outcome

The `missing/continuous` design on a logit scale. See `missing/README.md` for the
mechanisms, handling methods and the shared bug fixes.

| | |
|---|---|
| array | **9,900 jobs** |
| results | `../results/missing/binary/scenario_<k>/<n>/<type>/<prop>/<mechanism>/<method>/` |

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

## Running it

```bash
qsub missing/binary/jobscripts/bin_miss_1.sh        # 1-9900
Rscript missing/binary/bin_miss_check.R
qsub missing/binary/jobscripts/bin_miss_collect.sh
qsub missing/binary/jobscripts/bin_miss_metrics.sh
```

## Status

**Everything under `../results/missing/binary/` is superseded.** All three fixes
change the generated data or the target, so the whole study re-runs. Bug F
affects `dr_superlearner` here as elsewhere, but that is moot given the DGM
changes.

## Known issue found while profiling (not yet fixed)

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

**Root cause**, read from `R/cate_models.R`, not yet fixed there:
`pretest_superlearner()` (`R/cate_models.R:366`) drops any candidate algorithm
that errors, warns, or returns all-`NA` on a 2-fold inner CV, and returns
whatever survives - which can be `character(0)` if every candidate fails on a
given fold (exactly what the "Removed libraries" list above shows: all six
did). Both callers feed that result straight into `SuperLearner(..., SL.library
= <possibly empty>)` with no guard against the empty case:
- `nuisance_sl()` (`R/cate_models.R:306`), for `Y.hat`/`W.hat` per fold
- `stage_2_sl()` (`R/cate_models.R:408`), for the pseudo-outcome regression -
  this is the one that crashed above (`PRETEST_STAGE2 = TRUE` means it always
  uses the pretested library now, per bug F in `R/README.md`)

This is a different failure mode from bug F, which was about *which* library
variant gets used once pretesting leaves at least one survivor. This is what
happens when pretesting leaves zero. Worth a letter of its own once someone
fixes it - candidate "bug H" in `R/README.md`'s ledger (bug G is already
taken).

**Confirmed trigger:** `scenario = 1` (no `X3`/`X4`/`X5`, the covariate-poorest
scenario in `binary_missing`), `mechanism = MAR`, `run = 1`, `n_folds = 10`.
Not yet checked against other scenarios/mechanisms/runs in this study, or
against the other studies that call `pretest_superlearner()` (`binary/`,
`continuous/`, `crossfitting/`, `case_study/`) - flagged here because this is
where it surfaced during profiling, not because it is known to be
missing/binary-specific.
