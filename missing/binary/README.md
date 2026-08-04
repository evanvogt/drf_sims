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

## Status

**Everything under `../results/missing/binary/` is superseded.** All three fixes
change the generated data or the target, so the whole study re-runs. Bug F
affects `dr_superlearner` here as elsewhere, but that is moot given the DGM
changes.
