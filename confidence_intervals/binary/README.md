# Confidence intervals — binary outcome

The `confidence_intervals/continuous` design on a logit scale. See
`confidence_intervals/README.md` for the method.

| | |
|---|---|
| array | **20,000 jobs** (`bin_ci_1.sh` 1–10000, `bin_ci_2.sh` 10001–20000) |
| results | `../results/confidence_intervals/binary/scenario_<k>/<n>/<CI_sf>/` |

## ⚠ Bug A — this study ran on the continuous coefficient table

`bin_ci_dgms.R` declared its scenario parameters as
`b0 = c(0.4, 0.2, 0.3, ...)`, `b1 = -0.05`, `b2 = c(2, 2, ...)`, plus the
`s2`/`s4`/`s5` variance parameters — **the continuous table** — and then fed them
through `plogis` to generate a binary outcome. The main binary study uses
`b0 = -0.4`, `b1 = 0.5`, `b2 = 0.5`.

With `b2 = 2` on a logit scale the outcome model saturates: most of the
covariate range is pushed into the flat tails of `plogis`, so the true CATEs are
compressed towards zero and coverage is being measured against a target the
design never intended.

Scenario 10's treatment effect and the `needs_X3` flags came across from the
continuous file too.

Unlike `missing/binary`, the rest of this file was correctly converted — `bW`
uses `power.prop.test` and truth is computed with `plogis`. The coefficients were
the only thing wrong.

### The fix

`bin_ci_dgms.R` now names the scenario set `"binary_ci"`, which resolves through
`LEGACY_BIN_CI_PARAMS` in `R/dgm_scenarios.R`:

- `TRUE` → the continuous table (what produced the existing results)
- `FALSE` → the binary table (the default now)

With the flag `FALSE` this study's DGM is identical to `binary/`, which is the
point: the CI results become comparable with the main binary study.

## Status

**Everything under `../results/confidence_intervals/binary/` is superseded** —
different data, so every number moves.

`confidence_intervals/optimal_sf/bin_ci_sf_analysis.R` sources this same DGM, so
its 2,000 jobs re-run too. Together that is about 22,000 array jobs, the largest
single item in the re-run bill.

```bash
qsub confidence_intervals/binary/jobscripts/bin_ci_1.sh
qsub confidence_intervals/binary/jobscripts/bin_ci_2.sh
qsub confidence_intervals/optimal_sf/jobscripts/bin_ci_sf_1.sh
```
