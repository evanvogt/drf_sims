# Calibrating the bootstrap sample fraction

The half-sample bootstrap takes a `sample.fraction` for the forests it refits on
each half sample. Too small and those forests are noisy; too large and the halves
stop being independent enough for the roots to behave. This folder picks a value
rather than guessing one.

| | |
|---|---|
| array | **2,000 jobs** per outcome |
| results | `../results/confidence_intervals/<outcome>/sf_calibration/` |

## Method

`find_optimal_sf()` (in `R/bootstrap_ci.R`) treats the study's own point
estimates as a plug-in truth:

1. take the residuals of the pseudo-outcomes around `tau.hat`
2. resample them within each fold column, preserving the NA pattern that makes
   column `k` valid for fold `k`
3. re-estimate `tau` from the resampled pseudo-outcomes
4. build the bootstrap interval at each candidate `sample.fraction`
5. record how often it covers the original `tau.hat`

The value whose mean coverage is closest to `1 - alpha` wins.

Because the target is a plug-in rather than the real CATE, this calibrates the
bootstrap's *self-consistency*, not its true coverage. The `CI_sf` sweep in
`confidence_intervals/{continuous,binary}/` is what measures actual coverage
against the known truth.

## Files

| file | role |
|---|---|
| `cts_ci_sf_calibration.R` | a shim onto `R/bootstrap_ci.R` |
| `cts_ci_sf_analysis.R` | continuous entry point |
| `bin_ci_sf_analysis.R` | binary entry point |

The `cts_` prefix on `cts_ci_sf_calibration.R` is historical — `find_optimal_sf`
is outcome-agnostic and the binary analysis has always sourced it. Now that it
lives in `R/bootstrap_ci.R` the misleading filename only survives as a shim.

## Status

The continuous calibration is unaffected.

**`bin_ci_sf_analysis.R` re-runs**: it sources
`confidence_intervals/binary/bin_ci_dgms.R`, which carried bug A, so its results
were produced under the continuous coefficient table.
