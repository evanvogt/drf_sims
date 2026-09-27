# Confidence intervals — binary outcome

The `confidence_intervals/continuous` design with a binary outcome, on the
`binary/` study's DGM: the treatment effect is on the risk-difference scale (see
`binary/README.md`, "Outcome model and `bW` calibration"). See
`confidence_intervals/README.md` for the method.

| | |
|---|---|
| array | **20,000 jobs** (`bin_ci_1.sh` 1–10000, `bin_ci_2.sh` 10001–20000), plus **32,000** taking scenarios 1–4 to 500 runs (`bin_ci_extra_1.sh`–`bin_ci_extra_4.sh`, 20001–52000) |
| results | `../results/confidence_intervals/binary/scenario_<k>/<n>/<CI_sf>/` |

## ⚠ Bug A — fixed

This study initially ran on the continuous coefficient table instead of the binary table, invalidating results against the main binary study. Now fixed: the study always uses the correct binary coefficient table.

## Status

**Archive the old results first** - `R/archive_old_results.R` (root `README.md`, Status, step 0). They predate the current DGM and use the pre-2026-09-26 scenario numbers, so running into that tree would mix old and new results.

**Everything under `../results/confidence_intervals/binary/` is superseded** —
different data, so every number moves. Bug P (the binary `bW` calibration and
the modifier signs in scenarios 2, 5, 6 and 7) and then the risk-difference DGM
(see `binary/README.md`) change the data again, so run the re-run below on the
risk-difference code, not before it. With the effect on the risk-difference
scale, the query grid's true CATE no longer depends on where X1 and X2 are held
(`GRID_REFERENCE_VALUE`).

`confidence_intervals/optimal_sf/bin_ci_sf_analysis.R` sources this same DGM, so
its 2,000 jobs re-run too. Together that is about 22,000 array jobs, the largest
single item in the re-run bill.

```bash
qsub sample_size/confidence_intervals/binary/jobscripts/bin_ci_1.sh
qsub sample_size/confidence_intervals/binary/jobscripts/bin_ci_2.sh
for k in 1 2 3 4; do qsub sample_size/confidence_intervals/binary/jobscripts/bin_ci_extra_$k.sh; done
qsub sample_size/confidence_intervals/optimal_sf/jobscripts/bin_ci_sf_1.sh
```
