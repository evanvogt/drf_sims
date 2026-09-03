# Confidence intervals — binary outcome

The `confidence_intervals/continuous` design on a logit scale. See
`confidence_intervals/README.md` for the method.

| | |
|---|---|
| array | **20,000 jobs** (`bin_ci_1.sh` 1–10000, `bin_ci_2.sh` 10001–20000) |
| results | `../results/confidence_intervals/binary/scenario_<k>/<n>/<CI_sf>/` |

## ⚠ Bug A — fixed

This study initially ran on the continuous coefficient table instead of the binary table, invalidating results against the main binary study. Now fixed: the study always uses the correct binary coefficient table.

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
