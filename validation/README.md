# Validation

If a subgroup, a CATE ranking of covariates, or the spread of estimated
treatment effects looks real partway through a trial, does it still look real
once more data has accrued? These studies fit CATE estimators on the first
`interim_prop` of a trial and check whether what they find there — top/bottom
responder subgroups, the variance of the estimated CATEs, and the ranking of
covariates by treatment-effect variable importance (TE-VIM) — replicates on
the remaining chunk.

This is a robustness check for CATE-based subgroup discovery, not an
estimator comparison — it exists nowhere else in the repo.

| folder | |
|---|---|
| `continuous/` | continuous outcome |

Only a continuous outcome is validated today. A `binary/` or
`competing_risk/` sibling would slot in the same way `missing/binary/` and
`confidence_intervals/binary/` sit alongside their `continuous/` — same
`<prefix>_val_*.R` file split, its own `jobscripts/`, its own README.

See `continuous/README.md` for the design, estimators, file roles and status.
