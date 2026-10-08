# Validation

If a subgroup, a CATE ranking of covariates, or the spread of estimated
treatment effects looks real partway through a trial, does it still look real
once more data has accrued? These studies fit CATE estimators on the first
`interim_prop` of a trial and check whether what they find there — top/bottom
responder subgroups, the variance of the estimated CATEs, the ranking of
covariates by variable importance (TE-VIM and surrogate TreeSHAP), and whether
the single most important covariate still interacts with treatment —
replicates on the remaining chunk.

This is a robustness check for CATE-based subgroup discovery, not an
estimator comparison — it exists nowhere else in the repo.

| folder / file | |
|---|---|
| `continuous/` | continuous outcome |
| `binary/` | binary outcome, risk-difference scale |
| `val_common.R` | what both arms share — see below |

Both arms run the same design: scenario 2 on the correlated (rho = 0.5) set of
their outcome, one trial of 1000 split at eleven interim points, 100 runs, the
same three estimators and the same four chunk comparisons. A
`competing_risk/` sibling would slot in the same way: its own
`<prefix>_val_*.R` file split, `jobscripts/` and README, sourcing
`val_common.R`.

## Shared code

`val_common.R` holds everything the arms do identically:

- `split_trial()` and `chunk_folds()` — the one-trial split, and the DR
  SuperLearner's folds per chunk
- `fit_val_methods()` — the three estimators and both importance measures on
  one chunk
- the TE-VIM, surrogate-TreeSHAP and interaction-test helpers
- `chunk_validations()` — the four chunk comparisons

Each arm's `<prefix>_val_models.R` sources it and fixes the outcome family
(`gaussian()` / `binomial()`). Both arms' `<prefix>_val_analysis.R` call
`chunk_validations(robust = TRUE)`, so every interaction test uses HC3 standard
errors — see `continuous/README.md` and `binary/README.md` for why each needs
them. This code moved out of `continuous/` unchanged when the binary arm was
added (2026-10-08); the continuous arm's comparisons were checked identical
before and after, while it still used classical standard errors. It switched
to HC3 later the same day.

See each arm's README for its design, file roles and status.
