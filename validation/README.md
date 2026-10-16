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
estimator comparison — it exists nowhere else in the repo. Full design, both
arms, in `ADEMP.md`.

| folder / file | |
|---|---|
| `continuous/` | continuous outcome |
| `binary/` | binary outcome, risk-difference scale |
| `val_common.R` | what both arms share — see below |

Both arms run the same design: scenarios 1–4 (null, simple, complex,
non-linear) on the correlated (rho = 0.5) set of their outcome, one trial of
1000 split at eleven interim points, 100 runs, the same three estimators and
the same four chunk comparisons — 4400 array jobs per arm. In the null every
significant interaction test is a type I error, so it is the size reference for
every replication rate in the other three. A
`competing_risk/` sibling would slot in the same way: its own
`<prefix>_val_*.R` file split, `jobscripts/` and README, sourcing
`val_common.R`.

## Shared code

`val_common.R` holds everything the arms do identically:

- `split_trial()` and `chunk_folds()` — the one-trial split, and the DR
  SuperLearner's folds per chunk
- `fit_val_methods()` — the three estimators and both importance measures on
  one chunk
- the TE-VIM, surrogate-TreeSHAP and interaction-test helpers, including the
  two data-driven top-covariate tests (`data_driven_split()`,
  `projection_score()`), and the reference contrasts each test is read against
  (`tau_slope()`, `tau_slope_adj()`)
- `chunk_validations()` — the four chunk comparisons, for one kind of standard
  error
- `chunk_validations_se()` — `chunk_validations()` with classical and with HC3
  standard errors, stacked under an `se` column

Each arm's `<prefix>_val_models.R` sources it and fixes the outcome family
(`gaussian()` / `binomial()`). Both arms' `<prefix>_val_analysis.R` call
`chunk_validations_se()`, so every interaction test is saved with both
classical and HC3 standard errors. HC3 is the primary kind — see
`continuous/README.md` and `binary/README.md` for why each arm needs it.

Every interaction test also saves its `W:v` estimate and
standard error (`est_<test>`, `se_<test>`), and two reference contrasts:
`ref1_<test>`, the same contrast in chunk 1's estimated CATE (what chunk 1
claimed), and `true2_<test>`, the same contrast in chunk 2's true CATE (what
the test estimates). The signs say whether a replication went the right way.
See `ADEMP.md`, Methods, and either arm's README.

Each arm's report (`<prefix>_val_results.qmd`) renders one scenario at a time,
through a Quarto `scenario` parameter — see each arm's README.

See each arm's README for its design and file roles.
