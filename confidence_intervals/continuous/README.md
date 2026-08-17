# Confidence intervals — continuous outcome

The reference CI study. See `confidence_intervals/README.md` for the half-sample
bootstrap and what the metrics mean.

| | |
|---|---|
| array | **20,000 jobs** (`cts_ci_1.sh` 1–10000, `cts_ci_2.sh` 10001–20000) |
| results | `../results/confidence_intervals/continuous/scenario_<k>/<n>/<CI_sf>/` |
| metrics | `ci_cts_metrics.RDS` — coverage and width, no point metrics |

## Files

| file | role |
|---|---|
| `cts_ci_config.R` | the grid, including the `CI_sf` sweep |
| `cts_ci_dgms.R` | names the `continuous` scenario set — the same one `continuous/` uses |
| `cts_ci_models.R` | `profile = "ci"`, bootstrap on, SuperLearner and HTE tests off |
| `cts_ci_analysis.R` | array entry point |
| `cts_ci_results.R`, `cts_ci_results.Rmd`, `cts_ci_results.qmd` | summaries |

`cts_ci_dgms.R` used to carry the note *"just a straight copy of the cts dgm for
now"*, and it was — byte for byte. It now names the shared scenario set instead,
so the two studies cannot drift apart.

## Gotchas

**`CI_boot = 200` bootstrap draws per run, on top of 10 folds.** This study is
the most expensive in the repo per job; the `CI_sf` sweep multiplies it by ten.

**The nuisance object saved here is smaller than in the other studies.** The
`ci` profile does not compute the row-mean summaries (`po`, `Y.hat`, `W.hat`),
because only the BLP and independence tests consume them and this study runs
neither. Anything reading `nuisances_rf$po` from these files will get `NULL`.

**`dr_random_forest` is built inline** rather than through
`run_dr_random_forest()`, so it carries no BLP or independence test. Consistent
with the rest of the `ci` profile.

## Status

Unaffected by the bug ledger. Bug A is binary-only; bug F needs an `sl_lib`,
which this study never passes. **No re-run needed** — this is one of the two
studies that must stay bit-identical, and the regression harness treats it as
such.
