# Correlated covariates — confidence interval studies

`confidence_intervals/`'s three studies, run on `correlated/`'s copula DGM.
They ask whether the half-sample bootstrap's simultaneous bands, and the
data-driven choice of their `sample.fraction`, still behave when prognosis and
effect modification move together. The method, estimators and metrics are the
parent studies' (`../../confidence_intervals/README.md`). The DGM and pairing
are `correlated/`'s (`../README.md`).

| folder | parent | prefix |
|---|---|---|
| `continuous/` | `confidence_intervals/continuous/` | `cts_corr_ci_*` |
| `binary/` | `confidence_intervals/binary/` | `bin_corr_ci_*` |
| `optimal_sf/` | `confidence_intervals/optimal_sf/` | `cts_corr_ci_sf_*`, `bin_corr_ci_sf_*` |

## Design

| | `continuous/`, `binary/` | `optimal_sf/` (cts, bin) |
|---|---|---|
| scenario sets | `continuous_corr_<ρ>`, `binary_corr_<ρ>` | same |
| scenarios | 1–4 | 1–4 |
| ρ | 0, 0.5 (`CORR_RHOS`) | same |
| n | 500, 1000 | 500, 1000 |
| `CI_sf` | 0.05 to 0.5 by 0.05 (the parent's full sweep) | picked per run from the same values |
| runs | 100 | 100 |
| array per outcome | 16,000: rows 1–8000 are ρ = 0, 8001–16000 ρ = 0.5 | 1,600: rows 1–800 are ρ = 0, 801–1600 ρ = 0.5 |
| jobscripts | `*_1.sh` (1–10000), `*_2.sh` (10001–16000) — the PBS array cap, not the ρ boundary | `*_1.sh` (1–1600) |
| results | `../results/correlated/confidence_intervals/<outcome>/rho_<ρ>/scenario_<k>/<n>/<CI_sf>/` | `.../<outcome>/sf_calibration/rho_<ρ>/scenario_<k>/<n>/` |

Everything else is the parent's: `CI_boot = 200`, `alpha = 0.05` and
`profile = "ci"` for the CI studies. The optimal_sf calibration uses 50
resamples × `CI_boot = 100` per candidate, a final band at `CI_boot = 200`,
and `alpha = 0.1`.

## DGM and pairing

The scenario sets are `correlated/`'s (`R/dgm_scenarios.R`, "CORRELATED
SAMPLE-SIZE SETS"). τ(x) and m0(x) are the parent tables' rows 1–4, coefficient
for coefficient, so the binary set is the same table `binary_ci` resolves to.
Only the covariates' joint distribution changes.

Each run is seeded by its run index alone. So run *r* is the same dataset at
every `CI_sf`, and at ρ = 0 and 0.5 it is paired exactly as in `correlated/`:
the same W, raw latent normals, continuous noise and categorical pair.
Compare ρ with paired per-run differences (MCSE = sd(diff)/√runs). The ρ = 0
arm is **not** paired with the parent CI studies, because the copula draws in
a different order, so compare it to them with the unpaired MCSE.

## The query grid

The CI studies also score each band over the fixed covariate grid from
`build_query_grid()`, which is unchanged here. It crosses X4 and X5 over
[−2, 2] independently, with every other covariate at 0. τ at those points is
still the truth at ρ = 0.5. But (X4 = 2, X5 = −2) is about a 4-SD point in the
correlated latent rather than a 2-SD one, so the forests are extrapolating
there. A ρ = 0.5 fall in grid coverage can be extrapolation rather than
miscalibration. **The per-unit band is the primary ρ comparison.** The
bias-eliminated grid coverage helps separate the two: it scores the band
against the across-run mean estimate, so bias drops out.

## optimal_sf

`find_optimal_sf()` calibrates against the run's own estimate. It cannot see
that estimate's bias, and at ρ = 0.5 the prognostic function tracks τ, which
is exactly where bias would grow. The question this arm answers is whether
the calibration still delivers its nominal 0.90 against the *true* τ. Its
metrics record each run's pick and the calibration's own coverage at the pick
(`plugin_coverage`) alongside the true-τ coverage, so the gap between the two
is visible.

Neither parent optimal_sf study has a metrics script.
`*_corr_ci_sf_metrics.R` files each run's single arm under `dr_random_forest`
so that `compute_metrics()` and `grid_be_reference()` (`R/metrics.R`) run
unchanged. It could serve the parents as written, with only the config path
changed.

## Files

| | |
|---|---|
| `*_config.R` | the grid (`rho` is a path column, prefix `rho_`) |
| `*_analysis.R` | the parent's analysis script. DGM from `../<outcome>/*_corr_dgms.R`, models from the parent's `*_ci_models.R`, `set = corr_set(outcome, rho)` for the query grid, writes to `combo_dir()` |
| `*_check.R`, `*_collect.R`, `*_metrics.R` | as the parent's, outputs carry a `rho` column. optimal_sf's two studies share `optimal_sf/jobscripts/`, so their todo lists are `failed_cts_ids.txt` / `failed_bin_ids.txt` |
| `continuous/cts_corr_ci_results.qmd`, `binary/bin_corr_ci_results.qmd` | per outcome: coverage and length by ρ across the sweep, paired ρ differences, interval kinds (with the grid caveat), the best `CI_sf` by ρ, the optimal_sf pick and its true-τ coverage, a ρ = 0 sanity check against the parent CI study, and a headline table. The optimal_sf and sanity sections are skipped if their inputs are not on disk |

## Running it

On the cluster, run `Rscript make_log_dirs.R` from the repo root once; it
creates every new `logs_*/` directory. Then submit from each jobscripts folder,
because the scripts `cd "${PBS_O_WORKDIR}/.."`:

```bash
cd sample_size/correlated/confidence_intervals/continuous/jobscripts
qsub cts_corr_ci_1.sh              # 1-10000
qsub cts_corr_ci_2.sh              # 10001-16000
Rscript ../cts_corr_ci_check.R     # writes failed_ids.txt, points cts_corr_ci_rerun.sh at it
qsub cts_corr_ci_rerun.sh          # only if the check found failures
qsub cts_corr_ci_collect.sh
qsub cts_corr_ci_metrics.sh

cd ../../optimal_sf/jobscripts
qsub cts_corr_ci_sf_1.sh           # 1-1600
Rscript ../cts_corr_ci_sf_check.R  # writes failed_cts_ids.txt
qsub cts_corr_ci_sf_rerun.sh
qsub cts_corr_ci_sf_collect.sh
qsub cts_corr_ci_sf_metrics.sh
```

Binary is the same with `bin_`. Render a results qmd after both its CI
metrics and its optimal_sf metrics exist.

Resources are hand-set. The CI jobs take the parent's (30 min, 2 cores,
3 GB). The optimal_sf jobs take 8 h, 2 cores and 7 GB. About half the parent
optimal_sf runs outgrew 6 h / 5 GB, and the parent binary rerun went to
10 GB, so watch the binary memory.

Local smoke test, from the study folder: `Rscript cts_corr_ci_analysis.R 3`
(scenario 3, n = 500, `CI_sf` = 0.05, run 1, ρ = 0) and `8003` (the same run
at ρ = 0.5). For optimal_sf, rows `3` and `803` are the same pair, but a full
calibration takes hours. Shrink `n_sim`, `CI_boot` and `sf_grid` in a scratch
copy to check the wiring.

## Status

First run (2026-10-04), with nothing to archive. All four studies are
`first_run` in `R/study_registry.R` and skipped by `R/archive_old_results.R`.
