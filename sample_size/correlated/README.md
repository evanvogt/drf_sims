# Correlated covariates — sample size studies

`continuous/` and `binary/` draw every covariate independently. These two
studies rerun their scenarios 1–4 with correlated covariates, so the
estimators have to separate prognosis from effect modification when the two
move together. The estimators, metrics and HTE tests are exactly the parent
studies'.

| | continuous | binary |
|---|---|---|
| folder | `continuous/` (`cts_corr_*`) | `binary/` (`bin_corr_*`) |
| scenario sets | `continuous_corr_0`, `continuous_corr_0.5` | `binary_corr_0`, `binary_corr_0.5` |
| scenarios | 1–4 (null, simple, complex, non-linear) | same |
| ρ | 0, 0.5 (`CORR_RHOS`) | same |
| n | 100, 250, 500, 1000 | same |
| runs | 500 per (ρ, scenario, n) | same |
| array | 16,000: rows 1–8000 are ρ = 0, rows 8001–16000 are ρ = 0.5; the jobscripts split at the PBS cap instead, `*_corr_1.sh` 1–10000 and `*_corr_2.sh` 10001–16000 | same |
| results | `../results/correlated/continuous/rho_<ρ>/scenario_<k>/<n>/res_sim_<run>.RDS` | `.../correlated/binary/...` |

`confidence_intervals/` holds the CI studies' counterparts on the same DGM:
the continuous and binary `CI_sf` sweeps and both optimal_sf calibrations, at
n ∈ {500, 1000}. See `confidence_intervals/README.md`.

## DGM

The scenario tables are in `R/dgm_scenarios.R` ("CORRELATED SAMPLE-SIZE
SETS"). Each is the parent table's rows 1–4 plus a `rho` column, which sends
the generator down the copula branch already used by `missing/`
(`correlated_covariates()`):

- A latent `Z ~ N(0, R)` is drawn over X1–X5 and X01–X03, with R
  exchangeable at ρ.
- X1 and X3 are thresholded at their prevalences (0.4, 0.7). X2, X4 and X5
  are the latent columns, with SD 1. So every marginal is the parent study's.
- X01–X03 are correlated with X1–X5 too: they are noise for the outcome and
  the CATE, but partial proxies for the modifiers. X04 and X05 (the
  categorical pair) are independent.
- At ρ = 0.5 the observed correlations are about 0.5 between continuous
  covariates and 0.3–0.4 for pairs involving X1 or X3.

**τ(x) and m0(x) are the parent studies' functions, coefficient for
coefficient.** Only the covariates' joint distribution changes, and with it:

- bW and the ATE. The calibration integrates over the correlated latent:
  E[g] moves in scenario 3 (X4·X5), and the planned SD / control event rate
  moves through Cov(X1, X2).
- SD(τ) in scenario 3, and how strongly the prognostic m0 tracks τ.

From `Rscript R/calibration_report.R`, plus a large-sample simulation for
SD(τ) and cor(m0, τ):

| | ρ = 0 | ρ = 0.5 |
|---|---|---|
| **continuous** | | |
| true ATE, n = 100 / 1000 | −0.65 / −0.20 | −0.60 / −0.19 |
| realised power, scenarios 2 / 3 / 4 | 0.65–0.67 / 0.61–0.63 / 0.76–0.77 | 0.65–0.66 / 0.55–0.56 / 0.75–0.76 |
| SD τ, scenario 3 | 1.16 | 1.36 |
| cor(m0, τ), scenarios 2 / 3 / 4 | 0 / 0 / 0 | −0.44 / 0.39 / 0.01 |
| **binary** | | |
| true RD, n = 100 / 1000 | −0.248 / −0.085 | −0.247 / −0.085 |
| power | 0.80 | 0.80 |
| cor(m0, τ), scenarios 2 / 3 / 4 | 0 / 0 / 0 | 0.39 / −0.34 / 0.06 |
| treated-risk floor, scenario 3, n = 100 | 0.011 | **0.007** |

At ρ = 0 bW equals the parent studies' at every n.

**Binary scenario 3 floor.** `binary_corr_*` keeps `RD_SCALE`, so that τ(x)
is `binary/`'s. At ρ = 0.5, E[g] rises in scenario 3, bW falls, and the
lowest treated risk over the covariate support at n = 100 is 0.007. That is
inside [0, 1] but below `RD_EPS` = 0.01. It is accepted on purpose: a scale
of 0.049 would restore the floor at the cost of a slightly smaller HTE than
`binary/`. `bin_verify_hte.R` check 7 holds this one cell to > 0 and every
other cell to `RD_EPS`.

## Pairing, and how to compare ρ

Each run is seeded by its run index alone (`setup_rng_stream(run)`), and the
copula draws W, then one n × 8 block of raw normals, then the noise, then the
categorical pair. So run *r* at ρ = 0 and at ρ = 0.5 shares:

- W;
- the raw normals, multiplied by a different Cholesky factor;
- the continuous outcome's noise;
- X04 and X05.

Binary Y is coupled only through `rbinom`'s shared uniforms; about 91–99% of
outcomes agree across ρ.

- **ρ = 0.5 vs ρ = 0:** use paired per-run differences, with
  MCSE = sd(diff)/√runs. `corr_results.qmd` does this.
- **Within one ρ:** as in the parent studies, compare estimators against
  `dr_oracle` on the same run, and compute MCSEs per scenario, never pooled.
  Scenarios share seeds, so their errors are correlated (see
  `binary/README.md`, "Monte Carlo checks on the ATE bias").
- **ρ = 0 vs `continuous/` / `binary/`:** same distribution, different draws.
  The copula draws in a different order, so these runs are **not** paired with
  the parent studies. Use the unpaired MCSE. `corr_results.qmd`'s last section
  uses this as a sanity check.

## Files

| | |
|---|---|
| `*_corr_config.R` | the grid (`rho` is a path column, prefix `rho_`) |
| `*_corr_dgms.R` | the generator and oracle wrappers, with a `rho` argument |
| `*_corr_analysis.R` | the parent's analysis script; reuses `../continuous/cts_models.R` / `../binary/bin_models.R`, and writes to `combo_dir()` |
| `*_corr_check.R`, `*_corr_collect.R`, `*_corr_metrics.R` | as the parent's; outputs `*_corr_all.RDS`, `*_corr_metrics.RDS` and `*_corr_true_cate_tests.RDS`, all with a `rho` column |
| `corr_results.qmd` | both outcomes: levels by ρ, paired differences, tests, DR-oracle-relative ATE bias, and the ρ = 0 sanity check |
| `continuous/cts_corr_results.qmd`, `binary/bin_corr_results.qmd` | one outcome each, in the shape of the parent's `*_results.qmd`: every metric by ρ, the paired ρ = 0.5 − ρ = 0 differences, true-CATE tests, NA tables, the ρ = 0 sanity check and a headline table |

## Running it

On the cluster, run `Rscript make_log_dirs.R` from the repo root once; it
creates `logs_1/`, `logs_2/` and `logs_rerun/`. Then submit from the
jobscripts folder, because the scripts `cd "${PBS_O_WORKDIR}/.."`:

```bash
cd sample_size/correlated/continuous/jobscripts
qsub cts_corr_1.sh            # 1-10000 (rho = 0 is 1-8000)
qsub cts_corr_2.sh            # 10001-16000, all rho = 0.5
Rscript ../cts_corr_check.R   # writes failed_ids.txt, points cts_corr_rerun.sh at it
qsub cts_corr_rerun.sh        # only if the check found failures
qsub cts_corr_collect.sh
qsub cts_corr_metrics.sh
```

Binary is the same with `bin_corr`. Local smoke test, from the study folder:
`Rscript cts_corr_analysis.R 3` (scenario 3, n = 100, run 1, ρ = 0) and
`Rscript cts_corr_analysis.R 8003` (the same run at ρ = 0.5).

## Status

First run (2026-10-01). Nothing to archive:
`R/archive_old_results.R` skips both studies.
