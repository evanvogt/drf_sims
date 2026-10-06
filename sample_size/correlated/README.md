# Correlated covariates — sample size studies

`continuous/` and `binary/` draw every covariate independently. These two
studies rerun their scenarios with correlated covariates, so the estimators
have to separate prognosis from effect modification when the two move
together. The estimators, metrics and HTE tests are exactly the parent
studies'. Scenarios 1–4 are the main text's; 5–10 (added 2026-10-04) are the
appendix's and, like the parent studies' 5–10, run to 100.

| | continuous | binary |
|---|---|---|
| folder | `continuous/` (`cts_corr_*`) | `binary/` (`bin_corr_*`) |
| scenario sets | `continuous_corr_0`, `continuous_corr_0.5` | `binary_corr_0`, `binary_corr_0.5` |
| scenarios | 1–4 (null, simple, complex, non-linear), 5–10 | same |
| ρ | 0, 0.5 (`CORR_RHOS`) | same |
| n | 100, 250, 500, 1000 | same |
| runs | 500 per (ρ, scenario, n) for scenarios 1–4, 100 for 5–10 | same |
| array | 20,800. Scenarios 1–4 are rows 1–16000: 1–8000 are ρ = 0, 8001–16000 ρ = 0.5; the jobscripts split at the PBS cap instead, `*_corr_1.sh` 1–10000 and `*_corr_2.sh` 10001–16000. Scenarios 5–10 are rows 16001–20800: 16001–18400 are ρ = 0, 18401–20800 ρ = 0.5, all in `*_corr_extra.sh` | same |
| results | `../results/correlated/continuous/rho_<ρ>/scenario_<k>/<n>/res_sim_<run>.RDS` | `.../correlated/binary/...` |

`confidence_intervals/` holds the CI studies' counterparts on the same DGM:
the continuous and binary `CI_sf` sweeps and both optimal_sf calibrations, at
n ∈ {500, 1000}. See `confidence_intervals/README.md`.

## DGM

The scenario tables are in `R/dgm_scenarios.R` ("CORRELATED SAMPLE-SIZE
SETS"). Each is the parent table (all ten rows) plus a `rho` column, which sends
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
coefficient**, with one exception: binary scenario 9 (see below). Only the
covariates' joint distribution changes, and with it:

- bW and the ATE. The calibration integrates over the correlated latent:
  E[g] moves wherever g multiplies two modifiers (scenarios 3 and 9, X4·X5;
  7 and 8, X3·X4), and the planned SD / control event rate moves through
  Cov(X1, X2).
- SD(τ) in those scenarios, and how strongly the prognostic m0 tracks τ.

From `Rscript R/calibration_report.R` and, for SD(τ) and cor(m0, τ), all ten
scenarios, `Rscript corr_truth_summary.R`:

| | ρ = 0 | ρ = 0.5 |
|---|---|---|
| **continuous** | | |
| true ATE, n = 100 / 1000 | −0.65 / −0.20 | −0.60 / −0.19 |
| realised power, scenarios 2 / 3 / 4 | 0.65–0.67 / 0.61–0.63 / 0.76–0.77 | 0.65–0.66 / 0.55–0.56 / 0.75–0.76 |
| SD τ, scenario 3 | 1.16 | 1.36 |
| cor(m0, τ), scenarios 2 / 3 / 4 | 0 / 0 / 0 | −0.43 / 0.39 / 0.01 |
| **binary** | | |
| true RD, n = 100 / 1000 | −0.248 / −0.085 | −0.247 / −0.085 |
| power | 0.80 | 0.80 |
| cor(m0, τ), scenarios 2 / 3 / 4 | 0 / 0 / 0 | 0.39 / −0.34 / 0.06 |
| treated-risk floor, scenario 3, n = 100 | 0.011 | **0.007** |

The table covers scenarios 1–4. For all ten, `calibration_report.R` prints bW,
the true ATE, power and (binary) the treated-risk bounds. Continuous realised
power at ρ = 0.5 runs 0.62–0.80 in scenarios 5–10; scenario 8's drops most
(0.68 → 0.62–0.63).

At ρ = 0 bW equals the parent studies' at every n, in every scenario.

**Binary floors at ρ = 0.5.** `binary_corr_*` takes `RD_SCALE_CORR`, which is
`RD_SCALE` everywhere except scenario 9, so τ(x) is `binary/`'s wherever the
bounds allow it. At ρ = 0.5, E[g] rises in scenarios 3, 8 and 9, bW falls,
and the lowest treated risk over the covariate support at n = 100 falls:

| scenario | floor under `RD_SCALE` | handling |
|---|---|---|
| 3 | 0.007 | accepted: inside [0, 1], below `RD_EPS` = 0.01. 0.049 would restore it |
| 8 | 0.003 | accepted, the same way. 0.127 would restore it |
| 9 | −0.004 | the generator would stop. `RD_SCALE_CORR[9]` = 0.139 (from 0.164), the largest scale, floored to 3 dp, with the floor at `RD_EPS`, at both ρ so the arms stay paired. Floor 0.011 at ρ = 0.5, 0.023 at ρ = 0 |

So binary scenario 9's HTE is about 15% smaller than `binary/`'s, and its
ρ = 0 arm is not `binary/` scenario 9's distribution (bW still is: E[g] = 0
there). `R/bin_verify_hte.R` check 7 re-derives 0.139, holds scenarios 3 and 8
at ρ = 0.5 to > 0, and holds every other cell to `RD_EPS`.

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
  uses this as a sanity check, for all ten scenarios. Binary scenario 9 is
  expected to differ, because of its own scale.

## Files

| | |
|---|---|
| `corr_truth_summary.R` | SD(τ) (exact) and cor(m0, τ) (one seeded draw of 10^6) for every outcome × ρ × scenario; needs no simulation output |
| `*_corr_config.R` | the grid (`rho` is a path column, prefix `rho_`) |
| `*_corr_dgms.R` | the generator and oracle wrappers, with a `rho` argument |
| `*_corr_models.R` | the parent's models wrapper, moved here unchanged when the parent was retired (2026-10-06) |
| `*_corr_analysis.R` | the parent's analysis script; sources `*_corr_models.R` and writes to `combo_dir()` |
| `*_corr_check.R`, `*_corr_collect.R`, `*_corr_metrics.R` | as the parent's; outputs `*_corr_all.RDS`, `*_corr_metrics.RDS` and `*_corr_true_cate_tests.RDS`, all with a `rho` column |
| `corr_results.qmd` | both outcomes: levels by ρ, paired differences, tests, DR-oracle-relative ATE bias, and the ρ = 0 sanity check |
| `continuous/cts_corr_results.qmd`, `binary/bin_corr_results.qmd` | one outcome each, in the shape of the parent's `*_results.qmd`: every metric by ρ, the paired ρ = 0.5 − ρ = 0 differences, true-CATE tests, NA tables, the ρ = 0 sanity check and a headline table |

## Running it

On the cluster, run `Rscript make_log_dirs.R` from the repo root once; it
creates `logs_1/`, `logs_2/`, `logs_extra/` and `logs_rerun/`. Then submit
from the jobscripts folder, because the scripts `cd "${PBS_O_WORKDIR}/.."`:

```bash
cd sample_size/correlated/continuous/jobscripts
qsub cts_corr_1.sh            # 1-10000 (rho = 0 is 1-8000)
qsub cts_corr_2.sh            # 10001-16000, all rho = 0.5
qsub cts_corr_extra.sh        # 16001-20800, scenarios 5-10 (rho = 0 is 16001-18400)
Rscript ../cts_corr_check.R   # writes failed_ids.txt, points cts_corr_rerun.sh at it
qsub cts_corr_rerun.sh        # only if the check found failures
qsub cts_corr_collect.sh
qsub cts_corr_metrics.sh
```

The check covers the whole grid, so run it after `*_corr_extra.sh` finishes.
Otherwise it lists every scenario 5–10 row as failed.

Binary is the same with `bin_corr`. Local smoke test, from the study folder:
`Rscript cts_corr_analysis.R 3` (scenario 3, n = 100, run 1, ρ = 0),
`Rscript cts_corr_analysis.R 8003` (the same run at ρ = 0.5),
`Rscript cts_corr_analysis.R 16003` (scenario 7, n = 100, run 1, ρ = 0) and
`Rscript bin_corr_analysis.R 18405` (binary scenario 9, n = 100, run 1,
ρ = 0.5 — the cell with its own scale).

## Status

Scenarios 1–4: first run (2026-10-01). Scenarios 5–10: added 2026-10-04, not
yet run. Nothing to archive: `R/archive_old_results.R` skips both studies.
