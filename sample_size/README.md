# Sample-size studies

How well do doubly-robust and forest-based CATE estimators recover a
heterogeneous treatment effect as the trial gets smaller, for continuous and
binary outcomes? And do their confidence bands still cover?

Every study here uses one data-generating mechanism: ten treatment-effect
scenarios, with the covariates drawn from a Gaussian copula at latent
correlation ρ ∈ {0, 0.5}. ρ = 0.5 is the primary setting, where prognosis and
effect modification move together. ρ = 0 is the independent-covariate
comparison, paired run for run with ρ = 0.5. The full design is in
[`ADEMP.md`](ADEMP.md).

| study | folder | prefix | question | grid |
|---|---|---|---|---|
| continuous | `continuous/` | `cts_corr` | CATE estimation and HTE tests | scenarios 1–10 × n 100, 250, 500, 1000 × ρ; 500 runs for scenarios 1–4, 100 for 5–10 (20,800 rows) |
| binary | `binary/` | `bin_corr` | the same, binary outcome (risk-difference scale) | the same |
| CI, continuous | `confidence_intervals/continuous/` | `cts_corr_ci` | coverage of half-sample bootstrap bands across `CI_sf` | scenarios 1–4 × n 500, 1000 × ρ × 10 `CI_sf`; 100 runs (16,000 rows) |
| CI, binary | `confidence_intervals/binary/` | `bin_corr_ci` | the same, binary outcome | the same |
| optimal_sf | `confidence_intervals/optimal_sf/` | `cts_corr_ci_sf` | a data-driven choice of `CI_sf`, continuous outcome only | scenarios 1–4 × n 500, 1000 × ρ; 100 runs (1,600 rows) |
| SNR check | `snr/` | `bin_snr` | does the DR-learner catch the causal forest with more signal? | see its README |

Scenarios 1–4 (null, simple, complex, non-linear) are the main text's; 5–10 are
the appendix's. Results are written beside the repository, under
`../results/correlated/<outcome>/`, `../results/correlated/confidence_intervals/<outcome>/`
and `.../continuous/sf_calibration/` for optimal_sf.

## Checking the DGM

None of these needs simulation output:

```bash
Rscript R/calibration_report.R        # bW, true ATE, power and (binary) treated-risk bounds, every scenario and n
Rscript R/bin_verify_hte.R            # re-derives RD_SCALE / RD_SCALE_CORR and checks the binary floors
Rscript sample_size/corr_truth_summary.R   # SD(τ) and cor(m0, τ) for every outcome, ρ and scenario
```

After any change to the scenario tables, run the calibration report on the
cluster too: R 4.3.2 there could round a borderline bW differently.

## Files

Each study folder holds these files (the CI and optimal_sf studies use the
continuous and binary studies' `*_dgms.R`, and optimal_sf the CI studies'
`*_models.R`):

| file | |
|---|---|
| `<prefix>_config.R` | the parameter grid and results path, declared once. The array index is a row of `study$grid`: never filter or reorder it |
| `<prefix>_dgms.R` | the generator and oracle wrappers, with a `rho` argument |
| `<prefix>_models.R` | this study's call to `cate_methods()` |
| `<prefix>_analysis.R` | the array entry point: one grid row per index |
| `<prefix>_check.R` | lists missing runs in `jobscripts/failed_ids.txt` and points `<prefix>_rerun.sh` at it |
| `<prefix>_collect.R`, `<prefix>_metrics.R` | write `<prefix>_all.RDS` and `<prefix>_metrics.RDS` (and, for continuous and binary, `<prefix>_true_cate_tests.RDS`), all with a `rho` column |

In `sample_size/` itself:

| file | |
|---|---|
| `corr_results.qmd` | both outcomes: levels by ρ, paired differences, tests, ATE bias against `dr_oracle` |
| `corr_outcome_results.qmd` | one outcome (Quarto `outcome` param): every metric by ρ, paired differences, the true-CATE tests, NA tables, a headline table |
| `confidence_intervals/corr_ci_results.qmd` | the CI studies, one outcome per render (`outcome` param); the optimal_sf section, `corr_ci_results_sf.qmd`, is included for the continuous outcome only |
| `corr_truth_summary.R` | SD(τ) and cor(m0, τ) for every outcome, ρ and scenario; needs no simulation output |
| `corr_sl_failures.R` | summarises the SuperLearner arm's dropped learners and failed fits |
| `corr_blp_beta2.R` | pulls the saved BLP β2 estimates back out of the collected runs |

## Running a study

On the cluster, run `Rscript make_log_dirs.R` from the repo root once. Then
submit from the study's `jobscripts/` folder, because the scripts
`cd "${PBS_O_WORKDIR}/.."`:

```bash
cd sample_size/continuous/jobscripts
qsub cts_corr_1.sh            # rows 1-10000 (rho = 0 is 1-8000)
qsub cts_corr_2.sh            # rows 10001-16000, all rho = 0.5
qsub cts_corr_extra.sh        # rows 16001-20800, scenarios 5-10
Rscript ../../../study.R check sample_size/continuous/cts_corr_config.R   # run after all three finish
qsub cts_corr_rerun.sh        # only if the check found failures
qsub cts_corr_collect.sh
qsub cts_corr_metrics.sh
```

Binary is the same with `bin_corr`. The CI studies have `*_1.sh` (1–10000) and
`*_2.sh` (10001–16000); optimal_sf has one, `cts_corr_ci_sf_1.sh` (1–1600), and
its todo list is `failed_cts_ids.txt`. The jobscripts split at the PBS array
cap, not at the ρ boundary.

Local smoke tests, from the study folder:

- `Rscript cts_corr_analysis.R 3` (scenario 3, n = 100, run 1, ρ = 0) and
  `8003` (the same run at ρ = 0.5);
- `Rscript bin_corr_analysis.R 18405` (binary scenario 9, n = 100, run 1,
  ρ = 0.5: the cell with its own scale);
- `Rscript cts_corr_ci_analysis.R 3` (scenario 3, n = 500, `CI_sf` = 0.05,
  run 1, ρ = 0).

A full optimal_sf calibration takes hours; shrink `n_sim`, `CI_boot` and
`sf_grid` in a scratch copy to check the wiring.

**Resources** are set by hand. The continuous and binary jobs take one core,
2 GB and 30 minutes; the CI jobs two cores, 3 GB and 30 minutes; optimal_sf two
cores, 7 GB and 8 hours. Size for the n = 1000 cell, since one `#PBS -l` line
covers every n. Keep `workers × grf_threads` (the trailing arguments on the
`Rscript` line) within `ncpus`, and change the two together.
