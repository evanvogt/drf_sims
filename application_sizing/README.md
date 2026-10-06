# Application sizing

How long, and how much memory, does applying the CATE methods to two applied
datasets take? Mock datasets with the applied data's shape are run on the
cluster for **timing and memory only**; nothing here estimates performance,
and no estimate is saved.

The applied data itself never comes here. These datasets copy its size,
covariate types, missingness and outcome structure, with arm sizes rounded to
the nearest 50 and domains named D1-D5.

| file | |
|---|---|
| `as_schema.R` | the two schemas: a platform trial analysed as pairwise comparisons within domains (2 binary outcomes, a competing events outcome at 4 horizons) and a two-arm trial (3 continuous outcomes) |
| `as_mock.R` | `mock_platform_domain(d)` (every arm of a domain), `platform_comparison(dom, id)` (one pairwise comparison out of it), `mock_two_arm()`, and `as_cate_data()` / `as_surv_data()` |
| `as_mock_check.R` | checks the mock datasets against the schema; re-derives the missingness calibration |
| `as_config.R` | the timing grid (201 rows) and the analysis settings it times |
| `as_analysis.R` | times one grid row |
| `as_check.R`, `as_collect.R` | the usual todo list / collection; `as_collect.R` writes `as_fits.csv` and `as_jobs.csv` |
| `as_summary.R` | the totals: hours per trial and missing-data method, CI and sf-calibration cost per outcome, peak memory |

## What is timed

The analysis as planned, without the oracle-type arms:

- **Platform**, per pairwise comparison (9): `cate_methods(binomial())` on 6
  binary outcomes (the two binary outcomes and discharge by each horizon) and
  `all_cate_surv_models()` at the 4 horizons, RMTL for discharge only, with
  `csf_sh`, the pseudo-value CF / DR-RF / DR-SL arms and RSF DR in every
  crossfitting variant. Randomisation month is a propensity-only covariate
  (`X_ps`). Missing covariates are handled within domain, by arm, with
  `bin1`, `time` and `status` as the outcome predictors.
- **Two-arm**: `cate_methods(gaussian())` on 3 continuous outcomes.

Each missing-data method's total is built from rows: its imputation (timed
once per domain) plus analysis rows on data of the matching shape - see the
header of `as_config.R`. Multiple imputation is its imputation plus 50 x the
analysis of one completed dataset, rather than 50 analyses run here.

CIs are timed on one outcome per dataset (`ci` rows), and the sf calibration
one sample.fraction per row (`sf` rows; the full calibration is the sum over
the 10), so both come out as a cost per outcome.

## The mock data

- **Covariates** - Gaussian copula, exchangeable 0.3; a third of the
  continuous covariates log-normal; categoricals as dummies with the reference
  dropped. Platform datasets also carry other-domain allocation dummies, and
  `rand_month` separately as the propensity-only covariate (`ps_X`).
- **Allocation** - platform: starts equal and drifts monthly so each arm ends
  with its schema share; two-arm: 1:1 with a binary strata factor.
- **Missingness** - MCAR, correlated within the continuous and within the
  binary block. Platform: per-covariate rates in tiers, block correlation
  calibrated (`miss_rho`) to about 22% complete cases. Two-arm: about 82%.
- **Outcomes** - platform: death (cause 2) at about 33%, most recorded at day
  90; discharge (cause 1) otherwise; no censoring. Two-arm: three continuous
  outcomes driven by a shared, complete baseline value.

## Running

Local wiring check (minutes: 400 patients, 3 folds, a two-learner library,
results to a temporary directory), from this folder:

```
AS_SMOKE=1 Rscript as_analysis.R 33     # bash; PowerShell: $env:AS_SMOKE=1
Rscript as_mock_check.R
```

Cluster:

```
cd application_sizing/jobscripts
qsub as_1.sh                 # rows 1-201, 8 cores (4 workers x 2 grf threads), 32 GB, 24 h
Rscript ../as_check.R        # writes failed_ids.txt, updates as_rerun.sh
qsub as_rerun.sh
qsub as_collect.sh           # as_collect.R then as_summary.R
```

Resources are generous on purpose: memory is one of the things being
measured, and the laptop that will run the real analysis has 15 GB.

## Carrying the times over to the laptop

The real analysis runs locally (12 cores, 15 GB), with the same
`workers` / `grf_threads` pair as here (4 / 2) so the parallel layout matches.
The laptop's cores are not the cluster's, so time one of the smaller rows on
both - e.g. row 38 (`D3_3v1`, `binary`, `completed`, about 450 patients) - and
pass laptop / cluster seconds to `Rscript as_summary.R <ratio>`.

Memory: the R-side peak (`peak_r_mb`) sums the main process's and each
worker's high-water mark, an upper bound. `as_collect.R` also reads the PBS
peak from the logs in `jobscripts/logs/` if its job-summary footer has a
"(peak)" figure; check one log first, the footer format is assumed.

The half-sample bootstrap and `find_optimal_sf()` have no thread argument:
their forests take grf's default (every core the machine shows) inside each
future worker, so the `ci` and `sf` rows oversubscribe the job's 8 cores, as
they would the laptop's.
