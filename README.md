# drf_sims

Simulation study evaluating CATE (Conditional Average Treatment Effect) estimation methods across outcome types and data settings. Part of WP2 of a PhD at Imperial College London on causal machine learning methods in survival and competing-risks settings.

## Overview

This project benchmarks doubly-robust and forest-based CATE estimators under a range of data-generating processes (DGPs) covering continuous, binary, and time-to-event outcomes. Simulations are designed to run on an HPC cluster and cover scenarios with varying heterogeneous treatment effect (HTE) structures, missing covariates, and competing risks.

## Repository structure

```
drf_sims/
├── utils.R                        # Shared utilities (RNG streams, cross-fitting helpers)
│
├── continuous/                    # Continuous outcome simulations
│   ├── cts_dgms.R                 # DGM: 10 HTE scenarios
│   ├── cts_models.R               # CATE estimation functions
│   ├── cts_analysis.R             # HPC entry point (reads commandArgs)
│   ├── cts_collect.R              # Collects distributed results
│   ├── cts_metrics.R              # Computes performance metrics
│   └── results_cts.R              # Results summaries
│
├── binary/                        # Binary outcome simulations (same structure)
│   ├── bin_dgms.R
│   ├── bin_models.R
│   ├── bin_analysis.R
│   ├── bin_collect.R
│   └── bin_metrics.R
│
├── competing_risk/                # Competing risks simulations
│   ├── surv_dgm.R                 # DGM: 7 scenarios, Weibull joint hazards
│   ├── surv_models.R              # CATE estimation functions
│   ├── surv_analysis.R            # HPC entry point
│   ├── surv_collect.R             # Collects distributed results
│   └── surv_metrics.R             # Computes performance metrics
│
├── missing/
│   ├── continuous/                # Missing covariates, continuous outcome
│   └── binary/                    # Missing covariates, binary outcome
│
├── confidence_intervals/
│   ├── continuous/                # CI estimation for continuous outcome CATEs
│   ├── binary/                    # CI estimation for binary outcome CATEs
│   └── optimal_sf/                # Optimal sample-fraction calibration
│
├── case_study/                    # Factorial platform trial case study
│   ├── cs_surv_models.R           # Multi-arm competing risks CATE models
│   └── cs_bin_models.R            # Multi-arm binary CATE models
│
├── validation/                    # DGM and CATE validation scripts
├── results_processing/            # Thesis figure generation
│   └── thesis_figures/
│
└── renv/                          # renv lockfile for reproducibility
```

## Methods

| Method | Description |
|---|---|
| Causal Forest | `grf::causal_forest` with cross-fitted nuisances |
| DR-RF | Doubly-robust learner, random forest second stage |
| DR-SL | Doubly-robust learner, SuperLearner nuisances and second stage |
| DR-Oracle | DR learner with true nuisance functions |
| DR-Semi-Oracle | DR learner with known propensity (0.5), estimated outcome model |

For competing risks, additional methods are available:

| Method | Description |
|---|---|
| IPW-CF | Causal forest on IPW-transformed RMST outcomes (event 1, event 2, composite) |
| CSF-cs | Causal survival forest using cause-specific hazards |
| CSF-sh | Causal survival forest using subdistribution hazards |
| Pseudo-CF | Causal forest on pseudo-values (RMTL1, RMTL2, RMSTc) |
| Pseudo-DR | DR learner on jackknife pseudo-values |

Nuisance functions use double cross-fitting (10 folds, pairwise fold combinations). SuperLearner libraries are adapted to sample size.

## Simulation scenarios

### Continuous and binary outcomes (10 scenarios)
Scenarios vary the structure of the CATE function:
1. No HTE (ATE only)
2. Simple HTE — binary variable
3. Simple HTE — continuous variable
4. Two HTE variables (additive)
5. Continuous × binary interaction
6. Single effects + interaction
7. Continuous × continuous interaction
8. Single effects + different interaction
9. Cosine HTE
10. Exponential HTE

### Competing risks (7 scenarios)
Event of interest (EOI) and competing event (CE) generated from Weibull distributions via joint hazards (Beyersmann 2009). Scenarios vary ATE/HTE presence across EOI and CE:
1. ATE on EOI only
2. ATE on CE only
3. HTE on EOI, no ATE on CE
4. HTE on EOI, ATE on CE
5. HTE on CE, no ATE on EOI
6. HTE on CE, ATE on EOI
7. HTE on both events

Sample sizes: n ∈ {100, 250, 500, 1000}, 100 simulation runs per parameter combination.

## Performance metrics

- **Point estimation:** bias, ATE bias, RMSE, MAE, correlation, Spearman correlation, sign accuracy
- **HTE detection:** BLP test p-value (`GenericML`), independence test p-value (`coin`)

## Running the simulations

Scripts are designed for array job submission on an HPC cluster. The job index `i` selects a row from the full parameter grid:

```bash
Rscript continuous/cts_analysis.R $SGE_TASK_ID
```

Each script sources its own DGM and model files, sets up reproducible parallel RNG streams, fits all CATE methods, and saves results to `../results/<outcome_type>/`.

Results are collected with the corresponding `*_collect.R` script and metrics computed with `*_metrics.R`.

## Dependencies

Package versions are managed with `renv`. To restore the environment:

```r
renv::restore()
```

Key packages: `grf`, `SuperLearner`, `GenericML`, `pseudo`, `coin`, `furrr`, `future`, `ranger`, `glmnet`, `gam`, `dplyr`, `here`.
