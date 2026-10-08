#!/bin/bash
# The nuisance-arm pass: adds cv_shared and holdout to every completed run of
# the main study, writing to results/model_evaluation_strategies. See
# me_strategies.R's header for why this is a second pass rather than a rerun.
#
# PLACEHOLDER resources, for the same reason me_1.sh's are - check the first
# subjobs' resources_used before trusting them. The cost is not me_1.sh's: no
# data generation and no candidate fitting, but
# cv_shared re-runs the full 10-fold nuisance pipeline (about the cost of the
# old `cv` arm on its own) and holdout adds 10 more fits per target on small
# blocks - where H2O's per-call JVM overhead, not model size, is what
# dominates. Budget on the order of the original job's nuisance half.
#
# -J is 1-1200 and NOT the count of completed runs: the array index has to keep
# meaning the same row of the same 1200-row grid as me_1.sh. Runs that have no
# source file exit 0 with a message.
#
# 10 cores all go to n_cores (XGBoost nthread / H2O nthreads) - there are no
# candidate workers here, and the fold loops run one after another, so expect
# less speed-up than me_1.sh gets. mem covers the H2O JVM (max heap 10G) and
# one R session.
#
# %N throttle: same constraint as me_1.sh - each concurrent task starts its
# own H2O JVM cluster with a 10G heap.
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=10:ompthreads=10:mem=16gb
#PBS -J 1-1200%10
#PBS -N me_strat
#PBS -o logs_strat/
#PBS -e logs_strat/


module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

# Trailing arg is n_cores (XGBoost nthread / H2O nthreads). There is no
# `workers` knob here - nothing in this pass fits candidate models.
Rscript me_strategies.R "$PBS_ARRAY_INDEX" 10
