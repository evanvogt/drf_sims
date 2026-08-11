#!/bin/bash
# Deliberately over-provisioned so nothing fails for the wrong reason. ncpus=8
# covers the largest inner configuration in the sweep (workers=4 x
# grf_threads=2). ompthreads=8 matches ncpus so parallelly::availableCores()
# sees all 8 and plan(multisession, workers = 4) doesn't fail; cts_miss_profile.R
# still sets OMP_NUM_THREADS from the sweep grid before spawning its
# multisession workers, so per-worker thread counts are controlled there, not
# by ompthreads.
#
# #PBS -J 1-72, not the 9900-row production grid: cts_miss_profile.R uses a
# two-level outer(method x scenario x mechanism x run)/inner(workers x
# grf_threads) grid, not a flat one - one array index covers the OUTER 72-row
# grid and loops all 6 inner cells itself, reusing the one data generation +
# missingness handling pass across them (see cts_miss_profile.R's header).
# Walltime is longer than continuous/cts_profile.sh's for the same reason: each
# task does 6 sequential model-fitting cells plus, for missforest/multiple_imputation,
# a slow serial imputation step (missForest's iterative random forests /
# mice(m = 50, method = "rf")) before any of those cells even start.
#
# Memory and CPU are measured inside R by syrup, so there is no background
# qstat sampler here and nothing depends on the scheduler.
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=8:ompthreads=8:mem=32gb
#PBS -J 1-72
#PBS -N cts_miss_profile
#PBS -o logs_profile/
#PBS -e logs_profile/


module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

Rscript cts_miss_profile.R "$PBS_ARRAY_INDEX"
