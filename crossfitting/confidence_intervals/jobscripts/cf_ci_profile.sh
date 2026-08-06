#!/bin/bash
# Deliberately over-provisioned so nothing fails for the wrong reason. ncpus=8
# covers the largest configuration in the sweep (workers=4 x grf_threads=2).
# ompthreads=8 matches ncpus so parallelly::availableCores() sees all 8 and
# plan(multisession, workers = 4) doesn't fail; cf_ci_profile.R still sets
# OMP_NUM_THREADS from the sweep grid before spawning its multisession workers,
# so per-worker thread counts are controlled there, not by ompthreads.
#
# Memory and CPU are measured inside R by syrup, so there is no background
# qstat sampler here and nothing depends on the scheduler.
#PBS -l walltime=08:00:00
#PBS -l select=1:ncpus=8:ompthreads=8:mem=16gb
#PBS -J 1-48
#PBS -N cf_ci_profile
#PBS -o logs_ci_profile/
#PBS -e logs_ci_profile/


module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

Rscript cf_ci_profile.R "$PBS_ARRAY_INDEX"
