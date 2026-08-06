#!/bin/bash
# Runs cf_diagnose_sampler.R and cf_diagnose_multisession.R in exactly the
# environment where cf_profile.sh has failed/warned - same modules, same conda
# env, same working directory, same node type. Single job, under a minute each.
# Read the "=== reading ===" block at the end of each script's output in
# logs_diagnose/ for the verdict.
#
# ncpus=8:ompthreads=8 matches cf_profile.sh so the contention is comparable.
#PBS -l walltime=00:20:00
#PBS -l select=1:ncpus=8:ompthreads=8:mem=32gb
#PBS -N cf_diagnose
#PBS -o logs_diagnose/
#PBS -e logs_diagnose/


module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

Rscript cf_diagnose_sampler.R
Rscript cf_diagnose_multisession.R
