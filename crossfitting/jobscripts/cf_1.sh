#!/bin/bash
# PLACEHOLDER resources - run cf_profile.sh then cf_profile_summary.R, which
# rewrites the two -l lines below from measured timings and memory.
#PBS -l walltime=04:00:00
#PBS -l select=1:ncpus=2:ompthreads=2:mem=10gb
#PBS -J 1-2000%190
#PBS -N cf_1
#PBS -o logs_1/
#PBS -e logs_1/


module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

# Run R script with parameters
Rscript cf_analysis.R "$PBS_ARRAY_INDEX"
