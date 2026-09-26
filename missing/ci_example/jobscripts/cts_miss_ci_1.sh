#!/bin/bash
#PBS -l walltime=08:00:00
#PBS -l select=1:ncpus=1:ompthreads=1:mem=10gb
#PBS -J 1-500%380
#PBS -N cts_miss_ci
#PBS -o logs_1/
#PBS -e logs_1/

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

# Run R script with parameters. Trailing args are workers[1]/workers[2]/
# grf_threads, set by hand to match ncpus above - change them together.
Rscript cts_miss_ci_analysis.R "$PBS_ARRAY_INDEX" 1 1 1
