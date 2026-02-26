#!/bin/bash
#PBS -l walltime=24:00:00  
#PBS -l select=1:ncpus=10:ompthreads=10:mem=20gb
#PBS -J 1-500%43
#PBS -N cts_miss_ci
#PBS -o logs/
#PBS -e logs/

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

# Run R script with parameters
Rscript cts_miss_ci_analysis.R "$PBS_ARRAY_INDEX"
