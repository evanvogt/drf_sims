#!/bin/bash
#PBS -l walltime=01:00:00  
#PBS -l select=1:ncpus=3:ompthreads=3:mem=5gb
#PBS -J 1-10000%120
#PBS -N ci_cts_1
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
Rscript cts_ci_analysis.R "$PBS_ARRAY_INDEX"