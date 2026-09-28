#!/bin/bash
#PBS -l walltime=06:00:00
#PBS -l select=1:ncpus=2:ompthreads=2:mem=5gb
#PBS -J 1-800%100
#PBS -N ci_sf_1
#PBS -o logs_cts/
#PBS -j oe

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

# Run R script with parameters
Rscript cts_ci_sf_analysis.R "$PBS_ARRAY_INDEX"
