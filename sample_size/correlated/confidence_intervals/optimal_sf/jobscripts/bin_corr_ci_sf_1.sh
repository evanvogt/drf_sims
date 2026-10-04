#!/bin/bash
#PBS -l walltime=08:00:00
#PBS -l select=1:ncpus=2:ompthreads=2:mem=7gb
#PBS -J 1-1600%100
#PBS -N ci_sf_corr_bin_1
#PBS -o logs_bin/
#PBS -j oe

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

# Run R script with parameters - rows 1-800 are rho = 0, 801-1600 rho = 0.5
Rscript bin_corr_ci_sf_analysis.R "$PBS_ARRAY_INDEX"
