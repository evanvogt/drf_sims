#!/bin/bash

#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=3:ompthreads=2:mem=10gb
#PBS -J 1-6%100
#PBS -N ci_sf_bin_rerun
#PBS -o logs_bin_rerun/
#PBS -e logs_bin_rerun/

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

jobid=$(sed -n "${PBS_ARRAY_INDEX}p" "${PBS_O_WORKDIR}/failed_bin_ids.txt")

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

# Run R script with parameters
Rscript bin_ci_sf_analysis.R "$jobid"
