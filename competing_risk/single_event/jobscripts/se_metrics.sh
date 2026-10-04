#!/bin/bash
# Reads the same se_all.RDS as se_collect.sh, so the same memory. Single core:
# se_metrics.R does no parallel work.
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=1:ompthreads=1:mem=16gb
#PBS -N se_metrics
#PBS -j oe


module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to the script directory
cd "${PBS_O_WORKDIR}/.."

Rscript se_metrics.R
