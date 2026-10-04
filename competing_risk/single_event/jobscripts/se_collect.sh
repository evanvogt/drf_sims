#!/bin/bash
# se_all.RDS carries every run's data, truth and nuisances; 3,000 runs of 10
# single-target arms is well under competing_risk's 14,000 x 16 x 3, hence
# less than surv_collect.sh's 32gb.
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=2:ompthreads=2:mem=16gb
#PBS -N se_collect
#PBS -j oe


module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to the script directory
cd "${PBS_O_WORKDIR}/.."

Rscript se_collect.R
