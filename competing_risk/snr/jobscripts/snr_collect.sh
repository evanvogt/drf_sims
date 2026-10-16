#!/bin/bash
# Reads every per-run file, data and nuisances included (the n = 4000 ones are
# the bulk), so mem as surv_collect.sh.
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=2:ompthreads=2:mem=32gb
#PBS -N snr_collect
#PBS -j oe


module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to the script directory
cd "${PBS_O_WORKDIR}/.."

Rscript ../../study.R collect competing_risk/snr/snr_config.R
