#!/bin/bash
# Reads the oracle add-on's per-run files: one tau vector per estimand, no data
# or nuisances, so far smaller than snr_collect.sh's.
#PBS -l walltime=00:30:00
#PBS -l select=1:ncpus=2:ompthreads=2:mem=8gb
#PBS -N snr_oracle_collect
#PBS -j oe


module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to the script directory
cd "${PBS_O_WORKDIR}/.."

Rscript ../../study.R collect competing_risk/snr/snr_oracle_config.R
