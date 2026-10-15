#!/bin/bash
# Reads bin_snr_all.RDS, as bin_snr_collect.sh writes it. Single-threaded.
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=1:ompthreads=1:mem=16gb
#PBS -N bin_snr_metrics
#PBS -j oe


module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to the script directory
cd "${PBS_O_WORKDIR}/.."

Rscript bin_snr_metrics.R
