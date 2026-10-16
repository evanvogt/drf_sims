#!/bin/bash
# Reads every per-run file, data and nuisances included.
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=2:ompthreads=2:mem=16gb
#PBS -N bin_snr_collect
#PBS -j oe


module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to the script directory
cd "${PBS_O_WORKDIR}/.."

Rscript ../../study.R collect sample_size/snr/bin_snr_config.R
