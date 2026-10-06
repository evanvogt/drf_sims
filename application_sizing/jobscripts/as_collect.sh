#!/bin/bash
#PBS -l walltime=00:30:00
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -N as_collect
#PBS -j oe

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

# Run R scripts
Rscript as_collect.R
Rscript as_summary.R
