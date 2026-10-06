#!/bin/bash
# Reruns only the array indices listed in failed_ids.txt by as_check.R.
# -J and the resource request are rewritten by check_failed().
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=8:ompthreads=8:mem=32gb
#PBS -J 1-1%50
#PBS -N as_rerun
#PBS -o logs_rerun/
#PBS -j oe

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

jobid=$(sed -n "${PBS_ARRAY_INDEX}p" "${PBS_O_WORKDIR}/failed_ids.txt")

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

# Run R script with parameters
Rscript as_analysis.R "$jobid" 4 2
