#!/bin/bash
# Reruns only the array indices listed in failed_ids.txt by cf_check.R.
# Set -J to 1-<number of lines in failed_ids.txt> before submitting.
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=2:ompthreads=2:mem=10gb
#PBS -J 1-23
#PBS -N cf_rerun
#PBS -o logs_rerun/
#PBS -e logs_rerun/


module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

jobid=$(sed -n "${PBS_ARRAY_INDEX}p" "${PBS_O_WORKDIR}/failed_ids.txt")

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

Rscript cf_analysis.R "$jobid"
