#!/bin/bash
# Reruns only the array indices listed in failed_ids.txt by study.R check.
# -J and the resource request are rewritten by check_failed(); the values here
# are what it computes from bin_val_1.sh, so the first check leaves them alone.
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=1:ompthreads=1:mem=5gb
#PBS -J 1-4400
#PBS -N bin_val_rerun
#PBS -j oe
#PBS -o logs_rerun/

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

jobid=$(sed -n "${PBS_ARRAY_INDEX}p" "${PBS_O_WORKDIR}/failed_ids.txt")

# Navigate to the script directory
cd "${PBS_O_WORKDIR}/.."

# Run R script - jobid is a row number of study$grid (bin_val_config.R)
Rscript bin_val_analysis.R "$jobid" 1 1
