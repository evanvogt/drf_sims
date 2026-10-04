#!/bin/bash
# Reruns only the array indices listed in failed_ids.txt by cts_corr_ci_check.R.
# -J and the resource request are rewritten by check_failed(); the values here
# are what it computes from cts_corr_ci_1.sh.
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=3:ompthreads=2:mem=4gb
#PBS -J 1-2%100
#PBS -N ci_corr_cts_rerun
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
Rscript cts_corr_ci_analysis.R "$jobid"
