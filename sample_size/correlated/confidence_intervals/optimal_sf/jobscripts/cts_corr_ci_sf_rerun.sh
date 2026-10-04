#!/bin/bash
# Reruns only the array indices listed in failed_cts_ids.txt by cts_corr_ci_sf_check.R.
# -J and the resource request are rewritten by check_failed(); the values here
# are what it computes from cts_corr_ci_sf_1.sh.
#
# failed_cts_ids.txt, not failed_ids.txt: this jobscripts directory serves both
# optimal_sf studies, so the bin and cts todo lists are named apart.
#PBS -l walltime=16:00:00
#PBS -l select=1:ncpus=3:ompthreads=2:mem=9gb
#PBS -J 1-2%100
#PBS -N ci_sf_corr_cts_rerun
#PBS -o logs_cts_rerun/
#PBS -j oe

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

jobid=$(sed -n "${PBS_ARRAY_INDEX}p" "${PBS_O_WORKDIR}/failed_cts_ids.txt")

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

# Run R script with parameters
Rscript cts_corr_ci_sf_analysis.R "$jobid"
