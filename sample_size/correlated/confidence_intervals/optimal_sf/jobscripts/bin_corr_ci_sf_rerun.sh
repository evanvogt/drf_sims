#!/bin/bash
# Reruns only the array indices listed in failed_bin_ids.txt by bin_corr_ci_sf_check.R.
# -J and the resource request are rewritten by check_failed(); the values here
# are what it computes from bin_corr_ci_sf_1.sh.
#
# failed_bin_ids.txt, not failed_ids.txt: this jobscripts directory serves both
# optimal_sf studies, so the bin and cts todo lists are named apart.
#PBS -l walltime=08:00:00
#PBS -l select=1:ncpus=4:ompthreads=4:mem=9gb
#PBS -J 1-722%50
#PBS -N ci_sf_corr_bin_rerun
#PBS -o logs_bin_rerun/
#PBS -j oe

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

jobid=$(sed -n "${PBS_ARRAY_INDEX}p" "${PBS_O_WORKDIR}/failed_bin_ids.txt")

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

# Run R script with parameters
Rscript bin_corr_ci_sf_analysis.R "$jobid" 4
