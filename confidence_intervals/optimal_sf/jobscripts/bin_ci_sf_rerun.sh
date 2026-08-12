#!/bin/bash
# Reruns only the array indices listed in failed_bin_ids.txt by bin_ci_sf_check.R.
# -J and the resource request are rewritten by check_failed(); the values here
# are what it computes from bin_ci_sf_1.sh, so the first check leaves them alone.
#
# failed_bin_ids.txt, not failed_ids.txt: this jobscripts directory serves both
# optimal_sf studies, so the bin and cts todo lists have to be named apart or
# each check would clobber the other's. Every other study uses failed_ids.txt.
#
# Needs logs_bin_rerun/ to exist on the cluster - PBS aborts a job whose output
# directory is missing, and *logs*/ is gitignored so it is never checked out.
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=3:ompthreads=2:mem=6gb
#PBS -J 1-1%100
#PBS -N ci_sf_bin_rerun
#PBS -o logs_bin_rerun/
#PBS -e logs_bin_rerun/

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

jobid=$(sed -n "${PBS_ARRAY_INDEX}p" "${PBS_O_WORKDIR}/failed_bin_ids.txt")

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

# Run R script with parameters
Rscript bin_ci_sf_analysis.R "$jobid"
