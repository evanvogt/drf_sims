#!/bin/bash
# Reruns only the array indices listed in failed_ids_split.txt - the split
# tree's counterpart of me_strat_rerun.sh. Nothing writes that file
# automatically (me_check.R checks the split tree with write = FALSE, for the
# reason given in its header); write it by hand, or with grid_indices() as in
# README.md's "Rerunning part of a derived tree", and set -J to its line count. Keep the %N at 10 or
# below, as in me_split.sh.
#
# Resources and trailing args as me_split.sh's - change them together.
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=10:ompthreads=10:mem=24gb
#PBS -J 1-400%10
#PBS -N me_split_rerun
#PBS -o logs_split_rerun/
#PBS -e logs_split_rerun/


module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

jobid=$(sed -n "${PBS_ARRAY_INDEX}p" "${PBS_O_WORKDIR}/failed_ids_split.txt")

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

Rscript me_split.R "$jobid" 10 10
