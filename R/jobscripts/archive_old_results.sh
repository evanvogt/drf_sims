#!/bin/bash
# Archives every study's pre-2026-09-26 results tree, one .tar per study - see
# the header of R/archive_old_results.R. Submit from the repo root, after a dry
# run there (Rscript R/archive_old_results.R) has shown what it will do:
#   qsub R/jobscripts/archive_old_results.sh
# Packing and listing ~100k files is too much for a login node. The walltime
# is generous: tar is disk-bound, and a tree that doesn't finish is left intact.
#PBS -l walltime=08:00:00
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -N archive_old_results
#PBS -j oe

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

cd "${PBS_O_WORKDIR}/.."

Rscript archive_old_results.R --apply
