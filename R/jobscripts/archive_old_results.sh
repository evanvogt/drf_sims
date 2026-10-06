#!/bin/bash
# Archives results trees, one .tar per tree - see the header of
# R/archive_old_results.R. Like every jobscript here it cds to
# ${PBS_O_WORKDIR}/.., so submit from R/jobscripts/, after a dry run from the
# repo root (Rscript R/archive_old_results.R [args]) has shown what it will do:
#   cd R/jobscripts
#   qsub archive_old_results.sh                     # every study's pre-2026-09-26 tree
#   qsub -v ARCHIVE_ARGS="--label retired_2026-10 --trees continuous binary confidence_intervals" archive_old_results.sh
# ARCHIVE_ARGS is passed through to the R script; keep it comma-free (qsub -v
# splits on commas), which is why --trees is space-separated.
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

# unquoted on purpose: ARCHIVE_ARGS is several words
Rscript archive_old_results.R --apply ${ARCHIVE_ARGS}
