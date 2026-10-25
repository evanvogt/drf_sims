#!/bin/bash
# The case-study learners (me_cs_learners.R): scores DR-RF, the causal forest
# and DR-SL by the DR score in every simulated trial. Nothing is fit - it reads
# each strategies-tree run and the sample-size study's matching run - so this
# is one job, not an array. Pure I/O over ~2,400 files.
#
# Run after the strategies tree's me_metrics.sh (its SL1 check reads
# me_strat_metrics.RDS), and before
# results_processing/thesis_tables/me_cs_tables.R.
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=2:ompthreads=2:mem=16gb
#PBS -N me_cs_learners


module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to the script directory
cd "${PBS_O_WORKDIR}/.."

Rscript me_cs_learners.R
