#!/bin/bash
# The 80:20 single-split arm. n = 500 and 1000 only, 50 runs each, 4
# scenarios - 400 array indices, NOT 600. See me_config.R's study_split.
#
# Cheaper per replicate than me_1.sh despite refitting the candidates: each
# candidate is crossfit within the 80% once, and the nuisance is a single
# whole-set fit on the 20% rather than a 10-fold pipeline over all n.
# 10 cores, ompthreads=1 and mem for the same reasons as me_1.sh: the
# candidates' 10-fold crossfit runs in one round. Walltime is still a
# PLACEHOLDER - check the first subjobs' resources_used before trusting it.
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=10:ompthreads=1:mem=24gb
#PBS -J 1-400%4
#PBS -N me_split
#PBS -o logs_split/
#PBS -e logs_split/


module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

# Trailing args are workers (future plan for the candidate crossfit) and
# n_cores (XGBoost nthread / H2O nthreads).
Rscript me_split.R "$PBS_ARRAY_INDEX" 10 10
