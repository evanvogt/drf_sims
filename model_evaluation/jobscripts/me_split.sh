#!/bin/bash
# The 80:20 single-split arm. n = 500 and 1000 only, 30 runs each, 4
# scenarios - 240 array indices, NOT 360. See me_config.R's study_split.
#
# Cheaper per replicate than me_1.sh despite refitting the candidates: each
# candidate is crossfit within the 80% once, and the nuisance is a single
# whole-set fit on the 20% rather than a 10-fold pipeline over all n.
# Still a PLACEHOLDER - profile before trusting it.
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=4:ompthreads=2:mem=10gb
#PBS -J 1-240%4
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
Rscript me_split.R "$PBS_ARRAY_INDEX" 2 2
