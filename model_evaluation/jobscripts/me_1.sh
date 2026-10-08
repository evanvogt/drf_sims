#!/bin/bash
# The main study: 4 scenarios x 3 n x 50 runs = 600 grid rows.
#
# 10 cores (2026-10-08): the candidate phase crossfits each of the 9
# candidates over 10 folds in parallel (`workers`), so 10 workers run every
# fold loop in one round - 5 would take two, 8 still two. The nuisance phase
# (`n_cores`, XGBoost nthread / H2O nthreads) never overlaps it, so it gets the
# same 10. ompthreads=1 stops the 10 worker sessions each inheriting a
# multi-threaded BLAS; XGBoost and H2O set their own threads from n_cores.
#
# mem covers the H2O JVM (max heap 10G) plus the 10 worker R sessions, which
# stay alive through the nuisance phase, plus the main session.
#
# Walltime is still a PLACEHOLDER - check the first subjobs' resources_used
# (qstat -fx <jobid> | grep resources_used) before trusting it.
#
# %4: each concurrent task starts its own H2O JVM, and more at once makes the
# H2O calls fail. See README.md.
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=10:ompthreads=1:mem=24gb
#PBS -J 1-600%4
#PBS -N me_1
#PBS -o logs_1/
#PBS -e logs_1/


module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

# Run R script with parameters. Trailing args are workers/n_cores, set by
# hand to match ncpus above - change them together.
Rscript me_analysis.R "$PBS_ARRAY_INDEX" 10 10
