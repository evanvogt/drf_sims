#!/bin/bash
# mem=20gb, not continuous's 15gb: per-run objects here are heavier (9
# candidate models + 2 nuisance pipelines x 2 fold regimes of prediction
# data, vs. continuous's 5 CATE-model results). Measured at ~0.9MB/run via
# object.size() against ~0.18MB/run for a comparable continuous object -
# heavier per file, but this study has far fewer files (360 vs. 4000), so
# the aggregate is comparable order of magnitude. See README.md.
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=2:ompthreads=2:mem=20gb
#PBS -N me_collect


module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to the script directory
cd "${PBS_O_WORKDIR}/.."

Rscript me_collect.R
