#PBS -l walltime=00:30:00
#PBS -l select=1:ncpus=2:ompthreads=2:mem=15gb
#PBS -N bin_val_metrics

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to the script directory
cd "${PBS_O_WORKDIR}/.."

Rscript bin_val_metrics.R
