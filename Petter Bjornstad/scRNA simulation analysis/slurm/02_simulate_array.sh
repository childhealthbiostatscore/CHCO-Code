#!/bin/bash
#SBATCH --job-name=scD3_sim_analyze
#SBATCH --output=logs/02_sim_%A_%a.out
#SBATCH --error=logs/02_sim_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G             # simulate + run all 3 methods in memory
#SBATCH --time=02:00:00
#SBATCH --partition=compute
#SBATCH --array=1-21600%200   # 21,600 total tasks; max 200 running at once
                               # Update upper bound from param_grid if needed

# All I/O is via S3 (bucket: scrna).
# param_grid.rds and reference model are read from S3.
# Stats output files (~1 MB per task) are written to S3.
SCRIPT_DIR="$(cd "$(dirname "$0")/.." && pwd)/R"
N_CORES=4

mkdir -p logs

module purge
module load R/4.3.0            # adjust to your HPC R module

echo "── Task ${SLURM_ARRAY_TASK_ID} / ${SLURM_ARRAY_TASK_MAX} ──"

Rscript "${SCRIPT_DIR}/02_simulate_analyze.R"   \
    --array_id     "${SLURM_ARRAY_TASK_ID}"      \
    --n_cores      "${N_CORES}"                  \
    --nebula_method "LN"

echo "── Task ${SLURM_ARRAY_TASK_ID} done: $(date) ──"
