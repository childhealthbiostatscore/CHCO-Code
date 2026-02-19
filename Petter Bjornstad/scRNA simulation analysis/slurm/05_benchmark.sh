#!/bin/bash
#SBATCH --job-name=scD3_bench
#SBATCH --output=logs/05_benchmark_%j.out
#SBATCH --error=logs/05_benchmark_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --partition=compute

SCRIPT_DIR="$(cd "$(dirname "$0")/.." && pwd)/R"
RESULTS_DIR="$(cd "$(dirname "$0")/.." && pwd)/results"
N_CORES=32

mkdir -p "${RESULTS_DIR}/benchmark" logs

module purge
module load R/4.3.0

Rscript "${SCRIPT_DIR}/05_benchmark.R"                              \
    --param_grid   "${RESULTS_DIR}/param_grid/param_grid.rds"       \
    --stats_root   "${RESULTS_DIR}/stats"                           \
    --out_dir      "${RESULTS_DIR}/benchmark"                       \
    --n_cores      "${N_CORES}"

echo "Step 05 done: $(date)"
