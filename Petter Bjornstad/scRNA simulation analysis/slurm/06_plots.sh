#!/bin/bash
#SBATCH --job-name=scD3_plots
#SBATCH --output=logs/06_plots_%j.out
#SBATCH --error=logs/06_plots_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --partition=compute

SCRIPT_DIR="$(cd "$(dirname "$0")/.." && pwd)/R"
RESULTS_DIR="$(cd "$(dirname "$0")/.." && pwd)/results"

mkdir -p "${RESULTS_DIR}/plots" logs

module purge
module load R/4.3.0

Rscript "${SCRIPT_DIR}/06_plots.R" \
    --benchmark_avg "${RESULTS_DIR}/benchmark/benchmark_avg.rds" \
    --benchmark_sum "${RESULTS_DIR}/benchmark/benchmark_summary.rds" \
    --out_dir       "${RESULTS_DIR}/plots"

echo "Step 06 done: $(date)"
