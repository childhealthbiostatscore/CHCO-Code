#!/bin/bash
#SBATCH --job-name=scD3_param_grid
#SBATCH --output=logs/01_param_grid_%j.out
#SBATCH --error=logs/01_param_grid_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --partition=compute

SCRIPT_DIR="$(cd "$(dirname "$0")/.." && pwd)/R"
RESULTS_DIR="$(cd "$(dirname "$0")/.." && pwd)/results"
N_REPS=50

mkdir -p "${RESULTS_DIR}/param_grid" logs

module purge
module load R/4.3.0

Rscript "${SCRIPT_DIR}/01_parameter_grid.R" \
    --effect_summary "${RESULTS_DIR}/reference/effect_size_summary.rds" \
    --out_dir        "${RESULTS_DIR}/param_grid" \
    --n_reps         "${N_REPS}"

echo "Step 01 done: $(date)"
# Print grid size for next steps
Rscript -e "pg <- readRDS('${RESULTS_DIR}/param_grid/param_grid.rds'); cat('Grid rows:', nrow(pg), '\n')"
