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

# All I/O is via S3 (bucket: scrna).
# effect_size_summary.rds is read from S3; param_grid.rds/.csv are written to S3.
SCRIPT_DIR="$(cd "$(dirname "$0")/.." && pwd)/R"
N_REPS=50

mkdir -p logs

module purge
module load R/4.3.0

Rscript "${SCRIPT_DIR}/01_parameter_grid.R" \
    --n_reps "${N_REPS}"

echo "Step 01 done: $(date)"
