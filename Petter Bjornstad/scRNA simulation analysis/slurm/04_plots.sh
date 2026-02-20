#!/bin/bash
#SBATCH --job-name=scD3_plots
#SBATCH --output=logs/04_plots_%j.out
#SBATCH --error=logs/04_plots_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --partition=compute

# All I/O is via S3 (bucket: scrna).
# Reads benchmark_avg.rds from S3; saves PDF figures to S3.
SCRIPT_DIR="$(cd "$(dirname "$0")/.." && pwd)/R"

mkdir -p logs

module purge
module load R/4.3.0

Rscript "${SCRIPT_DIR}/04_plots.R"

echo "Step 06 done: $(date)"
