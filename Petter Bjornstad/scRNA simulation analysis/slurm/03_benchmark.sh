#!/bin/bash
#SBATCH --job-name=scD3_bench
#SBATCH --output=logs/03_benchmark_%j.out
#SBATCH --error=logs/03_benchmark_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --partition=compute

# All I/O is via S3 (bucket: scrna).
# Reads param_grid and per-task stats from S3; writes benchmark outputs to S3.
SCRIPT_DIR="$(cd "$(dirname "$0")/.." && pwd)/R"
N_CORES=32

mkdir -p logs

module purge
module load R/4.3.0

Rscript "${SCRIPT_DIR}/03_benchmark.R" \
    --n_cores "${N_CORES}"

echo "Step 05 done: $(date)"
