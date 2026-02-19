#!/bin/bash
#SBATCH --job-name=scD3_fit_ref
#SBATCH --output=logs/00_fit_reference_%j.out
#SBATCH --error=logs/00_fit_reference_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --partition=compute          # adjust to your HPC partition

# ─── Configuration ─────────────────────────────────────────────────────────
# attempt_so is loaded from S3 automatically by the R script (bucket: attempt).
# Set the desired cell type here.
CELL_TYPE="PT"                              # <-- UPDATE to desired cell type
SCRIPT_DIR="$(cd "$(dirname "$0")/.." && pwd)/R"
N_CORES=32

mkdir -p logs

module purge
module load R/4.3.0                         # adjust to your HPC module

Rscript "${SCRIPT_DIR}/00_fit_reference.R" \
    --cell_type "${CELL_TYPE}"             \
    --n_cores   "${N_CORES}"

echo "Step 00 done: $(date)"
