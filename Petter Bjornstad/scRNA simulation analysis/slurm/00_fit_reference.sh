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
# Edit these paths before submitting
SEURAT_PATH="/path/to/attempt_so.rds"      # <-- UPDATE THIS
CELL_TYPE="PT"                              # <-- UPDATE to desired cell type
SCRIPT_DIR="$(cd "$(dirname "$0")/.." && pwd)/R"
RESULTS_DIR="$(cd "$(dirname "$0")/.." && pwd)/results"
N_CORES=32

mkdir -p "${RESULTS_DIR}/reference" logs

module purge
module load R/4.3.0                         # adjust to your HPC module

Rscript "${SCRIPT_DIR}/00_fit_reference.R" \
    --seurat_path "${SEURAT_PATH}"          \
    --cell_type   "${CELL_TYPE}"            \
    --out_dir     "${RESULTS_DIR}/reference" \
    --n_cores     "${N_CORES}"

echo "Step 00 done: $(date)"
