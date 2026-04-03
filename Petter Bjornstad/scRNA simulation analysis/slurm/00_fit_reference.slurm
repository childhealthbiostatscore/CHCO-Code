#!/bin/bash
#SBATCH --job-name=scD3_fit_ref
#SBATCH --output="/mmfs1/gscratch/togo/yejichoi/project_logs/scD3_logs/00_fit_reference/output/00_fit_reference_%j.out"
#SBATCH --error="/mmfs1/gscratch/togo/yejichoi/project_logs/scD3_logs/00_fit_reference/error/00_fit_reference_%j.err"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=12:00:00
#SBATCH --partition=cpu-g2          # Adjust to your partition
#SBATCH --account=togo              # Adjust to your account

# ─── Configuration ─────────────────────────────────────────────────────────
CELL_TYPE="PT"   # <-- update if needed
N_CORES=8

# Project directories (quoted because of spaces)
BASE_DIR="/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad/scRNA simulation analysis"
SCRIPT_DIR="${BASE_DIR}/R"

# Load Apptainer module
module load apptainer  # or module load singularity

# Container path
CONTAINER="/mmfs1/gscratch/togo/YC_scRNA.sif"

# Move to base dir
cd "${BASE_DIR}"

echo "Working directory: $(pwd)"
echo "Script: ${SCRIPT_DIR}/00_fit_reference.R"
echo "CELL_TYPE: ${CELL_TYPE}"
echo "N_CORES: ${N_CORES}"
echo "Job started at: $(date)"

# Run inside container
apptainer exec --bind /mmfs1 "${CONTAINER}" \
  Rscript "${SCRIPT_DIR}/00_fit_reference.R" \
    --cell_type "${CELL_TYPE}" \
    --n_cores   "${N_CORES}"

echo "Step 00 done: $(date)"