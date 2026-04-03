#!/bin/bash
#SBATCH --job-name=scD3_volcano
#SBATCH --output=logs/05_volcano_%j.out
#SBATCH --error=logs/05_volcano_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --partition=cpu-g2          # Adjust to your partition
#SBATCH --account=togo              # Adjust to your account

# Run a single simulation + volcano plots for all 5 DE methods.
# All I/O is via S3 (bucket: scrna).
#
# USAGE:
#   sbatch 05_volcano.sh
#
# To customize scenario parameters, edit the Rscript arguments below.

BASE_DIR="/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad/scRNA simulation analysis"
SCRIPT_DIR="${BASE_DIR}/R"

mkdir -p "${BASE_DIR}/logs/output" "${BASE_DIR}/logs/error"

module load apptainer

CONTAINER="/mmfs1/gscratch/togo/YC_scRNA.sif"

cd "${BASE_DIR}"

echo "Working directory: $(pwd)"
echo "Script: ${SCRIPT_DIR}/05_volcano_single_sim.R"
echo "Job started at: $(date)"

apptainer exec --bind /mmfs1 "${CONTAINER}" \
  Rscript "${SCRIPT_DIR}/05_volcano_single_sim.R" \
    --prop_de       0.15 \
    --lfc_label     "med" \
    --indiv_var     "no" \
    --cells_label   "low" \
    --corr_cells    0.1 \
    --n_subjects    10 \
    --seed          42 \
    --n_cores       8 \
    --mast_max_cells 5000 \
    --output_name   "volcano_single_sim.pdf"

echo "Step 05 done: $(date)"
