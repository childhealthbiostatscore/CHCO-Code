#!/bin/bash
#SBATCH --job-name=scD3_plots
#SBATCH --output=logs/04_plots_%j.out
#SBATCH --error=logs/04_plots_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --partition=cpu-g2          # Adjust to your partition
#SBATCH --account=togo              # Adjust to your account

# All I/O is via S3 (bucket: scrna).
# Reads benchmark_avg.rds from S3; saves PDF figures to S3.
BASE_DIR="/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad/scRNA simulation analysis"
SCRIPT_DIR="${BASE_DIR}/R"

mkdir -p "${BASE_DIR}/logs/output" "${BASE_DIR}/logs/error"

# Load Apptainer module (matches your known-working pattern)
module load apptainer  # or module load singularity

# Container path
CONTAINER="/mmfs1/gscratch/togo/YC_scRNA.sif"

# Move to base dir (matches your known-working pattern)
cd "${BASE_DIR}"

echo "Working directory: $(pwd)"
echo "Script: ${SCRIPT_DIR}/04_plots.R"
echo "Job started at: $(date)"

apptainer exec --bind /mmfs1 "${CONTAINER}" \
  Rscript "${SCRIPT_DIR}/04_plots.R"

echo "Step 04 done: $(date)"
