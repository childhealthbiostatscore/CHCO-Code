#!/bin/bash
#SBATCH --job-name=scD3_param_grid
#SBATCH --output="/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad/scRNA simulation analysis/logs/output/01_param_grid_%j.out"
#SBATCH --error="/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad/scRNA simulation analysis/logs/error/01_param_grid_%j.err"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --partition=cpu-g2          # Adjust to your partition
#SBATCH --account=togo              # Adjust to your account

# All I/O is via S3 (bucket: scrna).
# effect_size_summary.rds is read from S3; param_grid.rds/.csv are written to S3.
BASE_DIR="/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad/scRNA simulation analysis"
SCRIPT_DIR="${BASE_DIR}/R"
N_REPS=50

mkdir -p "${BASE_DIR}/logs/output" "${BASE_DIR}/logs/error"

# Load Apptainer module (matches your known-working pattern)
module load apptainer  # or module load singularity

# Container path
CONTAINER="/mmfs1/gscratch/togo/YC_scRNA.sif"

# Move to base dir (matches your known-working pattern)
cd "${BASE_DIR}"

echo "Working directory: $(pwd)"
echo "Script: ${SCRIPT_DIR}/01_parameter_grid.R"
echo "Job started at: $(date)"

# Run inside container
apptainer exec --bind /mmfs1 "${CONTAINER}" \
  Rscript "${SCRIPT_DIR}/01_parameter_grid.R" \
      --n_reps "${N_REPS}"

echo "Step 01 done: $(date)"