#!/bin/bash
#SBATCH --job-name=scD3_fit_ref
#SBATCH --time=12:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=32
#SBATCH --partition=cpu-g2          # Adjust to your partition
#SBATCH --account=togo              # Adjust to your account
#SBATCH --output="/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad/scRNA simulation analysis/logs/00_fit_reference__%j.out"
#SBATCH --error="/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad/scRNA simulation analysis/logs/00_fit_reference__%j.err"

# Load Apptainer module (if needed on your system)
module load apptainer  # or module load singularity

# Set the container path
CONTAINER="/mmfs1/gscratch/togo/YC_scRNA.sif"

# Set the base directory for your project - PROPERLY QUOTED
BASE_DIR="/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad/Renal HEIRitage/nebula_parallel"

# Change to base directory
cd "${BASE_DIR}"

# Get the job parameters from config file
CONFIG_FILE="${BASE_DIR}/config/job_config.txt"
JOB_PARAMS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${CONFIG_FILE}")

# Parse analysis type and cell type
# Parse run name, cell type, and cell type variable
RUN_NAME=$(echo $JOB_PARAMS | cut -d' ' -f1)
CELL_TYPE=$(echo $JOB_PARAMS | cut -d' ' -f2)
CELL_TYPE_VAR=$(echo $JOB_PARAMS | cut -d' ' -f3)

echo "Working directory: $(pwd)"
echo "Config file: ${CONFIG_FILE}"
echo "Starting job ${SLURM_ARRAY_TASK_ID}: ${RUN_NAME} for ${CELL_TYPE} (variable: ${CELL_TYPE_VAR})"
echo "Job started at: $(date)"

# Run the R script using Apptainer
apptainer exec --bind /mmfs1 ${CONTAINER} Rscript "${BASE_DIR}/scripts/run_nebula_single.R" "${RUN_NAME}" "${CELL_TYPE}" "${CELL_TYPE_VAR}"
echo "Job completed at: $(date)"