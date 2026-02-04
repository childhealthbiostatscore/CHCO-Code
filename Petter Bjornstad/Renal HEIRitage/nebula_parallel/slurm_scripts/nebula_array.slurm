#!/bin/bash
#SBATCH --job-name=nebula_array
#SBATCH --array=1-18%10             # Run 168 jobs, max 7 concurrent
#SBATCH --time=03:00:00             # 3 hour per job (30min + buffer)
#SBATCH --mem=200G                  # Memory per job
#SBATCH --cpus-per-task=25          # CPUs (nebula uses 15 cores)
#SBATCH --partition=cpu-g2          # Adjust to your partition
#SBATCH --account=togo              # Adjust to your account
#SBATCH --output="/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad/Renal HEIRitage/nebula_parallel/logs/output/nebula_%A_%a.out"
#SBATCH --error="/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad/Renal HEIRitage/nebula_parallel/logs/error/nebula_%A_%a.err"

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