#!/bin/bash
#SBATCH --job-name=t1d_adip_nebula
#SBATCH --array=5001-5546%100           # up to run 5546
#SBATCH --time=03:00:00             # 3 hours per job (with buffer)
#SBATCH --mem=200G                  # Memory per job
#SBATCH --cpus-per-task=25          # CPUs (nebula uses 15 cores)
#SBATCH --partition=ckpt          # Adjust to your partition
#SBATCH --account=togo              # Adjust to your account
#SBATCH --output="/mmfs1/gscratch/togo/yejichoi/project_logs/t1d_adiposity_logs/output/nebula_%A_%a.out"
#SBATCH --error="/mmfs1/gscratch/togo/yejichoi/project_logs/t1d_adiposity_logs/error/nebula_%A_%a.err"

# Load Apptainer module
module load apptainer

# Set the container path
CONTAINER="/mmfs1/gscratch/togo/YC_scRNA.sif"

# Set the base directory for your project
BASE_DIR="/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad/T1D Adiposity/nebula"

# Change to base directory
cd "${BASE_DIR}"

# Get the job parameters from config file
CONFIG_FILE="${BASE_DIR}/config/job_config.txt"
JOB_PARAMS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${CONFIG_FILE}")

# Parse run name, cell type, and cell type variable (tab-delimited)
RUN_NAME=$(echo "$JOB_PARAMS" | cut -d$'\t' -f1)
CELL_TYPE=$(echo "$JOB_PARAMS" | cut -d$'\t' -f2)
CELL_TYPE_VAR=$(echo "$JOB_PARAMS" | cut -d$'\t' -f3)

echo "============================================================="
echo "SLURM Job ID: ${SLURM_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Working directory: $(pwd)"
echo "Config file: ${CONFIG_FILE}"
echo "Analysis: ${RUN_NAME}"
echo "Cell type: ${CELL_TYPE}"
echo "Cell type variable: ${CELL_TYPE_VAR}"
echo "Job started at: $(date)"
echo "============================================================="

# Run the R script using Apptainer
apptainer exec --bind /mmfs1 ${CONTAINER} Rscript "${BASE_DIR}/scripts/run_nebula_single.R" "${RUN_NAME}" "${CELL_TYPE}" "${CELL_TYPE_VAR}"

EXIT_CODE=$?
echo "============================================================="
echo "Job completed at: $(date)"
echo "Exit code: ${EXIT_CODE}"
echo "============================================================="

exit ${EXIT_CODE}
