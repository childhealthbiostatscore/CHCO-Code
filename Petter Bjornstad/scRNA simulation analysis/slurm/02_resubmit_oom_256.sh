#!/bin/bash
#SBATCH --job-name=scD3_oom256
#SBATCH --time=03:00:00
#SBATCH --mem=256G
#SBATCH --cpus-per-task=4
#SBATCH --partition=ckpt
#SBATCH --account=togo
#SBATCH --output="/mmfs1/gscratch/togo/yejichoi/project_logs/scD3_logs/02_sim/output/02_sim_oom256_%A_%a.out"
#SBATCH --error="/mmfs1/gscratch/togo/yejichoi/project_logs/scD3_logs/02_sim/error/02_sim_oom256_%A_%a.err"
#SBATCH --nodes=1
#SBATCH --ntasks=1

BASE_DIR="/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad/scRNA simulation analysis"
SCRIPT_DIR="${BASE_DIR}/R"
N_CORES=4

# Failure log
FAIL_LOG="/mmfs1/gscratch/togo/yejichoi/project_logs/scD3_logs/02_sim/failed_oom256_${SLURM_ARRAY_JOB_ID}.txt"

# File containing array IDs to rerun (one per line)
ID_FILE="mmfs1/gscratch/togo/yejichoi/project_logs/scD3_logs/02_sim/failed_arrays_33794735.txt"

# Get the PARAM_ID for this task from the ID file
# SLURM_ARRAY_TASK_ID is the line number (1-indexed)
PARAM_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${ID_FILE}")

if [ -z "${PARAM_ID}" ]; then
    echo "ERROR: No ID found at line ${SLURM_ARRAY_TASK_ID} in ${ID_FILE}"
    exit 1
fi

module load apptainer

CONTAINER="/mmfs1/gscratch/togo/YC_scRNA.sif"

cd "${BASE_DIR}"
echo "Working directory: $(pwd)"
echo "Resubmission (256GB): SLURM task ${SLURM_ARRAY_TASK_ID} → PARAM_ID ${PARAM_ID}"
echo "Job started at: $(date)"

apptainer exec --bind /mmfs1 "${CONTAINER}" \
  Rscript "${SCRIPT_DIR}/02_simulate_analyze.R"   \
      --array_id     "${PARAM_ID}"                 \
      --n_cores      "${N_CORES}"                  \
      --nebula_method "LN"

if [ $? -ne 0 ]; then
    echo "${PARAM_ID}" >> "${FAIL_LOG}"
    echo "── Task ${PARAM_ID} FAILED at: $(date) ──"
else
    echo "── Task ${PARAM_ID} done: $(date) ──"
fi