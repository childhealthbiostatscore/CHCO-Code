#!/bin/bash
#SBATCH --job-name=scD3_sim_analyze
#SBATCH --array=1-8100%50
#SBATCH --time=01:00:00
#SBATCH --mem=40G             # simulate + run all 3 methods in memory
#SBATCH --cpus-per-task=4
#SBATCH --partition=cpu-g2          # Adjust to your partition
#SBATCH --account=togo              # Adjust to your account
#SBATCH --output="/mmfs1/gscratch/togo/yejichoi/project_logs/scD3_logs/02_sim/output/02_sim_%j.out"
#SBATCH --error="/mmfs1/gscratch/togo/yejichoi/project_logs/scD3_logs/02_sim/error/02_sim_%j.err"
#SBATCH --nodes=1
#SBATCH --ntasks=1

# All I/O is via S3 (bucket: scrna).
# param_grid.rds and reference model are read from S3.
# Stats output files (~1 MB per task) are written to S3.
BASE_DIR="/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad/scRNA simulation analysis"
SCRIPT_DIR="${BASE_DIR}/R"
N_CORES=4

# Failure log: one file per job, collects all failed array IDs
FAIL_LOG="/mmfs1/gscratch/togo/yejichoi/project_logs/scD3_logs/02_sim/failed_arrays_${SLURM_ARRAY_JOB_ID}.txt"

# Load Apptainer module (matches your known-working pattern)
module load apptainer  # or module load singularity

# Container path
CONTAINER="/mmfs1/gscratch/togo/YC_scRNA.sif"

# Move to base dir (matches your known-working pattern)
cd "${BASE_DIR}"
echo "Working directory: $(pwd)"
echo "Script: ${SCRIPT_DIR}/02_simulate_analyze.R"
echo "Job started at: $(date)"

echo "── Task ${SLURM_ARRAY_TASK_ID} / ${SLURM_ARRAY_TASK_MAX} ──"

# Run inside container
apptainer exec --bind /mmfs1 "${CONTAINER}" \
  Rscript "${SCRIPT_DIR}/02_simulate_analyze.R"   \
      --array_id     "${SLURM_ARRAY_TASK_ID}"      \
      --n_cores      "${N_CORES}"                  \
      --nebula_method "LN"

# Log failed array IDs for easy resubmission
if [ $? -ne 0 ]; then
    echo "${SLURM_ARRAY_TASK_ID}" >> "${FAIL_LOG}"
    echo "── Task ${SLURM_ARRAY_TASK_ID} FAILED at: $(date) ──"
else
    echo "── Task ${SLURM_ARRAY_TASK_ID} done: $(date) ──"
fi
