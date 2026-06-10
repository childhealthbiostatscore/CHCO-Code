#!/bin/bash
#SBATCH --job-name=triad_nebula
#SBATCH --array=1-84%10
#SBATCH --time=24:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=6
#SBATCH --partition=cpu-2
#SBATCH --account=togo
#SBATCH --output="/mmfs1/gscratch/togo/leidholt/project_logs/triad_nebula/output/nebula_%A_%a.out"
#SBATCH --error="/mmfs1/gscratch/togo/leidholt/project_logs/triad_nebula/error/nebula_%A_%a.err"

module purge
module load apptainer

CONTAINER="/mmfs1/gscratch/togo/yzhangtufts_r_scrnaseq.sif"

BASE_DIR="/mmfs1/gscratch/togo/leidholt/CHCO-Code/Petter Bjornstad/T1D-T2D-multiOmics/nebula"

cd "${BASE_DIR}"

CONFIG_FILE="${BASE_DIR}/config/job_config.txt"

JOB_PARAMS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${CONFIG_FILE}")

RUN_NAME=$(echo "$JOB_PARAMS" | cut -d$'\t' -f1)
CELL_TYPE=$(echo "$JOB_PARAMS" | cut -d$'\t' -f2)


echo "============================================================="
echo "SLURM Job ID: ${SLURM_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Working directory: $(pwd)"
echo "Config file: ${CONFIG_FILE}"
echo "Run name / contrast: ${RUN_NAME}"
echo "Cell type: ${CELL_TYPE}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Job started at: $(date)"
echo "============================================================="

apptainer exec \
  --cleanenv \
  --bind /mmfs1:/mmfs1 \
  "${CONTAINER}" \
  Rscript "${BASE_DIR}/scripts/run_nebula.R" \
  "${RUN_NAME}" \
  "${CELL_TYPE}" 

EXIT_CODE=$?

echo "============================================================="
echo "Job completed at: $(date)"
echo "Exit code: ${EXIT_CODE}"
echo "============================================================="

exit ${EXIT_CODE}
