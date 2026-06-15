#!/bin/bash
#SBATCH --job-name=triad_nebula
#SBATCH --array=1-84%6
#SBATCH --time=72:00:00
#SBATCH --mem=100
#SBATCH --cpus-per-task=30
#SBATCH --partition=cpu-g2
#SBATCH --account=togo
#SBATCH --output="/mmfs1/gscratch/togo/leidholt/project_logs/triad_nebula.5/output/nebula_%A_%a.out"
#SBATCH --error="/mmfs1/gscratch/togo/leidholt/project_logs/triad_nebula.5/error/nebula_%A_%a.err"

module purge
module load apptainer

CONTAINER="/mmfs1/gscratch/togo/yzhangtufts_r_scrnaseq.sif"

BASE_DIR="/mmfs1/gscratch/togo/leidholt/CHCO-Code/Savanah Leidholt/T1D-T2D-multiOmics/nebula"

cd "${BASE_DIR}"

CONFIG_FILE="${BASE_DIR}/nebula_array.tsv"

JOB_PARAMS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${CONFIG_FILE}")

CELL_TYPE=$(echo "$JOB_PARAMS" | awk '{print $1}')
CONTRAST_NAME=$(echo "$JOB_PARAMS" | awk '{print $2}')


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
  Rscript "${BASE_DIR}/run_nebula_5.R" \
  "${CELL_TYPE}" \
  "${CONTRAST_NAME}" 

EXIT_CODE=$?

echo "============================================================="
echo "Job completed at: $(date)"
echo "Exit code: ${EXIT_CODE}"
echo "============================================================="

exit ${EXIT_CODE}
