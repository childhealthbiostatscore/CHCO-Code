#!/bin/bash
#SBATCH --job-name=prep_nebula
#SBATCH --array=1-14%2
#SBATCH --time=24:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=30
#SBATCH --partition=ckpt
#SBATCH --account=togo
#SBATCH --output="/mmfs1/gscratch/togo/leidholt/project_logs/triad_nebula.5/prep_output/prep_%A_%a.out"
#SBATCH --error="/mmfs1/gscratch/togo/leidholt/project_logs/triad_nebula.5/prep_error/prep_%A_%a.err"

module purge
module load apptainer

CONTAINER="/mmfs1/gscratch/togo/yzhangtufts_r_scrnaseq.sif"
BASE_DIR="/mmfs1/gscratch/togo/leidholt/CHCO-Code/Savanah Leidholt/T1D-T2D-multiOmics/nebula"
CONFIG_FILE="${BASE_DIR}/prep_celltype_array.tsv"

mkdir -p /mmfs1/gscratch/togo/leidholt/project_logs/triad_nebula.5/prep_output
mkdir -p /mmfs1/gscratch/togo/leidholt/project_logs/triad_nebula.5/prep_error

cd "${BASE_DIR}"

CELL_TYPE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${CONFIG_FILE}" | awk '{print $1}')

echo "============================================================="
echo "SLURM Job ID: ${SLURM_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Cell type: ${CELL_TYPE}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Job started at: $(date)"
echo "============================================================="

apptainer exec \
  --cleanenv \
  --bind /mmfs1:/mmfs1 \
  "${CONTAINER}" \
  Rscript "${BASE_DIR}/prep_nebula_celltype.R" \
  "${CELL_TYPE}"

EXIT_CODE=$?

echo "============================================================="
echo "Job completed at: $(date)"
echo "Exit code: ${EXIT_CODE}"
echo "============================================================="

exit ${EXIT_CODE}
