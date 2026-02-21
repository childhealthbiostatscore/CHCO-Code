#!/bin/bash
#SBATCH --job-name=scD3_sim_analyze
#SBATCH --array=1-5000%20
#SBATCH --array=5001-10000%20
#SBATCH --array=10001-15000%20
#SBATCH --array=15001-21600%20
#SBATCH --time=02:00:00
#SBATCH --mem=32G             # simulate + run all 3 methods in memory
#SBATCH --cpus-per-task=4
#SBATCH --partition=cpu-g2          # Adjust to your partition
#SBATCH --account=togo              # Adjust to your account
#SBATCH --output="/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad/scRNA simulation analysis/logs/output/02_sim_%j.out"
#SBATCH --error="/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad/scRNA simulation analysis/logs/error/02_sim_%j.err"
#SBATCH --nodes=1
#SBATCH --ntasks=1

# All I/O is via S3 (bucket: scrna).
# param_grid.rds and reference model are read from S3.
# Stats output files (~1 MB per task) are written to S3.
BASE_DIR="/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad/scRNA simulation analysis"
SCRIPT_DIR="${BASE_DIR}/R"
N_CORES=4

mkdir -p "${BASE_DIR}/logs/output" "${BASE_DIR}/logs/error"

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

echo "── Task ${SLURM_ARRAY_TASK_ID} done: $(date) ──"
