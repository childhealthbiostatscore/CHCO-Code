#!/bin/bash
#SBATCH --job-name=scD3_bench
#SBATCH --output="/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad/scRNA simulation analysis/logs/output/03_benchmark_%j.out""
#SBATCH --error="/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad/scRNA simulation analysis/logs/error/logs/03_benchmark_%j.err""
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --partition=cpu-g2          # Adjust to your partition
#SBATCH --account=togo              # Adjust to your account

# All I/O is via S3 (bucket: scrna).
# Reads param_grid and per-task stats from S3; writes benchmark outputs to S3.
BASE_DIR="/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad/scRNA simulation analysis"
SCRIPT_DIR="${BASE_DIR}/R"
N_CORES=32

mkdir -p "${BASE_DIR}/logs/output" "${BASE_DIR}/logs/error"

# Load Apptainer module (matches your known-working pattern)
module load apptainer  # or module load singularity

# Container path
CONTAINER="/mmfs1/gscratch/togo/YC_scRNA.sif"

# Move to base dir (matches your known-working pattern)
cd "${BASE_DIR}"

echo "Working directory: $(pwd)"
echo "Script: ${SCRIPT_DIR}/03_benchmark.R"
echo "Job started at: $(date)"

apptainer exec --bind /mmfs1 "${CONTAINER}" \
  Rscript "${SCRIPT_DIR}/03_benchmark.R" \
      --n_cores "${N_CORES}"

echo "Step 03 done: $(date)"
