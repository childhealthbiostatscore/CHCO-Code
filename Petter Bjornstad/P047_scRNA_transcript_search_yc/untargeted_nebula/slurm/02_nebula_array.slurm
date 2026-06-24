#!/bin/bash
#SBATCH --job-name=neb_run
#SBATCH --account=togo
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
# NOTE: --array=1-N is supplied by 01_split.slurm at submit time.
# To run a chunk manually:  sbatch --array=1-50 slurm/02_nebula_array.slurm

# ---------------------------------------------------------------------------
# Stage 2: one NEBULA run per array task (= one line of jobs.tsv).
# ---------------------------------------------------------------------------
set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad/P047_scRNA_transcript_search_yc/untargeted_nebula}"

# module load <your R module>
RSCRIPT="${RSCRIPT:-Rscript}"

cd "$PROJECT_DIR"
echo "[nebula] task ${SLURM_ARRAY_TASK_ID} on $(hostname) $(date)"
$RSCRIPT R/01_run_nebula.R "${SLURM_ARRAY_TASK_ID}"
echo "[nebula] task ${SLURM_ARRAY_TASK_ID} done $(date)"
