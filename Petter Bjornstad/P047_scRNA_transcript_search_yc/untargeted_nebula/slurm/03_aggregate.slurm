#!/bin/bash
#SBATCH --job-name=neb_agg
#SBATCH --account=togo
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# ---------------------------------------------------------------------------
# Stage 3: combine partials into one CSV per contrast.
# ---------------------------------------------------------------------------
set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad/P047_scRNA_transcript_search_yc/untargeted_nebula}"

# module load <your R module>
RSCRIPT="${RSCRIPT:-Rscript}"

cd "$PROJECT_DIR"
echo "[aggregate] starting $(date)"
$RSCRIPT R/02_aggregate.R
echo "[aggregate] done $(date)"
