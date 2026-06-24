#!/bin/bash
#SBATCH --job-name=neb_split
#SBATCH --account=togo
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=240G
#SBATCH --time=08:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# ---------------------------------------------------------------------------
# Stage 1: split each dataset into per-cell-type objects + build the manifest,
# then auto-submit the NEBULA array (sized from the manifest) and the
# aggregation job (dependent on the array).
# ---------------------------------------------------------------------------
set -euo pipefail

# EDIT THESE TWO PATHS for your checkout / scratch location:
PROJECT_DIR="/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad/P047_scRNA_transcript_search_yc/untargeted_nebula"
WORK_DIR="/gscratch/togo/yejichoi/scrna_untargeted"

# --- R environment (adjust to your Hyak setup) -----------------------------
# module load <your R module>            # e.g. module load r/4.3.2
# or: source ~/miniconda3/bin/activate <your-env>
RSCRIPT="${RSCRIPT:-Rscript}"

cd "$PROJECT_DIR"
mkdir -p "$WORK_DIR"

echo "[split] starting $(date)"
$RSCRIPT R/00_split_celltypes.R
echo "[split] done $(date)"

# number of runs = data lines in jobs.tsv (minus header)
NJOBS=$(($(wc -l < "$WORK_DIR/jobs.tsv") - 1))
echo "[split] manifest has $NJOBS runs"

# Submit the NEBULA array in chunks so we stay under SLURM MaxArraySize, and
# throttle concurrent tasks with %THROTTLE. Each chunk gets a row offset.
MAXARR="${MAXARR:-1000}"     # max array indices per submission
THROTTLE="${THROTTLE:-60}"   # max simultaneously running array tasks per chunk

DEPS=""
START=1
while [ "$START" -le "$NJOBS" ]; do
  END=$(( START + MAXARR - 1 ))
  [ "$END" -gt "$NJOBS" ] && END="$NJOBS"
  SIZE=$(( END - START + 1 ))
  OFFSET=$(( START - 1 ))
  JID=$(sbatch --parsable --array=1-"${SIZE}"%"${THROTTLE}" \
    --export=ALL,PROJECT_DIR="$PROJECT_DIR",WORK_DIR="$WORK_DIR",NEB_ROW_OFFSET="$OFFSET" \
    slurm/02_nebula_array.slurm)
  echo "[split] submitted NEBULA chunk rows $START-$END as job $JID (offset $OFFSET)"
  DEPS="${DEPS}:${JID}"
  START=$(( END + 1 ))
done

# aggregation depends on ALL chunks finishing (afterany: still combine whatever
# partials succeeded even if some runs fail)
AGG_JID=$(sbatch --parsable --dependency=afterany"${DEPS}" \
  --export=ALL,PROJECT_DIR="$PROJECT_DIR",WORK_DIR="$WORK_DIR" \
  slurm/03_aggregate.slurm)
echo "[split] submitted aggregation job $AGG_JID (after${DEPS})"
