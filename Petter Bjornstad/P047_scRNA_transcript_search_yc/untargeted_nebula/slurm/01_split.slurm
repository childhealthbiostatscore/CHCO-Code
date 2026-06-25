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

# EDIT THESE TWO PATHS for your checkout / scratch location.
# WORK_DIR must match config.R -> global$work_dir (00_split writes jobs.tsv
# there and this script reads it back to size the array).
PROJECT_DIR="/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad/P047_scRNA_transcript_search_yc/untargeted_nebula"
WORK_DIR="/mmfs1/gscratch/togo/yejichoi/project_logs/untargeted_nebula_logs"

# --- R via Apptainer container ---------------------------------------------
module load apptainer
CONTAINER="${CONTAINER:-/mmfs1/gscratch/togo/YC_scRNA.sif}"
# The container's site-library is read-only, so togolab/togolab.scrna install
# into this writable user library on /gscratch instead. It is passed in as
# NEB_RLIB (setup.R puts it first on .libPaths). PRE-INSTALL ONCE on a login
# node (has internet) so jobs only load — use an EXPLICIT lib= (don't rely on
# env vars, which the container's Renviron overrides):
#   module load apptainer
#   apptainer exec --bind /gscratch,/mmfs1 /mmfs1/gscratch/togo/YC_scRNA.sif Rscript -e '
#     lib <- "/mmfs1/gscratch/togo/yejichoi/R_libs/YC_scRNA";
#     dir.create(lib, recursive=TRUE, showWarnings=FALSE); .libPaths(c(lib,.libPaths()));
#     remotes::install_github(c("uwmdi-togo/togolab","uwmdi-togo/togolab.scrna"), lib=lib, upgrade="never")'
NEB_RLIB="${NEB_RLIB:-/mmfs1/gscratch/togo/yejichoi/R_libs/YC_scRNA}"
mkdir -p "$NEB_RLIB"
# bind the gscratch + mmfs1 trees and pass the user lib into the container.
# RSCRIPT is used unquoted below, so it expands to the full apptainer command.
RSCRIPT="apptainer exec --bind /gscratch,/mmfs1 --env NEB_RLIB=${NEB_RLIB} ${CONTAINER} Rscript"

cd "$PROJECT_DIR"
mkdir -p "$WORK_DIR"

# Set SKIP_SPLIT=1 to reuse the existing splits + jobs.tsv (e.g. after fixing a
# downstream bug) and jump straight to submitting the array.
if [ "${SKIP_SPLIT:-0}" = "1" ]; then
  echo "[split] SKIP_SPLIT=1 -> reusing existing $WORK_DIR/jobs.tsv"
else
  echo "[split] starting $(date)"
  $RSCRIPT R/00_split_celltypes.R
  echo "[split] done $(date)"
fi

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
