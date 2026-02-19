#!/bin/bash
################################################################################
# run_pipeline.sh
#
# Master orchestration script.  Submits all SLURM jobs with proper dependency
# chains.  Steps run in order:
#
#   00_fit_reference  (once; ~4-12 h)
#   01_param_grid     (depends on 00; < 1 min)
#   02_simulate_array (array; depends on 01; simulate + analyze combined)
#   05_benchmark      (depends on 02; ~2 h)
#   06_plots          (depends on 05; ~1 h)
#
# NOTE: The old separate 03_nebula and 04_pseudobulk array jobs are retired.
#       All analysis now happens inside 02_simulate_analyze.R.
#
# USAGE:
#   # First, submit (or manually run) Step 00, then:
#   bash slurm/run_pipeline.sh
################################################################################

set -euo pipefail

SLURM_DIR="$(cd "$(dirname "$0")" && pwd)"
RESULTS_DIR="$(cd "${SLURM_DIR}/.." && pwd)/results"
LOG_DIR="$(cd "${SLURM_DIR}/.." && pwd)/logs"
mkdir -p "${LOG_DIR}" "${RESULTS_DIR}"

# Determine grid size (reads from file if available, else uses default)
get_grid_size() {
  if [ -f "${RESULTS_DIR}/param_grid/param_grid.rds" ]; then
    Rscript -e "
      pg <- readRDS('${RESULTS_DIR}/param_grid/param_grid.rds')
      cat(nrow(pg))
    " 2>/dev/null
  else
    echo "21600"   # default: 432 unique scenarios x 50 reps
  fi
}

echo "════════════════════════════════════════"
echo " scDesign3 Simulation Pipeline"
echo " $(date)"
echo "════════════════════════════════════════"

# ── Step 00: Fit reference model ─────────────────────────────────────────────
# Uncomment to submit automatically; otherwise run manually first:
# JOB00=$(sbatch --parsable "${SLURM_DIR}/00_fit_reference.sh")
# echo "Submitted step 00 (fit reference): job ${JOB00}"
# DEPEND_01="--dependency=afterok:${JOB00}"
DEPEND_01=""   # set to "--dependency=afterok:${JOB00}" if submitting 00 above

# ── Step 01: Build parameter grid ────────────────────────────────────────────
JOB01=$(sbatch --parsable ${DEPEND_01} "${SLURM_DIR}/01_param_grid.sh")
echo "Submitted step 01 (param_grid):    job ${JOB01}"

GRID_SIZE=$(get_grid_size)
echo "Parameter grid size:               ${GRID_SIZE} rows"

# ── Step 02: Simulate + analyze (combined array) ──────────────────────────────
JOB02=$(sbatch --parsable \
    --dependency="afterok:${JOB01}" \
    --array="1-${GRID_SIZE}%200" \
    "${SLURM_DIR}/02_simulate_array.sh")
echo "Submitted step 02 (sim+analyze):   array job ${JOB02}[1-${GRID_SIZE}]"

# ── Step 05: Benchmark ────────────────────────────────────────────────────────
JOB05=$(sbatch --parsable \
    --dependency="afterok:${JOB02}" \
    "${SLURM_DIR}/05_benchmark.sh")
echo "Submitted step 05 (benchmark):     job ${JOB05}"

# ── Step 06: Plots ────────────────────────────────────────────────────────────
JOB06=$(sbatch --parsable \
    --dependency="afterok:${JOB05}" \
    "${SLURM_DIR}/06_plots.sh")
echo "Submitted step 06 (plots):         job ${JOB06}"

echo ""
echo "════════════════════════════════════════"
echo " All jobs submitted!"
echo " Monitor:  squeue -u \$USER"
echo " Logs:     ${LOG_DIR}"
echo "════════════════════════════════════════"
