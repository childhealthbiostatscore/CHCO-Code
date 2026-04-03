#!/bin/bash
#SBATCH --partition=cpu-g2
#SBATCH --account=togo
################################################################################
# launch_oom_resubmit.sh
#
# Helper script to submit the OOM rerun jobs.
# Run from the project base directory on Hyak.
################################################################################

BASE_DIR="/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad/scRNA simulation analysis"

# Count lines in each ID file
N_128=$(wc -l < "mmfs1/gscratch/togo/yejichoi/project_logs/scD3_logs/02_sim/failed_arrays_33794733.txt")
N_256=$(wc -l < "mmfs1/gscratch/togo/yejichoi/project_logs/scD3_logs/02_sim/failed_arrays_33794735.txt")

echo "OOM resubmission plan:"
echo "  128GB tier: ${N_128} jobs"
echo "  256GB tier: ${N_256} jobs"
echo ""

# Submit 128GB tier (throttle to 50 concurrent)
echo "Submitting 128GB tier..."
sbatch --array=1-${N_128}%50 "${BASE_DIR}/resubmit_oom_128.sh"

# Submit 256GB tier (throttle to 50 concurrent — larger memory footprint)
echo "Submitting 256GB tier..."
sbatch --array=1-${N_256}%50 "${BASE_DIR}/resubmit_oom_256.sh"

echo ""
echo "Done. Monitor with: squeue -u \$USER"