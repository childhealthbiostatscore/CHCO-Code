#!/bin/bash
# ---------------------------------------------------------------------------
# Submit the whole pipeline. The split job auto-submits the NEBULA array and
# the aggregation job, so this is the only command you need to run.
#
#   bash slurm/run_pipeline.sh
#
# Stages: 01_split  ->  02_nebula_array (auto, sized from manifest)  ->  03_aggregate
# ---------------------------------------------------------------------------
set -euo pipefail
cd "$(dirname "$0")/.."
JID=$(sbatch --parsable slurm/01_split.slurm)
echo "Submitted split job $JID. It will chain the NEBULA array + aggregation."
echo "Track with:  squeue --me   |   tail -f neb_*_${JID}.out"
