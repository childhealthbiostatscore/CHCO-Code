#!/bin/bash
#SBATCH --job-name=nebula_array
#SBATCH --output=logs/nebula_%A_%a.out
#SBATCH --error=logs/nebula_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=20
#SBATCH --partition=compute
#SBATCH --array=0-111%20  # (8 analysis x 14 cell types), max 20 running at once

# Load singularity/apptainer module
module load apptainer

# Path to Singularity container
CONTAINER="/mmfs1/gscratch/togo/yzhangtufts_r_scrnaseq.sif"

CELL_TYPES=("PT" "TAL" "PC" "IC" "DTL_ATL" "DCT_CNT" "EC" "Immune" "VSMC_P_FIB" "POD" "MC" "PEC" "Schwann" "Other")

# Analysis types from script
ANALYSIS_TYPES=("dkd100_vs_nondkd" "dkd100_vs_hc" "nondkd100_vs_hc" 
                "dkd30_vs_nondkd" "dkd30_vs_hc" "nondkd30_vs_hc"
                "glpn_vs_hc" "glpy_vs_glpn")

# Calculate which analysis and cell type based on array task ID
N_CELLS=${#CELL_TYPES[@]}
N_ANALYSES=${#ANALYSIS_TYPES[@]}

ANALYSIS_IDX=$((SLURM_ARRAY_TASK_ID / N_CELLS))
CELL_IDX=$((SLURM_ARRAY_TASK_ID % N_CELLS))

ANALYSIS=${ANALYSIS_TYPES[$ANALYSIS_IDX]}
CELL=${CELL_TYPES[$CELL_IDX]}

echo "=========================================="
echo "Job Array ID: $SLURM_ARRAY_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID (out of 112 total)"
echo "Analysis type: $ANALYSIS"
echo "Cell type: $CELL"
echo "Start time: $(date)"
echo "Node: $(hostname)"
echo "=========================================="

# Change to script directory
cd "/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad/Renal HEIRitage/"

# Run the parallelized R script
Rscript 04_scRNA_nebula.R ${ANALYSIS} ${CELL}

echo "=========================================="
echo "End time: $(date)"
echo "Exit code: $?"
echo "Job completed for $ANALYSIS - $CELL"
echo "=========================================="
