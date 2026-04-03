#!/bin/bash
# =============================================================================
# Setup and submit T1D Adiposity NEBULA jobs
# =============================================================================
# Run this script on the cluster to set up directories, generate config,
# and submit the SLURM array.
# =============================================================================

BASE_DIR="/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad/T1D Adiposity/nebula"

echo "Setting up T1D Adiposity NEBULA pipeline..."
echo "Base directory: ${BASE_DIR}"

# Create directory structure
mkdir -p "${BASE_DIR}/scripts"
mkdir -p "${BASE_DIR}/config"
mkdir -p "/mmfs1/gscratch/togo/yejichoi/project_logs/t1d_adiposity_logs/output"
mkdir -p "/mmfs1/gscratch/togo/yejichoi/project_logs/t1d_adiposity_logs/error"

echo "Directory structure created."

# Copy scripts (assumes files are in current directory)
cp run_nebula_single.R "${BASE_DIR}/scripts/"
cp generate_job_config.R "${BASE_DIR}/scripts/"
cp nebula_array.slurm "${BASE_DIR}/"
cp save_celltype_subsets.R "${BASE_DIR}/scripts/"

echo "Scripts copied."

# Step 1: Save cell type subsets (run this first, takes ~30min)
echo ""
echo "=== STEP 1: Save cell type subsets ==="
echo "This saves per-celltype Seurat objects to S3 so array jobs load faster."
echo "Run this interactively or via a single SLURM job:"
echo ""
echo "  # Option A: Interactive (on a compute node with enough memory)"
echo "  srun --mem=250G --cpus-per-task=4 --time=02:00:00 --partition=cpu-g2 --account=togo --pty bash"
echo "  module load apptainer"
echo "  apptainer exec --bind /mmfs1 /mmfs1/gscratch/togo/YC_scRNA.sif Rscript \"${BASE_DIR}/scripts/save_celltype_subsets.R\""
echo ""
echo "  # Option B: Single SLURM job"
echo "  sbatch --mem=250G --cpus-per-task=4 --time=02:00:00 --partition=cpu-g2 --account=togo \\"
echo "    --wrap='module load apptainer && apptainer exec --bind /mmfs1 /mmfs1/gscratch/togo/YC_scRNA.sif Rscript \"${BASE_DIR}/scripts/save_celltype_subsets.R\"'"
echo ""

# Step 2: Generate job config
echo "=== STEP 2: Generate job config ==="
echo "Run from the base directory:"
echo ""
echo "  cd \"${BASE_DIR}\""
echo "  module load apptainer"
echo "  apptainer exec --bind /mmfs1 /mmfs1/gscratch/togo/YC_scRNA.sif Rscript \"${BASE_DIR}/scripts/generate_job_config.R\""
echo ""

# Step 3: Submit array
echo "=== STEP 3: Submit SLURM array ==="
echo "After the config is generated, check the number of lines:"
echo ""
echo "  wc -l \"${BASE_DIR}/config/job_config.txt\""
echo ""
echo "Then update --array in nebula_array.slurm if needed, and submit:"
echo ""
echo "  cd \"${BASE_DIR}\""
echo "  sbatch nebula_array.slurm"
echo ""

echo "=== MONITORING ==="
echo "  squeue -u yejichoi                     # Check running jobs"
echo "  sacct -j <JOBID> --format=JobID,State   # Check completion"
echo "  ls ${BASE_DIR}/logs/error/ | head        # Check for errors"
echo ""
echo "  # Find failed jobs:"
echo "  grep -l 'ERROR\\|SKIPPING' \"${BASE_DIR}/logs/output/\"*.out"
echo ""
echo "  # Count completed jobs:"
echo "  grep -l 'COMPLETE' \"${BASE_DIR}/logs/output/\"*.out | wc -l"
