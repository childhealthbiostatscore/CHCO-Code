#!/bin/bash

# Create directory structure
mkdir -p nebula_parallel/{scripts,slurm_scripts,logs/{output,error},config}

# Check if config file exists
if [ ! -f "config/job_config.txt" ]; then
echo "Error: config/job_config.txt not found!"
exit 1
fi

# Count total jobs
TOTAL_JOBS=$(wc -l < config/job_config.txt)
echo "Submitting ${TOTAL_JOBS} jobs..."

# Submit the array job
sbatch slurm_scripts/nebula_array.slurm

# Monitor script (optional)
echo ""
echo "To monitor your jobs, use:"
echo "  squeue -u $USER"
echo "  sacct -j <job_id> --format=JobID,JobName,State,ExitCode,Elapsed"