#!/bin/bash

#SBATCH --job-name=Kidney_scRNA
#SBATCH --mail-type=END
#SBATCH --mail-user=tvigers@uw.edu
#SBATCH --account=togo
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem=200G
#SBATCH --time=01-00:00:00 # Max runtime in DD-HH:MM:SS format.
#SBATCH --export=all

# Your programs to run.
cd /gscratch/togo/tvigers
s3cmd sync "s3://uwmdi/Hailey Hampson" ./Data/UWMDI
apptainer exec --bind /gscratch ubuntu.sif Rscript GitHub/CHCO-Code/Petter\ Bjornstad/Kidney\ scRNA/Kidney\ scRNA/neg_binom_mixed_model.R