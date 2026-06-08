#!/bin/bash
#SBATCH --job-name=nebula
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/nebula_%j.out
#SBATCH --error=logs/nebula_%j.err

mkdir -p logs

module purge
module load apptainer

CONTAINER=/mmfs1/gscratch/togo/yzhangtufts_r_scrnaseq.sif

apptainer exec \
  --cleanenv \
  --bind /mmfs1:/mmfs1 \
  --env CELLTYPE=${CELLTYPE} \
  --env CONTRAST=${CONTRAST} \
  ${CONTAINER} \
  Rscript run_nebula.R ${CELLTYPE} ${CONTRAST}