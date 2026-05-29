#!/bin/bash
# Define input and output directories
INPUT_DIR="/home/tvigers/Documents/Data/UWMDI/kidney_organoids/data_raw/organoids"
OUTPUT_DIR="/home/tvigers/Documents/Data/UWMDI/kidney_organoids/data_clean/organoids/qc"
S3_BUCKET="s3://scrna/Kidney organoids"
# Create the output directory if it doesn't already exist
mkdir -p "$OUTPUT_DIR"
# Run FastQC
fastqc -t 4 -o "$OUTPUT_DIR" "$INPUT_DIR"/*.fastq-*.gz
echo "FastQC analysis complete!"
# Run MultiQC
echo "Starting MultiQC..."
multiqc "$OUTPUT_DIR" -o "$OUTPUT_DIR" --no-data-dir
# Sync to S3
echo "Syncing results to S3..."
s3togo sync "$OUTPUT_DIR/" "$S3_BUCKET/QC/" -v --check-md5 --delete-removed
echo "QC complete!"