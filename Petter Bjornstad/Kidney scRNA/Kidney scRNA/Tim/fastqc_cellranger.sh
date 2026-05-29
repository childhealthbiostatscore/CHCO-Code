#!/bin/bash
#-------------------------------------------------------------------------------
# QC
#-------------------------------------------------------------------------------
# Define input and output directories
alias s3togo='s3cmd -c ~/.s3cfg_togo'
INPUT_DIR="/home/tvigers/Documents/Data/UWMDI/kidney_organoids/data_raw/organoids"
OUTPUT_DIR="/home/tvigers/Documents/Data/UWMDI/kidney_organoids/data_clean/organoids"
S3_BUCKET="s3://scrna/Kidney organoids"
# Sync from S3 - faster if you only run once
s3togo() { s3cmd -c ~/.s3cfg_togo "$@"; }
export -f s3togo
mkdir -p "${INPUT_DIR}"
s3togo sync "${S3_BUCKET}/" "${INPUT_DIR}/" -v --check-md5 --exclude="data_clean/*"
# Create the output directory if it doesn't already exist
mkdir -p "${OUTPUT_DIR}/qc"
# Run FastQC
fastqc -t 4 -o "${OUTPUT_DIR}/qc" "${INPUT_DIR}"/*.fastq-*.gz
echo "FastQC analysis complete!"
# Run MultiQC
echo "Starting MultiQC..."
multiqc "${OUTPUT_DIR}/qc" -o "${OUTPUT_DIR}/qc" --no-data-dir
# Sync to S3
echo "Syncing results to S3..."
s3togo sync "${OUTPUT_DIR}/qc" "$S3_BUCKET/data_clean/qc/" -v --check-md5 --delete-removed
echo "QC complete!"
#-------------------------------------------------------------------------------
# Cell Ranger
#-------------------------------------------------------------------------------
# Define directories and parameters
TRANSCRIPTOME="/home/tvigers/Documents/References/refdata-gex-GRCh38-2024-A"
# Annoyingly, these FASTQ files have a weird number suffix (e.g. .fastq-030.gz).
# Copy and rename them into a temporary directory
TEMP_FASTQ_DIR="/home/tvigers/Documents/Data/UWMDI/kidney_organoids/data_clean/temp_fastqs"
mkdir -p "$TEMP_FASTQ_DIR"
for filepath in "$INPUT_DIR"/*.fastq-*.gz; do
    # Check if files actually exist to avoid errors if the folder is empty
    if [[ ! -e "$filepath" ]]; then
        echo "No files matching *.fastq-*.gz found in $INPUT_DIR"
        break
    fi
    # Extract just the filename (e.g., drops the /home/.../ folder path)
    filename=$(basename "$filepath")
    # Strip the weird ending and append the standard .fastq.gz
    new_name="${filename%.fastq*}.fastq.gz"
    # Copy the file to the temp directory with the new name
    cp "$filepath" "$TEMP_FASTQ_DIR/$new_name"
done
# Make output folder
mkdir -p "${OUTPUT_DIR}/cellranger"
cd "${OUTPUT_DIR}/cellranger" || exit
# Get sample names
SAMPLES=$(ls -1 "$TEMP_FASTQ_DIR"/*.fastq.gz | awk -F'/' '{print $NF}' | sed 's/_S[0-9]\+.*//' | sort -u)
echo "Found the following samples to process:"
echo "$SAMPLES"
# Loop through samples and run pipeline
for SAMPLE in $SAMPLES; do
    echo "Starting cellranger for sample: $SAMPLE"
    # Run the count pipeline
    # The --id flag dictates the name of the output folder for this sample
    cellranger count \
        --id="${SAMPLE}" \
        --transcriptome="$TRANSCRIPTOME" \
        --fastqs="$TEMP_FASTQ_DIR" \
        --sample="${SAMPLE}" \
        --localcores=16 \
        --localmem=64 \
        --create-bam false \
        --nosecondary \
        --disable-ui
    # Send results to S3
      s3togo sync --follow-symlinks --check-md5 \
        "${OUTPUT_DIR}/cellranger/${SAMPLE}/outs/" \
        "${S3_BUCKET}/data_clean/organoids/samples/${SAMPLE}/cellranger/outs/"
    echo "Finished processing $SAMPLE"
done
echo "All samples processed successfully!"
# Delete temp files
rm -rf "${TEMP_FASTQ_DIR}"