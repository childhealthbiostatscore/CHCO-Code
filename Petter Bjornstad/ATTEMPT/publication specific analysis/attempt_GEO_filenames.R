library(digest)

# Define the main folder path
user <- Sys.info()[["user"]]

if (user == "choiyej") {
  root_path <- "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/"
  git_path <- "/Users/choiyej/GitHub/CHCO-Code/Petter Bjornstad/"
} else if (user == "pylell") {
  root_path <- "/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive/"
  git_path <- "/Users/pylell/Documents/GitHub/CHCO-Code/Petter Bjornstad/"
} else {
  stop("Unknown user: please specify root path for this user.")
}

## Create an S3 client
keys <- fromJSON("/Users/choiyej/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Bjornstad Pyle Lab/keys.json")

Sys.setenv(
  "AWS_ACCESS_KEY_ID" = keys$MY_ACCESS_KEY,
  "AWS_SECRET_ACCESS_KEY" = keys$MY_SECRET_KEY,
  "AWS_DEFAULT_REGION" = "",
  "AWS_REGION" = "",
  "AWS_S3_ENDPOINT" = "s3.kopah.uw.edu"
)

########################################################################
############################## RDS -> TSV ##############################
########################################################################

# Load the combined RDS file
# rds_file <- "/Volumes/Peds Endo/Petter Bjornstad/scRNA/data_clean/seurat_data_CRC.RDS"
# seurat_object <-  s3readRDS(object = "cleaned_data/attempt_clean_so.rds", bucket = "attempt", region = "")

# Load necessary libraries
library(Seurat)
library(Matrix)
library(digest)

# Ensure the output directory exists
output_dir <- file.path(root_path, "ATTEMPT/Data Raw/ATTEMPT_ProcessedData")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Extract participant names (assuming participant IDs are stored in the 'orig.ident' metadata column)
participants <- unique(seurat_object@meta.data$orig.ident)

# Initialize a data frame to store participant IDs, file names, and checksums
file_info <- data.frame(
  Participant_ID = character(),
  MTX_File = character(),
  MTX_MD5 = character(),
  Barcodes_File = character(),
  Barcodes_MD5 = character(),
  Features_File = character(),
  Features_MD5 = character(),
  stringsAsFactors = FALSE
)

# Function to save data in MTX format and barcodes/features in TSV format, and compute MD5 checksums
save_data <- function(counts_matrix, barcodes, features, participant_id, output_dir) {
  # Replace hyphens with underscores in participant ID
  participant_id_fixed <- gsub("-", "_", participant_id)
  
  # Create a new folder for this participant
  participant_folder <- file.path(output_dir, participant_id_fixed)
  dir.create(participant_folder, showWarnings = FALSE)
  
  # Set up file paths with "_processed" suffix
  mtx_path <- file.path(participant_folder, paste0(participant_id_fixed, "_matrix_processed.mtx"))
  barcodes_path <- file.path(participant_folder, paste0(participant_id_fixed, "_barcodes_processed.tsv"))
  features_path <- file.path(participant_folder, paste0(participant_id_fixed, "_features_processed.tsv"))
  
  # Save counts matrix as MTX (Matrix Market format)
  Matrix::writeMM(counts_matrix, file = mtx_path)
  
  # Save barcodes (cell names)
  write.table(barcodes, file = barcodes_path, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # Save features (gene names)
  write.table(features, file = features_path, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # Calculate MD5 checksums
  mtx_md5 <- digest(mtx_path, algo = "md5", file = TRUE)
  barcodes_md5 <- digest(barcodes_path, algo = "md5", file = TRUE)
  features_md5 <- digest(features_path, algo = "md5", file = TRUE)
  
  # Return the file names (without paths) and MD5 checksums for recording
  return(list(
    mtx_file = basename(mtx_path), mtx_md5 = mtx_md5,
    barcodes_file = basename(barcodes_path), barcodes_md5 = barcodes_md5,
    features_file = basename(features_path), features_md5 = features_md5
  ))
}

# Loop through each participant and save their data
for (participant_id in participants) {
  # Subset the Seurat object for the specific participant
  participant_subset <- subset(seurat_object, subset = orig.ident == participant_id)
  
  # Extract the counts matrix for the participant
  counts_matrix <- GetAssayData(participant_subset, slot = "counts")
  
  # Extract barcodes (cell names) and features (gene names)
  barcodes <- colnames(counts_matrix)
  features <- rownames(counts_matrix)
  
  # Save the counts matrix, barcodes, and features, and calculate MD5 checksums
  file_paths <- save_data(counts_matrix, barcodes, features, participant_id, output_dir)
  
  # Add participant ID, file names, and MD5 checksums to the data frame
  file_info <- rbind(file_info, data.frame(
    Participant_ID = gsub("-", "_", participant_id),
    MTX_File = file_paths$mtx_file,
    MTX_MD5 = file_paths$mtx_md5,
    Barcodes_File = file_paths$barcodes_file,
    Barcodes_MD5 = file_paths$barcodes_md5,
    Features_File = file_paths$features_file,
    Features_MD5 = file_paths$features_md5,
    stringsAsFactors = FALSE
  ))
}

# Save the file info data frame to a CSV file, including the MD5 checksums
write.csv(file_info, file = file.path(output_dir, "processed_participant_file_info_with_md5.csv"), row.names = FALSE)

# Optional: Print a message when done
cat("All participant data, barcodes, and features saved to MTX and TSV formats, and MD5 checksums compiled into CSV.")

#####################################################
# Pull meta data
#####################################################
md5_csv <- read.csv(file.path(output_dir, "processed_participant_file_info_with_md5.csv"))
seurat_object_meta <- seurat_object@meta.data

md5_csv_meta <- seurat_object_meta %>%
  mutate(Participant_ID = gsub("-", "_", orig.ident)) %>%
  dplyr::select(Participant_ID, visit, treatment) %>%
  right_join(md5_csv) %>%
  distinct(Participant_ID, .keep_all = T)

write.csv(md5_csv_meta, file = file.path(output_dir, "processed_participant_file_info_with_md5_meta.csv"), row.names = FALSE)
