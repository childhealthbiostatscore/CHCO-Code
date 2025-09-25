# GTF to BED Converter Function
# Converts GENCODE GTF.gz files to BED format for specified feature types

suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(dplyr))

#' Convert GTF to BED format for specific feature types
#' 
#' @param gtf_url URL to the GTF.gz file from GENCODE
#' @param feature_type Character vector of feature types to extract (e.g., "gene", "transcript", "start_codon")
#' @param output_dir Directory to save BED files (default: current directory)
#' @param download_first Logical, whether to download GTF file first (default: TRUE)
#' 
#' @return List of data frames in BED format, one for each feature type
#' @export

gtf_to_bed_converter <- function(gtf_url, 
                                 feature_type = c("gene", "transcript", "start_codon"),
                                 output_dir = ".",
                                 download_first = TRUE) {
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Download GTF file if requested
  if (download_first) {
    cat("Downloading GTF file...\n")
    local_gtf <- file.path(output_dir, "gencode.gtf.gz")
    download.file(gtf_url, local_gtf, method = "auto")
    gtf_file <- local_gtf
  } else {
    gtf_file <- gtf_url
  }
  
  # Read GTF file
  cat("Reading GTF file...\n")
  gtf_data <- import(gtf_file, format = "gtf")
  
  # Initialize results list
  bed_results <- list()
  
  # Process each feature type
  for (feat in feature_type) {
    cat(paste("Processing feature type:", feat, "\n"))
    
    # Filter for specific feature type
    feature_data <- gtf_data[gtf_data$type == feat]
    
    if (length(feature_data) == 0) {
      warning(paste("No features found for type:", feat))
      next
    }
    
    my_names = feature_data$transcript_id
    
    if(feat == "gene"){my_names = feature_data$gene_id}
    if(feat == "exon"){my_names = feature_data$exon_id}
    
    # Convert to BED format (0-based coordinates)
    bed_df <- data.frame(
      chrom = as.character(seqnames(feature_data)),
      chromStart = start(feature_data) - 1,  # Convert to 0-based
      chromEnd = end(feature_data),
      name = my_names,
      score = ifelse(!is.na(feature_data$score), feature_data$score, 0),
      strand = as.character(strand(feature_data))
    )

    # Store result
    bed_results[[feat]] <- bed_df
    
    # Save to file
    output_file <- file.path(output_dir, paste0("gencode_", feat, ".bed"))
    write.table(bed_df, output_file, 
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    cat(paste("Saved", nrow(bed_df), "features to", output_file, "\n"))
  }
  
  cat("Conversion completed!\n")
  return(bed_results)
}

# Alternative function for custom GTF parsing (if rtracklayer has issues)
gtf_to_bed_manual <- function(gtf_url, feature_type = c("gene", "transcript", "start_codon"), 
                              output_dir = ".") {
  
  library(readr)
  library(stringr)
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Download and read GTF
  cat("Downloading and reading GTF file...\n")
  local_gtf <- file.path(output_dir, "temp_gencode.gtf.gz")
  download.file(gtf_url, local_gtf)
  
  # Read GTF file
  gtf_lines <- read_lines(local_gtf)
  gtf_lines <- gtf_lines[!grepl("^#", gtf_lines)]  # Remove comment lines
  
  bed_results <- list()
  
  for (feat in feature_type) {
    cat(paste("Processing feature type:", feat, "\n"))
    
    # Filter lines for feature type
    feat_lines <- gtf_lines[grepl(paste0("\t", feat, "\t"), gtf_lines)]
    
    if (length(feat_lines) == 0) {
      warning(paste("No features found for type:", feat))
      next
    }
    
    # Parse GTF lines
    parsed_data <- do.call(rbind, lapply(feat_lines, function(line) {
      fields <- str_split(line, "\t")[[1]]
      
      # Extract attributes
      attributes <- fields[9]
      gene_name <- str_extract(attributes, 'gene_name "([^"]+)"')
      gene_name <- ifelse(is.na(gene_name), "", str_replace(gene_name, 'gene_name "([^"]+)"', "\\1"))
      
      gene_id <- str_extract(attributes, 'gene_id "([^"]+)"')
      gene_id <- ifelse(is.na(gene_id), "", str_replace(gene_id, 'gene_id "([^"]+)"', "\\1"))
      
      # Extract appropriate ID based on feature type
      transcript_id <- str_extract(attributes, 'transcript_id "([^"]+)"')
      transcript_id <- ifelse(is.na(transcript_id), "", str_replace(transcript_id, 'transcript_id "([^"]+)"', "\\1"))
      
      feature_id <- switch(feat,
                           "gene" = gene_id,
                           "transcript" = transcript_id,
                           "start_codon" = transcript_id,
                           "stop_codon" = transcript_id,
                           "CDS" = transcript_id
      )
      
      data.frame(
        chrom = fields[1],
        chromStart = as.numeric(fields[4]) - 1,  # Convert to 0-based
        chromEnd = as.numeric(fields[5]),
        name = feature_id
      )
    }))
    
    bed_results[[feat]] <- parsed_data
    
    # Save to file
    output_file <- file.path(output_dir, paste0("gencode_", feat, "_manual.bed"))
    write.table(parsed_data, output_file, 
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    cat(paste("Saved", nrow(parsed_data), "features to", output_file, "\n"))
  }
  
  # Clean up temporary file
  unlink(local_gtf)
  
  return(bed_results)
}

# Example usage:
# URL for GENCODE human GTF (update version as needed)
# gencode_url <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz"
# 
# # Convert specific feature types
# bed_files <- gtf_to_bed_converter(
#   gtf_url = gencode_url,
#   feature_type = c("gene", "transcript", "start_codon"),
#   output_dir = "./bed_output"
# )
# 
# # Access results
# genes_bed <- bed_files$gene
# transcripts_bed <- bed_files$transcript
# start_codons_bed <- bed_files$start_codon