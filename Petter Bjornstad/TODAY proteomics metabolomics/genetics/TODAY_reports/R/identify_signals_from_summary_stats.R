library(data.table)

#' Identify unique GWAS signals and lead variants
#' 
#' @param rds_file Path to RDS file containing GWAS summary statistics
#' @param p_threshold P-value threshold for significance (default: 5e-8)
#' @param clump_distance Distance in base pairs for clumping variants (default: 1MB)
#' @param min_p_diff Minimum p-value difference to consider variants independent (default: 100)
#' 
#' @return Data frame with columns: CHROM, POS, ID, REF, ALT, P, SIGNAL_ID
#' 
#' @examples
#' signals <- get_gwas_signals("gwas_results.rds")
#' signals <- get_gwas_signals("gwas_results.rds", p_threshold = 1e-6, clump_distance = 500000)
get_gwas_signals <- function(rds_file, 
                             p_threshold = 5e-8, 
                             clump_distance = 1000000,
                             min_p_diff = 100) {
  
  # Load and validate data
  if (!file.exists(rds_file)) {
    stop("RDS file not found: ", rds_file)
  }
  
  gwas_data <- readRDS(rds_file)
  
  # Check required columns
  required_cols <- c("CHROM", "POS", "ID", "REF", "ALT", "P")
  missing_cols <- setdiff(required_cols, names(gwas_data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Convert to data.table for efficient processing
  dt0 <- as.data.table(gwas_data %>% filter(is.na(P) == FALSE))
  
  # Filter for significant variants and remove missing p-values
  dt <- dt0[P <= p_threshold]
  
  if (nrow(dt) == 0) {
    return(head(as.data.frame(dt0) %>% 
                  filter(P <= min(P)) %>% 
                  mutate(SIGNAL_ID = "SIGNAL_000", N_VARIANTS = 0, START = POS, END = POS), n = 1))
  }
  
  # Ensure chromosome is properly formatted
  dt[, CHROM := as.character(CHROM)]
  dt[, POS := as.numeric(POS)]
  dt[, P := as.numeric(P)]
  
  # Sort by chromosome and position
  setorder(dt, CHROM, POS, P)
  
  # Initialize results
  signals <- data.table()
  signal_counter <- 1
  
  # Process each chromosome separately
  chroms <- unique(dt$CHROM)
  
  for (chrom in chroms) {
    chrom_data <- dt[CHROM == chrom]
    
    # Continue until all variants are processed
    while (nrow(chrom_data) > 0) {
      
      # Find the most significant variant (lead variant)
      lead_idx <- which.min(chrom_data$P)
      lead_variant <- chrom_data[lead_idx]
      
      # Define the clumping window around the lead variant
      window_start <- lead_variant$POS - clump_distance
      window_end <- lead_variant$POS + clump_distance
      
      # Find all variants in the clumping window
      clumped_variants <- chrom_data[POS >= window_start & POS <= window_end]
      
      # Check if there are other very significant variants that might be independent
      # This helps identify cases where there might be multiple independent signals
      other_very_sig <- clumped_variants[P < lead_variant$P / min_p_diff & 
                                           abs(POS - lead_variant$POS) > clump_distance/4]
      
      if (nrow(other_very_sig) > 0) {
        # If there are other very significant variants, use a smaller clumping window
        window_start <- lead_variant$POS - clump_distance/2
        window_end <- lead_variant$POS + clump_distance/2
        clumped_variants <- chrom_data[POS >= window_start & POS <= window_end]
      }
      
      # Add signal ID and store the lead variant
      lead_variant[, `:=`(SIGNAL_ID = paste0("SIGNAL_", sprintf("%03d", signal_counter)), 
                          N_VARIANTS = nrow(clumped_variants),
                          START = window_start, 
                          END = window_end)] 
      
      signals <- rbind(signals, lead_variant, fill = TRUE)
      
      # Remove all clumped variants from further consideration
      chrom_data <- chrom_data[!(POS >= window_start & POS <= window_end)]
      
      signal_counter <- signal_counter + 1
    }
  }
  
  # Sort results by chromosome and position
  setorder(signals, CHROM, POS)
  
  # Convert back to regular data frame
  result <- as.data.frame(signals)
  
  return(result)
}

#' Helper function to summarize GWAS signals
#' 
#' @param signals Output from get_gwas_signals()
#' 
#' @return Summary statistics
summarize_signals <- function(signals) {
  if (nrow(signals) == 0) {
    return("No signals found")
  }
  
  summary_stats <- list(
    total_signals = nrow(signals),
    chromosomes = length(unique(signals$CHROM)),
    min_p_value = min(signals$P),
    signals_per_chrom = table(signals$CHROM),
    top_signals = signals[order(signals$P)[1:min(5, nrow(signals))], 
                          c("CHROM", "POS", "ID", "P", "SIGNAL_ID")]
  )
  
  return(summary_stats)
}

#' Export signals to various formats
#' 
#' @param signals Output from get_gwas_signals()
#' @param output_prefix Prefix for output files
#' 
export_signals <- function(signals, output_prefix = "gwas_signals") {
  
  # Save as RDS
  saveRDS(signals, paste0(output_prefix, ".rds"))
  
  # Save as tab-delimited text
  write.table(signals, paste0(output_prefix, ".txt"), 
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Save BED format for genomic coordinates
  bed_format <- data.frame(
    chrom = paste0("chr", signals$CHROM),
    chromStart = signals$POS - 1,  # BED is 0-based
    chromEnd = signals$POS,
    name = signals$ID,
    score = -log10(signals$P),
    strand = "."
  )
  
  write.table(bed_format, paste0(output_prefix, ".bed"), 
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  cat("Results exported to:", 
      paste0(output_prefix, c(".rds", ".txt", ".bed")), "\n")
}

# Example usage:
# signals <- get_gwas_signals("my_gwas_results.rds")
# summary <- summarize_signals(signals)
# export_signals(signals, "my_gwas_signals")