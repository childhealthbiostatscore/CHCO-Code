# R GATK CRAM to VCF Pipeline for Large Cohorts (>100 samples)
# Optimized for GenomicsDB and parallel processing

# Load required packages
library(R6)
library(parallel)
library(data.table)
library(stringr)
library(readr)
library(dplyr)

# Define null coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x

#' GATK CRAM to VCF Pipeline Class
#' 
#' @description
#' R wrapper for GATK-based variant calling pipeline optimized for large cohorts
#' 
#' @examples
#' pipeline <- GATKCramPipeline$new(
#'   cram_dir = "/path/to/cram/files",
#'   output_dir = "/path/to/output", 
#'   reference_genome = "/path/to/reference.fa"
#' )
#' pipeline$run_pipeline()
GATKCramPipeline <- R6Class("GATKCramPipeline",
                            
                            public = list(
                              
                              # Configuration
                              cram_dir = NULL,
                              output_dir = NULL,
                              reference_genome = NULL,
                              gatk_path = "gatk",
                              samtools_path = "samtools",
                              bcftools_path = "bcftools",
                              
                              # Parameters
                              threads = NULL,
                              memory = "8g",
                              batch_size = 50,  # For GenomicsDB batching
                              intervals_file = NULL,
                              
                              # Quality thresholds
                              min_base_quality = 20,
                              min_mapping_quality = 20,
                              min_depth = 10,
                              max_depth = 500,
                              
                              # Internal
                              sample_info = NULL,
                              log_file = NULL,
                              
                              #' Initialize pipeline
                              #' 
                              #' @param cram_dir Directory containing CRAM files
                              #' @param output_dir Output directory for results
                              #' @param reference_genome Path to reference genome FASTA
                              #' @param gatk_path Path to GATK executable
                              #' @param threads Number of threads (default: detectCores())
                              #' @param memory Memory allocation for GATK (default: "8g")
                              #' @param batch_size Batch size for GenomicsDB (default: 50)
                              #' @param auto_create_intervals Automatically create intervals file (default: TRUE)
                              initialize = function(cram_dir, output_dir, reference_genome, 
                                                    gatk_path = "gatk", threads = NULL, 
                                                    memory = "8g", batch_size = 50, auto_create_intervals = TRUE) {
                                
                                self$cram_dir <- normalizePath(cram_dir, mustWork = TRUE)
                                self$output_dir <- normalizePath(output_dir, mustWork = FALSE)
                                self$reference_genome <- normalizePath(reference_genome, mustWork = TRUE)
                                self$gatk_path <- gatk_path
                                self$threads <- threads %||% parallel::detectCores()
                                self$memory <- memory
                                self$batch_size <- batch_size
                                
                                # Create output directories
                                private$create_directories()
                                
                                # Setup logging
                                self$log_file <- file.path(self$output_dir, "pipeline.log")
                                private$log(paste("GATK CRAM Pipeline initialized at", Sys.time()))
                                private$log(paste("CRAM directory:", self$cram_dir))
                                private$log(paste("Output directory:", self$output_dir))
                                private$log(paste("Reference genome:", self$reference_genome))
                                private$log(paste("Threads:", self$threads))
                                private$log(paste("Memory:", self$memory))
                                private$log(paste("Batch size:", self$batch_size))
                                
                                # Create intervals file if needed
                                if (auto_create_intervals) {
                                  private$log("Auto-creating intervals file...")
                                  self$create_intervals_file(chromosome_only = TRUE)
                                }
                              },
                              
                              #' Create intervals file from reference genome
                              #' 
                              #' @param chromosome_only If TRUE, only include main chromosomes
                              create_intervals_file = function(chromosome_only = TRUE) {
                                
                                private$log("Creating intervals file from reference genome...")
                                
                                intervals_file <- file.path(self$output_dir, "intervals.list")
                                
                                # Try to get chromosomes from sequence dictionary first
                                dict_file <- sub("\\.(fasta|fa)$", ".dict", self$reference_genome)
                                
                                if (file.exists(dict_file)) {
                                  
                                  tryCatch({
                                    dict_lines <- readLines(dict_file)
                                    seq_lines <- dict_lines[grepl("^@SQ", dict_lines)]
                                    
                                    if (length(seq_lines) > 0) {
                                      # Extract sequence names and lengths
                                      seq_names <- gsub(".*SN:([^\\s]+).*", "\\1", seq_lines)
                                      seq_lengths <- as.numeric(gsub(".*LN:([0-9]+).*", "\\1", seq_lines))
                                      
                                      if (chromosome_only) {
                                        # Filter for main chromosomes only
                                        main_chrom_idx <- grepl("^(chr)?([0-9]+|X|Y|MT|M)$", seq_names)
                                        seq_names <- seq_names[main_chrom_idx]
                                        seq_lengths <- seq_lengths[main_chrom_idx]
                                      }
                                      
                                      # Sort chromosomes properly
                                      if (any(grepl("^chr", seq_names))) {
                                        # Has "chr" prefix
                                        chr_nums <- ifelse(grepl("^chr[0-9]+$", seq_names), 
                                                           as.numeric(gsub("chr", "", seq_names)), 
                                                           99)
                                        chr_order <- order(chr_nums)
                                      } else {
                                        # No "chr" prefix
                                        chr_nums <- ifelse(grepl("^[0-9]+$", seq_names),
                                                           as.numeric(seq_names),
                                                           99)
                                        chr_order <- order(chr_nums)
                                      }
                                      
                                      seq_names <- seq_names[chr_order]
                                      
                                      # Write intervals file
                                      writeLines(seq_names, intervals_file)
                                      
                                      private$log(paste("Created intervals file with", length(seq_names), "chromosomes"))
                                      private$log(paste("Intervals file:", intervals_file))
                                      
                                      self$intervals_file <- intervals_file
                                      return(intervals_file)
                                    }
                                  }, error = function(e) {
                                    private$log(paste("Error reading sequence dictionary:", e$message))
                                  })
                                }
                                
                                # Fallback: create standard human chromosomes
                                private$log("Using standard human chromosome intervals")
                                
                                # Check if reference uses "chr" prefix by looking at fasta index
                                fai_file <- paste0(self$reference_genome, ".fai")
                                has_chr_prefix <- FALSE
                                
                                if (file.exists(fai_file)) {
                                  fai_first_line <- readLines(fai_file, n = 1)
                                  has_chr_prefix <- grepl("^chr", fai_first_line)
                                }
                                
                                if (has_chr_prefix) {
                                  chromosomes <- c(paste0("chr", 1:22), "chrX", "chrY")
                                } else {
                                  chromosomes <- c(1:22, "X", "Y")
                                }
                                
                                writeLines(chromosomes, intervals_file)
                                
                                private$log(paste("Created standard intervals file with", length(chromosomes), "chromosomes"))
                                
                                self$intervals_file <- intervals_file
                                return(intervals_file)
                              },
                              
                              #' Find and catalog CRAM files with better existing file detection
                              #' 
                              #' @param bam_dir Directory containing existing BAM files (optional)
                              #' @param gvcf_dir Directory containing existing GVCF files (optional)
                              #' @return data.frame with CRAM file information
                              find_cram_files = function(bam_dir = NULL, gvcf_dir = NULL) {
                                
                                private$log("Finding CRAM files...")
                                
                                # Find all CRAM files (recursive search)
                                cram_files <- list.files(self$cram_dir, pattern = "\\.cram$", 
                                                         full.names = TRUE, recursive = TRUE)
                                
                                if (length(cram_files) == 0) {
                                  stop("No CRAM files found in ", self$cram_dir)
                                }
                                
                                # Create sample information
                                sample_ids <- tools::file_path_sans_ext(basename(cram_files))
                                
                                self$sample_info <- data.frame(
                                  sample_id = sample_ids,
                                  cram_file = cram_files,
                                  bam_file = file.path(self$output_dir, "bam", paste0(sample_ids, ".bam")),
                                  gvcf_file = file.path(self$output_dir, "gvcf", paste0(sample_ids, ".g.vcf.gz")),
                                  stringsAsFactors = FALSE
                                )
                                
                                private$log(paste("Found", nrow(self$sample_info), "CRAM files"))
                                
                                # Check for existing files with better detection
                                private$detect_existing_files(bam_dir, gvcf_dir)
                                
                                return(self$sample_info)
                              },
                              
                              #' Manually set file paths for existing BAM and GVCF files
                              #' 
                              #' @param bam_dir Directory containing BAM files
                              #' @param gvcf_dir Directory containing GVCF files
                              #' @param bam_pattern Pattern for BAM file names (default: "\\.bam$")
                              #' @param gvcf_pattern Pattern for GVCF file names (default: "\\.g?\\.vcf\\.gz$")
                              set_existing_file_paths = function(bam_dir = NULL, gvcf_dir = NULL, 
                                                                 bam_pattern = "\\.bam$", gvcf_pattern = "\\.g?\\.vcf\\.gz$") {
                                
                                private$log("Manually setting existing file paths...")
                                
                                if (is.null(self$sample_info)) {
                                  stop("Run find_cram_files() first")
                                }
                                
                                # Reset file existence flags
                                self$sample_info$bam_exists <- FALSE
                                self$sample_info$gvcf_exists <- FALSE
                                
                                # Find BAM files
                                if (!is.null(bam_dir) && dir.exists(bam_dir)) {
                                  private$log(paste("Searching for BAM files in:", bam_dir))
                                  
                                  bam_files <- list.files(bam_dir, pattern = bam_pattern, full.names = TRUE)
                                  private$log(paste("Found", length(bam_files), "BAM files"))
                                  
                                  if (length(bam_files) > 0) {
                                    # Show first few BAM files for debugging
                                    private$log(paste("Example BAM files:", paste(head(basename(bam_files), 3), collapse = ", ")))
                                    
                                    # Try to match BAM files to samples
                                    private$match_files_to_samples(bam_files, "bam")
                                  }
                                }
                                
                                # Find GVCF files
                                if (!is.null(gvcf_dir) && dir.exists(gvcf_dir)) {
                                  private$log(paste("Searching for GVCF files in:", gvcf_dir))
                                  
                                  gvcf_files <- list.files(gvcf_dir, pattern = gvcf_pattern, full.names = TRUE)
                                  private$log(paste("Found", length(gvcf_files), "GVCF files"))
                                  
                                  if (length(gvcf_files) > 0) {
                                    # Show first few GVCF files for debugging
                                    private$log(paste("Example GVCF files:", paste(head(basename(gvcf_files), 3), collapse = ", ")))
                                    
                                    # Try to match GVCF files to samples
                                    private$match_files_to_samples(gvcf_files, "gvcf")
                                  }
                                }
                                
                                private$log(paste("Final count - BAM files:", sum(self$sample_info$bam_exists)))
                                private$log(paste("Final count - GVCF files:", sum(self$sample_info$gvcf_exists)))
                                
                                return(invisible(self))
                              },
                              #' 
                              #' @param force_rerun Force reprocessing of existing BAM files
                              convert_cram_to_bam = function(force_rerun = FALSE) {
                                
                                private$log("Converting CRAM files to BAM format...")
                                
                                if (is.null(self$sample_info)) {
                                  self$find_cram_files()
                                }
                                
                                # Determine which files need processing
                                if (force_rerun) {
                                  to_process <- self$sample_info
                                } else {
                                  to_process <- self$sample_info[!self$sample_info$bam_exists, ]
                                }
                                
                                if (nrow(to_process) == 0) {
                                  private$log("All BAM files already exist. Use force_rerun=TRUE to reprocess.")
                                  return(invisible(self))
                                }
                                
                                private$log(paste("Processing", nrow(to_process), "CRAM files"))
                                
                                # Parallel processing function
                                convert_single <- function(i) {
                                  row <- to_process[i, ]
                                  
                                  tryCatch({
                                    # Convert CRAM to BAM
                                    cmd <- sprintf("%s view -@ %d -T %s -b -o %s %s",
                                                   self$samtools_path, self$threads, 
                                                   self$reference_genome, row$bam_file, row$cram_file)
                                    
                                    result <- system(cmd, intern = FALSE)
                                    if (result != 0) {
                                      stop("samtools view failed")
                                    }
                                    
                                    # Index BAM
                                    cmd <- sprintf("%s index -@ %d %s", 
                                                   self$samtools_path, self$threads, row$bam_file)
                                    
                                    result <- system(cmd, intern = FALSE)
                                    if (result != 0) {
                                      stop("samtools index failed")
                                    }
                                    
                                    return(list(sample_id = row$sample_id, status = "success"))
                                    
                                  }, error = function(e) {
                                    return(list(sample_id = row$sample_id, status = paste("error:", e$message)))
                                  })
                                }
                                
                                # Use parallel processing
                                cl <- makeCluster(min(self$threads, nrow(to_process)))
                                clusterExport(cl, c("to_process", "self"), envir = environment())
                                
                                results <- parLapply(cl, 1:nrow(to_process), convert_single)
                                stopCluster(cl)
                                
                                # Process results
                                success_count <- sum(sapply(results, function(x) x$status == "success"))
                                failed_samples <- sapply(results[sapply(results, function(x) x$status != "success")], 
                                                         function(x) x$sample_id)
                                
                                private$log(paste("Successfully converted", success_count, "files"))
                                if (length(failed_samples) > 0) {
                                  private$log(paste("Failed samples:", paste(failed_samples, collapse = ", ")))
                                }
                                
                                # Update sample info
                                self$sample_info$bam_exists <- file.exists(self$sample_info$bam_file)
                                
                                return(invisible(self))
                              },
                              
                              #' Call variants using GATK HaplotypeCaller in parallel
                              #' 
                              #' @param force_rerun Force reprocessing of existing GVCF files
                              call_variants_per_sample = function(force_rerun = FALSE) {
                                
                                private$log("Calling variants per sample with GATK HaplotypeCaller...")
                                
                                # Ensure BAM files exist
                                if (!all(self$sample_info$bam_exists)) {
                                  private$log("Some BAM files missing. Converting CRAM files first...")
                                  self$convert_cram_to_bam(force_rerun = FALSE)
                                }
                                
                                # Determine which files need processing
                                if (force_rerun) {
                                  to_process <- self$sample_info[self$sample_info$bam_exists, ]
                                } else {
                                  to_process <- self$sample_info[self$sample_info$bam_exists & !self$sample_info$gvcf_exists, ]
                                }
                                
                                if (nrow(to_process) == 0) {
                                  private$log("All GVCF files already exist. Use force_rerun=TRUE to reprocess.")
                                  return(invisible(self))
                                }
                                
                                private$log(paste("Processing", nrow(to_process), "samples for variant calling"))
                                
                                # Get intervals argument
                                intervals_arg <- if (!is.null(self$intervals_file)) {
                                  paste("-L", self$intervals_file)
                                } else {
                                  private$get_intervals_argument()
                                }
                                
                                # Parallel processing function
                                call_variants_single <- function(i) {
                                  row <- to_process[i, ]
                                  
                                  tryCatch({
                                    # Build GATK HaplotypeCaller command (prevent MNPs)
                                    cmd <- sprintf(
                                      "%s --java-options '-Xmx%s' HaplotypeCaller -R %s -I %s -O %s -ERC GVCF --native-pair-hmm-threads 1 --max-mnp-distance 0 --min-base-quality-score %d --min-mapping-quality-score %d --tmp-dir %s %s",
                                      self$gatk_path,
                                      self$memory,
                                      self$reference_genome,
                                      row$bam_file,
                                      row$gvcf_file,
                                      self$min_base_quality,
                                      self$min_mapping_quality,
                                      file.path(self$output_dir, "temp"),
                                      intervals_arg
                                    )
                                    
                                    result <- system(cmd, intern = FALSE)
                                    if (result != 0) {
                                      stop("GATK HaplotypeCaller failed")
                                    }
                                    
                                    # Index GVCF
                                    cmd <- sprintf("%s IndexFeatureFile -I %s", self$gatk_path, row$gvcf_file)
                                    result <- system(cmd, intern = FALSE)
                                    if (result != 0) {
                                      stop("GATK IndexFeatureFile failed")
                                    }
                                    
                                    return(list(sample_id = row$sample_id, status = "success"))
                                    
                                  }, error = function(e) {
                                    return(list(sample_id = row$sample_id, status = paste("error:", e$message)))
                                  })
                                }
                                
                                # Use parallel processing with smaller cluster size for memory management
                                max_parallel <- min(8, self$threads, nrow(to_process))  # Limit parallel jobs for memory
                                cl <- makeCluster(max_parallel)
                                clusterExport(cl, c("to_process", "self", "intervals_arg"), envir = environment())
                                
                                results <- parLapply(cl, 1:nrow(to_process), call_variants_single)
                                stopCluster(cl)
                                
                                # Process results
                                success_count <- sum(sapply(results, function(x) x$status == "success"))
                                failed_samples <- sapply(results[sapply(results, function(x) x$status != "success")], 
                                                         function(x) x$sample_id)
                                
                                private$log(paste("Successfully called variants for", success_count, "samples"))
                                if (length(failed_samples) > 0) {
                                  private$log(paste("Failed samples:", paste(failed_samples, collapse = ", ")))
                                }
                                
                                # Update sample info
                                self$sample_info$gvcf_exists <- file.exists(self$sample_info$gvcf_file)
                                
                                return(invisible(self))
                              },
                              
                              #' Clean GVCFs to remove MNPs for GenomicsDB compatibility
                              #' 
                              #' @param force_rerun Force reprocessing of existing cleaned files
                              clean_gvcfs_for_genomicsdb = function(force_rerun = FALSE) {
                                
                                private$log("Cleaning GVCFs to remove MNPs for GenomicsDB compatibility...")
                                
                                available_samples <- self$sample_info[self$sample_info$gvcf_exists, ]
                                
                                if (nrow(available_samples) == 0) {
                                  stop("No GVCF files found. Run call_variants_per_sample() first.")
                                }
                                
                                # Add cleaned GVCF paths
                                available_samples$cleaned_gvcf <- sub("\\.g\\.vcf\\.gz$", ".cleaned.g.vcf.gz", 
                                                                      available_samples$gvcf_file)
                                
                                # Determine which files need cleaning
                                if (force_rerun) {
                                  to_clean <- available_samples
                                } else {
                                  to_clean <- available_samples[!file.exists(available_samples$cleaned_gvcf), ]
                                }
                                
                                if (nrow(to_clean) == 0) {
                                  private$log("All cleaned GVCF files already exist.")
                                  # Update sample_info with cleaned paths
                                  self$sample_info$gvcf_file[self$sample_info$gvcf_exists] <- available_samples$cleaned_gvcf
                                  return(invisible(self))
                                }
                                
                                private$log(paste("Cleaning", nrow(to_clean), "GVCF files"))
                                
                                # Clean GVCFs in parallel
                                clean_single_gvcf <- function(i) {
                                  row <- to_clean[i, ]
                                  
                                  tryCatch({
                                    # Remove MNPs using bcftools
                                    cmd <- sprintf(
                                      "%s view -v ^mnps -O z -o %s %s",
                                      self$bcftools_path,
                                      row$cleaned_gvcf,
                                      row$gvcf_file
                                    )
                                    
                                    result <- system(cmd, intern = FALSE)
                                    if (result != 0) {
                                      stop("bcftools MNP removal failed")
                                    }
                                    
                                    # Index cleaned GVCF
                                    cmd <- sprintf("%s index %s", self$bcftools_path, row$cleaned_gvcf)
                                    result <- system(cmd, intern = FALSE)
                                    if (result != 0) {
                                      stop("bcftools index failed")
                                    }
                                    
                                    return(list(sample_id = row$sample_id, status = "success"))
                                    
                                  }, error = function(e) {
                                    return(list(sample_id = row$sample_id, status = paste("error:", e$message)))
                                  })
                                }
                                
                                # Use parallel processing
                                cl <- makeCluster(min(8, self$threads))
                                clusterExport(cl, c("to_clean", "self"), envir = environment())
                                
                                results <- parLapply(cl, 1:nrow(to_clean), clean_single_gvcf)
                                stopCluster(cl)
                                
                                # Process results
                                success_count <- sum(sapply(results, function(x) x$status == "success"))
                                
                                private$log(paste("Successfully cleaned", success_count, "GVCF files"))
                                
                                # Update sample_info to use cleaned files
                                self$sample_info$gvcf_file[self$sample_info$gvcf_exists] <- available_samples$cleaned_gvcf
                                
                                return(invisible(self))
                              },
                              
                              #' Emergency fix: merge existing batch VCF files using bcftools
                              #' This is a workaround for the sample batching issue
                              fix_batch_vcf_merge = function() {
                                
                                private$log("Attempting to fix batch VCF merging issue...")
                                
                                vcf_dir <- file.path(self$output_dir, "vcf")
                                batch_vcfs <- list.files(vcf_dir, pattern = "^batch_[0-9]+_variants\\.vcf\\.gz$", full.names = TRUE)
                                
                                if (length(batch_vcfs) == 0) {
                                  stop("No batch VCF files found to merge")
                                }
                                
                                private$log(paste("Found", length(batch_vcfs), "batch VCF files to merge"))
                                
                                # Sort batch files numerically
                                batch_numbers <- as.numeric(gsub(".*batch_([0-9]+)_variants\\.vcf\\.gz$", "\\1", basename(batch_vcfs)))
                                batch_vcfs <- batch_vcfs[order(batch_numbers)]
                                
                                # Use bcftools merge instead of GATK MergeVcfs
                                # bcftools merge can handle different samples across VCFs
                                final_vcf <- file.path(vcf_dir, "raw_variants_merged.vcf.gz")
                                
                                # First, ensure all VCFs are indexed
                                private$log("Indexing batch VCF files...")
                                for (vcf_file in batch_vcfs) {
                                  if (!file.exists(paste0(vcf_file, ".tbi"))) {
                                    cmd <- sprintf("%s index %s", self$bcftools_path, vcf_file)
                                    system(cmd, intern = FALSE)
                                  }
                                }
                                
                                # Merge using bcftools
                                private$log("Merging VCF files with bcftools...")
                                vcf_list <- paste(batch_vcfs, collapse = " ")
                                
                                cmd <- sprintf(
                                  "%s merge -O z -o %s --threads %d %s",
                                  self$bcftools_path, final_vcf, self$threads, vcf_list
                                )
                                
                                result <- system(cmd, intern = FALSE)
                                if (result != 0) {
                                  stop("bcftools merge failed")
                                }
                                
                                # Index the merged VCF
                                cmd <- sprintf("%s index %s", self$bcftools_path, final_vcf)
                                system(cmd, intern = FALSE)
                                
                                private$log(paste("Successfully merged batch VCFs into:", final_vcf))
                                
                                return(final_vcf)
                              },
                              joint_genotype = function() {
                                
                                private$log("Performing joint genotyping with GenomicsDB...")
                                
                                # Get samples with existing GVCFs (should be cleaned at this point)
                                available_samples <- self$sample_info[self$sample_info$gvcf_exists, ]
                                
                                if (nrow(available_samples) == 0) {
                                  stop("No GVCF files found. Run call_variants_per_sample() first.")
                                }
                                
                                private$log(paste("Joint genotyping", nrow(available_samples), "samples"))
                                
                                # Create sample map for GenomicsDB
                                sample_map_file <- file.path(self$output_dir, "sample_map.txt")
                                write.table(
                                  data.frame(sample_id = available_samples$sample_id, 
                                             gvcf_path = available_samples$gvcf_file),
                                  file = sample_map_file,
                                  sep = "\t",
                                  quote = FALSE,
                                  row.names = FALSE,
                                  col.names = FALSE
                                )
                                
                                # For large cohorts, process in batches using GenomicsDB
                                if (nrow(available_samples) > self$batch_size) {
                                  return(private$joint_genotype_batched(available_samples, sample_map_file))
                                } else {
                                  return(private$joint_genotype_single(sample_map_file))
                                }
                              },
                              
                              #' Filter variants using GATK best practices
                              #' 
                              #' @param raw_vcf Path to raw VCF file
                              filter_variants = function(raw_vcf) {
                                
                                private$log("Filtering variants using GATK best practices...")
                                
                                vcf_dir <- file.path(self$output_dir, "vcf")
                                
                                # Extract SNPs
                                snps_vcf <- file.path(vcf_dir, "raw_snps.vcf.gz")
                                cmd <- sprintf("%s SelectVariants -R %s -V %s -select-type SNP -O %s",
                                               self$gatk_path, self$reference_genome, raw_vcf, snps_vcf)
                                
                                result <- system(cmd, intern = FALSE)
                                if (result != 0) stop("SNP extraction failed")
                                
                                # Extract INDELs  
                                indels_vcf <- file.path(vcf_dir, "raw_indels.vcf.gz")
                                cmd <- sprintf("%s SelectVariants -R %s -V %s -select-type INDEL -O %s",
                                               self$gatk_path, self$reference_genome, raw_vcf, indels_vcf)
                                
                                result <- system(cmd, intern = FALSE)
                                if (result != 0) stop("INDEL extraction failed")
                                
                                # Filter SNPs
                                filtered_snps_vcf <- file.path(vcf_dir, "filtered_snps.vcf.gz")
                                cmd <- sprintf(
                                  '%s VariantFiltration -R %s -V %s -O %s --filter-expression "QD < 2.0" --filter-name "QD2" --filter-expression "QUAL < 30.0" --filter-name "QUAL30" --filter-expression "SOR > 3.0" --filter-name "SOR3" --filter-expression "FS > 60.0" --filter-name "FS60" --filter-expression "MQ < 40.0" --filter-name "MQ40" --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"',
                                  self$gatk_path, self$reference_genome, snps_vcf, filtered_snps_vcf
                                )
                                
                                result <- system(cmd, intern = FALSE)
                                if (result != 0) stop("SNP filtering failed")
                                
                                # Filter INDELs
                                filtered_indels_vcf <- file.path(vcf_dir, "filtered_indels.vcf.gz")
                                cmd <- sprintf(
                                  '%s VariantFiltration -R %s -V %s -O %s --filter-expression "QD < 2.0" --filter-name "QD2" --filter-expression "QUAL < 30.0" --filter-name "QUAL30" --filter-expression "FS > 200.0" --filter-name "FS200" --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20"',
                                  self$gatk_path, self$reference_genome, indels_vcf, filtered_indels_vcf
                                )
                                
                                result <- system(cmd, intern = FALSE)
                                if (result != 0) stop("INDEL filtering failed")
                                
                                # Merge filtered variants
                                final_vcf <- file.path(vcf_dir, "filtered_variants.vcf.gz")
                                
                                # Create sequence dictionary if it doesn't exist
                                dict_file <- paste0(tools::file_path_sans_ext(self$reference_genome), ".dict")
                                if (!file.exists(dict_file)) {
                                  private$log("Creating sequence dictionary...")
                                  cmd <- sprintf("%s CreateSequenceDictionary -R %s -O %s",
                                                 self$gatk_path, self$reference_genome, dict_file)
                                  system(cmd, intern = FALSE)
                                }
                                
                                cmd <- sprintf("%s MergeVcfs -I %s -I %s -O %s",
                                               self$gatk_path, filtered_snps_vcf, filtered_indels_vcf, final_vcf)
                                
                                result <- system(cmd, intern = FALSE)
                                if (result != 0) stop("VCF merging failed")
                                
                                # Create PASS-only version
                                pass_vcf <- file.path(vcf_dir, "pass_variants.vcf.gz")
                                cmd <- sprintf("%s SelectVariants -R %s -V %s -O %s --exclude-filtered",
                                               self$gatk_path, self$reference_genome, final_vcf, pass_vcf)
                                
                                result <- system(cmd)
                                if (result == 0) {
                                  private$log(paste("PASS-only VCF created:", pass_vcf))
                                  return(pass_vcf)
                                } else {
                                  private$log("Warning: Could not create PASS-only VCF")
                                  return(final_vcf)
                                }
                              },
                              
                              #' Fix existing batch VCF files that were created with sample-based batching
                              #' 
                              #' @param force_rerun Force reprocessing even if merged file exists
                              fix_batch_vcf_issue = function(force_rerun = FALSE) {
                                
                                private$log("Checking for batch VCF merging issues...")
                                
                                vcf_dir <- file.path(self$output_dir, "vcf")
                                
                                # Look for batch VCF files
                                batch_vcfs <- list.files(vcf_dir, pattern = "^batch_[0-9]+_variants\\.vcf\\.gz$", full.names = TRUE)
                                merged_vcf <- file.path(vcf_dir, "raw_variants_merged.vcf.gz")
                                
                                if (length(batch_vcfs) == 0) {
                                  private$log("No batch VCF files found to merge")
                                  return(NULL)
                                }
                                
                                if (file.exists(merged_vcf) && !force_rerun) {
                                  private$log("Merged VCF already exists. Use force_rerun=TRUE to recreate.")
                                  return(merged_vcf)
                                }
                                
                                private$log(paste("Found", length(batch_vcfs), "batch VCF files that need merging"))
                                
                                # Use the private method to fix the batch VCFs
                                merged_vcf <- private$fix_existing_batch_vcfs(batch_vcfs)
                                
                                return(merged_vcf)
                              },
                              
                              #' Generate comprehensive statistics and QC plots
                              generate_statistics = function(vcf_files) {
                                
                                private$log("Generating statistics and QC plots...")
                                
                                stats_dir <- file.path(self$output_dir, "stats")
                                
                                # Generate bcftools stats for each VCF
                                for (vcf_file in vcf_files) {
                                  if (file.exists(vcf_file)) {
                                    vcf_name <- tools::file_path_sans_ext(basename(vcf_file))
                                    stats_file <- file.path(stats_dir, paste0(vcf_name, "_stats.txt"))
                                    
                                    cmd <- sprintf("%s stats %s > %s", self$bcftools_path, vcf_file, stats_file)
                                    system(cmd, intern = FALSE)
                                  }
                                }
                                
                                # Generate pipeline summary
                                private$generate_pipeline_summary(vcf_files)
                                
                                # Generate QC plots
                                private$generate_qc_plots(vcf_files)
                                
                                return(invisible(self))
                              },
                              
                              #' Run the complete pipeline
                              #' 
                              #' @param force_rerun Force reprocessing of existing files
                              #' @param skip_bam_conversion Skip BAM conversion if BAM files exist
                              #' @param skip_variant_calling Skip variant calling if GVCF files exist
                              run_pipeline = function(force_rerun = FALSE, skip_bam_conversion = FALSE, 
                                                      skip_variant_calling = FALSE) {
                                
                                private$log("Starting complete GATK CRAM to VCF pipeline...")
                                
                                start_time <- Sys.time()
                                
                                tryCatch({
                                  
                                  # Step 1: Find CRAM files
                                  private$log("Step 1: Finding CRAM files...")
                                  self$find_cram_files()
                                  
                                  # Step 2: Convert CRAM to BAM
                                  if (!skip_bam_conversion) {
                                    private$log("Step 2: Converting CRAM to BAM...")
                                    self$convert_cram_to_bam(force_rerun = force_rerun)
                                  } else {
                                    private$log("Step 2: Skipping BAM conversion")
                                  }
                                  
                                  # Step 3: Call variants per sample
                                  if (!skip_variant_calling) {
                                    private$log("Step 3: Calling variants per sample...")
                                    self$call_variants_per_sample(force_rerun = force_rerun)
                                  } else {
                                    private$log("Step 3: Skipping variant calling")
                                  }
                                  
                                  # Step 4: Clean GVCFs (remove MNPs for GenomicsDB compatibility)
                                  private$log("Step 4: Cleaning GVCFs for GenomicsDB...")
                                  self$clean_gvcfs_for_genomicsdb(force_rerun = force_rerun)
                                  
                                  # Step 5: Joint genotyping (will auto-detect and fix batch VCF issues)
                                  private$log("Step 5: Joint genotyping...")
                                  raw_vcf <- self$joint_genotype()
                                  
                                  # Step 5.5: Check if we need to fix batch VCF merging issues
                                  if (is.null(raw_vcf) || !file.exists(raw_vcf)) {
                                    private$log("Joint genotyping may have failed, checking for batch VCF fix...")
                                    raw_vcf <- self$fix_batch_vcf_issue(force_rerun = force_rerun)
                                    
                                    if (is.null(raw_vcf) || !file.exists(raw_vcf)) {
                                      stop("Joint genotyping failed and no batch VCF files found to merge")
                                    }
                                  }
                                  
                                  # Step 6: Variant filtering
                                  private$log("Step 6: Filtering variants...")
                                  final_vcf <- self$filter_variants(raw_vcf)
                                  
                                  # Step 7: Generate statistics
                                  private$log("Step 7: Generating statistics...")
                                  self$generate_statistics(c(raw_vcf, final_vcf))
                                  
                                  end_time <- Sys.time()
                                  runtime <- difftime(end_time, start_time, units = "hours")
                                  
                                  private$log(paste("Pipeline completed successfully in", round(runtime, 2), "hours"))
                                  private$log(paste("Final VCF:", final_vcf))
                                  
                                  return(final_vcf)
                                  
                                }, error = function(e) {
                                  private$log(paste("Pipeline failed:", e$message))
                                  stop(e)
                                })
                              }
                            ),
                            
                            private = list(
                              
                              # Detect existing BAM and GVCF files with flexible matching
                              detect_existing_files = function(bam_dir = NULL, gvcf_dir = NULL) {
                                
                                private$log("Detecting existing BAM and GVCF files...")
                                
                                # Default directories
                                bam_search_dirs <- c(
                                  file.path(self$output_dir, "bam"),
                                  bam_dir,
                                  self$output_dir,
                                  dirname(self$sample_info$cram_file[1])  # Same dir as CRAM files
                                )
                                bam_search_dirs <- unique(bam_search_dirs[!is.null(bam_search_dirs)])
                                
                                gvcf_search_dirs <- c(
                                  file.path(self$output_dir, "gvcf"),
                                  gvcf_dir,
                                  self$output_dir,
                                  dirname(self$sample_info$cram_file[1])
                                )
                                gvcf_search_dirs <- unique(gvcf_search_dirs[!is.null(gvcf_search_dirs)])
                                
                                # Initialize as FALSE
                                self$sample_info$bam_exists <- FALSE
                                self$sample_info$gvcf_exists <- FALSE
                                
                                # Find BAM files
                                for (search_dir in bam_search_dirs) {
                                  if (dir.exists(search_dir)) {
                                    private$log(paste("Searching for BAM files in:", search_dir))
                                    
                                    # Find all BAM files in directory
                                    bam_files <- list.files(search_dir, pattern = "\\.bam$", full.names = TRUE)
                                    
                                    if (length(bam_files) > 0) {
                                      bam_basenames <- tools::file_path_sans_ext(basename(bam_files))
                                      
                                      # Match sample IDs to BAM files
                                      for (i in 1:nrow(self$sample_info)) {
                                        sample_id <- self$sample_info$sample_id[i]
                                        
                                        # Try exact match first
                                        exact_match <- match(sample_id, bam_basenames)
                                        if (!is.na(exact_match)) {
                                          self$sample_info$bam_file[i] <- bam_files[exact_match]
                                          self$sample_info$bam_exists[i] <- TRUE
                                          next
                                        }
                                        
                                        # Try partial matches
                                        partial_matches <- grep(sample_id, bam_basenames, fixed = TRUE)
                                        if (length(partial_matches) > 0) {
                                          self$sample_info$bam_file[i] <- bam_files[partial_matches[1]]
                                          self$sample_info$bam_exists[i] <- TRUE
                                          next
                                        }
                                        
                                        # Try reverse partial match (BAM name in sample ID)
                                        reverse_matches <- grep(bam_basenames, sample_id, fixed = TRUE)
                                        if (length(reverse_matches) > 0) {
                                          self$sample_info$bam_file[i] <- bam_files[reverse_matches[1]]
                                          self$sample_info$bam_exists[i] <- TRUE
                                        }
                                      }
                                    }
                                  }
                                }
                                
                                # Find GVCF files
                                for (search_dir in gvcf_search_dirs) {
                                  if (dir.exists(search_dir)) {
                                    private$log(paste("Searching for GVCF files in:", search_dir))
                                    
                                    # Find all GVCF files in directory
                                    gvcf_files <- list.files(search_dir, pattern = "\\.g?\\.vcf\\.gz$", full.names = TRUE)
                                    
                                    if (length(gvcf_files) > 0) {
                                      # Extract sample names from GVCF files (remove .g.vcf.gz or .vcf.gz)
                                      gvcf_basenames <- gsub("\\.(g\\.)?vcf\\.gz$", "", basename(gvcf_files))
                                      
                                      # Match sample IDs to GVCF files
                                      for (i in 1:nrow(self$sample_info)) {
                                        sample_id <- self$sample_info$sample_id[i]
                                        
                                        # Try exact match first
                                        exact_match <- match(sample_id, gvcf_basenames)
                                        if (!is.na(exact_match)) {
                                          self$sample_info$gvcf_file[i] <- gvcf_files[exact_match]
                                          self$sample_info$gvcf_exists[i] <- TRUE
                                          next
                                        }
                                        
                                        # Try partial matches
                                        partial_matches <- grep(sample_id, gvcf_basenames, fixed = TRUE)
                                        if (length(partial_matches) > 0) {
                                          self$sample_info$gvcf_file[i] <- gvcf_files[partial_matches[1]]
                                          self$sample_info$gvcf_exists[i] <- TRUE
                                          next
                                        }
                                        
                                        # Try reverse partial match
                                        reverse_matches <- grep(gvcf_basenames, sample_id, fixed = TRUE)
                                        if (length(reverse_matches) > 0) {
                                          self$sample_info$gvcf_file[i] <- gvcf_files[reverse_matches[1]]
                                          self$sample_info$gvcf_exists[i] <- TRUE
                                        }
                                      }
                                    }
                                  }
                                }
                                
                                private$log(paste("Found existing BAM files:", sum(self$sample_info$bam_exists)))
                                private$log(paste("Found existing GVCF files:", sum(self$sample_info$gvcf_exists)))
                                
                                # Debug information
                                private$debug_file_detection()
                              },
                              
                              # Debug file detection issues
                              debug_file_detection = function() {
                                
                                private$log("=== FILE DETECTION DEBUG ===")
                                
                                # Show first few samples and their file status
                                debug_samples <- head(self$sample_info, 5)
                                
                                for (i in 1:nrow(debug_samples)) {
                                  sample_id <- debug_samples$sample_id[i]
                                  private$log(paste("Sample:", sample_id))
                                  private$log(paste("  CRAM:", debug_samples$cram_file[i], "- exists:", file.exists(debug_samples$cram_file[i])))
                                  private$log(paste("  BAM: ", debug_samples$bam_file[i], "- exists:", debug_samples$bam_exists[i]))
                                  private$log(paste("  GVCF:", debug_samples$gvcf_file[i], "- exists:", debug_samples$gvcf_exists[i]))
                                }
                                
                                # List actual files in directories
                                bam_dir <- file.path(self$output_dir, "bam")
                                gvcf_dir <- file.path(self$output_dir, "gvcf")
                                
                                if (dir.exists(bam_dir)) {
                                  actual_bams <- list.files(bam_dir, pattern = "\\.bam$")
                                  private$log(paste("Actual BAM files in", bam_dir, ":", length(actual_bams)))
                                  if (length(actual_bams) > 0) {
                                    private$log(paste("  Examples:", paste(head(actual_bams, 3), collapse = ", ")))
                                  }
                                }
                                
                                if (dir.exists(gvcf_dir)) {
                                  actual_gvcfs <- list.files(gvcf_dir, pattern = "\\.vcf\\.gz$")
                                  private$log(paste("Actual GVCF files in", gvcf_dir, ":", length(actual_gvcfs)))
                                  if (length(actual_gvcfs) > 0) {
                                    private$log(paste("  Examples:", paste(head(actual_gvcfs, 3), collapse = ", ")))
                                  }
                                }
                                
                                private$log("=== END DEBUG ===")
                              },
                              create_directories = function() {
                                dirs <- c("bam", "gvcf", "vcf", "stats", "temp", "genomicsdb")
                                
                                for (dir in dirs) {
                                  dir_path <- file.path(self$output_dir, dir)
                                  if (!dir.exists(dir_path)) {
                                    dir.create(dir_path, recursive = TRUE)
                                  }
                                }
                              },
                              
                              # Logging function
                              log = function(message) {
                                timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
                                log_message <- paste0("[", timestamp, "] ", message)
                                
                                # Print to console
                                cat(log_message, "\n")
                                
                                # Write to log file
                                if (!is.null(self$log_file)) {
                                  cat(log_message, "\n", file = self$log_file, append = TRUE)
                                }
                              },
                              
                              # Get intervals argument for GATK commands
                              get_intervals_argument = function() {
                                
                                # Option 1: Try to find an existing intervals file
                                possible_intervals <- c(
                                  sub("\\.fasta$", ".bed", self$reference_genome),
                                  sub("\\.fasta$", ".intervals", self$reference_genome), 
                                  sub("\\.fa$", ".bed", self$reference_genome),
                                  sub("\\.fa$", ".intervals", self$reference_genome)
                                )
                                
                                for (intervals_file in possible_intervals) {
                                  if (file.exists(intervals_file)) {
                                    private$log(paste("Using intervals file:", intervals_file))
                                    return(paste("-L", intervals_file))
                                  }
                                }
                                
                                # Option 2: Create intervals from sequence dictionary
                                dict_file <- sub("\\.(fasta|fa)$", ".dict", self$reference_genome)
                                
                                if (file.exists(dict_file)) {
                                  private$log("Creating intervals from sequence dictionary...")
                                  intervals <- private$create_intervals_from_dict(dict_file)
                                  if (!is.null(intervals)) {
                                    return(intervals)
                                  }
                                }
                                
                                # Option 3: Use standard human chromosomes (fallback)
                                private$log("Using standard human chromosome intervals (fallback)")
                                
                                # Check if chromosomes have chr prefix
                                fai_file <- paste0(self$reference_genome, ".fai")
                                if (file.exists(fai_file)) {
                                  fai_first_line <- readLines(fai_file, n = 1)
                                  if (grepl("^chr", fai_first_line)) {
                                    chromosomes <- c(paste0("chr", 1:22), "chrX", "chrY")
                                  } else {
                                    chromosomes <- c(1:22, "X", "Y")
                                  }
                                } else {
                                  chromosomes <- c(paste0("chr", 1:22), "chrX", "chrY")
                                }
                                
                                interval_args <- paste("-L", chromosomes[1:5], collapse = " ")  # Start with first 5 chromosomes
                                
                                return(interval_args)
                              },
                              
                              # Create intervals from sequence dictionary
                              create_intervals_from_dict = function(dict_file) {
                                
                                tryCatch({
                                  # Read sequence dictionary
                                  dict_lines <- readLines(dict_file)
                                  seq_lines <- dict_lines[grepl("^@SQ", dict_lines)]
                                  
                                  if (length(seq_lines) == 0) {
                                    return(NULL)
                                  }
                                  
                                  # Extract sequence names
                                  seq_names <- gsub(".*SN:([^\\s]+).*", "\\1", seq_lines)
                                  
                                  # Filter for main chromosomes (avoid contigs/scaffolds for efficiency)
                                  main_chromosomes <- seq_names[grepl("^(chr)?([0-9]+|X|Y|MT|M)$", seq_names)]
                                  
                                  if (length(main_chromosomes) == 0) {
                                    # If no standard chromosomes found, use first 25 sequences
                                    main_chromosomes <- head(seq_names, 25)
                                  }
                                  
                                  private$log(paste("Found", length(main_chromosomes), "main chromosomes"))
                                  private$log(paste("Using chromosomes:", paste(head(main_chromosomes, 10), collapse = ", "), 
                                                    if(length(main_chromosomes) > 10) "..." else ""))
                                  
                                  # Create interval arguments
                                  interval_args <- paste("-L", head(main_chromosomes, 10), collapse = " ")  # Limit to first 10
                                  
                                  return(interval_args)
                                  
                                }, error = function(e) {
                                  private$log(paste("Error reading sequence dictionary:", e$message))
                                  return(NULL)
                                })
                              },
                              
                              # Joint genotyping for large cohorts using chromosome-based batching
                              joint_genotype_batched = function(available_samples, sample_map_file) {
                                
                                private$log(paste("Using chromosome-based batching for", nrow(available_samples), "samples"))
                                
                                # Check if there are existing batch VCF files that need fixing
                                vcf_dir <- file.path(self$output_dir, "vcf")
                                existing_batch_vcfs <- list.files(vcf_dir, pattern = "^batch_[0-9]+_variants\\.vcf\\.gz$", full.names = TRUE)
                                
                                if (length(existing_batch_vcfs) > 0) {
                                  private$log(paste("Found", length(existing_batch_vcfs), "existing batch VCF files"))
                                  private$log("These were created with sample-based batching and need to be merged using bcftools")
                                  
                                  # Fix the existing batch VCF files using bcftools merge
                                  return(private$fix_existing_batch_vcfs(existing_batch_vcfs))
                                }
                                
                                # Otherwise, proceed with chromosome-based approach for new processing
                                # Instead of sample batching, we'll use chromosome batching for very large cohorts
                                # This allows all samples to be processed together but reduces memory per chromosome
                                
                                # Get chromosomes from intervals file or create standard ones
                                if (!is.null(self$intervals_file) && file.exists(self$intervals_file)) {
                                  chromosomes <- readLines(self$intervals_file)
                                } else {
                                  # Check if reference uses chr prefix
                                  fai_file <- paste0(self$reference_genome, ".fai")
                                  if (file.exists(fai_file)) {
                                    fai_first_line <- readLines(fai_file, n = 1)
                                    has_chr_prefix <- grepl("^chr", fai_first_line)
                                  } else {
                                    has_chr_prefix <- TRUE  # Default assumption
                                  }
                                  
                                  if (has_chr_prefix) {
                                    chromosomes <- c(paste0("chr", 1:22), "chrX", "chrY")
                                  } else {
                                    chromosomes <- c(1:22, "X", "Y")
                                  }
                                }
                                
                                private$log(paste("Processing", length(chromosomes), "chromosomes"))
                                
                                # For very large cohorts (>500 samples), process by chromosome
                                # For moderate cohorts, try single GenomicsDB first
                                if (nrow(available_samples) > 500) {
                                  return(private$process_by_chromosome(available_samples, sample_map_file, chromosomes))
                                } else {
                                  # Try single GenomicsDB approach first
                                  tryCatch({
                                    return(private$joint_genotype_single(sample_map_file))
                                  }, error = function(e) {
                                    private$log(paste("Single GenomicsDB failed:", e$message))
                                    private$log("Falling back to chromosome-based processing...")
                                    return(private$process_by_chromosome(available_samples, sample_map_file, chromosomes))
                                  })
                                }
                              },
                              
                              # Fix existing batch VCF files using bcftools merge
                              fix_existing_batch_vcfs = function(batch_vcfs) {
                                
                                private$log("Fixing existing batch VCF files with bcftools merge...")
                                
                                vcf_dir <- dirname(batch_vcfs[1])
                                
                                # Sort batch files numerically
                                batch_numbers <- as.numeric(gsub(".*batch_([0-9]+)_variants\\.vcf\\.gz$", "\\1", basename(batch_vcfs)))
                                batch_vcfs <- batch_vcfs[order(batch_numbers)]
                                
                                private$log(paste("Found batch VCF files:", paste(basename(batch_vcfs), collapse = ", ")))
                                
                                # Check and create indices for all batch VCFs
                                private$log("Checking/creating VCF indices...")
                                for (vcf_file in batch_vcfs) {
                                  index_file <- paste0(vcf_file, ".tbi")
                                  if (!file.exists(index_file)) {
                                    private$log(paste("Creating index for", basename(vcf_file)))
                                    cmd <- sprintf("%s index %s", self$bcftools_path, vcf_file)
                                    result <- system(cmd, intern = FALSE)
                                    if (result != 0) {
                                      stop("Failed to index ", vcf_file)
                                    }
                                  }
                                }
                                
                                # Use bcftools merge to combine different samples from batch VCFs
                                final_vcf <- file.path(vcf_dir, "raw_variants_merged.vcf.gz")
                                
                                private$log("Merging batch VCF files using bcftools...")
                                vcf_list <- paste(batch_vcfs, collapse = " ")
                                
                                # Build regions argument if intervals file is specified
                                regions_arg <- ""
                                if (!is.null(self$intervals_file) && file.exists(self$intervals_file)) {
                                  private$log(paste("Restricting merge to regions in:", self$intervals_file))
                                  
                                  # Check if it's a BED file or interval list
                                  if (grepl("\\.(bed|BED)$", self$intervals_file)) {
                                    regions_arg <- paste("--regions-file", self$intervals_file)
                                  } else {
                                    # Assume it's an interval list, convert to regions format
                                    intervals <- readLines(self$intervals_file)
                                    # Remove empty lines and comments
                                    intervals <- intervals[!grepl("^$|^#", intervals)]
                                    if (length(intervals) > 0) {
                                      # Join intervals with comma for bcftools --regions
                                      regions_string <- paste(intervals, collapse = ",")
                                      regions_arg <- paste("--regions", shQuote(regions_string))
                                    }
                                  }
                                }
                                
                                cmd <- sprintf(
                                  "%s merge -O z -o %s --threads %d %s %s",
                                  self$bcftools_path, final_vcf, self$threads, regions_arg, vcf_list
                                )
                                
                                private$log(paste("Running command:", cmd))
                                result <- system(cmd, intern = FALSE)
                                
                                if (result != 0) {
                                  stop("bcftools merge failed")
                                }
                                
                                # Index the merged VCF
                                private$log("Indexing merged VCF...")
                                cmd <- sprintf("%s index -t %s", self$bcftools_path, final_vcf)
                                system(cmd, intern = FALSE)
                                
                                # Get statistics
                                cmd <- sprintf("%s query -l %s | wc -l", self$bcftools_path, final_vcf)
                                n_samples <- as.numeric(system(cmd, intern = TRUE))
                                
                                cmd <- sprintf("%s view -H %s | wc -l", self$bcftools_path, final_vcf)
                                n_variants <- as.numeric(system(cmd, intern = TRUE))
                                
                                private$log(paste("Successfully merged batch VCFs!"))
                                private$log(paste("Final VCF:", final_vcf))
                                private$log(paste("Samples:", n_samples))
                                private$log(paste("Variants:", n_variants))
                                
                                if (!is.null(self$intervals_file) && file.exists(self$intervals_file)) {
                                  private$log("Variants restricted to specified intervals only")
                                }
                                
                                # Optional: Clean up original batch files after successful merge
                                # Uncomment the next few lines if you want to remove the batch files
                                # private$log("Cleaning up original batch VCF files...")
                                # for (vcf_file in batch_vcfs) {
                                #   file.remove(vcf_file)
                                #   file.remove(paste0(vcf_file, ".tbi"))
                                # }
                                
                                return(final_vcf)
                              },
                              
                              # Process samples by chromosome for very large cohorts
                              process_by_chromosome = function(available_samples, sample_map_file, chromosomes) {
                                
                                private$log("Processing by chromosome for large cohort...")
                                
                                chr_vcfs <- character(length(chromosomes))
                                
                                for (i in seq_along(chromosomes)) {
                                  chr <- chromosomes[i]
                                  private$log(paste("Processing chromosome", chr, "(", i, "of", length(chromosomes), ")"))
                                  
                                  # Create GenomicsDB for this chromosome with ALL samples
                                  genomicsdb_path <- file.path(self$output_dir, "genomicsdb", paste0("chr_", chr))
                                  if (dir.exists(genomicsdb_path)) {
                                    unlink(genomicsdb_path, recursive = TRUE)
                                  }
                                  
                                  cmd <- sprintf(
                                    "%s --java-options '-Xmx%s' GenomicsDBImport --sample-name-map %s --genomicsdb-workspace-path %s -R %s --tmp-dir %s -L %s",
                                    self$gatk_path, self$memory, sample_map_file, genomicsdb_path, self$reference_genome,
                                    file.path(self$output_dir, "temp"), chr
                                  )
                                  
                                  result <- system(cmd, intern = FALSE)
                                  if (result != 0) {
                                    stop(paste("GenomicsDBImport failed for chromosome", chr))
                                  }
                                  
                                  # Joint genotyping for this chromosome with ALL samples
                                  chr_vcf <- file.path(self$output_dir, "vcf", paste0("chr_", chr, "_variants.vcf.gz"))
                                  cmd <- sprintf(
                                    "%s --java-options '-Xmx%s' GenotypeGVCFs -R %s -V gendb://%s -O %s --tmp-dir %s",
                                    self$gatk_path, self$memory, self$reference_genome, genomicsdb_path, chr_vcf,
                                    file.path(self$output_dir, "temp")
                                  )
                                  
                                  result <- system(cmd, intern = FALSE)
                                  if (result != 0) {
                                    stop(paste("GenotypeGVCFs failed for chromosome", chr))
                                  }
                                  
                                  chr_vcfs[i] <- chr_vcf
                                  private$log(paste("Completed chromosome", chr))
                                }
                                
                                # Merge chromosome VCFs (this is the correct use of MergeVcfs - same samples, different variants)
                                private$log("Merging chromosome VCFs...")
                                
                                final_vcf <- file.path(self$output_dir, "vcf", "raw_variants.vcf.gz")
                                input_args <- paste("-I", chr_vcfs, collapse = " ")
                                
                                cmd <- sprintf("%s MergeVcfs %s -O %s", self$gatk_path, input_args, final_vcf)
                                
                                result <- system(cmd, intern = FALSE)
                                if (result != 0) {
                                  stop("Chromosome VCF merging failed")
                                }
                                
                                # Clean up chromosome VCFs
                                file.remove(chr_vcfs)
                                
                                return(final_vcf)
                              },
                              
                              # Joint genotyping for smaller cohorts
                              joint_genotype_single = function(sample_map_file) {
                                
                                private$log("Using single GenomicsDB approach")
                                
                                # Create GenomicsDB
                                genomicsdb_path <- file.path(self$output_dir, "genomicsdb", "cohort")
                                if (dir.exists(genomicsdb_path)) {
                                  unlink(genomicsdb_path, recursive = TRUE)
                                }
                                
                                # Get intervals
                                intervals_arg <- if (!is.null(self$intervals_file)) {
                                  paste("-L", self$intervals_file)
                                } else {
                                  private$get_intervals_argument()
                                }
                                
                                cmd <- sprintf(
                                  "%s --java-options '-Xmx%s' GenomicsDBImport --sample-name-map %s --genomicsdb-workspace-path %s -R %s --tmp-dir %s %s",
                                  self$gatk_path, self$memory, sample_map_file, genomicsdb_path, self$reference_genome,
                                  file.path(self$output_dir, "temp"), intervals_arg
                                )
                                
                                result <- system(cmd, intern = FALSE)
                                if (result != 0) stop("GenomicsDBImport failed")
                                
                                # Joint genotyping
                                raw_vcf <- file.path(self$output_dir, "vcf", "raw_variants.vcf.gz")
                                cmd <- sprintf(
                                  "%s --java-options '-Xmx%s' GenotypeGVCFs -R %s -V gendb://%s -O %s --tmp-dir %s",
                                  self$gatk_path, self$memory, self$reference_genome, genomicsdb_path, raw_vcf,
                                  file.path(self$output_dir, "temp")
                                )
                                
                                result <- system(cmd, intern = FALSE)
                                if (result != 0) stop("GenotypeGVCFs failed")
                                
                                return(raw_vcf)
                              },
                              
                              # Generate pipeline summary
                              generate_pipeline_summary = function(vcf_files) {
                                
                                summary_file <- file.path(self$output_dir, "stats", "pipeline_summary.txt")
                                
                                # Get variant counts
                                variant_counts <- list()
                                sample_counts <- list()
                                
                                for (vcf_file in vcf_files) {
                                  if (file.exists(vcf_file)) {
                                    vcf_name <- tools::file_path_sans_ext(basename(vcf_file))
                                    
                                    # Count variants
                                    cmd <- sprintf("%s view -H %s | wc -l", self$bcftools_path, vcf_file)
                                    result <- system(cmd, intern = TRUE)
                                    variant_counts[[vcf_name]] <- as.numeric(result)
                                    
                                    # Count samples
                                    cmd <- sprintf("%s query -l %s | wc -l", self$bcftools_path, vcf_file)
                                    result <- system(cmd, intern = TRUE)
                                    sample_counts[[vcf_name]] <- as.numeric(result)
                                  }
                                }
                                
                                # Write summary
                                summary_text <- c(
                                  "GATK CRAM to VCF Pipeline Summary",
                                  "================================",
                                  paste("Date:", Sys.time()),
                                  paste("CRAM directory:", self$cram_dir),
                                  paste("Output directory:", self$output_dir),
                                  paste("Reference genome:", self$reference_genome),
                                  paste("GATK path:", self$gatk_path),
                                  paste("Threads:", self$threads),
                                  paste("Memory:", self$memory),
                                  paste("Batch size:", self$batch_size),
                                  "",
                                  "Sample Information:",
                                  paste("Total CRAM files found:", nrow(self$sample_info)),
                                  paste("BAM files created:", sum(self$sample_info$bam_exists)),
                                  paste("GVCF files created:", sum(self$sample_info$gvcf_exists)),
                                  "",
                                  "Variant Counts:"
                                )
                                
                                for (vcf_name in names(variant_counts)) {
                                  summary_text <- c(summary_text, 
                                                    paste(vcf_name, ":", variant_counts[[vcf_name]], "variants,", 
                                                          sample_counts[[vcf_name]], "samples"))
                                }
                                
                                writeLines(summary_text, summary_file)
                                private$log(paste("Pipeline summary written to:", summary_file))
                              },
                              
                              # Generate QC plots
                              generate_qc_plots = function(vcf_files) {
                                
                                if (!requireNamespace("ggplot2", quietly = TRUE)) {
                                  private$log("ggplot2 not available, skipping QC plots")
                                  return(invisible(NULL))
                                }
                                
                                private$log("Generating QC plots...")
                                
                                # Sample completion plot
                                if (!is.null(self$sample_info)) {
                                  p1 <- ggplot2::ggplot(data.frame(
                                    Stage = c("CRAM", "BAM", "GVCF"),
                                    Count = c(nrow(self$sample_info), 
                                              sum(self$sample_info$bam_exists),
                                              sum(self$sample_info$gvcf_exists))
                                  ), ggplot2::aes(x = Stage, y = Count)) +
                                    ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
                                    ggplot2::labs(title = "Sample Processing Progress", 
                                                  y = "Number of Samples", x = "Processing Stage") +
                                    ggplot2::theme_minimal()
                                  
                                  ggplot2::ggsave(file.path(self$output_dir, "stats", "sample_progress.png"), 
                                                  p1, width = 8, height = 6)
                                }
                                
                                private$log("QC plots saved to stats directory")
                              }
                            )
)

#' Create and run GATK CRAM to VCF pipeline
#' 
#' @param cram_dir Directory containing CRAM files
#' @param output_dir Output directory
#' @param reference_genome Reference genome FASTA file
#' @param gatk_path Path to GATK executable (default: "gatk")
#' @param threads Number of threads (default: detectCores())
#' @param memory Memory for GATK (default: "8g")
#' @param batch_size Batch size for large cohorts (default: 50)
#' @param intervals_file Optional BED file for specific regions
#' @param force_rerun Force reprocessing of existing files
#' @param skip_bam_conversion Skip BAM conversion step
#' @param skip_variant_calling Skip variant calling step
#' 
#' @return Path to final filtered VCF file
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' # Basic usage
#' final_vcf <- run_gatk_cram_pipeline(
#'   cram_dir = "/path/to/cram/files",
#'   output_dir = "/path/to/output",
#'   reference_genome = "/path/to/reference.fa"
#' )
#' 
#' # Advanced usage with custom parameters
#' final_vcf <- run_gatk_cram_pipeline(
#'   cram_dir = "/path/to/cram/files",
#'   output_dir = "/path/to/output", 
#'   reference_genome = "/path/to/reference.fa",
#'   threads = 16,
#'   memory = "32g",
#'   batch_size = 100,
#'   intervals_file = "/path/to/regions.bed"
#' )
#' }
run_gatk_cram_pipeline <- function(cram_dir, output_dir, reference_genome,
                                   gatk_path = "gatk", threads = NULL, 
                                   memory = "8g", batch_size = 50,
                                   intervals_file = NULL, force_rerun = FALSE,
                                   skip_bam_conversion = FALSE, 
                                   skip_variant_calling = FALSE) {
  
  # Load required packages
  required_packages <- c("R6", "parallel", "data.table", "dplyr")
  
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Required package", pkg, "is not installed. Install with: install.packages('", pkg, "')", sep = ""))
    }
  }
  
  # Create pipeline object
  pipeline <- GATKCramPipeline$new(
    cram_dir = cram_dir,
    output_dir = output_dir,
    reference_genome = reference_genome,
    gatk_path = gatk_path,
    threads = threads,
    memory = memory,
    batch_size = batch_size
  )
  
  # Set intervals file if provided
  if (!is.null(intervals_file)) {
    pipeline$intervals_file <- intervals_file
  }
  
  # Run pipeline
  final_vcf <- pipeline$run_pipeline(
    force_rerun = force_rerun,
    skip_bam_conversion = skip_bam_conversion,
    skip_variant_calling = skip_variant_calling
  )
  
  return(final_vcf)
}

#' Test GATK pipeline with example data
#' 
#' @export
test_gatk_pipeline <- function() {
  
  cat("Testing GATK CRAM to VCF Pipeline\n")
  cat("=================================\n")
  
  # Check required packages
  required_packages <- c("R6", "parallel", "data.table", "dplyr", "stringr", "readr")
  missing_packages <- required_packages[!sapply(required_packages, function(pkg) requireNamespace(pkg, quietly = TRUE))]
  
  if (length(missing_packages) > 0) {
    cat("Missing required packages:", paste(missing_packages, collapse = ", "), "\n")
    cat("Install with: install.packages(c(", paste0("'", missing_packages, "'", collapse = ", "), "))\n")
    return(FALSE)
  }
  
  # Check required software
  software_checks <- list(
    "GATK" = system("gatk --version", ignore.stdout = TRUE, ignore.stderr = TRUE) == 0,
    "samtools" = system("samtools --version", ignore.stdout = TRUE, ignore.stderr = TRUE) == 0,
    "bcftools" = system("bcftools --version", ignore.stdout = TRUE, ignore.stderr = TRUE) == 0
  )
  
  cat("Software availability:\n")
  for (software in names(software_checks)) {
    status <- if (software_checks[[software]]) "✓ Available" else "✗ Missing"
    cat(sprintf("  %s: %s\n", software, status))
  }
  
  all_available <- all(unlist(software_checks))
  
  if (all_available) {
    cat("\nAll dependencies satisfied! Pipeline is ready to use.\n")
    cat("\nExample usage:\n")
    cat('final_vcf <- run_gatk_cram_pipeline(\n')
    cat('  cram_dir = "/path/to/cram/files",\n')
    cat('  output_dir = "/path/to/output",\n')
    cat('  reference_genome = "/path/to/reference.fa",\n')
    cat('  threads = 16,\n')
    cat('  memory = "32g",\n')
    cat('  batch_size = 100\n')
    cat(')\n')
  } else {
    cat("\nPlease install missing software before running the pipeline.\n")
  }
  
  return(all_available)
}

#' Scan directories to find existing BAM and GVCF files
#' 
#' @param output_dir Directory to scan
#' @param show_details Show detailed file listings
#' 
#' @export
scan_existing_files <- function(output_dir, show_details = TRUE) {
  
  cat("Scanning for existing files in:", output_dir, "\n")
  cat("=" , rep("=", nchar(output_dir) + 30), "\n")
  
  # Scan for BAM files
  bam_dirs_to_check <- c(
    file.path(output_dir, "bam"),
    output_dir,
    file.path(output_dir, "alignments"),
    file.path(output_dir, "processed")
  )
  
  bam_files_found <- c()
  for (bam_dir in bam_dirs_to_check) {
    if (dir.exists(bam_dir)) {
      bams <- list.files(bam_dir, pattern = "\\.bam$", full.names = TRUE)
      if (length(bams) > 0) {
        cat("BAM files in", bam_dir, ":", length(bams), "\n")
        if (show_details) {
          cat("  Examples:", paste(head(basename(bams), 3), collapse = ", "), "\n")
        }
        bam_files_found <- c(bam_files_found, bams)
      }
    }
  }
  
  # Scan for GVCF files
  gvcf_dirs_to_check <- c(
    file.path(output_dir, "gvcf"),
    output_dir,
    file.path(output_dir, "variants"),
    file.path(output_dir, "vcf")
  )
  
  gvcf_files_found <- c()
  for (gvcf_dir in gvcf_dirs_to_check) {
    if (dir.exists(gvcf_dir)) {
      gvcfs <- list.files(gvcf_dir, pattern = "\\.g?\\.vcf\\.gz$", full.names = TRUE)
      if (length(gvcfs) > 0) {
        cat("GVCF files in", gvcf_dir, ":", length(gvcfs), "\n")
        if (show_details) {
          cat("  Examples:", paste(head(basename(gvcfs), 3), collapse = ", "), "\n")
        }
        gvcf_files_found <- c(gvcf_files_found, gvcfs)
      }
    }
  }
  
  cat("\nSummary:\n")
  cat("Total BAM files found:", length(unique(bam_files_found)), "\n")
  cat("Total GVCF files found:", length(unique(gvcf_files_found)), "\n")
  
  return(list(
    bam_files = unique(bam_files_found),
    gvcf_files = unique(gvcf_files_found)
  ))
}


#' 
#' @param sample_map_file Path to sample map file
#' @param output_dir Output directory
#' @param bcftools_path Path to bcftools executable
#' 
#' @export
fix_mnp_gvcfs <- function(sample_map_file, output_dir, bcftools_path = "bcftools") {
  
  cat("Fixing MNP issue in existing GVCF files...\n")
  
  # Read sample map to get GVCF files
  sample_map <- read.table(sample_map_file, sep = "\t", stringsAsFactors = FALSE)
  colnames(sample_map) <- c("sample_id", "gvcf_file")
  
  # Create cleaned versions
  sample_map$cleaned_gvcf <- sub("\\.g\\.vcf\\.gz$", ".cleaned.g.vcf.gz", sample_map$gvcf_file)
  
  # Function to clean single GVCF
  clean_gvcf <- function(i) {
    row <- sample_map[i, ]
    
    cat("Cleaning GVCF for sample:", row$sample_id, "\n")
    
    tryCatch({
      # Remove MNPs using bcftools
      cmd <- sprintf("%s view -v ^mnps -O z -o %s %s", 
                     bcftools_path, row$cleaned_gvcf, row$gvcf_file)
      
      result <- system(cmd, intern = FALSE)
      if (result != 0) {
        stop("bcftools failed")
      }
      
      # Index cleaned GVCF
      cmd <- sprintf("%s index %s", bcftools_path, row$cleaned_gvcf)
      result <- system(cmd, intern = FALSE)
      if (result != 0) {
        stop("bcftools index failed")
      }
      
      return(TRUE)
      
    }, error = function(e) {
      cat("Error cleaning", row$sample_id, ":", e$message, "\n")
      return(FALSE)
    })
  }
  
  # Process in parallel
  cl <- makeCluster(min(8, detectCores()))
  clusterExport(cl, c("sample_map", "bcftools_path"), envir = environment())
  
  results <- parLapply(cl, 1:nrow(sample_map), clean_gvcf)
  stopCluster(cl)
  
  success_count <- sum(unlist(results))
  cat("Successfully cleaned", success_count, "out of", nrow(sample_map), "GVCF files\n")
  
  # Create new sample map with cleaned files
  cleaned_sample_map_file <- file.path(output_dir, "cleaned_sample_map.txt")
  write.table(
    data.frame(sample_id = sample_map$sample_id, gvcf_path = sample_map$cleaned_gvcf),
    file = cleaned_sample_map_file,
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
  )
  
  cat("Created cleaned sample map:", cleaned_sample_map_file, "\n")
  return(cleaned_sample_map_file)
}