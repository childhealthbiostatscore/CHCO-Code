# R GATK CRAM to VCF Pipeline for Large Cohorts (>100 samples)
# Optimized for GenomicsDB and parallel processing

suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

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
GATKCramPipeline <- R6::R6Class("GATKCramPipeline",
  
  public = list(
    
    # Configuration
    cram_dir = NULL,
    output_dir = NULL,
    reference_genome = NULL,
    gatk_path     = "gatk",
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
    initialize = function(cram_dir, output_dir, reference_genome, 
                         gatk_path = "gatk", 
                         threads = NULL, 
                         memory = "8g", 
                         batch_size = 50) {
      
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
    },
    
    #' Find and catalog CRAM files
    #' 
    #' @return data.frame with CRAM file information
    find_cram_files = function() {
      
      private$log("Finding CRAM files...")
      
      # Find all CRAM files
      #cram_files <- list.files(self$cram_dir, pattern = "\\.cram$", full.names = TRUE, recursive = FALSE)
      cram_files <- list.files(self$cram_dir, pattern = "\\.cram$", full.names = TRUE, recursive = TRUE)[1:3]
      
      if (length(cram_files) == 0) {
        cram_files = readLines("data/TODAY_WHOLE_EXOMES.list_files.txt")
        #stop("No CRAM files found in ", self$cram_dir)
      }
      
      # Create sample information
      self$sample_info <- data.frame(
        sample_id = tools::file_path_sans_ext(basename(cram_files)),
        cram_file = cram_files,
        bam_file = file.path(self$output_dir, "bam", 
                            paste0(tools::file_path_sans_ext(basename(cram_files)), ".bam")),
        gvcf_file = file.path(self$output_dir, "gvcf",
                             paste0(tools::file_path_sans_ext(basename(cram_files)), ".g.vcf.gz")),
        stringsAsFactors = FALSE
      )
      
      private$log(paste("Found", nrow(self$sample_info), "CRAM files"))
      
      # Check for existing files
      self$sample_info$bam_exists <- file.exists(self$sample_info$bam_file)
      self$sample_info$gvcf_exists <- file.exists(self$sample_info$gvcf_file)
      
      private$log(paste("Existing BAM files:", sum(self$sample_info$bam_exists)))
      private$log(paste("Existing GVCF files:", sum(self$sample_info$gvcf_exists)))
      
      return(self$sample_info)
    },
    
    #' Convert CRAM files to BAM format in parallel
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
          
          print(cmd)
          
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
      
      # Parallel processing function
      call_variants_single <- function(i) {
        row <- to_process[i, ]

        tryCatch({
          # Build GATK HaplotypeCaller command
          cmd <- sprintf(
            #"%s HaplotypeCaller -R %s -I %s -O %s -ERC GVCF --native-pair-hmm-threads 1 --max-mnp-distance 2 --min-base-quality-score %d --min-mapping-quality-score %d --tmp-dir %s %s --java-options '-Xmx%s'",
            "%s HaplotypeCaller -R %s -I %s -O %s -ERC GVCF --native-pair-hmm-threads 1 --max-mnp-distance 2 --min-base-quality-score %d --minimum-mapping-quality %d --tmp-dir %s %s --java-options '-Xmx%s'",
            self$gatk_path,
            self$reference_genome,
            row$bam_file,
            row$gvcf_file,
            self$min_base_quality,
            self$min_mapping_quality,
            file.path(self$output_dir, "temp"),
            if (!is.null(self$intervals_file)) paste("-L", self$intervals_file) else "",
            self$memory
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
      clusterExport(cl, c("to_process", "self"), envir = environment())
      
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
    
    #' Perform joint genotyping using GenomicsDB (optimized for large cohorts)
    joint_genotype = function() {
      
      private$log("Performing joint genotyping with GenomicsDB...")
      
      # Get samples with existing GVCFs
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
        
        # Step 4: Joint genotyping
        private$log("Step 4: Joint genotyping...")
        raw_vcf <- self$joint_genotype()
        
        # Step 5: Variant filtering
        private$log("Step 5: Filtering variants...")
        final_vcf <- self$filter_variants(raw_vcf)
        
        # Step 6: Generate statistics
        private$log("Step 6: Generating statistics...")
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
    
    # Create output directory structure
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
    
    # Joint genotyping for large cohorts using batched approach
    joint_genotype_batched = function(available_samples, sample_map_file) {
      
      private$log(paste("Using batched approach for", nrow(available_samples), "samples"))
      
      # Split samples into batches
      n_batches <- ceiling(nrow(available_samples) / self$batch_size)
      batch_vcfs <- character(n_batches)
      
      for (batch_i in 1:n_batches) {
        private$log(paste("Processing batch", batch_i, "of", n_batches))
        
        # Get batch samples
        start_idx <- (batch_i - 1) * self$batch_size + 1
        end_idx <- min(batch_i * self$batch_size, nrow(available_samples))
        batch_samples <- available_samples[start_idx:end_idx, ]
        
        # Create batch sample map
        batch_sample_map <- file.path(self$output_dir, paste0("sample_map_batch_", batch_i, ".txt"))
        write.table(
          data.frame(sample_id = batch_samples$sample_id, 
                    gvcf_path = batch_samples$gvcf_file),
          file = batch_sample_map,
          sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
        )
        
        # Create GenomicsDB for this batch
        genomicsdb_path <- file.path(self$output_dir, "genomicsdb", paste0("batch_", batch_i))
        if (dir.exists(genomicsdb_path)) {
          unlink(genomicsdb_path, recursive = TRUE)
        }
        
        cmd <- sprintf(
          "%s  --java-options '-Xmx%s' GenomicsDBImport --sample-name-map %s --genomicsdb-workspace-path %s -R %s --tmp-dir %s -L %s",
          self$gatk_path, self$memory, sample_map_file, genomicsdb_path, self$reference_genome,
          file.path(self$output_dir, "temp"),
          if (!is.null(self$intervals_file)) self$intervals_file else sub("fasta", "bed", self$reference_genome)
        )
        
        result <- system(cmd, intern = FALSE)
        if (result != 0) stop(paste("GenomicsDBImport failed for batch", batch_i))
        
        # Joint genotyping for this batch
        batch_vcf <- file.path(self$output_dir, "vcf", paste0("batch_", batch_i, "_variants.vcf.gz"))
        cmd <- sprintf(
          "%s GenotypeGVCFs -R %s -V gendb://%s -O %s --tmp-dir %s --java-options '-Xmx%s'",
          self$gatk_path, self$reference_genome, genomicsdb_path, batch_vcf,
          file.path(self$output_dir, "temp"), self$memory
        )
        
        result <- system(cmd, intern = FALSE)
        if (result != 0) stop(paste("GenotypeGVCFs failed for batch", batch_i))
        
        batch_vcfs[batch_i] <- batch_vcf
        private$log(paste("Completed batch", batch_i))
      }
      
      # Merge all batch VCFs
      if (n_batches > 1) {
        private$log("Merging batch VCFs...")
        
        final_vcf <- file.path(self$output_dir, "vcf", "raw_variants.vcf.gz")
        input_args <- paste("-I", batch_vcfs, collapse = " ")
        
        cmd <- sprintf("%s MergeVcfs %s -O %s", self$gatk_path, input_args, final_vcf)
        
        result <- system(cmd, intern = FALSE)
        if (result != 0) stop("VCF merging failed")
        
        # Clean up batch VCFs
        file.remove(batch_vcfs)
        
        return(final_vcf)
      } else {
        return(batch_vcfs[1])
      }
    },
    
    # Joint genotyping for smaller cohorts
    joint_genotype_single = function(sample_map_file) {
      
      private$log("Using single GenomicsDB approach")
      
      # Create GenomicsDB
      genomicsdb_path <- file.path(self$output_dir, "genomicsdb", "cohort")
      if (dir.exists(genomicsdb_path)) {
        unlink(genomicsdb_path, recursive = TRUE)
      }
      
      cmd <- sprintf(
        "%s  --java-options '-Xmx%s' GenomicsDBImport --sample-name-map %s --genomicsdb-workspace-path %s -R %s --tmp-dir %s -L %s",
        self$gatk_path, self$memory, sample_map_file, genomicsdb_path, self$reference_genome,
        file.path(self$output_dir, "temp"),
        if (!is.null(self$intervals_file)) self$intervals_file else sub("fasta", "bed", self$reference_genome)
      )

      result <- system(cmd, intern = FALSE)
      if (result != 0) stop("GenomicsDBImport failed")
      
      # Joint genotyping
      raw_vcf <- file.path(self$output_dir, "vcf", "raw_variants.vcf.gz")
      cmd <- sprintf(
        "%s GenotypeGVCFs -R %s -V gendb://%s -O %s --tmp-dir %s --java-options '-Xmx%s'",
        self$gatk_path, self$reference_genome, genomicsdb_path, raw_vcf,
        file.path(self$output_dir, "temp"), self$memory
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
        p1 <- ggplot(data.frame(
          Stage = c("CRAM", "BAM", "GVCF"),
          Count = c(nrow(self$sample_info), 
                   sum(self$sample_info$bam_exists),
                   sum(self$sample_info$gvcf_exists))
        ), aes(x = Stage, y = Count)) +
          geom_bar(stat = "identity", fill = "steelblue") +
          labs(title = "Sample Processing Progress", 
               y = "Number of Samples", x = "Processing Stage") +
          theme_minimal()
        
        ggsave(file.path(self$output_dir, "stats", "sample_progress.png"), 
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
                                  gatk_path = "gatk", 
                                  threads = NULL, 
                                  memory = "8g", 
                                  batch_size = 50,
                                  intervals_file = NULL, 
                                  force_rerun = FALSE,
                                  skip_bam_conversion = FALSE, 
                                  skip_variant_calling = FALSE) {
  
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

# Example usage and testing function
#' Test GATK pipeline with example data
#' 
#' @export
test_gatk_pipeline <- function(gatk_path, samtools_path, bcftools_path) {
  
  cat("Testing GATK CRAM to VCF Pipeline\n")
  cat("=================================\n")
  
  # Check required packages
  required_packages <- c("R6", "parallel", "data.table", "dplyr", "ggplot2")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  
  if (length(missing_packages) > 0) {
    cat("Missing required packages:", paste(missing_packages, collapse = ", "), "\n")
    cat("Install with: install.packages(c(", paste0("'", missing_packages, "'", collapse = ", "), "))\n")
    return(FALSE)
  }
  
  # Check required software
  software_checks <- list(
    "GATK"     = system(paste(gatk_path    , "--version"), ignore.stdout = TRUE, ignore.stderr = TRUE) == 0,
    "samtools" = system(paste(samtools_path, "--version"), ignore.stdout = TRUE, ignore.stderr = TRUE) == 0,
    "bcftools" = system(paste(bcftools_path, "--version"), ignore.stdout = TRUE, ignore.stderr = TRUE) == 0
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