suppressPackageStartupMessages(library(reticulate    ))
reticulate::use_python("/mmfs1/gscratch/scrubbed/mdanto/envs/bioinfo/bin/python") 
use_condaenv("bioinfo", required = TRUE)

Sys.setenv(PATH = paste("/mmfs1/gscratch/scrubbed/mdanto/envs/bioinfo/bin", Sys.getenv("PATH"), sep = ":"))

source("R/functions.R")

reference_genome = here("data/reference", "Homo_sapiens_assembly38.fasta")
cram_dir         = here("data/TODAY_WHOLE_EXOMES")
output_dir       = here("data/genotyping")
intervals_file   = here("data/reference/Homo_sapiens_assembly38.bed")

# Create a bed file with all the chromosome lengths
dict_to_bed <- function(dict_path, bed_path, chromosomes = NULL) {
  # Read dict file
  dict_lines <- readLines(dict_path) 
  
  # Keep only contig lines
  sq_lines <- grep("^@SQ", dict_lines, value = TRUE)
  
  # Extract contig name and length
  contigs <- sub(".*SN:([^ \t]+).*", "\\1", sq_lines)
  lengths <- as.integer(sub(".*LN:([0-9]+).*", "\\1", sq_lines))
  
  # Build data frame
  bed_df <- data.frame(
    chrom = contigs,
    start = 0,
    end   = lengths,
    stringsAsFactors = FALSE
  )
  
  # Filter if specific chromosomes are requested
  if (!is.null(chromosomes)) {
    bed_df <- bed_df[bed_df$chrom %in% chromosomes, ]
  }
  
  # Write BED
  write.table(
    bed_df,
    file = bed_path,
    sep = "\t",
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE
  )
  
  message("BED file saved to: ", bed_path)
}

# Full genome
dict_to_bed(dict_path   = here("data/reference", "Homo_sapiens_assembly38.dict"), 
            bed_path    = intervals_file,
            chromosomes = paste0("chr", c(1:22, "X", "Y"))
            )


# Source the pipeline code
source("R/r_gatk_cram_pipeline.r")

# Basic usage for large cohort (>100 samples)
final_vcf <- run_gatk_cram_pipeline(
  cram_dir                = cram_dir,
  output_dir              = output_dir,
  reference_genome        = reference_genome,
  intervals_file          = intervals_file,
  skip_bam_conversion     = TRUE,  # Set to TRUE if BAMs already exist
  skip_variant_calling    = TRUE,  # Set to TRUE if GVCFs already exist
  threads                 = 50,
  memory                  = "1000g",
  batch_size              = 20
)

# Divide VCF file into chromosomes and create plink files
source("/mmfs1/gscratch/togo/matteo/projects/software/R_scripts/vcf_chr_splitter.r")

outpur_dir_split = here("data/genotyping/vcf_split")

results <- split_vcf_by_chromosomes(
  vcf_file       = final_vcf,
  intervals_file = intervals_file,
  output_folder  = outpur_dir_split,
  #output_format  = "all",
  output_format  = "plink2",
  plink_memory   = 960000,   
  plink_threads  = 50 
)

# Impute variants using Beagle
source("/mmfs1/gscratch/togo/matteo/projects/software/R_scripts/imputation_wrapper.r")

outpur_dir_imputation = here("data/genotyping/imputation")
ref_panel_dir         = "/mmfs1/gscratch/togo/matteo/projects/data/1kgp_reference_imputation"
genetic_map_dir       = "/mmfs1/gscratch/togo/matteo/projects/data/genetic_maps"

results <- fast_imputation(
  vcf_dir          = here("data/genotyping/vcf_split/vcf/"),
  output_dir       = outpur_dir_imputation,
  ref_panel_dir    = ref_panel_dir,
  genetic_map_dir  = genetic_map_dir,
  beagle_jar       = "/mmfs1/gscratch/togo/matteo/projects/software/beagle.27Feb25.75f.jar",
  chromosomes      = 1:22,
  total_cores      = 50, 
  total_memory     = 1000
)

