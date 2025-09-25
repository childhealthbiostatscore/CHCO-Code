source("R/functions.R")
source("R/idat_to_vcf_wrapper.r")

# Detect what array type you have
array_info <- detect_array_type(here("data/array"))

# See the results
print(array_info)

# Get specific information
cat("Detected array:", array_info$array_type)
cat("Confidence:", array_info$confidence) 
cat("Manifest URL:", array_info$manifest_url)
cat("Manifest files to look for:", paste(array_info$manifest_files, collapse = ", "))

# See recommendations
cat(array_info$recommendations, sep = "\n")

# Run pipeline
idat_to_vcf(here("data/array"), 
            output_vcf    = here("output/variants.vcf.gz"), 
            manifest_file = here("data/infinium-global-screening-array-24-v1-0-c2-manifest-file-csv-build38/GSA-24v1-0_C2.csv"))


# Test how it worked
library(vcfR)

vcf_file <- here("output/variants.vcf.gz")
vcf      <- read.vcfR(vcf_file)

# Explore
vcf
vcf@gt[1:5, 1:5]   # First 5 genotypes
vcf@fix[1:5, ]     # First 5 fixed fields
