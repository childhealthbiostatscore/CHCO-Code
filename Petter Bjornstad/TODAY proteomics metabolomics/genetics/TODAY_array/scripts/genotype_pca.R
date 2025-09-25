suppressPackageStartupMessages(library(reticulate    ))
reticulate::use_python("/mmfs1/gscratch/scrubbed/mdanto/envs/bioinfo/bin/python") 
use_condaenv("bioinfo", required = TRUE)

Sys.setenv(PATH = paste("/mmfs1/gscratch/scrubbed/mdanto/envs/bioinfo/bin", Sys.getenv("PATH"), sep = ":"))

source("R/functions.R")

# Look at the VCF files
library(vcfR)
vcf_file <- here("output", "variants.vcf.gz")
vcf      <- read.vcfR(vcf_file)


# Run genotype PCA on the 1000 Genomes chromosome 21 as test
source("/mmfs1/gscratch/togo/matteo/projects/software/R_scripts/plink_pca_wrapper.r")

system("wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr21.filtered.shapeit2-duohmm-phased.vcf.gz")
system("wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr21.filtered.shapeit2-duohmm-phased.vcf.gz.tbi")

dir.create(here("output", "genotype_pca_1kgp_chr21"), showWarnings = FALSE)
genotype_pca <- plink_vcf_pca("CCDG_14151_B01_GRM_WGS_2020-08-05_chr21.filtered.shapeit2-duohmm-phased.vcf.gz", 
                              output_prefix   = here("output/genotype_pca_1kgp_chr21", "pca"),
                              plink_version   = "2.0",
                              ld_window_size  = 1000,
                              ld_r2_threshold = 0.1,
                              mind            = 0.05,
                              geno            = 0.01)

# Run genotype PCA on dataset
source("/mmfs1/gscratch/togo/matteo/projects/software/R_scripts/plink_pca_wrapper.r")

dir.create(here("output", "genotype_pca"), showWarnings = FALSE)
#genotype_pca <- plink_vcf_pca(here("output", "variants.vcf.gz"), output_prefix = here("output/genotype_pca", "pca"))
genotype_pca <- plink_vcf_pca(here("output", "variants.vcf.gz"), 
                              output_prefix   = here("output/genotype_pca", "pca"),
                              plink_version   = "2.0",
                              ld_window_size  = 1000,
                              ld_r2_threshold = 0.1,
                              mind            = 0.05,
                              geno            = 0.01)





