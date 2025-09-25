suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(reticulate    ))
suppressPackageStartupMessages(library(tibble    ))
suppressPackageStartupMessages(library(tidyr    ))

my_lib_path = "/mmfs1/gscratch/togo/matteo/software/R"
.libPaths(my_lib_path)

conda_env_path = "/mmfs1/gscratch/togo/matteo/software/conda/r_deps"
# tell the compiler where to find headers/libraries
Sys.setenv(
  #PYTHONHOME = "/gscratch/scrubbed/mdanto/envs/r_deps/bin/python", # OLD, maybe useless
  CPATH           = paste(paste(conda_env_path, "include", sep = "/"), Sys.getenv("CPATH"          ), sep = ":"),
  LIBRARY_PATH    = paste(paste(conda_env_path, "lib"    , sep = "/"), Sys.getenv("LIBRARY_PATH"   ), sep = ":"),
  LD_LIBRARY_PATH = paste(paste(conda_env_path, "lib"    , sep = "/"), Sys.getenv("LD_LIBRARY_PATH"), sep = ":"),
  PATH            = paste(paste(conda_env_path, "bin"    , sep = "/"), Sys.getenv("PATH"           ), sep = ":")
)

################################################################################
# Download Gencode data (same as eQTL catalog)
# Create a SH file and download via command line

sh_file = here("scripts/download_gencode.sh")

writeLines(c(paste("mkdir", here("data/gencode_v39")),
             paste("cd"   , here("data/gencode_v39")),
             "",
             paste("wget", "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gtf.gz"),
             paste("wget", "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.metadata.EntrezGene.gz")
), 
sh_file, 
sep = "\n")

################################################################################
# Convert Gencode gene file to BED

source(here("R/convert_gtf_to_bed.R"))

bed_files <- gtf_to_bed_converter(
  gtf_url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gtf.gz",
  feature_type = c("gene", "transcript", "start_codon", "stop_codon", "exon"),
  output_dir = here("data/gencode_v39/")
)

################################################################################
# Annotate VCF files: gencode and dbSNP

dir.create(here("data", "vcf_annotated"), showWarnings = FALSE)

annotate_snpeff = function(chrom, dbsnp_file)
{
  message(chrom)
  
  vcf_in   = paste0("/mmfs1/gscratch/togo/matteo/projects/TODAY_genetics/data/genotyping/imputation/chr", chrom, "_imputed.vcf.gz")
  vcf_out  = here("data", "vcf_annotated", paste0("chr", chrom, "_annotated.vcf.gz"))
  command1 = paste("snpEff", "ann", 
                   "-Xmx400g", "-Xms100g", 
                   "-dataDir", "/mmfs1/gscratch/togo/matteo/software/snpeff_data", "GRCh38_gencode39",
                   vcf_in,
                   "|", 
                   "bgzip", "-c",
                   ">", vcf_out)
  
  command2 = paste("bcftools", "index", vcf_out)

  system(command1, intern = TRUE)
  system(command2, intern = TRUE)

  print(vcf_out)
}

lapply(22:1, function(chrom){annotate_snpeff(chrom, dbsnp_file)})
  

################################################################################
# Get a correspondence between ID and rsID

dbsnp_file = "/mmfs1/gscratch/togo/matteo/projects/data/dbsnp/All_20180418.vcf.gz"
command    = paste("bcftools", "query", "-H",
                   #"-r", "1:1-1000000",
                   "-f", "'%CHROM\\t%POS\\t%REF\\t%ALT\\t%ID\\t%VC\\t%TOPMED\\n'",
                   dbsnp_file
                   #"| head -n 10000"
                   )

id2rsid = as.data.frame(fread(cmd = command, sep = "\t", header = TRUE, data.table = FALSE) %>% 
                          setNames(c("chrom", "pos", "ref", "alt", "rsid", "var_type", "topmed")) %>%
                          mutate(topmed = ifelse(topmed == ".", NA, sub("^[^,]*,", "", topmed)))  %>%
                          separate_rows(alt, topmed, sep = ",")) %>%
  mutate(topmed = ifelse(is.na(topmed) == FALSE & topmed == ".", NA, topmed),
         ID = paste(chrom, pos, ref, alt, sep = ":")
         ) %>%
  select(ID, rsid, var_type, topmed)


fwrite(id2rsid, here("data", "id2rsid.txt"), sep = "\t", col.names = TRUE)

