suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(dplyr))

################################################################################
# Download eQTL catalogue metadata

metadata_file = here("data", "metadata.tsv") # downloaded from "https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/data_tables/dataset_metadata.tsv"
metadata      = fread(metadata_file, sep = "\t", header = TRUE, data.table = FALSE) %>% 
  filter(quant_method == "ge" & condition_label %in% c("naive", "statin") & sample_size >= 100) %>%
  group_by(tissue_label, condition_label, study_type) %>%
  slice(which.max(sample_size)) %>%
  ungroup() %>%
  mutate(summary_stats_ftp     = paste("https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats", study_id, dataset_id, paste(dataset_id, "all.tsv.gz"    , sep = "."), sep = "/"),
         summary_stats_ftp_tbi = paste("https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats", study_id, dataset_id, paste(dataset_id, "all.tsv.gz.tbi", sep = "."), sep = "/")
         )

metadata = as.data.frame(metadata)

fwrite(metadata, metadata_file, sep = "\t", col.names = TRUE, row.names = FALSE)

################################################################################
# Download eQTL summary statistics
# Create a SH file and download via command line

sh_file = here("scripts/download_eqtl_catalogue.sh")

writeLines(c(paste("mkdir", here("data/eqtl_catalog_summary_stats")),
             paste("cd"   , here("data/eqtl_catalog_summary_stats")),
             "",
             paste("wget", metadata$summary_stats_ftp    ), 
             paste("wget", metadata$summary_stats_ftp_tbi),
             ), 
           sh_file, 
           sep = "\n")


