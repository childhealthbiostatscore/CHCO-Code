source("R/describe_sumstats_functions.R")

progress = readRDS(here("data/progress.rds"))
data     = readRDS("/mmfs1/gscratch/togo/matteo/projects/TODAY_qtl/data/data.rds")

################################################################################
# Find top SNPs: remember to rerun the qtl_report quarto app to update to all 7k analyses

sumstats_report_file = here("data/top_snps_sumstats.txt")
top_snps = get_all_sumstat_files(progress, sumstats_report_file)

fwrite(top_snps, sumstats_report_file, sep = "\t", col.names = TRUE, row.names = FALSE)

################################################################################
# Add information about proteins
outfile              = here("data", "top_snps_description.txt")
top_snps_description = merge(data$proteomics$analytes %>% 
                               select(AptName, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol) %>%
                               rename(phenotype = AptName),
                             top_snps, 
                             all.y = TRUE
                             ) %>% arrange(P)

fwrite(top_snps_description %>% arrange(P), outfile, sep = "\t", col.names = TRUE, row.names = FALSE)

################################################################################
