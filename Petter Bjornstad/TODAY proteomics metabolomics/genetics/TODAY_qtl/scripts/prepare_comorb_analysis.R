source("R/functions.R")

data       = readRDS(here("data", "data.rds"))
comorb     = data$clinical$comorb
outcome    = data$clinical$outcome
comorb     = comorb[,grepl("^DAYSTO", colnames(comorb)) == FALSE]
ids        = fread(here("data/to_qtl", "id_conversion.txt"), sep = "\t", header = TRUE, data.table = FALSE)
phenotypes = merge(outcome %>% select(releaseid, outcome),
                   comorb  
                   ) %>%
  rename(sample_id = releaseid) %>% 
  filter(sample_id %in% ids$phenotype_id)

phenotypes_bin = as.matrix(phenotypes[,2:ncol(phenotypes)])

phenotypes_bin[is.na(phenotypes_bin)] = 0

phenotypes_new           = as.data.frame(phenotypes_bin + 1)
phenotypes_new$sample_id = phenotypes$sample_id

fwrite(phenotypes[,colnames(phenotypes)], file = here("data/to_qtl", "phenotypes_outcomes.txt" ), sep = "\t", col.names = TRUE, row.names = FALSE)

################################################################################
# Add CKM information: copy from C:\Users\s0416\OneDrive - UW\Laura Pyle's files - Biostatistics Core Shared Drive\TODAY subaward\CKM\Results

library(fastDummies)

ckm = fread(here("data", "input/CKM stage summary.csv"), sep = ",", header = TRUE, data.table = FALSE, 
            select = c("RELEASEID", "fup_time", "CKM_syn_base", "CKM_max", "progress_CKM", 
                       "days_to_first_ckm_prog", "num_prog", "CKM_syn_first_progression", 
                       "CKM_syn_num_first_progression", "CKM_max_char")) %>% 
  rename(sample_id = RELEASEID) %>% 
  filter(sample_id %in% ids$phenotype_id)

ckm_bin            = dummy_cols(ckm %>% select(CKM_syn_base, CKM_syn_first_progression, CKM_max_char), remove_first_dummy = FALSE, remove_selected_columns = TRUE)
ckm_full           = cbind(ckm %>% select(sample_id, progress_CKM), ckm_bin)
colnames(ckm_full) = gsub("\\+", "_plus", gsub(" ", "_", colnames(ckm_full)))

fwrite(ckm_full, file = here("data/to_qtl", "phenotypes_ckm.txt" ), sep = "\t", col.names = TRUE, row.names = FALSE)

