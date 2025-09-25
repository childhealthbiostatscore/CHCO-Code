source("R/functions.R")

# Import files to the data folder
dir.create(here("data/input"), showWarnings = FALSE)

# For now, import all important files here manually

# Expected File Formats:
# Phenotype file: Sample ID in first column, phenotypes in subsequent columns
# ID conversion: Two columns (genotype_id, phenotype_id)
# PCA file: Sample ID + PC1, PC2, PC3, etc.
# Covariate file: Sample ID + your specified covariate columns


# Step 1: clean the ID conversion (remove the array ID)
id_conversion = read.xlsx(here("data/input", "TODAY Genetics ID Linkage File.xlsx"), sheet = 1, detectDates = TRUE) %>%
  setNames(c("phenotype_id", "genotype_id", "array_id")) %>%
  select(genotype_id, phenotype_id) %>%
  drop_na()

# Step 2: clean the genotype PCA data
genotype_pcs_file = here("data/input", "genotype_pca.txt")
suppressWarnings(file.link("/mmfs1/gscratch/togo/matteo/projects/TODAY_genetics/output/genotype_pca/target_projected.txt", genotype_pcs_file))
genotype_pcs = fread(here("data/input", "genotype_pca.txt"), sep = "\t", header = TRUE, data.table = FALSE) %>% 
  select(-FID) %>%
  rename(genotype_id = IID)

# Step 3: fix phenotypes

load(here("data/input", "analytes.Rdata"))
load(here("data/input", "soma.Rdata"    ))
load(here("data/input", "plasma.Rdata"))
load(here("data/input", "urine.Rdata" ))

proteomics_analytes = analytes
proteomics_soma     = soma   %>% filter(visnum == 0) %>% select(-visnum)
metabolomics_plasma = plasma %>% filter(visnum == 0) %>% select(-visnum)
metabolomics_urine  = urine  %>% filter(visnum == 0) %>% select(-visnum)

metabolomics_remove = c("current_label" ,"Freezerworks.ID", "Sample.Name", "site", "material_type", "bsi_id", "Date.Drawn", "location", "MASK.ID", "SAMPLE_ID")
proteomics_remove   = setdiff(colnames(proteomics_soma), c("releaseid", "visnum", proteomics_analytes$AptName))

metabolomics_manifest = read.xlsx(here("data/input", "ElghormliLaure_190032_Manifest.xlsx"), sheet = "ElghormliLaure_190032_Manifest", detectDates = TRUE)
proteomics_manifest   = fread    (here("data/input", "samplesheet_WUS-22-001_updated.txt" ), sep = "\t", header = TRUE, data.table = FALSE)

# Clean manifests

metabolomics_manifest = metabolomics_manifest %>%
  setNames(c("releaseid", "visnum", "sample_type", "sample_id", "mask_id", "date_drawn", "location", "box", "position"))

proteomics_manifest = proteomics_manifest %>%
  mutate(mask_id = Optional1) %>%
  mutate(
    releaseid = case_when(
      mask_id == "EDTA" ~ NA_character_,
      TRUE ~ str_extract(mask_id, "^.*(?=-[^-]+$)")
    ),
    visnum = case_when(
      mask_id == "EDTA" ~ NA_character_,
      TRUE ~ str_extract(mask_id, "[^-]+$")
    )
  )

phenotypes_proteomics = proteomics_soma %>%
  select(all_of(c("releaseid", proteomics_analytes$AptName))) %>%
  rename(sample_id = releaseid)

phenotypes_metabolomics_plasma = metabolomics_plasma %>% select(-all_of(metabolomics_remove))
phenotypes_metabolomics_plasma = phenotypes_metabolomics_plasma %>% 
  select(all_of(c("releaseid", colnames(phenotypes_metabolomics_plasma)[1:(ncol(phenotypes_metabolomics_plasma) - 1)]))) %>%
  rename(sample_id = releaseid)

phenotypes_metabolomics_urine = metabolomics_urine %>% select(-all_of(metabolomics_remove))
phenotypes_metabolomics_urine = phenotypes_metabolomics_urine %>% 
  select(all_of(c("releaseid", colnames(phenotypes_metabolomics_urine)[1:(ncol(phenotypes_metabolomics_urine) - 1)]))) %>%
  rename(sample_id = releaseid)

# Covariates:
# - Standard: age, sex, self-reported ancestry
# - additional: 
#   - BMI, SBP, DBP (from baseline); 
#   - LDL, HDL, VLDL, Chol, Glucose, HbA1c (from CBL); 
#   - CHOL, DFIB, FAT, GLUC, KCAL, MFA, OMEGA3X, PRO, TCHO, TSUGAR, PERCFAT, PERCCARB, PERCPRO, PERCSFA, PERCMFA, PERCPFA (from FFQ)

clinical_baseline = as.data.frame(fread(here("data/input/BASELINE.csv"), sep = ",", header = TRUE, data.table = FALSE) %>% 
                                    group_by(releaseid) %>%
                                    slice(which.min(days))) %>%
  select(releaseid, bmi, sbp, dbp)

clinical_cbl      = as.data.frame(fread(here("data/input/CBL.csv"), sep = ",", header = TRUE, data.table = FALSE) %>% 
                                    group_by(releaseid) %>%
                                    slice(which.min(abs(days)))) %>%
  select(releaseid, LDL, HDL, VLDL, Chol, Glucose, HbA1c)

clinical_ffq      = as.data.frame(fread(here("data/input/FFQ.csv"), sep = ",", header = TRUE, data.table = FALSE) %>% 
                                    group_by(releaseid) %>%
                                    slice(which.min(abs(days)))) %>% 
  select(releaseid, CHOL, DFIB, FAT, GLUC, KCAL, MFA, OMEGA3X, PRO, TCHO, TSUGAR, PERCFAT, PERCCARB, PERCPRO, PERCSFA, PERCMFA, PERCPFA)

clinical_pat      = as.data.frame(fread(here("data/input/PAT.csv"), sep = ",", header = TRUE, data.table = FALSE) %>% 
                                    group_by(releaseid) %>%
                                    slice(which.min(abs(days)))) %>%
  select(releaseid, sex, age, race) %>% 
  mutate(race = factor(race)) %>%
  tidyr::pivot_wider(
    names_from = race,
    values_from = race,
    names_prefix = "race",
    values_fill = 0,
    values_fn = ~1
  )

clinical_outcome = fread(here("data/input/PRIMOUT.csv"), sep = ",", header = TRUE, data.table = FALSE)

load(here("data/input", "comorb.Rdata"))
clinical_comorb = comorb


covariates = as.data.frame(reduce(list(clinical_pat, clinical_baseline, clinical_cbl, clinical_ffq), left_join, by = "releaseid") %>%
                             rename(sample_id = releaseid) %>%
                             mutate(sex = sex - 1) %>%
                             drop_na())

# Finalize the tables for QTL analysis
ids_all = Reduce(intersect, list((id_conversion %>% filter(genotype_id %in% genotype_pcs$genotype_id))[,"phenotype_id"], 
                                 phenotypes_proteomics$sample_id, 
                                 phenotypes_metabolomics_plasma$sample_id,
                                 phenotypes_metabolomics_urine$sample_id,
                                 covariates$sample_id))


dir.create(here("data/to_qtl"), showWarnings = FALSE)

phenotypes_metabolomics_plasma = phenotypes_metabolomics_plasma %>% rename_with(~ gsub("[()/]", "_", .x))
phenotypes_metabolomics_urine  = phenotypes_metabolomics_urine  %>% rename_with(~ gsub("[()/]", "_", .x))

genotype_pcs_clean = genotype_pcs %>% filter(genotype_id %in% (id_conversion %>% filter(phenotype_id %in% ids_all))[,"genotype_id"])

fwrite(id_conversion                  %>% filter(phenotype_id %in% ids_all), file = here("data/to_qtl", "id_conversion.txt"                 ), sep = "\t", col.names = TRUE, row.names = FALSE)
fwrite(phenotypes_proteomics          %>% filter(sample_id    %in% ids_all), file = here("data/to_qtl", "phenotypes_proteomics.txt"         ), sep = "\t", col.names = TRUE, row.names = FALSE)
fwrite(phenotypes_metabolomics_plasma %>% filter(sample_id    %in% ids_all), file = here("data/to_qtl", "phenotypes_metabolomics_plasma.txt"), sep = "\t", col.names = TRUE, row.names = FALSE)
fwrite(phenotypes_metabolomics_urine  %>% filter(sample_id    %in% ids_all), file = here("data/to_qtl", "phenotypes_metabolomics_urine.txt" ), sep = "\t", col.names = TRUE, row.names = FALSE)
fwrite(covariates                     %>% filter(sample_id    %in% ids_all), file = here("data/to_qtl", "covariates.txt"                    ), sep = "\t", col.names = TRUE, row.names = FALSE)
fwrite(genotype_pcs_clean                                                  , file = here("data/to_qtl", "genotype_pcs.txt"                  ), sep = "\t", col.names = TRUE, row.names = FALSE)

saveRDS(list(clinical = list(cbl      = clinical_cbl,
                             baseline = clinical_baseline,
                             ffq      = clinical_ffq,
                             comorb   = clinical_comorb,
                             outcome  = clinical_outcome
                             ),
             metabolomics = list(plasma = metabolomics_plasma %>% select(-all_of(metabolomics_remove)),
                                 urine  = metabolomics_urine  %>% select(-all_of(metabolomics_remove))
                                 ),
             proteomics = list(analytes = proteomics_analytes,
                               soma     = proteomics_soma %>% select(-all_of(proteomics_remove))
                               )
             ), 
        here("data", "data.rds"))


