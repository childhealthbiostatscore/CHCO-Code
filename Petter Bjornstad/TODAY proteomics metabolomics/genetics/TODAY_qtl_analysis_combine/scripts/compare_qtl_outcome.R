source("R/functions.R")

################################################################################
# Get all input data

sumstats_file                      = here("data/top_snps_description.txt")
covariates_file                    = here("data/covariates.txt")
genotype_pcs_file                  = here("data/genotype_pcs.txt")
id_conversion_file                 = here("data/id_conversion.txt")
phenotype_metabolomics_plasma_file = here("data/phenotypes_metabolomics_plasma.txt")
phenotype_metabolomics_urine_file  = here("data/phenotypes_metabolomics_urine.txt")
phenotype_proteomics_file          = here("data/phenotypes_proteomics.txt")
phenotype_outcomes_file            = here("data/phenotypes_outcomes.txt")
phenotype_ckm_file                 = here("data/phenotypes_ckm.txt")

suppressWarnings(file.link("/mmfs1/gscratch/togo/matteo/projects/TODAY_reports/data/top_snps_description.txt"             , sumstats_file                     ))
suppressWarnings(file.link("/mmfs1/gscratch/togo/matteo/projects/TODAY_qtl/data/to_qtl/covariates.txt"                    , covariates_file                   ))
suppressWarnings(file.link("/mmfs1/gscratch/togo/matteo/projects/TODAY_qtl/data/to_qtl/genotype_pcs.txt"                  , genotype_pcs_file                 ))
suppressWarnings(file.link("/mmfs1/gscratch/togo/matteo/projects/TODAY_qtl/data/to_qtl/id_conversion.txt"                 , id_conversion_file                ))
suppressWarnings(file.link("/mmfs1/gscratch/togo/matteo/projects/TODAY_qtl/data/to_qtl/phenotypes_metabolomics_plasma.txt", phenotype_metabolomics_plasma_file))
suppressWarnings(file.link("/mmfs1/gscratch/togo/matteo/projects/TODAY_qtl/data/to_qtl/phenotypes_metabolomics_urine.txt" , phenotype_metabolomics_urine_file ))
suppressWarnings(file.link("/mmfs1/gscratch/togo/matteo/projects/TODAY_qtl/data/to_qtl/phenotypes_proteomics.txt"         , phenotype_proteomics_file         ))
suppressWarnings(file.link("/mmfs1/gscratch/togo/matteo/projects/TODAY_qtl/data/to_qtl/phenotypes_outcomes.txt"           , phenotype_outcomes_file           ))
suppressWarnings(file.link("/mmfs1/gscratch/togo/matteo/projects/TODAY_qtl/data/to_qtl/phenotypes_ckm.txt"                , phenotype_ckm_file                ))

suppressWarnings(file.symlink("/mmfs1/gscratch/togo/matteo/projects/TODAY_genetics/data/genotyping/imputation/plink2/merged_imputed-merge.pgen", here("data/genotypes.pgen")))
suppressWarnings(file.symlink("/mmfs1/gscratch/togo/matteo/projects/TODAY_genetics/data/genotyping/imputation/plink2/merged_imputed-merge.pvar", here("data/genotypes.pvar")))
suppressWarnings(file.symlink("/mmfs1/gscratch/togo/matteo/projects/TODAY_genetics/data/genotyping/imputation/plink2/merged_imputed-merge.psam", here("data/genotypes.psam")))

sumstats                      =  fread(sumstats_file                      , sep = "\t", header = TRUE, data.table = FALSE)
covariates                    =  fread(covariates_file                    , sep = "\t", header = TRUE, data.table = FALSE)
genotype_pcs                  =  fread(genotype_pcs_file                  , sep = "\t", header = TRUE, data.table = FALSE)
id_conversion                 =  fread(id_conversion_file                 , sep = "\t", header = TRUE, data.table = FALSE)
phenotype_metabolomics_plasma =  fread(phenotype_metabolomics_plasma_file , sep = "\t", header = TRUE, data.table = FALSE)
phenotype_metabolomics_urine  =  fread(phenotype_metabolomics_urine_file  , sep = "\t", header = TRUE, data.table = FALSE)
phenotype_proteomics          =  fread(phenotype_proteomics_file          , sep = "\t", header = TRUE, data.table = FALSE)
phenotype_outcomes            =  fread(phenotype_outcomes_file            , sep = "\t", header = TRUE, data.table = FALSE)
phenotype_ckm                 =  fread(phenotype_ckm_file                 , sep = "\t", header = TRUE, data.table = FALSE)

phenotype_outcomes = merge(phenotype_outcomes, phenotype_ckm, by = "sample_id")

################################################################################
# Prepare data for PLINK
# Convert IDs to genotype IDs

covariates_to_qtl = merge(id_conversion, covariates, by.x = "phenotype_id", by.y = "sample_id")
covariates_to_qtl = merge(covariates_to_qtl, genotype_pcs) %>% 
  select(all_of(c("genotype_id", "age", "sex", paste0("race", 1:3), paste0("PC", 1:10)))) %>%
  mutate(across(paste0("PC", 1:10), scale)) %>%
  mutate(sex = sex + 1,
         race1 = race1 + 1,
         race2 = race2 + 1,
         race3 = race3 + 1
  ) %>%
  rename(IID = genotype_id)

phenotpye_outcomes_to_qtl = merge(id_conversion, phenotype_outcomes, by.x = "phenotype_id", by.y = "sample_id") %>%
  select(-phenotype_id) %>%
  rename(IID = genotype_id)


phenotypes_bin = as.matrix(phenotpye_outcomes_to_qtl[,2:ncol(phenotpye_outcomes_to_qtl)])

phenotypes_bin[is.na(phenotypes_bin)] = 0
phenotypes_bin                        = ifelse(phenotypes_bin >0, yes = 2, no = 1)

phenotypes_new     = as.data.frame(phenotypes_bin) %>% select(where(~ is.numeric(.) && sd(.) != 0))
phenotypes_new$IID = phenotpye_outcomes_to_qtl$IID
phenotypes_new     = phenotypes_new[,c("IID", setdiff(colnames(phenotypes_new), c("IID")))]

dir.create(here("data/to_qtl"), showWarnings = FALSE)

covariates_to_qtl_file = here("data/to_qtl", "covariates.txt")
phenotypes_to_qtl_file = here("data/to_qtl", "phenotypes.txt")

fwrite(covariates_to_qtl, covariates_to_qtl_file, sep = "\t", col.names = TRUE, row.names = FALSE)
fwrite(phenotypes_new   , phenotypes_to_qtl_file, sep = "\t", col.names = TRUE, row.names = FALSE)



################################################################################
# Run PLINK

dir.create(here("output/outcomes"), showWarnings = FALSE)

maf_threshold = 0.01
hwe_threshold = 1e-6
plink_output  = here("output/outcomes/sumstats")

done = sub("sumstats.", "", sub(".glm.logistic.hybrid", "", list.files(here("output/outcomes"))))
to_run = setdiff(colnames(phenotypes_new), c("IID", done))

command = paste("plink2",
                "--double-id",
                "--pfile", here("data/genotypes"),
                "--geno" , 0.05,
                "--mind" , 0.05,
                "--maf"  , maf_threshold,
                "--hwe"  , hwe_threshold,
                "--pheno", phenotypes_to_qtl_file,
                "--covar", covariates_to_qtl_file,
                "--glm"  , "firth-fallback", "hide-covar",
                "--pheno-name", paste(to_run, collapse = ","),
                "--out"  , plink_output,
                "--threads", 100,
                "--memory" , 1000000,
                ""
                )


if(length(to_run) > 0)
{
  system(command)
}

################################################################################

################################################################################
# for each phenotype and each outcome, create a model for: 
# 1. just the genotype
# 2. just the phenotype
# 3. both

sumstats_full = sumstats %>% arrange(P) %>% filter(P < 5e-8)
phenotypes    = unique(sumstats_full %>% select(phenotype, analysis))

#sumstats_full %>% filter(phenotype == "seq.3309.2")

phenotype_data = list(metabolomics_plasma.standard_covariates = phenotype_metabolomics_plasma,
                      metabolomics_urine.standard_covariates  = phenotype_metabolomics_urine ,
                      proteomics.standard_covariates          = phenotype_proteomics
                      )

phenotpye_outcomes_to_qtl <- phenotpye_outcomes_to_qtl %>% select(where(~ any(. != 0)))


tmp_folder   = here("output/tmp")

dir.create(here("output/models"), showWarnings = FALSE)
dir.create(tmp_folder           , showWarnings = FALSE)

source(here("R/test_qtl_vs_outcome.R"))

ii = 1
find_variants(ii, phenotypes, sumstats_full, id_conversion, covariates_to_qtl, phenotype_data, phenotpye_outcomes_to_qtl, tmp_folder)
#a = readRDS(here("output/models/proteomics.standard_covariates.seq.3309.2.rds"))

invisible(lapply(1:nrow(phenotypes), function(ii)
{
  try(find_variants(ii, phenotypes, sumstats_full, id_conversion, covariates_to_qtl, phenotype_data, phenotpye_outcomes_to_qtl, tmp_folder))
}))



#phenotypes %>% filter(phenotype %in% c("seq.10855.55", "seq.9216.100"))
