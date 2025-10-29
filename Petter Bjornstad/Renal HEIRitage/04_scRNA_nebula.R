library(aws.s3)
library(jsonlite)
library(biomaRt)

# specify user information
user <- Sys.info()[["user"]]

if (user == "choiyej") { # local version
  root_path <- "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive"
  git_path <- "/Users/choiyej/GitHub/CHCO-Code/Petter Bjornstad"
} else if (user == "rameshsh") { # hyak version
  root_path <- ""
  git_path <- "/mmfs1/gscratch/togo/rameshsh/CHCO-Code/Petter Bjornstad"
} else if (user == "yejichoi") { # hyak version
  root_path <- "/mmfs1/gscratch/togo/yejichoi/"
  git_path <- "/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad"
} else if (user == "pylell") {
  root_path <- "/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive"
  git_path <- "/Users/pylell/Documents/GitHub/CHCO-Code/Petter Bjornstad"
} else {
  stop("Unknown user: please specify root path for this user.")
}

# Set up environment for Kopah
keys <- fromJSON("/mmfs1/home/yejichoi/keys.json")

Sys.setenv(
  "AWS_ACCESS_KEY_ID" = keys$MY_ACCESS_KEY,
  "AWS_SECRET_ACCESS_KEY" = keys$MY_SECRET_KEY,
  "AWS_DEFAULT_REGION" = "",
  "AWS_REGION" = "",
  "AWS_S3_ENDPOINT" = "s3.kopah.uw.edu"
)

# Pull in necessary datasets
pb90_subset_clean <- s3readRDS(object = "data_clean/subset/pb90_ckd_analysis_subset.rds", bucket = "scrna", region = "")
rh_rh2_croc_improve_unique <- 
  s3readRDS(object = "data_clean/rh_rh2_croc_improve_unique.RDS", 
            bucket = "harmonized.dataset",
            region = "")

source(file.path(git_path, "Renal HEIRitage/RH_RH2_IMPROVE_scRNA_functions.R"))
source(file.path(git_path, "Renal HEIRitage/RH_RH2_IMPROVE_functions.R"))
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")


# DKD 100
# set up nebula for DKD vs. non-DKD analysis in T2D/OB only (eGFR <90 | acr_u >= 100) (7 DKD vs. 30 non_DKD)
for (cell in names(celltype_groups)) {
  gc()
  pb90_celltype <- subset(pb90_subset_clean, 
                          KPMP_celltype_general == cell & 
                            group != "Lean_Control" & 
                            !is.na(pb90_subset_clean@meta.data$dkd_group_100))
  print(cell)
  table(pb90_celltype@meta.data %>%
          distinct(record_id, .keep_all = T) %>%
          pull(dkd_group_100))
  
  nebula_res <- run_nebula_parallel(pb90_celltype,
                                    n_cores = 15,
                                    subject_id_col = "record_id",
                                    formula = ~ dkd_group_100,
                                    s3_bucket = "scrna",
                                    s3_key = paste0("Projects/CKD/RH_RH2/Results/nebula/DKD_vs_nonDKD_100/", cell, "/", cell, "_rh_rh2_imp_nebula_kpmp_dkd100.rds"),
                                    group = F)
  
  processed <- process_nebula_results(nebula_res$results,
                                      pval_col = "p_dkd_group_100DKD")
  
  # Get gene symbols
  gene_symbols <- processed$results$Gene
  
  # Query Ensembl
  gene_info <- getBM(
    attributes = c("hgnc_symbol", "description", "gene_biotype"),
    filters = "hgnc_symbol",
    values = gene_symbols,
    mart = mart
  ) %>%
    dplyr::rename(Gene = hgnc_symbol)
  
  # Join annotation
  annotated_df <- processed$results %>%
    left_join(gene_info, by = "Gene")
  
  # Assign to variable dynamically
  var_name <- paste0(tolower(cell), "_kpmp_dkd_nondkd_100")
  assign(var_name, annotated_df, envir = .GlobalEnv)
  
  s3saveRDS(annotated_df, object = paste0("Projects/CKD/RH_RH2/Results/nebula/DKD_vs_nonDKD_100/", cell, "/", cell, "_rh_rh2_imp_nebula_kpmp_dkd100_processed.rds"), 
            bucket = "scrna", region = "")
}

# set up nebula for DKD in T2D/OB vs. HC analysis (eGFR <90 | acr_u >= 100) (7 DKD vs. 12 HC)
for (cell in names(celltype_groups)) {
  gc()
  pb90_celltype <- subset(pb90_subset_clean, 
                          KPMP_celltype_general == cell & 
                            dkd_group_100_hc != "non_DKD" &
                            !is.na(pb90_subset_clean@meta.data$dkd_group_100_hc))
  
  pb90_celltype$dkd_group_100_hc <- factor(pb90_celltype$dkd_group_100_hc)
  pb90_celltype$dkd_group_100_hc <- relevel(pb90_celltype$dkd_group_100_hc, ref = "HC")
  levels(pb90_celltype$dkd_group_100_hc)
  
  print(cell)
  table(pb90_celltype@meta.data %>%
          distinct(record_id, .keep_all = T) %>%
          pull(dkd_group_100_hc))
  
  nebula_res <- run_nebula_parallel(pb90_celltype,
                                    n_cores = 15,
                                    subject_id_col = "record_id",
                                    formula = ~ dkd_group_100_hc,
                                    s3_bucket = "scrna",
                                    s3_key = paste0("Projects/CKD/RH_RH2/Results/nebula/DKD_100_vs_hc/", cell, "/", cell, "_rh_rh2_imp_nebula_kpmp_dkd100_hc.rds"),
                                    group = F)
  
  processed <- process_nebula_results(nebula_res$results,
                                      pval_col = "p_dkd_group_100_hcDKD")
  
  # Get gene symbols
  gene_symbols <- processed$results$Gene
  
  # Query Ensembl
  gene_info <- getBM(
    attributes = c("hgnc_symbol", "description", "gene_biotype"),
    filters = "hgnc_symbol",
    values = gene_symbols,
    mart = mart
  ) %>%
    dplyr::rename(Gene = hgnc_symbol)
  
  # Join annotation
  annotated_df <- processed$results %>%
    left_join(gene_info, by = "Gene")
  
  # Assign to variable dynamically
  var_name <- paste0(tolower(cell), "_kpmp_dkd_nondkd_100_hc")
  assign(var_name, annotated_df, envir = .GlobalEnv)
  
  s3saveRDS(annotated_df, object = paste0("Projects/CKD/RH_RH2/Results/nebula/DKD_100_vs_hc/", cell, "/", cell, "_rh_rh2_imp_nebula_kpmp_dkd100_hc_processed.rds"), 
            bucket = "scrna", region = "")
}

# set up nebula for non-DKD in T2D/OB vs. HC analysis (eGFR <90 | acr_u >= 100) (30 non-DKD vs. 12 HC)
for (cell in names(celltype_groups)) {
  gc()
  
  pb90_celltype <- subset(pb90_subset_clean, 
                          KPMP_celltype_general == cell & 
                            dkd_group_100_hc != "DKD" &
                            !is.na(pb90_subset_clean@meta.data$dkd_group_100_hc))
  
  pb90_celltype$dkd_group_100_hc <- factor(pb90_celltype$dkd_group_100_hc)
  pb90_celltype$dkd_group_100_hc <- relevel(pb90_celltype$dkd_group_100_hc, ref = "HC")
  levels(pb90_celltype$dkd_group_100_hc)
  
  print(cell)
  table(pb90_celltype@meta.data %>%
          distinct(record_id, .keep_all = T) %>%
          pull(dkd_group_100_hc))
  
  nebula_res <- run_nebula_parallel(pb90_celltype,
                                    n_cores = 15,
                                    subject_id_col = "record_id",
                                    formula = ~ dkd_group_100_hc,
                                    s3_bucket = "scrna",
                                    s3_key = paste0("Projects/CKD/RH_RH2/Results/nebula/nonDKD_100_vs_hc/", cell, "/", cell, "_rh_rh2_imp_nebula_kpmp_nondkd100_hc.rds"),
                                    group = F)
  
  processed <- process_nebula_results(nebula_res$results,
                                      pval_col = "p_dkd_group_100_hcnon_DKD")
  
  # Get gene symbols
  gene_symbols <- processed$results$Gene
  
  # Query Ensembl
  gene_info <- getBM(
    attributes = c("hgnc_symbol", "description", "gene_biotype"),
    filters = "hgnc_symbol",
    values = gene_symbols,
    mart = mart
  ) %>%
    dplyr::rename(Gene = hgnc_symbol)
  
  # Join annotation
  annotated_df <- processed$results %>%
    left_join(gene_info, by = "Gene")
  
  # Assign to variable dynamically
  var_name <- paste0(tolower(cell), "_kpmp_nondkd_100_hc")
  assign(var_name, annotated_df, envir = .GlobalEnv)
  
  s3saveRDS(annotated_df, object = paste0("Projects/CKD/RH_RH2/Results/nebula/nonDKD_100_vs_hc/", cell, "/", cell, "_rh_rh2_imp_nebula_kpmp_nondkd100_hc_processed.rds"), 
            bucket = "scrna", region = "")
}


# DKD 30

# set up nebula for DKD vs. non-DKD analysis in T2D/OB only (eGFR <90 | acr_u >= 30) (11 DKD vs. 26 non_DKD)
for (cell in names(celltype_groups)) {
  gc()
  
  pb90_celltype <- subset(pb90_subset_clean, 
                          KPMP_celltype_general == cell & 
                            group != "Lean_Control" & 
                            !is.na(pb90_subset_clean@meta.data$dkd_group_30))
  print(cell)
  table(pb90_celltype@meta.data %>%
          distinct(record_id, .keep_all = T) %>%
          pull(dkd_group_30))
  
  nebula_res <- run_nebula_parallel(pb90_celltype,
                                    n_cores = 15,
                                    subject_id_col = "record_id",
                                    formula = ~ dkd_group_30,
                                    s3_bucket = "scrna",
                                    s3_key = paste0("Projects/CKD/RH_RH2/Results/nebula/DKD_vs_nonDKD_30/", cell, "/", cell, "_rh_rh2_imp_nebula_kpmp_dkd30.rds"),
                                    group = F)
  
  processed <- process_nebula_results(nebula_res$results,
                                      pval_col = "p_dkd_group_30DKD")
  
  # Get gene symbols
  gene_symbols <- processed$results$Gene
  
  # Query Ensembl
  gene_info <- getBM(
    attributes = c("hgnc_symbol", "description", "gene_biotype"),
    filters = "hgnc_symbol",
    values = gene_symbols,
    mart = mart
  ) %>%
    dplyr::rename(Gene = hgnc_symbol)
  
  # Join annotation
  annotated_df <- processed$results %>%
    left_join(gene_info, by = "Gene")
  
  # Assign to variable dynamically
  var_name <- paste0(tolower(cell), "_kpmp_dkd_nondkd_30")
  assign(var_name, annotated_df, envir = .GlobalEnv)
  
  s3saveRDS(annotated_df, object = paste0("Projects/CKD/RH_RH2/Results/nebula/DKD_vs_nonDKD_30/", cell, "/", cell, "_rh_rh2_imp_nebula_kpmp_dkd30_processed.rds"), 
            bucket = "scrna", region = "")
}

# set up nebula for DKD in T2D/OB vs. HC analysis (eGFR <90 | acr_u >= 30) (11 DKD vs. 12 HC)
for (cell in names(celltype_groups)) {
  gc()
  
  pb90_celltype <- subset(pb90_subset_clean, 
                          KPMP_celltype_general == cell & 
                            dkd_group_30_hc != "non_DKD" &
                            !is.na(pb90_subset_clean@meta.data$dkd_group_30_hc))
  
  pb90_celltype$dkd_group_30_hc <- factor(pb90_celltype$dkd_group_30_hc)
  pb90_celltype$dkd_group_30_hc <- relevel(pb90_celltype$dkd_group_30_hc, ref = "HC")
  levels(pb90_celltype$dkd_group_30_hc)
  
  print(cell)
  table(pb90_celltype@meta.data %>%
          distinct(record_id, .keep_all = T) %>%
          pull(dkd_group_30_hc))
  
  nebula_res <- run_nebula_parallel(pb90_celltype,
                                    n_cores = 15,
                                    subject_id_col = "record_id",
                                    formula = ~ dkd_group_30_hc,
                                    s3_bucket = "scrna",
                                    s3_key = paste0("Projects/CKD/RH_RH2/Results/nebula/DKD_30_vs_hc/", cell, "/", cell, "_rh_rh2_imp_nebula_kpmp_dkd30_hc.rds"),
                                    group = F)
  
  processed <- process_nebula_results(nebula_res$results,
                                      pval_col = "p_dkd_group_30_hcDKD")
  
  # Get gene symbols
  gene_symbols <- processed$results$Gene
  
  # Query Ensembl
  gene_info <- getBM(
    attributes = c("hgnc_symbol", "description", "gene_biotype"),
    filters = "hgnc_symbol",
    values = gene_symbols,
    mart = mart
  ) %>%
    dplyr::rename(Gene = hgnc_symbol)
  
  # Join annotation
  annotated_df <- processed$results %>%
    left_join(gene_info, by = "Gene")
  
  # Assign to variable dynamically
  var_name <- paste0(tolower(cell), "_kpmp_dkd_nondkd_30_hc")
  assign(var_name, annotated_df, envir = .GlobalEnv)
  
  s3saveRDS(annotated_df, object = paste0("Projects/CKD/RH_RH2/Results/nebula/DKD_30_vs_hc/", cell, "/", cell, "_rh_rh2_imp_nebula_kpmp_dkd30_hc_processed.rds"), 
            bucket = "scrna", region = "")
}

# set up nebula for non-DKD in T2D/OB vs. HC analysis (eGFR <90 | acr_u >= 30) (26 non-DKD vs. 12 HC)
for (cell in names(celltype_groups)) {
  gc()
  
  pb90_celltype <- subset(pb90_subset_clean, 
                          KPMP_celltype_general == cell & 
                            dkd_group_30_hc != "DKD" &
                            !is.na(pb90_subset_clean@meta.data$dkd_group_30_hc))
  
  pb90_celltype$dkd_group_30_hc <- factor(pb90_celltype$dkd_group_30_hc)
  pb90_celltype$dkd_group_30_hc <- relevel(pb90_celltype$dkd_group_30_hc, ref = "HC")
  levels(pb90_celltype$dkd_group_30_hc)
  
  print(cell)
  table(pb90_celltype@meta.data %>%
          distinct(record_id, .keep_all = T) %>%
          pull(dkd_group_30_hc))
  
  nebula_res <- run_nebula_parallel(pb90_celltype,
                                    n_cores = 15,
                                    subject_id_col = "record_id",
                                    formula = ~ dkd_group_30_hc,
                                    s3_bucket = "scrna",
                                    s3_key = paste0("Projects/CKD/RH_RH2/Results/nebula/nonDKD_30_vs_hc/", cell, "/", cell, "_rh_rh2_imp_nebula_kpmp_nondkd30_hc.rds"),
                                    group = F)
  
  processed <- process_nebula_results(nebula_res$results,
                                      pval_col = "p_dkd_group_30_hcnon_DKD")
  
  # Get gene symbols
  gene_symbols <- processed$results$Gene
  
  # Query Ensembl
  gene_info <- getBM(
    attributes = c("hgnc_symbol", "description", "gene_biotype"),
    filters = "hgnc_symbol",
    values = gene_symbols,
    mart = mart
  ) %>%
    dplyr::rename(Gene = hgnc_symbol)
  
  # Join annotation
  annotated_df <- processed$results %>%
    left_join(gene_info, by = "Gene")
  
  # Assign to variable dynamically
  var_name <- paste0(tolower(cell), "_kpmp_nondkd_30_hc")
  assign(var_name, annotated_df, envir = .GlobalEnv)
  
  s3saveRDS(annotated_df, object = paste0("Projects/CKD/RH_RH2/Results/nebula/nonDKD_30_vs_hc/", cell, "/", cell, "_rh_rh2_imp_nebula_kpmp_nondkd30_hc_processed.rds"), 
            bucket = "scrna", region = "")
}

# GLP1- vs. HC

# set up nebula for GLP1- (T2D/OB) vs. Healthy Control analysis (24 GLP- vs. 12 HC)
for (cell in names(celltype_groups)) {
  gc()
  
  pb90_celltype <- subset(pb90_subset_clean, 
                          KPMP_celltype_general == cell & 
                            glp_t2dob != "GLP_Y" & 
                            !is.na(pb90_subset_clean@meta.data$glp_t2dob))
  
  pb90_celltype$glp_t2dob <- factor(pb90_celltype$glp_t2dob)
  pb90_celltype$glp_t2dob <- relevel(pb90_celltype$glp_t2dob, ref = "HC")
  levels(pb90_celltype$glp_t2dob)
  
  print(cell)
  table(pb90_celltype@meta.data %>%
          distinct(record_id, .keep_all = T) %>%
          pull(glp_t2dob))
  
  nebula_res <- run_nebula_parallel(pb90_celltype,
                                    n_cores = 15,
                                    subject_id_col = "record_id",
                                    formula = ~ glp_t2dob,
                                    s3_bucket = "scrna",
                                    s3_key = paste0("Projects/CKD/RH_RH2/Results/nebula/GLP_N_vs_HC/", cell, "/", cell, "_rh_rh2_imp_nebula_kpmp_glpn_hc.rds"),
                                    group = F)
  
  processed <- process_nebula_results(nebula_res$results,
                                      pval_col = "p_glp_t2dobGLP_N")
  
  # Get gene symbols
  gene_symbols <- processed$results$Gene
  
  # Query Ensembl
  gene_info <- getBM(
    attributes = c("hgnc_symbol", "description", "gene_biotype"),
    filters = "hgnc_symbol",
    values = gene_symbols,
    mart = mart
  ) %>%
    dplyr::rename(Gene = hgnc_symbol)
  
  # Join annotation
  annotated_df <- processed$results %>%
    left_join(gene_info, by = "Gene")
  
  # Assign to variable dynamically
  var_name <- paste0(tolower(cell), "_kpmp_glpn_hc")
  assign(var_name, annotated_df, envir = .GlobalEnv)
  
  s3saveRDS(annotated_df, object = paste0("Projects/CKD/RH_RH2/Results/nebula/GLP_N_vs_HC/", cell, "/", cell, "_rh_rh2_imp_nebula_kpmp_glpn_hc_processed.rds"), 
            bucket = "scrna", region = "")
}

# GLP- vs. GLP+
# set up nebula for GLP1+ (T2D/OB) vs. GLP1- (T2D/OB) analysis (24 GLP- vs. 13 GLP+)
for (cell in names(celltype_groups)) {
  gc()
  
  pb90_celltype <- subset(pb90_subset_clean, 
                          KPMP_celltype_general == cell & 
                            group != "Lean_Control" & 
                            !is.na(pb90_subset_clean@meta.data$glp_t2dob))
  
  pb90_celltype$glp_t2dob <- factor(pb90_celltype$glp_t2dob)
  pb90_celltype$glp_t2dob <- relevel(pb90_celltype$glp_t2dob, ref = "GLP_N")
  levels(pb90_celltype$glp_t2dob)
  
  print(cell)
  
  table(pb90_celltype@meta.data %>%
          distinct(record_id, .keep_all = T) %>%
          pull(glp_t2dob))
  
  nebula_res <- run_nebula_parallel(pb90_celltype,
                                    n_cores = 15,
                                    subject_id_col = "record_id",
                                    formula = ~ glp_t2dob,
                                    s3_bucket = "scrna",
                                    s3_key = paste0("Projects/CKD/RH_RH2/Results/nebula/GLP_Y_vs_GLP_N/", cell, "/", cell, "_rh_rh2_imp_nebula_kpmp_glpy_glpn.rds"),
                                    group = F) 
  
  processed <- process_nebula_results(nebula_res$results,
                                      pval_col = "p_glp_t2dobGLP_Y")
  
  # Get gene symbols
  gene_symbols <- processed$results$Gene
  
  # Query Ensembl
  gene_info <- getBM(
    attributes = c("hgnc_symbol", "description", "gene_biotype"),
    filters = "hgnc_symbol",
    values = gene_symbols,
    mart = mart
  ) %>%
    dplyr::rename(Gene = hgnc_symbol)
  
  # Join annotation
  annotated_df <- processed$results %>%
    left_join(gene_info, by = "Gene")
  
  # Assign to variable dynamically
  var_name <- paste0(tolower(cell), "_kpmp_glpy_glpn")
  assign(var_name, annotated_df, envir = .GlobalEnv)
  
  s3saveRDS(annotated_df, object = paste0("Projects/CKD/RH_RH2/Results/nebula/GLP_Y_vs_GLP_N/", cell, "/", cell, "_rh_rh2_imp_nebula_kpmp_glpy_glpn_processed.rds"), 
            bucket = "scrna", region = "")
}

