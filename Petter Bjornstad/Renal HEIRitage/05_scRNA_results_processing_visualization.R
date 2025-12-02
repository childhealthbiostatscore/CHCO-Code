library(aws.s3)
library(jsonlite)
library(biomaRt)

user <- Sys.info()[["user"]]

if (user == "choiyej") { # local version
  root_path <- "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive"
  git_path <- "/Users/choiyej/GitHub/CHCO-Code/Petter Bjornstad"
  keys <- fromJSON("/Users/choiyej/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Bjornstad Pyle Lab/keys.json")
} else if (user == "rameshsh") { # hyak version
  root_path <- ""
  git_path <- "/mmfs1/gscratch/togo/rameshsh/CHCO-Code/Petter Bjornstad"
} else if (user == "yejichoi") { # hyak version
  root_path <- "/mmfs1/gscratch/togo/yejichoi/"
  git_path <- "/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad"
  keys <- fromJSON("/mmfs1/home/yejichoi/keys.json")
} else if (user == "pylell") {
  root_path <- "/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive"
  git_path <- "/Users/pylell/Documents/GitHub/CHCO-Code/Petter Bjornstad"
  keys <- fromJSON("/mmfs1/home/pylell/keys.json")
} else {
  stop("Unknown user: please specify root path for this user.")
}

source(file.path(git_path, "Renal HEIRitage/RH_RH2_IMPROVE_functions.R"))
source(file.path(git_path, "Renal HEIRitage/RH_RH2_IMPROVE_scRNA_functions.R"))

# Pull in nebula results
folders <- c(
  
  # --- Base DKD vs nonDKD ---
  "DKD_vs_nonDKD_30" = "dkd30",
  "DKD_vs_nonDKD_100" = "dkd100",
  
  # --- DKD vs nonDKD with T2D ---
  "DKD_vs_nonDKD_30_t2d" = "dkd30_t2d",
  "DKD_vs_nonDKD_100_t2d" = "dkd100_t2d",
  
  # --- DKD vs HC (30 / 100) ---
  "DKD_30_vs_HC" = "dkd30_hc",
  "DKD_100_vs_HC" = "dkd100_hc",
  
  # --- DKD vs HC with T2D ---
  "DKD_30_t2d_vs_HC" = "dkd30_t2d_hc",
  "DKD_100_t2d_vs_HC" = "dkd100_t2d_hc",
  
  # --- nonDKD vs HC (30 / 100) ---
  "nonDKD_30_vs_HC" = "nondkd30_hc",
  "nonDKD_100_vs_HC" = "nondkd100_hc",
  
  # --- nonDKD vs HC with T2D ---
  "nonDKD_30_t2d_vs_HC" = "nondkd30_t2d_hc",
  "nonDKD_100_t2d_vs_HC" = "nondkd100_t2d_hc",
  
  # --- GLP (within DKD) ---
  "DKD_30_GLP_Y_vs_DKD_30_GLP_N" = "dkd_30_glpy_glpn",
  "DKD_100_GLP_Y_vs_DKD_100_GLP_N" = "dkd_100_glpy_glpn",
  
  # --- GLP (within nonDKD) ---
  "nonDKD_30_GLP_Y_vs_nonDKD_30_GLP_N" = "nondkd_30_glpy_glpn",
  "nonDKD_100_GLP_Y_vs_nonDKD_100_GLP_N" = "nondkd_100_glpy_glpn",
  
  # --- GLP N/Y vs HC (nonDKD) ---
  "nonDKD_30_GLP_N_vs_HC" = "nondkd30_glpn_hc",
  "nonDKD_30_GLP_Y_vs_HC" = "nondkd30_glpy_hc",
  "nonDKD_100_GLP_N_vs_HC" = "nondkd100_glpn_hc",
  "nonDKD_100_GLP_Y_vs_HC" = "nondkd100_glpy_hc",
  
  # --- GLP N/Y vs HC (DKD) ---
  "DKD_30_GLP_N_vs_HC" = "dkd30_glpn_hc",
  "DKD_30_GLP_Y_vs_HC" = "dkd30_glpy_hc",
  "DKD_100_GLP_N_vs_HC" = "dkd100_glpn_hc",
  "DKD_100_GLP_Y_vs_HC" = "dkd100_glpy_hc",
  
  # --- SGLT2i (within DKD) ---
  # "DKD_30_SGLT2i_Y_vs_DKD_30_SGLT2i_N" = "dkd30_sglt2iy_sglt2in",
  # "DKD_100_SGLT2i_Y_vs_DKD_100_SGLT2i_N" = "dkd100_sglt2iy_sglt2in",
  
  # --- SGLT2i (within nonDKD) ---
  # "nonDKD_30_SGLT2i_Y_vs_nonDKD_30_SGLT2i_N" = "nondkd30_sglt2iy_sglt2in",
  # "nonDKD_100_SGLT2i_Y_vs_nonDKD_100_SGLT2i_N" = "nondkd100_sglt2iy_sglt2in",
  
  # --- SGLT2i N/Y vs HC (nonDKD) ---
  "nonDKD_30_SGLT2i_N_vs_HC" = "nondkd30_sglt2in_hc",
  "nonDKD_30_SGLT2i_Y_vs_HC" = "nondkd30_sglt2iy_hc",
  "nonDKD_100_SGLT2i_N_vs_HC" = "nondkd100_sglt2in_hc",
  "nonDKD_100_SGLT2i_Y_vs_HC" = "nondkd100_sglt2iy_hc"
  
  # --- SGLT2i N/Y vs HC (DKD) ---
  # "DKD_30_SGLT2i_N_vs_HC" = "dkd30_sglt2in_hc",
  # "DKD_30_SGLT2i_Y_vs_HC" = "dkd30_sglt2iy_hc",
  # "DKD_100_SGLT2i_N_vs_HC" = "dkd100_sglt2in_hc",
  # "DKD_100_SGLT2i_Y_vs_HC" = "dkd100_sglt2iy_hc"
)

# Pull in nebula results
folders <- c(
  "nonDKD_30_GLP_Y_vs_HC" = "nondkd30_glpy_hc"
)

# Define common parameters
celltype_groups <- list(
  PT = c("PT-S1/S2", "PT-S3", "aPT"),
  aPT = "aPT",
  `PT-S1/S2` = "PT-S1/S2",
  `PT-S3` = "PT-S3",
  `PT-1` = "PT-1",
  `PT-2` = "PT-2",
  `PT-3` = "PT-3",
  `PT-4` = "PT-4",
  `PT-5` = "PT-5",
  TAL = c("C-TAL-1", "C-TAL-2", "aTAL", "dTAL"),
  PC = c("CCD-PC", "CNT-PC", "dCCD-PC", "M-PC", "tPC-IC"),
  IC = c("IC-A", "IC-B", "aIC"),
  DTL_ATL = c("DTL", "aDTL", "ATL"),
  DCT_CNT = c("DCT", "dDCT", "CNT"),
  EC = c("EC-AVR", "EC-GC", "EC-PTC", "EC-AEA", "EC-LYM", "EC/VSMC", "EC-A"),
  Immune = c("MAC", "MON", "cDC", "pDC", "CD4+ T", "CD8+ T", "B", "NK", "cycT"),
  VSMC_P_FIB = c("VSMC/P", "FIB"),
  POD = "POD",
  MC = "MC",
  PEC = "PEC",
  Schwann = "SchwannCells"
)

# Load all run metadata and save as one file
meta_df_compiled <- data.frame(analysis_type = NULL, 
                               celltype = NULL,
                               celltype_variable = NULL,
                               group_variable = NULL,
                               reference_level = NULL,
                               pval_column = NULL,
                               logfc_column = NULL,
                               n_samples_per_group = NULL,
                               n_cells_per_group = NULL, 
                               total_samples = NULL,
                               total_cells = NULL)

for (folder in names(folders)) {
  for (cell in names(celltype_groups)) {
    # Construct the S3 path
    s3_path <- paste0("Projects/CKD/RH_RH2/Results/nebula/", folder, "/", cell, "/", cell, "_rh_rh2_imp_nebula_kpmp_", folders[folder], "_metadata.csv")
    
    # Try to read the file with error handling
    meta_df <- tryCatch({
      s3read_using_region(FUN = read.csv, object = s3_path, bucket = "scrna", region = "")
    }, error = function(e) {
      # Print error message with specific cell and folder information
      message(paste0("ERROR: Failed to read file for folder '", folder, "' and cell type '", cell, "'"))
      message(paste0("  Path attempted: ", s3_path))
      message(paste0("  Error details: ", e$message))
      return(NULL)
    })
    
    # Skip to next iteration if file reading failed
    if (is.null(meta_df)) {
      message(paste0("Skipping folder '", folder, "' and cell type '", cell, "'\n"))
      next
    }
    
    meta_df_compiled <- rbind(meta_df_compiled, meta_df)
  }
}

write.csv(meta_df_compiled, file.path(root_path, paste0("Renal HERITAGE/Results/nebula/csv/run_metadata_compiled.csv")), 
          row.names = F, fileEncoding = "UTF-8")

# Load all data first and save as csv
for (folder in names(folders)) {
  for (cell in names(celltype_groups)) {
    # Construct the S3 path
    s3_path <- paste0("Projects/CKD/RH_RH2/Results/nebula/", folder, "/", cell, "/", cell, "_rh_rh2_imp_nebula_kpmp_", folders[folder], "_processed.rds")
    
    # Try to read the file with error handling
    processed_df <- tryCatch({
      s3readRDS(object = s3_path, bucket = "scrna", region = "")
    }, error = function(e) {
      # Print error message with specific cell and folder information
      message(paste0("ERROR: Failed to read file for folder '", folder, "' and cell type '", cell, "'"))
      message(paste0("  Path attempted: ", s3_path))
      message(paste0("  Error details: ", e$message))
      return(NULL)
    })
    
    # Skip to next iteration if file reading failed
    if (is.null(processed_df)) {
      message(paste0("Skipping folder '", folder, "' and cell type '", cell, "'\n"))
      next
    }
    
    # If successful, continue with processing
    var_name <- paste0(tolower(cell), "_", folders[folder])
    var_name <- gsub("/", "_", var_name)
    assign(var_name, processed_df, envir = .GlobalEnv)
    
    # Write CSV files
    write.csv(processed_df, 
              file.path(root_path, paste0("Renal HERITAGE/Results/nebula/csv/full/", var_name, "_kpmp.csv")), 
              row.names = F, fileEncoding = "UTF-8")
    
    write.csv(subset(processed_df, processed_df[[6]] < 0.05), 
              file.path(root_path, paste0("Renal HERITAGE/Results/nebula/csv/pval/", var_name, "_kpmp_pval.csv")), 
              row.names = F, fileEncoding = "UTF-8")
    
    write.csv(subset(processed_df, fdr < 0.05), 
              file.path(root_path, paste0("Renal HERITAGE/Results/nebula/csv/fdr/", var_name, "_kpmp_fdr.csv")), 
              row.names = F, fileEncoding = "UTF-8")
    
    # Optional: Print success message
    message(paste0("Successfully processed: ", var_name))
  }
}




# Define plot parameters for each comparison
plot_params <- list(
  # Existing parameters
  "dkd100" = list(
    fc_col = "logFC_dkd_group_100DKD",
    p_col = "p_dkd_group_100DKD",
    x_label = "logFC DKD vs. non-DKD\n(DKD: eGFR < 90 or UACR >= 100)",
    positive_text = "Positive with DKD",
    negative_text = "Negative with DKD",
    formula = "DKD"
  ),
  "dkd100_hc" = list(
    fc_col = "logFC_dkd_group_100_hcDKD",
    p_col = "p_dkd_group_100_hcDKD",
    x_label = "logFC DKD vs. HC\n(DKD: eGFR < 90 or UACR >= 100)",
    positive_text = "Positive with DKD",
    negative_text = "Negative with DKD",
    formula = "DKD"
  ),
  "nondkd100_hc" = list(
    fc_col = "logFC_dkd_group_100_hcnon_DKD",
    p_col = "p_dkd_group_100_hcnon_DKD",
    x_label = "logFC non-DKD vs. HC\n(non-DKD: eGFR >= 90 and UACR < 100)",
    positive_text = "Positive with non-DKD",
    negative_text = "Negative with non-DKD",
    formula = "NonDKD"
  ),
  "dkd30" = list(
    fc_col = "logFC_dkd_group_30DKD",
    p_col = "p_dkd_group_30DKD",
    x_label = "logFC DKD vs. non-DKD\n(DKD: eGFR < 90 or UACR >= 30)",
    positive_text = "Positive with DKD",
    negative_text = "Negative with DKD",
    formula = "DKD"
  ),
  "dkd30_hc" = list(
    fc_col = "logFC_dkd_group_30_hcDKD",
    p_col = "p_dkd_group_30_hcDKD",
    x_label = "logFC DKD vs. HC\n(DKD: eGFR < 90 or UACR >= 30)",
    positive_text = "Positive with DKD",
    negative_text = "Negative with DKD",
    formula = "DKD"
  ),
  "nondkd30_hc" = list(
    fc_col = "logFC_dkd_group_30_hcnon_DKD",
    p_col = "p_dkd_group_30_hcnon_DKD",
    x_label = "logFC non-DKD vs. HC\n(non-DKD: eGFR >= 90 and UACR < 30)",
    positive_text = "Positive with non-DKD",
    negative_text = "Negative with non-DKD",
    formula = "nonDKD"
  ),
  "glpn_hc" = list(
    fc_col = "logFC_glp_t2dobGLP_N",
    p_col = "p_glp_t2dobGLP_N",
    x_label = "logFC GLP1-RA No vs. HC",
    positive_text = "Positive with No GLP1-RA",
    negative_text = "Negative with No GLP1-RA",
    formula = "GLP_N"
  ),
  "glpy_glpn" = list(
    fc_col = "logFC_glp_t2dobGLP_Y",
    p_col = "p_glp_t2dobGLP_Y",
    x_label = "logFC GLP1-RA Yes vs. No",
    positive_text = "Positive with GLP1-RA",
    negative_text = "Negative with GLP1-RA",
    formula = "GLP_Y"
  ),
  "nondkd100_glpn_hc" = list(
    fc_col = "logFC_dkd_group_100_hcnon_DKD",
    p_col = "p_dkd_group_100_hcnon_DKD",
    x_label = "logFC non-DKD (GLP1-RA-) vs. HC\n(non-DKD: eGFR >= 90 and UACR < 100)",
    positive_text = "Positive with non-DKD\n(GLP-)",
    negative_text = "Negative with non-DKD\n(GLP-)",
    formula = "nonDKD_GLP_N"
  ),
  "nondkd100_glpy_hc" = list(
    fc_col = "logFC_dkd_group_100_glp_hcnon_DKD_GLP_Y",
    p_col = "p_dkd_group_100_glp_hcnon_DKD_GLP_Y",
    x_label = "logFC non-DKD (GLP1-RA+) vs. HC\n(non-DKD: eGFR >= 90 and UACR < 100)",
    positive_text = "Positive with non-DKD\n(GLP+)",
    negative_text = "Negative with non-DKD\n(GLP+)",
    formula = "nonDKD_GLP_Y"
  ),
  "nondkd30_glpn_hc" = list(
    fc_col = "logFC_dkd_group_30_hcnon_DKD",
    p_col = "p_dkd_group_30_hcnon_DKD",
    x_label = "logFC non-DKD (GLP1-RA-) vs. HC\n(non-DKD: eGFR >= 90 and UACR < 30)",
    positive_text = "Positive with non-DKD\n(GLP-)",
    negative_text = "Negative with non-DKD\n(GLP-)",
    formula = "nonDKD_GLP_N"
  ),
  "nondkd30_glpy_hc" = list(
    fc_col = "logFC_dkd_group_30_hcnon_DKD",
    p_col = "p_dkd_group_30_hcnon_DKD",
    x_label = "logFC non-DKD (GLP1-RA+) vs. HC\n(non-DKD: eGFR >= 90 and UACR < 30)",
    positive_text = "Positive with non-DKD\n(GLP+)",
    negative_text = "Negative with non-DKD\n(GLP+)",
    formula = "nonDKD_GLP_Y"
  ),
  
  # # T2D comparisons (inferred)
  "dkd100_t2d" = list(
    fc_col = "logFC_dkd_group_100DKD",
    p_col = "p_dkd_group_100DKD",
    x_label = "logFC DKD vs. non-DKD (T2D only)\n(DKD: eGFR < 90 or UACR >= 100)",
    positive_text = "Positive with DKD (T2D)",
    negative_text = "Negative with DKD (T2D)",
    formula = "DKD_T2D"
  ),
  "kpmp_dkd100_t2d_hc" = list(
    fc_col = "logFC_dkd_group_100_t2d_hcDKD",
    p_col = "p_dkd_group_100_t2d_hcDKD",
    x_label = "logFC DKD (T2D) vs. HC\n(DKD: eGFR < 90 or UACR >= 100)",
    positive_text = "Positive with DKD (T2D)",
    negative_text = "Negative with DKD (T2D)",
    formula = "DKD_T2D"
  ),
  "nondkd100_t2d_hc" = list(
    fc_col = "logFC_dkd_group_100_hcnon_DKD",
    p_col = "p_dkd_group_100_hcnon_DKD",
    x_label = "logFC non-DKD (T2D) vs. HC\n(non-DKD: eGFR >= 90 and UACR < 100)",
    positive_text = "Positive with non-DKD (T2D)",
    negative_text = "Negative with non-DKD (T2D)",
    formula = "nonDKD_T2D"
  ),
  "dkd30_t2d" = list(
    fc_col = "logFC_dkd_group_30DKD",
    p_col = "p_dkd_group_30DKD",
    x_label = "logFC DKD vs. non-DKD (T2D only)\n(DKD: eGFR < 90 or UACR >= 30)",
    positive_text = "Positive with DKD (T2D)",
    negative_text = "Negative with DKD (T2D)",
    formula = "DKD_T2D"
  ),
  "dkd30_t2d_hc" = list(
    fc_col = "logFC_dkd_group_30_t2d_hcDKD",
    p_col = "p_dkd_group_30_t2d_hcDKD",
    x_label = "logFC DKD (T2D) vs. HC\n(DKD: eGFR < 90 or UACR >= 30)",
    positive_text = "Positive with DKD (T2D)",
    negative_text = "Negative with DKD (T2D)",
    formula = "DKD_T2D"
  ),
  "nondkd30_t2d_hc" = list(
    fc_col = "logFC_dkd_group_30_hcnon_DKD",
    p_col = "p_dkd_group_30_hcnon_DKD",
    x_label = "logFC non-DKD (T2D) vs. HC\n(non-DKD: eGFR >= 90 and UACR < 30)",
    positive_text = "Positive with non-DKD (T2D)",
    negative_text = "Negative with non-DKD (T2D)",
    formula = "nonDKD_T2D"
  ),

  # SGLT2i comparisons (inferred)
  "sglt2in_hc" = list(
    fc_col = "logFC_sglt2i_t2dobSGLT2i_N",
    p_col = "p_sglt2i_t2dobSGLT2i_N",
    x_label = "logFC SGLT2i No vs. HC",
    positive_text = "Positive with T2D/OB\n(SGLT2i-)",
    negative_text = "Negative with T2D/OB\n(SGLT2i-)",
    formula = "SGLT2i_N"
  ),
  "sglt2iy_sglt2in" = list(
    fc_col = "logFC_sglt2i_t2dobSGLT2i_Y",
    p_col = "p_sglt2i_t2dobSGLT2i_Y",
    x_label = "logFC SGLT2i Yes vs. No",
    positive_text = "Positive with SGLT2i\n(in T2D/OB)",
    negative_text = "Negative with SGLT2i\n(in T2D/OB)",
    formula = "SGLT2i_Y"
  ),
  "nondkd100_sglt2in_hc" = list(
    fc_col = "logFC_dkd_group_100_hcnon_DKD",
    p_col = "p_dkd_group_100_hcnon_DKD",
    x_label = "logFC non-DKD (SGLT2i-) vs. HC\n(non-DKD: eGFR >= 90 and UACR < 100)",
    positive_text = "Positive with non-DKD\n(SGLT2i-)",
    negative_text = "Negative with non-DKD\n(SGLT2i-)",
    formula = "nonDKD_SGLT2i_N"
  ),
  "nondkd100_sglt2iy_hc" = list(
    fc_col = "logFC_dkd_group_100_hcnon_DKD",
    p_col = "p_dkd_group_100_hcnon_DKD",
    x_label = "logFC non-DKD (SGLT2i+) vs. HC\n(non-DKD: eGFR >= 90 and UACR < 100)",
    positive_text = "Positive with non-DKD\n(SGLT2i+)",
    negative_text = "Negative with non-DKD\n(SGLT2i+)",
    formula = "nonDKD_SGLT2i_Y"
  ),
  "nondkd30_sglt2in_hc" = list(
    fc_col = "logFC_dkd_group_30_sglt2i_hcnon_DKD",
    p_col = "p_dkd_group_30_sglt2i_hcnon_DKD",
    x_label = "logFC non-DKD (SGLT2i-) vs. HC\n(non-DKD: eGFR >= 90 and UACR < 30)",
    positive_text = "Positive with non-DKD\n(SGLT2i-)",
    negative_text = "Negative with non-DKD\n(SGLT2i-)",
    formula = "nonDKD_SGLT2i_N"
  ),
  "nondkd30_sglt2iy_hc" = list(
    fc_col = "logFC_dkd_group_30_hcnon_DKD",
    p_col = "p_dkd_group_30_hcnon_DKD",
    x_label = "logFC non-DKD (SGLT2i+) vs. HC\n(non-DKD: eGFR >= 90 and UACR < 30)",
    positive_text = "Positive with non-DKD\n(SGLT2i+)",
    negative_text = "Negative with non-DKD\n(SGLT2i+)",
    formula = "nonDKD_SGLT2i_Y"
  ),

  # Combo (GLP1-RA + SGLT2i) comparisons (inferred)
  "combon_hc" = list(
    fc_col = "logFC_combo_t2dobCombo_N",
    p_col = "p_combo_t2dobCombo_N",
    x_label = "logFC Combo No vs. HC\n(No GLP1-RA and SGLT2i)",
    positive_text = "Positive with No Combo",
    negative_text = "Negative with No Combo",
    formula = "Combo_N"
  ),
  "comboy_combon" = list(
    fc_col = "logFC_combo_t2dobCombo_Y",
    p_col = "p_combo_t2dobCombo_Y",
    x_label = "logFC Combo Yes vs. No\n(GLP1-RA and SGLT2i)",
    positive_text = "Positive with Combo",
    negative_text = "Negative with Combo",
    formula = "Combo_Y"
  ),
  "nondkd100_combon_hc" = list(
    fc_col = "logFC_dkd_group_100_combo_hcnon_DKD_Combo_N",
    p_col = "p_dkd_group_100_combo_hcnon_DKD_Combo_N",
    x_label = "logFC non-DKD (Combo-) vs. HC\n(non-DKD: eGFR >= 90 and UACR < 100)",
    positive_text = "Positive with non-DKD\n(Combo-)",
    negative_text = "Negative with non-DKD\n(Combo-)",
    formula = "nonDKD_Combo_N"
  ),
  "nondkd100_comboy_hc" = list(
    fc_col = "logFC_dkd_group_100_combo_hcnon_DKD_Combo_Y",
    p_col = "p_dkd_group_100_combo_hcnon_DKD_Combo_Y",
    x_label = "logFC non-DKD (Combo+) vs. HC\n(non-DKD: eGFR >= 90 and UACR < 100)",
    positive_text = "Positive with non-DKD\n(Combo+)",
    negative_text = "Negative with non-DKD\n(Combo+)",
    formula = "nonDKD_Combo_Y"
  ),
  "nondkd30_combon_hc" = list(
    fc_col = "logFC_dkd_group_30_combo_hcnon_DKD_Combo_N",
    p_col = "p_dkd_group_30_combo_hcnon_DKD_Combo_N",
    x_label = "logFC non-DKD (Combo-) vs. HC\n(non-DKD: eGFR >= 90 and UACR < 30)",
    positive_text = "Positive with non-DKD\n(Combo-)",
    negative_text = "Negative with non-DKD\n(Combo-)",
    formula = "nonDKD_Combo_N"
  ),
  "nondkd30_comboy_hc" = list(
    fc_col = "logFC_dkd_group_30_combo_hcnon_DKD_Combo_Y",
    p_col = "p_dkd_group_30_combo_hcnon_DKD_Combo_Y",
    x_label = "logFC non-DKD (Combo+) vs. HC\n(non-DKD: eGFR >= 90 and UACR < 30)",
    positive_text = "Positive with non-DKD\n(Combo+)",
    negative_text = "Negative with non-DKD\n(Combo+)",
    formula = "nonDKD_Combo_Y"
  ),    # --- DKD vs non-DKD (no SGLT2i) ---
  "dkd100_nosglt2i" = list(
    fc_col = "logFC_dkd_group_100DKD",
    p_col  = "p_dkd_group_100DKD",
    x_label = "logFC DKD vs non-DKD\n(DKD: eGFR < 90 or UACR ≥ 100)",
    positive_text = "Positive with DKD",
    negative_text = "Negative with DKD",
    formula = "DKD"
  ),
  "dkd30_nosglt2i" = list(
    fc_col = "logFC_dkd_group_30DKD",
    p_col  = "p_dkd_group_30DKD",
    x_label = "logFC DKD vs non-DKD\n(DKD: eGFR < 90 or UACR ≥ 30)",
    positive_text = "Positive with DKD",
    negative_text = "Negative with DKD",
    formula = "DKD"
  ),
  
  # --- GLP within DKD/non-DKD (ACR ≥ 100) ---
  "dkd_100_glpy_glpn" = list(
    fc_col = "logFC_glp_t2dobGLP_Y",
    p_col  = "p_glp_t2dobGLP_Y",
    x_label = "logFC GLP+ vs GLP– (DKD, ACR ≥ 100)",
    positive_text = "Positive with GLP+",
    negative_text = "Negative with GLP+",
    formula = "GLP_Y"
  ),
  "nondkd_100_glpy_glpn" = list(
    fc_col = "logFC_glp_t2dobGLP_Y",
    p_col  = "p_glp_t2dobGLP_Y",
    x_label = "logFC GLP+ vs GLP– (non-DKD, ACR ≥ 100)",
    positive_text = "Positive with GLP+",
    negative_text = "Negative with GLP+",
    formula = "GLP_Y"
  ),
  
  # --- GLP within DKD/non-DKD (ACR ≥ 30) ---
  "dkd_30_glpy_glpn" = list(
    fc_col = "logFC_glp_t2dobGLP_Y",
    p_col  = "p_glp_t2dobGLP_Y",
    x_label = "logFC GLP+ vs GLP– (DKD, ACR ≥ 30)",
    positive_text = "Positive with GLP+",
    negative_text = "Negative with GLP+",
    formula = "GLP_Y"
  ),
  "nondkd_30_glpy_glpn" = list(
    fc_col = "logFC_glp_t2dobGLP_Y",
    p_col  = "p_glp_t2dobGLP_Y",
    x_label = "logFC GLP+ vs GLP– (non-DKD, ACR ≥ 30)",
    positive_text = "Positive with GLP+",
    negative_text = "Negative with GLP+",
    formula = "GLP_Y"
  ),
  
  # --- SGLT2i within DKD (vs HC) ---
  "dkd30_sglt2in_hc" = list(
    fc_col = "logFC_sglt2i_t2dobSGLT2i_N",
    p_col  = "p_sglt2i_t2dobSGLT2i_N",
    x_label = "logFC SGLT2i– vs HC (DKD, ACR ≥ 30)",
    positive_text = "Positive with SGLT2i–",
    negative_text = "Negative with SGLT2i–",
    formula = "SGLT2i_N"
  ),
  "dkd30_sglt2iy_hc" = list(
    fc_col = "logFC_sglt2i_t2dobSGLT2i_Y",
    p_col  = "p_sglt2i_t2dobSGLT2i_Y",
    x_label = "logFC SGLT2i+ vs HC (DKD, ACR ≥ 30)",
    positive_text = "Positive with SGLT2i+",
    negative_text = "Negative with SGLT2i+",
    formula = "SGLT2i_Y"
  ),
  "dkd30_sglt2iy_sglt2in" = list(
    fc_col = "logFC_sglt2i_t2dobSGLT2i_Y",
    p_col  = "p_sglt2i_t2dobSGLT2i_Y",
    x_label = "logFC SGLT2i+ vs SGLT2i– (DKD, ACR ≥ 30)",
    positive_text = "Positive with SGLT2i+",
    negative_text = "Negative with SGLT2i+",
    formula = "SGLT2i_Y"
  ),
  "dkd100_sglt2iy_sglt2in" = list(
    fc_col = "logFC_sglt2i_t2dobSGLT2i_Y",
    p_col  = "p_sglt2i_t2dobSGLT2i_Y",
    x_label = "logFC SGLT2i+ vs SGLT2i– (DKD, ACR ≥ 100)",
    positive_text = "Positive with SGLT2i+",
    negative_text = "Negative with SGLT2i+",
    formula = "SGLT2i_Y"
  ),
  
  # --- SGLT2i within non-DKD (Y vs N) ---
  "nondkd30_sglt2iy_sglt2in" = list(
    fc_col = "logFC_sglt2i_t2dobSGLT2i_Y",
    p_col  = "p_sglt2i_t2dobSGLT2i_Y",
    x_label = "logFC SGLT2i+ vs SGLT2i– (non-DKD, ACR ≥ 30)",
    positive_text = "Positive with SGLT2i+",
    negative_text = "Negative with SGLT2i+",
    formula = "SGLT2i_Y"
  ),
  "nondkd100_sglt2iy_sglt2in" = list(
    fc_col = "logFC_sglt2i_t2dobSGLT2i_Y",
    p_col  = "p_sglt2i_t2dobSGLT2i_Y",
    x_label = "logFC SGLT2i+ vs SGLT2i– (non-DKD, ACR ≥ 100)",
    positive_text = "Positive with SGLT2i+",
    negative_text = "Negative with SGLT2i+",
    formula = "SGLT2i_Y"
  ),
  
  # --- SGLT2i within DKD (ACR ≥ 100) vs HC ---
  "dkd100_sglt2in_hc" = list(
    fc_col = "logFC_sglt2i_t2dobSGLT2i_N",
    p_col  = "p_sglt2i_t2dobSGLT2i_N",
    x_label = "logFC SGLT2i– vs HC (DKD, ACR ≥ 100)",
    positive_text = "Positive with SGLT2i–",
    negative_text = "Negative with SGLT2i–",
    formula = "SGLT2i_N"
  ),
  "dkd100_sglt2iy_hc" = list(
    fc_col = "logFC_sglt2i_t2dobSGLT2i_Y",
    p_col  = "p_sglt2i_t2dobSGLT2i_Y",
    x_label = "logFC SGLT2i+ vs HC (DKD, ACR ≥ 100)",
    positive_text = "Positive with SGLT2i+",
    negative_text = "Negative with SGLT2i+",
    formula = "SGLT2i_Y"
  ),
  
  # --- Combo global vs HC (Neither / SGLT2i_only / GLP_only / Both) ---
  "neither_hc" = list(
    fc_col = "logFC_combo_t2dobNeither",
    p_col  = "p_combo_t2dobNeither",
    x_label = "logFC Neither vs HC\n(No GLP1-RA or SGLT2i)",
    positive_text = "Positive with Neither",
    negative_text = "Negative with Neither",
    formula = "Neither"
  ),
  "sglt2ionly_hc" = list(
    fc_col = "logFC_combo_t2dobSGLT2i_only",
    p_col  = "p_combo_t2dobSGLT2i_only",
    x_label = "logFC SGLT2i only vs HC",
    positive_text = "Positive with SGLT2i only",
    negative_text = "Negative with SGLT2i only",
    formula = "SGLT2i_only"
  ),
  "glponly_hc" = list(
    fc_col = "logFC_combo_t2dobGLP_only",
    p_col  = "p_combo_t2dobGLP_only",
    x_label = "logFC GLP only vs HC",
    positive_text = "Positive with GLP only",
    negative_text = "Negative with GLP only",
    formula = "GLP_only"
  ),
  "both_hc" = list(
    fc_col = "logFC_combo_t2dobBoth",
    p_col  = "p_combo_t2dobBoth",
    x_label = "logFC Both vs HC\n(GLP1-RA + SGLT2i)",
    positive_text = "Positive with Both",
    negative_text = "Negative with Both",
    formula = "Both"
  ),
  
  # --- Combo vs Neither (no treatment) ---
  "sglt2ionly_neither" = list(
    fc_col = "logFC_combo_t2dobSGLT2i_only",
    p_col  = "p_combo_t2dobSGLT2i_only",
    x_label = "logFC SGLT2i only vs Neither",
    positive_text = "Positive with SGLT2i only",
    negative_text = "Negative with SGLT2i only",
    formula = "SGLT2i_only"
  ),
  "glponly_neither" = list(
    fc_col = "logFC_combo_t2dobGLP_only",
    p_col  = "p_combo_t2dobGLP_only",
    x_label = "logFC GLP only vs Neither",
    positive_text = "Positive with GLP only",
    negative_text = "Negative with GLP only",
    formula = "GLP_only"
  ),
  "both_neither" = list(
    fc_col = "logFC_combo_t2dobBoth",
    p_col  = "p_combo_t2dobBoth",
    x_label = "logFC Both vs Neither",
    positive_text = "Positive with Both",
    negative_text = "Negative with Both",
    formula = "Both"
  ),
  
  # --- Combo: Both vs single treatments ---
  "both_sglt2ionly" = list(
    fc_col = "logFC_combo_t2dobBoth",
    p_col  = "p_combo_t2dobBoth",
    x_label = "logFC Both vs SGLT2i only",
    positive_text = "Positive with Both",
    negative_text = "Negative with Both",
    formula = "Both"
  ),
  "both_glponly" = list(
    fc_col = "logFC_combo_t2dobBoth",
    p_col  = "p_combo_t2dobBoth",
    x_label = "logFC Both vs GLP only",
    positive_text = "Positive with Both",
    negative_text = "Negative with Both",
    formula = "Both"
  ),
  "sglt2ionly_glponly" = list(
    fc_col = "logFC_combo_t2dobSGLT2i_only",
    p_col  = "p_combo_t2dobSGLT2i_only",
    x_label = "logFC SGLT2i only vs GLP only",
    positive_text = "Positive with SGLT2i only",
    negative_text = "Negative with SGLT2i only",
    formula = "SGLT2i_only"
  ),
  
  # --- Combo within non-DKD (ACR ≥ 100) vs HC ---
  "nondkd100_neither_hc" = list(
    fc_col = "logFC_dkd_group_100_hcnon_DKD",
    p_col  = "p_dkd_group_100_hcnon_DKD",
    x_label = "logFC non-DKD vs HC\n(non-DKD: ACR < 100, subgroup: Neither)",
    positive_text = "Positive with non-DKD",
    negative_text = "Negative with non-DKD",
    formula = "nonDKD"
  ),
  "nondkd100_sglt2ionly_hc" = list(
    fc_col = "logFC_dkd_group_100_hcnon_DKD",
    p_col  = "p_dkd_group_100_hcnon_DKD",
    x_label = "logFC non-DKD vs HC\n(non-DKD: ACR < 100, subgroup: SGLT2i only)",
    positive_text = "Positive with non-DKD",
    negative_text = "Negative with non-DKD",
    formula = "nonDKD"
  ),
  "nondkd100_glponly_hc" = list(
    fc_col = "logFC_dkd_group_100_hcnon_DKD",
    p_col  = "p_dkd_group_100_hcnon_DKD",
    x_label = "logFC non-DKD vs HC\n(non-DKD: ACR < 100, subgroup: GLP only)",
    positive_text = "Positive with non-DKD",
    negative_text = "Negative with non-DKD",
    formula = "nonDKD"
  ),
  "nondkd100_both_hc" = list(
    fc_col = "logFC_dkd_group_100_hcnon_DKD",
    p_col  = "p_dkd_group_100_hcnon_DKD",
    x_label = "logFC non-DKD vs HC\n(non-DKD: ACR < 100, subgroup: Both)",
    positive_text = "Positive with non-DKD",
    negative_text = "Negative with non-DKD",
    formula = "nonDKD"
  ),
  
  # --- Combo within non-DKD (ACR ≥ 30) vs HC ---
  "nondkd30_neither_hc" = list(
    fc_col = "logFC_dkd_group_30_hcnon_DKD",
    p_col  = "p_dkd_group_30_hcnon_DKD",
    x_label = "logFC non-DKD vs HC\n(non-DKD: ACR < 30, subgroup: Neither)",
    positive_text = "Positive with non-DKD",
    negative_text = "Negative with non-DKD",
    formula = "nonDKD"
  ),
  "nondkd30_sglt2ionly_hc" = list(
    fc_col = "logFC_dkd_group_30_hcnon_DKD",
    p_col  = "p_dkd_group_30_hcnon_DKD",
    x_label = "logFC non-DKD vs HC\n(non-DKD: ACR < 30, subgroup: SGLT2i only)",
    positive_text = "Positive with non-DKD",
    negative_text = "Negative with non-DKD",
    formula = "nonDKD"
  ),
  "nondkd30_glponly_hc" = list(
    fc_col = "logFC_dkd_group_30_hcnon_DKD",
    p_col  = "p_dkd_group_30_hcnon_DKD",
    x_label = "logFC non-DKD vs HC\n(non-DKD: ACR < 30, subgroup: GLP only)",
    positive_text = "Positive with non-DKD",
    negative_text = "Negative with non-DKD",
    formula = "nonDKD"
  ),
  "nondkd30_both_hc" = list(
    fc_col = "logFC_dkd_group_30_hcnon_DKD",
    p_col  = "p_dkd_group_30_hcnon_DKD",
    x_label = "logFC non-DKD vs HC\n(non-DKD: ACR < 30, subgroup: Both)",
    positive_text = "Positive with non-DKD",
    negative_text = "Negative with non-DKD",
    formula = "nonDKD"
  ),
  "dkd30_glpn_hc" = list(
    fc_col = "logFC_dkd_group_30_hcDKD",
    p_col  = "p_dkd_group_30_hcDKD",
    x_label = "logFC DKD GLP1RA- vs HC\n(non-DKD: ACR < 30)",
    positive_text = "Positive with DKD, GLP-",
    negative_text = "Negative with DKD, GLP-",
    formula = "DKD, GLP-"
    ),
  "dkd30_glpy_hc" = list(
    fc_col = "logFC_dkd_group_30_hcDKD",
    p_col  = "p_dkd_group_30_hcDKD",
    x_label = "logFC DKD GLP1RA+ vs HC\n(non-DKD: ACR < 30)",
    positive_text = "Positive with DKD, GLP+",
    negative_text = "Negative with DKD, GLP+",
    formula = "DKD, GLP+"
  ),
  "dkd100_glpn_hc" = list(
    fc_col = "logFC_dkd_group_100_hcDKD",
    p_col  = "p_dkd_group_100_hcDKD",
    x_label = "logFC DKD GLP1RA- vs HC\n(non-DKD: ACR < 100)",
    positive_text = "Positive with DKD, GLP-",
    negative_text = "Negative with DKD, GLP-",
    formula = "DKD, GLP-"
  ),
  "dkd100_glpy_hc" = list(
    fc_col = "logFC_dkd_group_100_hcDKD",
    p_col  = "p_dkd_group_100_hcDKD",
    x_label = "logFC DKD GLP1RA+ vs HC\n(non-DKD: ACR < 100)",
    positive_text = "Positive with DKD, GLP+",
    negative_text = "Negative with DKD, GLP+",
    formula = "DKD, GLP+"
  )
)

# Loop through all combinations and create plots
for (comparison_name in names(folders)) {
  comparison_abbrev <- folders[comparison_name]
  params <- plot_params[[comparison_abbrev]]
  
  # Skip if parameters not defined
  if (is.null(params)) {
    message(paste0("Skipping ", comparison_name, ": parameters not defined"))
    next
  }
  
  for (cell_type in names(celltype_groups)) {
    # Get the data frame
    df_name <- paste0(tolower(gsub("/", "_", cell_type)), "_", comparison_abbrev)
    
    # Check if data frame exists
    if (!exists(df_name)) {
      message(paste0("Data frame ", df_name, " not found, skipping"))
      next
    }
    
    df <- get(df_name)
    
    # Plot with p-value
    tryCatch({
      plot_volcano(df,
                   fc = params$fc_col,
                   p_col = params$p_col,
                   title = NULL,
                   x_axis = params$x_label,
                   y_axis = "-log10(p.value)",
                   file_suffix = paste0("pval/", comparison_abbrev, "_", tolower(gsub("/", "_", cell_type))),
                   p_thresh = 0.05,
                   positive_text = params$positive_text,
                   negative_text = params$negative_text,
                   formula = params$formula,
                   cell_type = cell_type)
      
      message(paste0("Created p-value volcano plot for ", cell_type, " - ", comparison_name))
    }, error = function(e) {
      message(paste0("Error creating p-value plot for ", cell_type, " - ", comparison_name, ": ", e$message))
    })
    
    # Plot with FDR if FDR column exists
    if ("fdr" %in% colnames(df)) {
      tryCatch({
        plot_volcano(df,
                     fc = params$fc_col,
                     p_col = "fdr",
                     title = NULL,
                     x_axis = params$x_label,
                     y_axis = "-log10(FDR)",
                     file_suffix = paste0("fdr/", comparison_abbrev, "_", tolower(gsub("/", "_", cell_type)), "_fdr"),
                     p_thresh = 0.05,
                     positive_text = params$positive_text,
                     negative_text = params$negative_text,
                     formula = params$formula,
                     cell_type = cell_type)
        
        message(paste0("Created FDR volcano plot for ", cell_type, " - ", comparison_name))
      }, error = function(e) {
        message(paste0("Error creating FDR plot for ", cell_type, " - ", comparison_name, ": ", e$message))
      })
    }
  }
}

# Pathways using Hallmark

for (folder in folders) {
  for (cell in names(celltype_groups)) {
    # Create the variable name
    var_name <- paste0(tolower(gsub("/", "_", cell)), "_", folder)
    
    # Try to process with error handling
    tryCatch({
      # Check if the variable exists
      if (!exists(var_name)) {
        message(paste0("WARNING: Variable '", var_name, "' does not exist. Skipping..."))
        next
      }
      
      # Get the dataframe
      df <- get(var_name)
      
      # Run GSEA analysis
      gsea_list <- run_fgsea_analysis(results_annotated = df,
                                      stat_col = 2,
                                      gene_col = "Gene",
                                      nPermSimple = 100000,
                                      references = c("hallmark"))
      
      # Save the results
      output_path <- file.path(root_path, paste0("Renal HERITAGE/Results/GSEA/hallmark_", folder, "_", tolower(gsub("/", "_", cell)), "_gsea.RDS"))
      saveRDS(gsea_list, output_path)
      
      message(paste0("Successfully processed: ", folder, " - ", cell))
      
    }, error = function(e) {
      message(paste0("ERROR: Failed to process folder '", folder, "' and cell type '", cell, "'"))
      message(paste0("  Variable name: ", var_name))
      message(paste0("  Error details: ", e$message))
      message("  Skipping to next combination...\n")
    })
  }
}

# Pathway visualization
# Initialize a list to store all results
gsea_results <- list()

for (folder in folders) {
    for (cell in names(celltype_groups)) {
    file_path <- file.path(root_path, paste0("Renal HERITAGE/Results/GSEA/hallmark_", folder, "_", tolower(gsub("/", "_", cell)), "_gsea.RDS"))
    
    # Check if file exists and read it
    if (file.exists(file_path)) {
      # Create a unique name for this combination
      result_name <- paste0(folder, "_", tolower(gsub("/", "_", cell)))
      
      # Read the RDS file
      gsea_results[[result_name]] <- readRDS(file_path)
      message(paste0("Successfully loaded: ", result_name))
    } else {
      message(paste0("File not found: ", file_path))
    }
  }
}

for (folder in folders) {
  for (cell in names(celltype_groups)) {
    message(paste0(cell, ": ", folder))
    gsea_res_sub <- gsea_results[[paste0(folder, "_", tolower(gsub("/", "_", cell)))]]
    plot_gsea_results(gsea_res_sub, 
                      cell_name = cell)
    ggsave(file.path(root_path, paste0("Renal HERITAGE/Results/Figures/Pathways/hallmark_", folder, "_", tolower(gsub("/", "_", cell)), "_gsea.png")))
  }
}



########################################################################################################################################################################################################################

# look for overlaps and directions of overlaps
ec_nondkd30_glpn_hc$direction_nondkd_30_glpn_hc <- ifelse(ec_nondkd30_glpn_hc$logFC_dkd_group_30_hcnon_DKD > 0 & ec_nondkd30_glpn_hc$p_dkd_group_30_hcnon_DKD < 0.05, "+",
                                                          ifelse(ec_nondkd30_glpn_hc$logFC_dkd_group_30_hcnon_DKD < 0 & ec_nondkd30_glpn_hc$p_dkd_group_30_hcnon_DKD < 0.05, "-",
                                                                 "NS"))

ec_nondkd30_glpy_hc$direction_nondkd_30_glpy_hc <- ifelse(ec_nondkd30_glpy_hc$logFC_dkd_group_30_hcnon_DKD > 0 & ec_nondkd30_glpy_hc$p_dkd_group_30_hcnon_DKD < 0.05, "+",
                                                          ifelse(ec_nondkd30_glpy_hc$logFC_dkd_group_30_hcnon_DKD < 0 & ec_nondkd30_glpy_hc$p_dkd_group_30_hcnon_DKD < 0.05, "-",
                                                                 "NS"))

ec_nondkd30_glpn_hc_subset <- ec_nondkd30_glpn_hc %>%
  select(Gene, direction_nondkd_30_glpn_hc)

ec_nondkd30_glpy_hc_subset <- ec_nondkd30_glpy_hc %>%
  select(Gene, direction_nondkd_30_glpy_hc)

ec_genes <- unique(c(ec_nondkd30_glpn_hc_subset$Gene, ec_nondkd30_glpy_hc_subset$Gene))

ec_overlap <- data.frame(Gene = ec_genes) %>% 
  left_join(ec_nondkd30_glpn_hc_subset) %>%
  left_join(ec_nondkd30_glpy_hc_subset) %>%
  mutate(category = case_when(direction_nondkd_30_glpn_hc == "+" & direction_nondkd_30_glpy_hc == "-" ~ "Reversed",
                              direction_nondkd_30_glpn_hc == "-" & direction_nondkd_30_glpy_hc == "+" ~ "Reversed",
                              direction_nondkd_30_glpn_hc != "NS" & direction_nondkd_30_glpy_hc == "NS" ~ "Normalized",
                              direction_nondkd_30_glpn_hc == "NS" & direction_nondkd_30_glpy_hc != "NS" ~ "Treatment specific",
                              direction_nondkd_30_glpn_hc == "+" & direction_nondkd_30_glpy_hc == "+" ~ "Persistent",
                              direction_nondkd_30_glpn_hc == "-" & direction_nondkd_30_glpy_hc == "-" ~ "Persistent")) %>%
  filter(!is.na(category))

ec_overlap %>%
  dplyr::count(category) %>%
  ggplot(aes(x = factor(category, levels = c("Treatment specific", "Normalized", "Persistent", "Reversed")), y = n, fill = category)) +
  geom_col() + 
  geom_text(aes(label = n), vjust = 0) + 
  labs(x = NULL, y = "Count") +
  theme(panel.background = element_blank(),
        text = element_text(size = 15),
        axis.text.x = element_text(angle = 60, hjust = 1),
        legend.position = "none") +
  scale_fill_manual(values = c("#52796f", "#89c2d9",
                               "#e26d5c", "#3d5a80"))
# ggsave()

ec_combined_df <- data.frame(Gene = ec_genes) %>% 
  left_join(ec_nondkd30_glpn_hc, by = join_by(Gene)) %>%
  left_join(ec_nondkd30_glpy_hc, by = join_by(Gene)) %>%
  left_join(ec_overlap) %>%
  filter(!is.na(category))

ec_combined_df %>%
  filter(category %in% c("Reversed", "Persistent")) %>%
  ggplot(aes(x = logFC_dkd_group_30_hcnon_DKD.x, y = logFC_dkd_group_30_hcnon_DKD.y, color = category)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  theme(panel.background = element_blank()) +
  labs(x = "logFC (non-DKD/GLP-)",
       y = "logFC (non-DKD/GLP+)",
       color = NULL) +
  scale_color_manual(values = c("#89c2d9",
                                "#e26d5c"))

ec_combined_df %>%
  filter(category %nin% c("Reversed", "Persistent")) %>%
  ggplot(aes(x = logFC_dkd_group_30_hcnon_DKD.x, y = logFC_dkd_group_30_hcnon_DKD.y, color = category)) +
  geom_point(alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  theme(panel.background = element_blank()) +
  labs(x = "logFC (non-DKD/GLP-)",
       y = "logFC (non-DKD/GLP+)",
       color = NULL) +
  scale_color_manual(values = c("#52796f",
                                "#3d5a80"))

## volcano plots with these highlights
# - **VEGFA** (master angiogenic regulator): Downregulated in GLP+ (log2FC=-0.61, p=0.0028) - Treatment-specific effect suggesting reduced angiogenic signaling
# - **ICAM1** (adhesion molecule): Downregulated in GLP+ (log2FC=-0.38) - Suggests reduced inflammatory adhesion
# - **COL4A1** (basement membrane collagen): Upregulated in GLP- (log2FC=0.32) but **normalized in GLP+** - Suggests prevention of early fibrotic changes
# - **SLC2A4/GLUT4** (insulin-responsive glucose transporter): Downregulated in GLP+ (log2FC=-1.23) - Unexpected finding, may reflect reduced metabolic stress or altered insulin signaling
# - **NOX4** (NADPH oxidase): Upregulated in GLP+ (log2FC=0.96) - **CONCERNING** - increased oxidative stress marker
# - **CAT** (Catalase): Upregulated in GLP- but normalized in GLP+ - Suggests reduced oxidative stress response need

# c("VEGFA", "ICAM1", "COL4A1", "SLC2A4", "GLUT4", "NOX4", "CAT")
