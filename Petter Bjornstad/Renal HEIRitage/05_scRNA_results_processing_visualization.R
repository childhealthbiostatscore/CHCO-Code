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
folders <- c("DKD_vs_nonDKD_100" = "dkd100",
             "DKD_100_vs_hc" = "dkd100_hc",
             "nonDKD_100_vs_hc" = "nondkd100_hc",
             "DKD_vs_nonDKD_30" = "dkd30",
             "DKD_30_vs_hc" = "dkd30_hc",
             "nonDKD_30_vs_hc" = "nondkd30_hc",
             "GLP_N_vs_HC" = "glpn_hc",
             "GLP_Y_vs_GLP_N" = "glpy_glpn",
             "nonDKD_100_GLP_N_vs_HC" = "nondkd100_glpn_hc",
             "nonDKD_100_GLP_Y_vs_HC" = "nondkd100_glpy_hc",
             "nonDKD_30_GLP_N_vs_HC" = "nondkd30_glpn_hc",
             "nonDKD_30_GLP_Y_vs_HC" = "nondkd30_glpy_hc"
)

# Define common parameters
celltype_groups <- list(
  PT = c("PT-S1/S2", "PT-S3", "aPT"),
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

# Load all data first
for (folder in names(folders)) {
  for (cell in names(celltype_groups)){
    processed_df <- s3readRDS(object = paste0("Projects/CKD/RH_RH2/Results/nebula/", folder, "/", cell, "/", cell, "_rh_rh2_imp_nebula_kpmp_", folders[folder], "_processed.rds"),
                              bucket = "scrna", region = "")
    var_name <- paste0(tolower(cell), "_", folders[folder])
    assign(var_name, processed_df, envir = .GlobalEnv)
    
    write.csv(processed_df, file.path(root_path, paste0("Renal HERITAGE/Results/nebula/csv/full/", var_name, "_kpmp.csv")), row.names = F, fileEncoding = "UTF-8")
    
    write.csv(subset(processed_df, processed_df[[6]] < 0.05), file.path(root_path, paste0("Renal HERITAGE/Results/nebula/csv/pval/", var_name, "_kpmp_pval.csv")), row.names = F, fileEncoding = "UTF-8")
    
    write.csv(subset(processed_df, fdr < 0.05), file.path(root_path, paste0("Renal HERITAGE/Results/nebula/csv/fdr/", var_name, "_kpmp_fdr.csv")), row.names = F, fileEncoding = "UTF-8")
  }
}

# Define plot parameters for each comparison
plot_params <- list(
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
    p_col = "p_dkd_group_100_hcnonDKD",
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
    df_name <- paste0(tolower(cell_type), "_", comparison_abbrev)
    
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
                   file_suffix = paste0("pval/", comparison_abbrev, "_", tolower(cell_type)),
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
                     file_suffix = paste0("fdr/", comparison_abbrev, "_", tolower(cell_type), "_fdr"),
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


# Overlaps in models per cell type
pt_mod1 <- pt_nondkd100_hc %>%
  dplyr::select(Gene, 
                logFC = logFC_dkd_group_100_hcnon_DKD, 
                p.value = p_dkd_group_100_hcnon_DKD, 
                fdr) %>%
  mutate(model = "A_DKDn_HC") %>%
  distinct(Gene, .keep_all = T)

pt_mod2 <- pt_dkd100_hc %>%
  dplyr::select(Gene, 
                logFC = logFC_dkd_group_100_hcDKD, 
                p.value = p_dkd_group_100_hcDKD, 
                fdr) %>%
  mutate(model = "B_DKD_HC")%>%
  distinct(Gene, .keep_all = T)

pt_mod3 <- pt_dkd100 %>%
  dplyr::select(Gene,
                logFC = logFC_dkd_group_100DKD,
                p.value = p_dkd_group_100DKD, 
                fdr) %>%
  mutate(model = "C_DKDy_DKDn")%>%
  distinct(Gene, .keep_all = T)
  
pt_mod4 <- pt_glpn_hc %>%
  dplyr::select(Gene, 
                logFC = logFC_glp_t2dobGLP_N,
                p.value = p_glp_t2dobGLP_N,
                fdr) %>%
  mutate(model = "D_GLPn_HC")%>%
  distinct(Gene, .keep_all = T)

pt_mod5 <- pt_glpy_glpn %>%
  dplyr::select(Gene,
                logFC = logFC_glp_t2dobGLP_Y,
                p.value = p_glp_t2dobGLP_Y,
                fdr) %>%
  mutate(model = "E_GLPy_GLPn")%>%
  distinct(Gene, .keep_all = T)

pt_100_pvalue <- rbind(pt_mod1, pt_mod2, pt_mod3, pt_mod4, pt_mod5) %>%
  filter(p.value < 0.05) %>%
  pivot_wider(
    id_cols = Gene,
    names_from = model,
    values_from = c(logFC, p.value, fdr),
    names_glue = "{model}_{.value}"
  ) %>%
  mutate(overlaps = paste0(
    case_when(
      is.na(A_DKDn_HC_logFC) ~ "",
      A_DKDn_HC_logFC > 0 ~ "A+",
      A_DKDn_HC_logFC < 0 ~ "A-",
      TRUE ~ ""
    ),
    case_when(
      is.na(B_DKD_HC_logFC) ~ "",
      B_DKD_HC_logFC > 0 ~ "B+",
      B_DKD_HC_logFC < 0 ~ "B-",
      TRUE ~ ""
    ),
    case_when(
      is.na(C_DKDy_DKDn_logFC) ~ "",
      C_DKDy_DKDn_logFC > 0 ~ "C+",
      C_DKDy_DKDn_logFC < 0 ~ "C-",
      TRUE ~ ""
    ),
    case_when(
      is.na(D_GLPn_HC_logFC) ~ "",
      D_GLPn_HC_logFC > 0 ~ "D+",
      D_GLPn_HC_logFC < 0 ~ "D-",
      TRUE ~ ""
    ),
    case_when(
      is.na(E_GLPy_GLPn_logFC) ~ "",
      E_GLPy_GLPn_logFC > 0 ~ "E+",
      E_GLPy_GLPn_logFC < 0 ~ "E-",
      TRUE ~ ""
    )
  ),
  n_overlap = rowSums(
    !is.na(dplyr::select(cur_data(),
                         A_DKDn_HC_logFC, B_DKD_HC_logFC,
                         C_DKDy_DKDn_logFC,
                         D_GLPn_HC_logFC, E_GLPy_GLPn_logFC
    ))
  )
  )

pt_100_fdr <- rbind(pt_mod1, pt_mod2, pt_mod3, pt_mod4, pt_mod5) %>%
  filter(fdr < 0.05) %>%
  pivot_wider(
    id_cols = Gene,
    names_from = model,
    values_from = c(logFC, p.value, fdr),
    names_glue = "{model}_{.value}"
  ) %>%
  mutate(overlaps = paste0(
    case_when(
      is.na(A_DKDn_HC_logFC) ~ "",
      A_DKDn_HC_logFC > 0 ~ "A+",
      A_DKDn_HC_logFC < 0 ~ "A-",
      TRUE ~ ""
    ),
    case_when(
      is.na(B_DKD_HC_logFC) ~ "",
      B_DKD_HC_logFC > 0 ~ "B+",
      B_DKD_HC_logFC < 0 ~ "B-",
      TRUE ~ ""
    ),
    # case_when(
    #   is.na(C_DKDy_DKDn_logFC) ~ "",
    #   C_DKDy_DKDn_logFC > 0 ~ "C+",
    #   C_DKDy_DKDn_logFC < 0 ~ "C-",
    #   TRUE ~ ""
    # ),

    case_when(
      is.na(D_GLPn_HC_logFC) ~ "",
      D_GLPn_HC_logFC > 0 ~ "D+",
      D_GLPn_HC_logFC < 0 ~ "D-",
      TRUE ~ ""
    ),
    case_when(
      is.na(E_GLPy_GLPn_logFC) ~ "",
      E_GLPy_GLPn_logFC > 0 ~ "E+",
      E_GLPy_GLPn_logFC < 0 ~ "E-",
      TRUE ~ ""
    )
  ),
  n_overlap = rowSums(
    !is.na(dplyr::select(cur_data(),
                  A_DKDn_HC_logFC, B_DKD_HC_logFC,
                  D_GLPn_HC_logFC, E_GLPy_GLPn_logFC
    ))
  )
  )


pt_100_fdr <- make_overlap_table(
  modA = pt_nondkd100_hc %>% dplyr::select(Gene,
                                    logFC   = logFC_dkd_group_100_hcnon_DKD,
                                    p.value = p_dkd_group_100_hcnon_DKD,
                                    fdr),
  modB = pt_dkd100_hc %>% dplyr::select(Gene,
                                 logFC   = logFC_dkd_group_100_hcDKD,
                                 p.value = p_dkd_group_100_hcDKD,
                                 fdr),
  modC = pt_dkd100 %>% dplyr::select(Gene,
                              logFC   = logFC_dkd_group_100DKD,
                              p.value = p_dkd_group_100DKD,
                              fdr),
  modD = pt_glpn_hc %>% dplyr::select(Gene,
                               logFC   = logFC_glp_t2dobGLP_N,
                               p.value = p_glp_t2dobGLP_N,
                               fdr),
  modE = pt_glpy_glpn %>% dplyr::select(Gene,
                                 logFC   = logFC_glp_t2dobGLP_Y,
                                 p.value = p_glp_t2dobGLP_Y,
                                 fdr),
  sig = "fdr", alpha = 0.05)
