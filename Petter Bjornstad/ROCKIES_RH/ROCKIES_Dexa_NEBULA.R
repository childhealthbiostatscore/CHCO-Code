#DEXA Analysis
library(scran)
library(future)
library(future.apply)
library(tidyverse)
library(colorspace)
library(patchwork)
library(ggdendro)
library(cowplot)
library(ggpubr)
library(rstatix)
library(arsenal)
library(Biobase)
library(msigdbr)
library(kableExtra)
library(knitr)
library(REDCapR)
library(data.table)
library(emmeans)
library(NMF)
library(pheatmap)
library(UpSetR)
library(enrichR)
library(WriteXLS)
library(SAVER)
library(readxl)
library(limma)
library(edgeR)
library(BiocGenerics)
library(GSEABase)
library(slingshot)
library(SingleCellExperiment)
library(MAST)
library(muscat)
library(scater)
library(Seurat)
library(jsonlite)
library(dplyr)
library(glmmTMB)
library(reshape2)
library(broom.mixed)
library(nebula)
library(doParallel)




harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')


dat <- harmonized_data %>%
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, first(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, first(na.omit(.x)))),
                   .by = c(record_id, visit))


dat2 <- dat %>% dplyr::select(record_id, mrn, visit, epic_sglti2_1, starts_with('dexa_')) %>% 
  filter(!is.na(epic_sglti2_1))


dat2 <- dat2 %>% filter(!is.na(dexa_body_fat))






dex_var <- c("dexa_ag_ratio",  "dexa_body_fat",            "dexa_bone_mineral_density", "dexa_est_vat",             
             "dexa_fat_kg",               "dexa_lean_kg",              "dexa_lean_mass",            "dexa_trunk_kg",            
            "dexa_trunk_mass"  )




gene_list <- c(
      "LEPR", "ADIPOR1", "ADIPOR2",
      "PPARG", "PPARA", "PPARD",
      "SREBF1", "SREBF2",
      "STAT3", "JAK2", "SOCS3",
      "TNF", "IL6", "IL1B", "CCL2",
      "NFKB1", "RELA", "IKBKB",
      "TLR2", "TLR4", "MYD88",
      "NLRP3", "CASP1", "IL18",
      "FASN", "ACACA", "SCD",
      "CPT1A", "CPT2", "ACOX1",
      "LDLR", "SCARB1", "ABCA1",
      "APOE", "APOB", "LPL",
      "DGAT1", "DGAT2", "PLIN1",
      "TGFB1", "TGFB2", "SMAD2", "SMAD3",
      "COL1A1", "COL3A1", "COL4A1",
      "FN1", "CTGF", "ACTA2",
      "MMP2", "MMP9", "TIMP1",
      "NOX1", "NOX2", "NOX4",
      "SOD1", "SOD2", "CAT", "GPX1",
      "NRF2", "KEAP1", "NQO1",
      "HMOX1", "GCLC", "GCLM",
      "HMOX1", "GCLC", "GCLM",
      "IRS1", "IRS2", "INSR",
      "GLUT4", "GLUT2",
      "SOCS1", "SOCS3", "PTPN1",
      "PRKAA1", "PRKAA2",
      "HSPA5", "EIF2AK3", "ERN1", "ATF6", # UPR sensors
      "XBP1", "ATF4", "DDIT3",
      "EDEM1", "DNAJB9", "HYOU1"
    )



#cell-type specific analyses 


 
"PT" = c(
                        # Fatty acid uptake and metabolism
                        "CD36", "FABP1", "FABP3", "SLC27A2", "SLC27A4",
                        # Glucose/lipid interaction
                        "PPARGC1A", "FOXO1", "SIRT1", "AMPK",
                        # PT-specific transporters affected by lipids
                        "SLC5A2", "SLC9A3", "SLC34A1"
                      )
              "POD" = c(
                        # Podocyte lipotoxicity
                        "NPHS1", "NPHS2", "PODXL", "SYNPO", "CD2AP",
                        # Cholesterol efflux
                        "ABCA1", "ABCG1", "APOE",
                      "CHOP", "BIP", "XBP1s"

)
"Endothelial" = c(
  # Endothelial dysfunction
  "VCAM1", "ICAM1", "SELE", "SELP",
  # NO synthesis
  "NOS3", "GUCY1A3", "GUCY1B3",
  # Angiogenesis
  "VEGFA", "VEGFR2", "ANGPT2"
)
"Immune" = c(
  # Macrophage polarization
  "CD68", "CD86", "CD163", "CD206",
  # Cytokines
  "IL10", "IL12", "TGFB1",
  # Chemokines
  "CCL2", "CCL3", "CCL5", "CXCL10"
)
"TAL" = c(
  # Metabolic adaptation
  "PPARGC1A", "NRF1", "TFAM",
  # Transport regulation
  "SLC12A1", "KCNJ1", "ATP1A1")

"DCT" = c(
  # Calcium handling
  "SLC12A3", "TRPM6", "CALB1",
  # Metabolic genes
  "HIF1A", "EPAS1"
)







#LC vs. T2D (no SLGT2)



load('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/No_Med_line700.Rdata')

remove(so_kpmp_sc)

#dat_groups <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/ROCKIES_GroupAssignments.txt')
#dat_groups <- dat_groups %>% filter(group2 %in% c('Lean Control', 'T2D-No SGLTi2'))

#so_subset <- subset(so_subset, record_id == dat_groups$record_id)
test <- so_subset@meta.data

test <- test %>% dplyr::select(-starts_with('dexa_'))


test <- test %>% left_join(dat2, by='record_id')


so_subset@meta.data$dexa_ag_ratio <- test$dexa_ag_ratio
so_subset@meta.data$dexa_body_fat <- test$dexa_body_fat
so_subset@meta.data$dexa_bone_mineral_density <- test$dexa_bone_mineral_density
so_subset@meta.data$dexa_est_vat <- test$dexa_est_vat
so_subset@meta.data$dexa_fat_kg <- test$dexa_fat_kg
so_subset@meta.data$dexa_lean_kg <- test$dexa_lean_kg
so_subset@meta.data$dexa_lean_mass <- test$dexa_lean_mass
so_subset@meta.data$dexa_trunk_kg <- test$dexa_trunk_kg
so_subset@meta.data$dexa_trunk_mass <- test$dexa_trunk_mass


#Make sure exposure/independent/x variable or group variable is a factor variable
so_subset$group <- factor(so_subset$group)
#Make sure to set reference level
so_subset$group  <- relevel(so_subset$group ,ref="Lean_Control")

counts_path <- round(GetAssayData(so_subset, layer = "counts")) # load counts and round


T2D_LC_Dexa_Analysis <- function(data, dir.results, celltype, variable){
  if(celltype == 'All'){
    so_celltype <- so_subset 
  }else if(celltype %in% c('TAL', 'EC', 'POD', 'PT')){
    so_celltype <- subset(so_subset,celltype2==celltype)
    DefaultAssay(so_celltype) <- "RNA" 
  }else if(celltype != 'DCTall'){
    so_celltype <- subset(so_subset,KPMP_celltype==celltype)
    DefaultAssay(so_celltype) <- "RNA" 
  }else if(celltype == 'DCTall'){
    so_celltype <- subset(so_subset, DCT_celltype=='DCT')
  }
  
  
  
  
  
  
  nrow(so_celltype) #34 genes
  ncol(so_celltype) #13534 PT cells
  
  celltype2 <- str_replace_all(celltype,"/","_")
  celltype2 <- str_replace_all(celltype2,"-","_")
  
  #Make sure exposure/independent/x variable or group variable is a factor variable
  so_celltype$group <- factor(so_celltype$group)
  #Make sure to set reference level
  so_celltype$group  <- relevel(so_celltype$group ,ref="Lean_Control")
  
  
  counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round
  
  # With parallelization
  #TCA Cycle
  # List of genes
  genes_list <- tca_genes
  
  cl <- makeCluster(1)
  registerDoParallel(cl)
  test2 <- test %>% filter(record_id %in% unique(so_celltype@meta.data$record_id))
  t2d_count <- test2 %>% filter(group == 'Type_2_Diabetes') %>% nrow()
  lc_count <- test2 %>% filter(group == 'Lean_Control') %>% nrow()
  
  
  start_time <- Sys.time()
  
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- subset(so_celltype,features=g)@meta.data
      pred_gene <- model.matrix(~group, data = meta_gene)
      # library <- meta_gene$library_size
      library <- meta_gene$pooled_offset
      data_g_gene <- group_cell(count = count_gene, id = meta_gene$kit_id, pred = pred_gene,offset=library)
      
      if (is.null(data_g_gene)) {
        data_g_gene <- list(count = count_gene, id = meta_gene$kit_id, pred = pred_gene, offset = library)
      }
      
      #With offset
      result <- nebula(count = data_g_gene$count, id = data_g_gene$id, pred = data_g_gene$pred, ncore = 1, reml=T,model="NBLMM",output_re = T,covariance=T,offset=data_g_gene$library)
      
      list(gene = g, result = result)  # return both gene name and result
      
    }, error = function(e) {
      NULL
    })
  }
  
  stopCluster(cl)
  end_time <- Sys.time()
  print(end_time - start_time)
  
  # set the names of results based on gene names
  nebula_results_list <- Filter(Negate(is.null), nebula_results_list)  # remove NULLs first
  names(nebula_results_list) <- sapply(nebula_results_list, function(x) x$gene)  # set names
  nebula_results_list <- lapply(nebula_results_list, function(x) x$result)  # clean list back to just results
  
  PT_nebula_converged <- map_dfr(
    names(nebula_results_list),
    function(gene_name) {
      converged <- nebula_results_list[[gene_name]]$convergence
      df <- data.frame(Gene = gene_name,
                       Convergence_Code = converged)
      return(df)
    }
  )
  
  nebula_summaries <- map_dfr(
    names(nebula_results_list),
    function(gene_name) {
      df <- nebula_results_list[[gene_name]]$summary
      df <- df %>% mutate(Gene = gene_name)
      return(df)
    }
  )
  nonconverge_genes <- unique(PT_nebula_converged$Gene[which(PT_nebula_converged$Convergence_Code==-40)]) 
  
  #Make dataframe of final results
  full_results <- as.data.frame(nebula_summaries)
  #Calculate number of genes filtered out for low expression 
  low_exp <- length(tca_genes)-length(full_results$gene)
  #Filter out non-converging genes
  full_results <- full_results %>% 
    filter(!gene %in%  nonconverge_genes)
  #Calculate nonconvergence rate
  nebula_nonconverged_percent <- paste0(round((1-(length(tca_genes)-length(nonconverge_genes))/length(tca_genes))*100,3),"%")
  # nebula_nonconverged_percent <- (length(rownames(counts_path))-length(unique(full_results$gene)))/length(rownames(counts_path))
  # print(paste0(nebula_nonconverged_percent*100, "% failed to converge"))
  full_results <- full_results %>%
    mutate(fdr=p.adjust(`p_groupType_2_Diabetes`,method="fdr"))  
  # mutate(fdr3=p.adjust(PValue3,method="fdr"))
  full_results$PValue10 <- -log10(pmax(full_results$`p_groupType_2_Diabetes`, 1e-10))  # Avoid log(0)
  
  write.csv(full_results,fs::path(dir.results,paste0("NEBULA_TCA_cycle_",celltype2,"_cells_LC_T2D_NoMed_unadjusted_pooled_offset.csv")))
  
  
  
  
  
  
  
  
  
  
}
















#T2D SGLT2 Analysis 










