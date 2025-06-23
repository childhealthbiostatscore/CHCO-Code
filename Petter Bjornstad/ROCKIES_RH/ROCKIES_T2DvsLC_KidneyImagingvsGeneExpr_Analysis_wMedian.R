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



#Get files 


so_kpmp_sc <- readRDS("C:/Users/netio/Downloads/PB_90_RPCAFix_Metadata_LR_RM_kpmpV1labelled.rds")

load("C:/Users/netio/Downloads/TCA_genes.rda")
load("C:/Users/netio/Downloads/OxPhos_genes.rda")


so_kpmp_sc$kit_id[which(so_kpmp_sc$kit_id=="KI-0014643")] <- "KL-0014643"
so_kpmp_sc$kit_id[which(so_kpmp_sc$kit_id=="kl-0023998")] <- "KL-0023998"


#Get harmonized data
harmonized_data <- read.csv("C:/Users/netio/Documents/Harmonized_data/harmonized_dataset.csv", na = '')


dat <- harmonized_data %>% dplyr::select(-dob) %>% 
  arrange(screen_date) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, first(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, first(na.omit(.x)))),
                   .by = c(record_id, visit))


length(unique(dat$mrn))#467
length(unique(dat$record_id))#644

#Filter to renal heir/renal heritage/crocodile 
dat <- dat %>% 
  filter(grepl("RH",record_id) |grepl("RH2",record_id)|grepl("CRC",record_id))  #grepl("IT",record_id)
length(unique(dat$mrn))#173
length(unique(dat$record_id))#206

#Filter to baseline visits only (no post surgery)
dat <- dat %>% 
  filter(visit=="baseline") 
length(unique(dat$mrn))#173
length(unique(dat$record_id))#206

#Filter to T2D or LC only
dat <- dat %>% 
  filter(group=="Type 2 Diabetes" | group=="Lean Control")
length(unique(dat$mrn))#113
length(unique(dat$record_id))#130


#Filter to those NOT on SGLT2s
#dat <- dat %>% 
#  filter(epic_sglti2_1=="No")
length(unique(dat$mrn))#50
length(unique(dat$record_id))#52

#Filter to only those with a kit id and biospy in pb 90
dat <- dat %>% 
  filter(kit_id %in% so_kpmp_sc$kit_id)
length(unique(dat$mrn))#31
length(unique(dat$record_id))#32

# #Find coenrolled individuals
dat$record_id[which(duplicated(dat$mrn))] #"RH-23-T/"RH2-14-T"
dat$mrn[which(duplicated(dat$mrn))] #1664581
dat$kit_id[which(dat$record_id=="RH-23-T")] #KL-0014632
dat$kit_id[which(dat$record_id=="RH2-14-T")] #KL-0029535, Exclude their second biopsy. Maintian baseline no repeated measures only

dat <- dat %>%
  filter(record_id!="RH2-14-T")
length(unique(dat$mrn))#31
length(unique(dat$record_id))#31

#Final dataset = 31
table(dat$study,dat$group)
#               Lean Control Type 2 Diabetes
# CROCODILE                 13               0
# RENAL-HEIR                 0              11
# RENAL-HEIRitage            0               7

table(dat$group,dat$epic_sglti2_1)



medications <- readxl::read_xlsx("C:/Users/netio/Documents/Harmonized_data/Biopsies_w_mrn_Oct3.xlsx")
medications <- medications %>% dplyr::select(mrn, ends_with('_1'), -starts_with('ever_'))
names(medications) <- str_replace(names(medications), pattern = '_1', replacement = '')

RH2_55_meds <- medications %>% filter(mrn == 2204106)
RH2_55_meds

med_desc <- readxl::read_xlsx('C:/Users/netio/Documents/Harmonized_data/Biopsies_w_mrn_Oct3.xlsx', sheet = 2, col_names=F)
med_desc <- med_desc[seq(from = 1, to = 57, by =4),]
names(med_desc) <- c('Medication', 'Description')
med_desc$Medication <- str_replace(med_desc$Medication, pattern='_1', replacement = '')





length(unique(dat$mrn)) #91 LC and T2D at baseline
length(unique(dat$record_id)) #25 no med
length(unique(which(dat$group=="Lean Control"))) #13
length(unique(which(dat$group=="Type 2 Diabetes"))) #18



dat <- dat %>% mutate(group2 = ifelse(group == 'Lean Control', 'Lean Control', 
                                      ifelse(epic_sglti2_1 == 'Yes', 'T2D-SGLTi2', 'T2D-No SGLTi2')))



table1::table1(~age + sex + bmi+epic_mfm_1+epic_insulin_1+epic_sglti2_1+epic_glp1ra_1| group2,data=dat)
dat$record_id



ids <- c(dat$kit_id)
so_kpmp_sc <- subset(so_kpmp_sc, kit_id %in% ids)

meta_kidney_sc <-  so_kpmp_sc@meta.data
rownames(meta_kidney_sc) <- rownames(so_kpmp_sc@meta.data)

#Merge metadata from 83 participants at baseline into seurat object metadata
meta_kidney_sc <- meta_kidney_sc %>%
  left_join(dat,by="kit_id")
rownames(meta_kidney_sc) <- rownames(so_kpmp_sc@meta.data)

#Merge metadata back into seurat object
so_kpmp_sc <- AddMetaData(so_kpmp_sc, meta_kidney_sc)


so_subset <- so_kpmp_sc

rm(so_kpmp_sc)






fixed_data <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/T2D_rawimagingdata.txt')

# #Calculate K2 and F variables
fixed_data <- fixed_data %>%
  mutate(avg_c_k2 = (lc_k2+rc_k2)/2) %>%
  mutate(avg_m_k2 = (lm_k2+rm_k2)/2) %>%
  mutate(avg_c_f = (lc_f+rc_f)/2) %>%
  mutate(avg_m_f = (lm_f+rm_f)/2)
fixed_data <- fixed_data %>%
  rowwise() %>%
  mutate(avg_c_k2_f = (avg_c_k2/avg_c_f)) %>%
  mutate(avg_m_k2_f = (avg_m_k2/avg_m_f))

med_c_k2 <- median(fixed_data$avg_c_k2, na.rm=T)
med_m_k2 <- median(fixed_data$avg_m_k2, na.rm=T)
med_c_f <- median(fixed_data$avg_c_f, na.rm=T)
med_m_f <- median(fixed_data$avg_m_f, na.rm=T)
med_c_k2_f <- median(fixed_data$avg_c_k2_f, na.rm=T)
med_m_k2_f <- median(fixed_data$avg_m_k2_f, na.rm=T)


fixed_data <- fixed_data %>% 
  rowwise() %>% 
  mutate(avg_c_k2_med = ifelse(avg_c_k2 >= med_c_k2, 'Above Median', 'Below Median'),
         avg_m_k2_med = ifelse(avg_m_k2 >= med_m_k2, 'Above Median', 'Below Median'),
         avg_c_f_med = ifelse(avg_c_f >= med_c_f, 'Above Median', 'Below Median'),
         avg_m_f_med = ifelse(avg_m_f >= med_m_f, 'Above Median', 'Below Median'),
         avg_c_k2_f_med = ifelse(avg_c_k2_f >= med_c_k2_f, 'Above Median', 'Below Median'),
         avg_m_k2_f_med = ifelse(avg_m_k2_f >= med_m_k2_f, 'Above Median', 'Below Median'))



#Try to add this data from dat 
so_subset@meta.data <- so_subset@meta.data[, !colnames(so_subset@meta.data) %in% c('avg_c_k2', 'avg_m_k2', 'avg_c_f', 
                                                                                   'avg_m_f', 
                                                                                   'avg_c_k2_f', 'avg_m_k2_f',
                                                                                   'lc_k2', 'rc_k2', 'lm_k2', 
                                                                                   'rm_k2', 'lm_f', 'rm_f',
                                                                                   'avg_c_k2_med', 'avg_m_k2_med',
                                                                                   'avg_c_f_med', 'avg_m_f_med', 
                                                                                   'avg_c_k2_f_med', 'avg_m_k2_f_med')]

so_subset@meta.data <- so_subset@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  left_join(fixed_data, by = "mrn") %>%
  tibble::column_to_rownames("cell_id")


so_subset$avg_c_k2_med <- factor(so_subset$avg_c_k2_med)
so_subset$avg_c_k2_med <- relevel(so_subset$avg_c_k2_med,"Below Median")

so_subset$avg_m_k2_med <- factor(so_subset$avg_m_k2_med)
so_subset$avg_m_k2_med <- relevel(so_subset$avg_m_k2_med,"Below Median")

so_subset$avg_c_f_med <- factor(so_subset$avg_c_f_med)
so_subset$avg_c_f_med <- relevel(so_subset$avg_c_f_med,"Below Median")

so_subset$avg_m_f_med <- factor(so_subset$avg_m_f_med)
so_subset$avg_m_f_med <- relevel(so_subset$avg_m_f_med,"Below Median")

so_subset$avg_c_k2_f_med <- factor(so_subset$avg_c_k2_f_med)
so_subset$avg_c_k2_f_med <- relevel(so_subset$avg_c_k2_f_med,"Below Median")

so_subset$avg_m_k2_f_med <- factor(so_subset$avg_m_k2_f_med)
so_subset$avg_m_k2_f_med <- relevel(so_subset$avg_m_k2_f_med,"Below Median")


save.image('C:/Users/netio/Documents/UofW/Rockies/Increased_N/LC_vs_T2D/LC_vs_T2D_Line250_Analysis.RData')






kidneyimaging_analysis <- function(celltype, genes, gene_list_name = 'TCA', median = F, adjustment = NULL,
                                   dir.results, cl_number = 1, cpc = 0.005, 
                                   set_cutoff = F, logFC_thresh = 10, pvalue_thresh = 0.1){
  
  if(median == F){
    k2_vars <- c("avg_c_k2","avg_m_k2","avg_c_f","avg_m_f","avg_c_k2_f","avg_m_k2_f")
  }else{
    k2_vars <- c("avg_c_k2_med","avg_m_k2_med","avg_c_f_med","avg_m_f_med","avg_c_k2_f_med","avg_m_k2_f_med")
  }
  
  if(celltype == 'PT'){
    so_celltype <- subset(so_subset,celltype2 == celltype)
    cat('PT Cells')
  }else if(celltype == 'TAL'){
    so_celltype <- subset(so_subset, TAL_celltype == celltype)
    cat('TAL Cells')
  }else if(celltype == 'DCT'){
    so_celltype <- subset(so_subset, DCT_celltype == celltype)
    cat('DCT Cells')
  }else{
    so_celltype <- subset(so_subset, KPMP_celltype == celltype)
    cat('Other celltypes')
  }
  # so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
  k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
  so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
  DefaultAssay(so_celltype) <- "RNA" 
  
  nrow(so_celltype) #34 genes
  ncol(so_celltype) #4926 PT cells
  
  counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round
  
  celltype <- str_replace_all(celltype, pattern='/', replacement = '_')
  
  # With parallelization
  #TCA Cycle
  # List of genes
  genes_list <- genes
  
  total_results <- data.frame()
  for (exposure in k2_vars) {
    cl <- makeCluster(cl_number)
    registerDoParallel(cl)
    
    start_time <- Sys.time()
    nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
      tryCatch({
        count_gene <- counts_path[g, , drop = FALSE]
        meta_gene <- subset(so_celltype,features=g)@meta.data
        
        if(!is.null(adjustment)){
          tmp.formula <- as.formula(paste0('~', exposure, '+', adjustment))
        }else{
          tmp.formula <- as.formula(paste0('~', exposure))
        }
        
        pred.formula <- as.formula(tmp.formula)
        pred_gene <- model.matrix(pred.formula, data = meta_gene)
        # library <- meta_gene$library_size
        library <- meta_gene$pooled_offset
        data_g_gene <- group_cell(count = count_gene, id = meta_gene$kit_id, pred = pred_gene,offset=library)
        
        if (is.null(data_g_gene)) {
          data_g_gene <- list(count = count_gene, id = meta_gene$kit_id, pred = pred_gene, offset = library)
        }
        
        #With offset
        result <- nebula(count = data_g_gene$count, id = data_g_gene$id, pred = data_g_gene$pred, 
                         ncore = 1, reml=T,model="NBLMM",output_re = T,covariance=T,
                         offset=data_g_gene$library, cpc= cpc)
        
        list(gene = g, result = result)  # return both gene name and result
        
      }, error = function(e) {
        list(gene = g, summary = NA, overdispersion = NA, convergence = NA, 
             algorithm = NA, covariance = NA, random_effect = NA)
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
        # Safely extract convergence code
        converged <- tryCatch({
          conv <- nebula_results_list[[gene_name]]$convergence
          if (is.null(conv) || length(conv) == 0) NA else conv
        }, error = function(e) NA)
        
        data.frame(Gene = gene_name,
                   Convergence_Code = converged)
      }
    )
    
    nebula_summaries <- map_dfr(
      names(nebula_results_list),
      function(gene_name) {
        tryCatch({
          # Check if the result exists and has summary info
          result <- nebula_results_list[[gene_name]]
          
          if (is.null(result) || is.null(result$summary)) {
            # If no summary info, return NULL (will be filtered out by map_dfr)
            return(NULL)
          } else {
            df <- result$summary
            
            # Check if summary is empty or not a data.frame
            if (is.null(df) || nrow(df) == 0 || !is.data.frame(df)) {
              return(NULL)
            } else {
              df <- df %>% mutate(Gene = gene_name)
              return(df)
            }
          }
        }, error = function(e) {
          # If any error occurs, return NULL (will be filtered out)
          cat("Error processing gene", gene_name, ":", e$message, "\n")
          return(NULL)
        })
      }
    )
    
    
    
    
    nonconverge_genes <- unique(PT_nebula_converged$Gene[which(PT_nebula_converged$Convergence_Code==-40)]) 
    
    #Make dataframe of final results
    full_results <- as.data.frame(nebula_summaries)
    
    if(nrow(full_results) == 0){
      next
    }
    
    if(is.null(adjustment)){
      colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept",
                                  "p_value","gene_id","gene","Gene")
    }else{
      colnames(full_results) <- c("logFC_Intercept","logFC",'logFC_SGLTi2', "se_Intercept",
                                  "se",'se_SGLTi2', "p_Intercept","p_value",'p_value_SGLTi2', "gene_id","gene",
                                  "Gene")
    }
    
    
    
    full_results$Variable <- exposure
    #Calculate number of genes filtered out for low expression 
    full_results$low_exp <- length(tca_genes)-length(full_results$gene)
    #Filter out non-converging genes
    full_results <- full_results %>% 
      filter(!gene %in%  nonconverge_genes)
    #Calculate nonconvergence rate
    full_results$nebula_nonconverged_percent <- paste0(round((1-(length(tca_genes)-length(nonconverge_genes))/length(tca_genes))*100,3),"%")
    # nebula_nonconverged_percent <- (length(rownames(counts_path))-length(unique(full_results$gene)))/length(rownames(counts_path))
    # print(paste0(nebula_nonconverged_percent*100, "% failed to converge"))
    full_results$fdr <- p.adjust(full_results$p_value,method="fdr")  
    # mutate(fdr3=p.adjust(PValue3,method="fdr"))
    full_results$PValue10 <- -log10(pmax(full_results[,6], 1e-10))
    
    
    
    total_results <- rbind(total_results,full_results)
  }
  
  if(is.null(adjustment)){
    file.name <- paste0(dir.results, 'NEBULA_', gene_list_name, 
                        '_median', median, '_', celltype, '_cells_PET_unadjusted_pooled_offset_T2D.txt')
  }else{
    file.name <- paste0(dir.results, 'NEBULA_', gene_list_name, 
                        '_median', median, '_', celltype, '_cells_PET_adjusted_pooled_offset_T2D.txt')
  }
  
  
  write.table(total_results,file.name, row.names=F, quote=F, sep='\t')
  
  cat(paste0(nrow(full_results), ' genes passed QC and were analysed in NEBULA \n'))
  # total_results <- read.csv(fs::path(dir.results,"NEBULA_TCA_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))
  
  # Define significance stars
  total_results <- total_results %>%
    mutate(signif = case_when(
      fdr < 0.01 ~ "**",
      fdr < 0.05 ~ "*",
      TRUE ~ ""
    ))
  
  # Select only the needed columns and rename LogFC for clarity
  if(median == T){
    subtitle1 <- paste0(celltype, ' (Median), ')
  }else{
    subtitle1 <- paste0(celltype, ', ')
  }
  
  if(!is.null(adjustment)){
    subtitle2 <- paste0('with adjustment for ', adjustment)
  }else{
    subtitle2 <- paste0('with no adjustment')
  }
  subtitle <- paste0(subtitle1, subtitle2)
  
  
  if(set_cutoff == TRUE){
    rem_index <- which(abs(total_results$logFC) > logFC_thresh & total_results$p_value > pvalue_thresh)
    total_results$logFC[rem_index] <- NA
  }
  
  
  heatmap_data <- total_results %>%
    dplyr::select(Gene, Variable, logFC, signif)
  
  if(median == TRUE){
    custom_order <- c("avg_c_k2_med", "avg_m_k2_med", "avg_c_f_med", "avg_m_f_med", "avg_c_k2_f_med", "avg_m_k2_f_med")
  }else{
    custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
  }
  
  heatmap_data$Variable <- factor(heatmap_data$Variable, levels = custom_order)
  custom_labels <- c("Average Cortical K2","Average Medulla K2","Average Cortical F","Average Medulla F",
                     "Average Cortical K2/F","Average Medulla K2/F")
  
  heat_map_p <- ggplot(heatmap_data, aes(x = Variable, y = Gene, fill = logFC)) +
    geom_tile(color = "grey90") +
    geom_text(aes(label = signif), size = 3, color = "black") +
    scale_fill_gradient2(
      low = "#264653", mid = "white", high = "darkred",
      midpoint = 0,
      name = "LogFC"
    ) +
    theme_minimal() +
    labs(title = paste0(gene_list_name, " Genes vs. PET Variables (T2D)"),
         subtitle = subtitle,
         x = "Exposure",
         y = "Gene") +
    scale_x_discrete(labels = setNames(custom_labels, custom_order))+
    theme(
      text = element_text(face="bold"),
      axis.text.y = element_text(size = 6),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.ticks.y = element_blank(),
      axis.title.x = element_blank(),
      panel.grid = element_blank()
    )
  
  # custom_colors <- c("#264653", "#2a9d8f", "#e9c46a", "#f4a261", "#e76f51", "darkred")
  
  if(is.null(adjustment)){
    file.name <- paste0(dir.results, 'NEBULA_', gene_list_name, 
                        'NEBULA_','median', median, '_', celltype, '_PET_unadjusted_pooled_offset_T2D.png')
  }else{
    file.name <- paste0(dir.results, 'NEBULA_', gene_list_name, 
                        'NEBULA_', 'median', median, '_', celltype, '_PET_adjusted_pooled_offset_T2D.png')  
  }
  
  
  # Ensure clean graphics device state
  while(dev.cur() > 1) dev.off()
  
  # Remove existing file if it exists
  if (file.exists(file.name)) {
    file.remove(file.name)
  }
  
  # Create PNG with error handling
  tryCatch({
    png(file.name, 
        width = 1500, height = 2000, res = 300)
    print(heat_map_p)
  }, finally = {
    dev.off()
  })
  
  # Verify file was created
  if(file.exists(file.name)) {
    cat("PNG file successfully created:", file.name, "\n")
  } else {
    cat("Failed to create PNG file:", file.name, "\n")
  }
  
}





















