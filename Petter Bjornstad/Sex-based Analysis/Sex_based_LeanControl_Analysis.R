#### Sex-Specific Analyses in Lean Controls




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
#library(table1)
library(clusterProfiler)
library('org.Hs.eg.db')







load('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/No_Med_line700.Rdata')

so_subset <- so_kpmp_sc
remove(so_kpmp_sc)

#dat_groups <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/ROCKIES_GroupAssignments.txt')
#dat_groups <- dat_groups %>% filter(group2 %in% c('Lean Control', 'T2D-No SGLTi2'))

#so_subset <- subset(so_subset, record_id == dat_groups$record_id)
test <- so_subset@meta.data %>% dplyr::select(record_id, group) %>% filter(!duplicated(record_id))








dir.results <- 'C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/LeanControl_Only/'


so_subset <- subset(so_subset, subset = record_id != 'CRC-55')
so_subset <- subset(so_subset, subset = group == 'Lean_Control')
test <- so_subset@meta.data %>% dplyr::select(record_id, group) %>% filter(!duplicated(record_id))






harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

dat <- harmonized_data %>% 
  dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(
    across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
    across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm=T))),
    .by = c(record_id, visit)
  )





## demographics tables 

dat_small <- dat %>% filter(record_id %in% test$record_id)



library(gt)
library(gtsummary)

desc_table1_fixed <- dat_small %>%
  select(age, sex, race_ethnicity, bmi, hba1c, study) %>%
  tbl_summary(
    by = sex,
    type = list(
      age ~ "continuous",
      bmi ~ "continuous", 
      hba1c ~ "continuous",
      race_ethnicity ~ "categorical",
      study ~ "categorical"
      
    ),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      age ~ 1,
      bmi ~ 1,
      hba1c ~ 2,
      all_categorical() ~ c(0, 1)
    ),
    label = list(
      age ~ "Age, years",
      race_ethnicity ~ "Race/Ethnicity",
      bmi ~ "BMI, kg/m²",
      hba1c ~ "HbA1c, %",
      study ~ "Study"
    ),
    missing_text = "Missing"
  ) %>%
  add_p(test = list(
    all_continuous() ~ "t.test"
    # Skip categorical p-values if they cause issues
  )) %>%
  add_overall(col_label = "**Overall**\nN = {N}") %>%
  modify_header(label ~ "**Characteristic**") %>%
  modify_spanning_header(all_stat_cols() ~ "**Group**") %>%
  modify_footnote(all_stat_cols() ~ "Mean (SD) for continuous variables; n (%) for categorical variables")

# Save version with epic
desc_table1_fixed %>%
  as_gt() %>%
  tab_options(
    table.font.size = 11,
    heading.title.font.size = 14,
    column_labels.font.size = 12
  ) %>%
  gtsave(paste0(dir.results, "LC_demographics.png"), 
         vwidth = 1200, vheight = 800)












#### Celltype Analyses






#function
LC_NEBULA_Analysis <- function(so_subset, dir.results, celltype, genes){
  if(celltype == 'All'){
    so_celltype <- so_subset 
  }else if(celltype %in% c('TAL', 'EC', 'POD', 'PT')){
    so_celltype <- subset(so_subset,celltype2==celltype)
    DefaultAssay(so_celltype) <- "RNA" 
  }else if(celltype == 'IC'){
    so_celltype <- subset(so_subset, KPMP_celltype %in% c(
      "cDC",
      "cycT",
      "CD4+ T",
      "CD8+ T",
      "NK",
      "B",
      "MON",
      "MAC",
      "MC"))
  }else if(celltype == 'DCTall'){
    so_celltype <- subset(so_subset, DCT_celltype=='DCT')
  }
  
  
  
  
  full_analysis <- FindVariableFeatures(so_celltype, selection.method = "vst", nfeatures = 2000)
  hvgs <- VariableFeatures(full_analysis)
  so_celltype <- subset(so_celltype, features = hvgs)

  
  celltype2 <- str_replace_all(celltype,"/","_")
  celltype2 <- str_replace_all(celltype2,"-","_")
  

  
  print(paste0(celltype2, ' is running.'))
  #FOR LOOP OVER VARIABLES 
 
    counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round
    count_gene <- counts_path
    meta_gene <- subset(so_celltype)@meta.data
    
    
    
    complete_idx <- complete.cases(meta_gene$sex)
    cat("Cells with complete data:", sum(complete_idx), "\n")
    
    # Step 3: Filter all your data to only include cells with complete predictor data
    meta_gene <- meta_gene[complete_idx, ]
    count_gene <- count_gene[, complete_idx]  # Note: subsetting columns for cells
    
    
    num_cells <- nrow(meta_gene)
    num_part <- unique(meta_gene$record_id) %>% length()
    
    tmp_df <- meta_gene %>% dplyr::select(record_id, group, sex) %>% filter(!duplicated(record_id))
    
    num_male <- tmp_df %>% filter(sex == 'Male') %>% nrow()
    num_female <- tmp_df %>% filter(sex == 'Female') %>% nrow()
    
    # Step 4: Create prediction matrix from the complete data
    pred_gene <- model.matrix(~sex, data = meta_gene)
    
    
    # library <- meta_gene$library_size
    library <- meta_gene$pooled_offset
    data_g_gene <- group_cell(count = count_gene, id = meta_gene$kit_id, pred = pred_gene,offset=library)
    
    if (is.null(data_g_gene)) {
      data_g_gene <- list(count = count_gene, id = meta_gene$kit_id, pred = pred_gene, offset = library)
    }
    
    #With offset
    result <- nebula(count = data_g_gene$count, id = data_g_gene$id, 
                     pred = data_g_gene$pred, ncore = 1, reml=T,model="NBLMM",output_re = T,covariance=T,offset=data_g_gene$library)
    
    
    
    #Make dataframe of final results
    full_results <- as.data.frame(result)
    #Calculate number of genes filtered out for low expression 
    
    
    full_results$num_cells <- num_cells
    full_results$num_male <- num_male
    full_results$num_female <- num_female
    
    write.table(full_results,paste0(dir.results,"NEBULA_", 
                                    celltype2, "_cells__LC_pooledoffset.csv"),
                row.names=F, quote=F, sep=',')

 
  print(paste0(celltype2, ' is done.'))
  
  
}


celltypes_vec <- c('All', 
                   'PT', 
                   'TAL', 
                   'EC', 
                   'DCTall', 
                   'IC')

for(i  in 1:length(celltypes_vec)){
  LC_NEBULA_Analysis(so_subset = so_subset, dir.results = dir.results, celltype = celltypes_vec[i], genes = gene_list)
  print(celltypes_vec[i])
}







### Full Analysis (all Genes)



#function
LC_NEBULA_Analysis_full <- function(so_subset, dir.results, celltype, genes){
  if(celltype == 'All'){
    so_celltype <- so_subset 
  }else if(celltype %in% c('TAL', 'EC', 'POD', 'PT')){
    so_celltype <- subset(so_subset,celltype2==celltype)
    DefaultAssay(so_celltype) <- "RNA" 
  }else if(celltype == 'IC'){
    so_celltype <- subset(so_subset, KPMP_celltype %in% c(
      "cDC",
      "cycT",
      "CD4+ T",
      "CD8+ T",
      "NK",
      "B",
      "MON",
      "MAC",
      "MC"))
  }else if(celltype == 'DCTall'){
    so_celltype <- subset(so_subset, DCT_celltype=='DCT')
  }
  
  
  
  
 # full_analysis <- FindVariableFeatures(so_celltype, selection.method = "vst", nfeatures = 2000)
#  hvgs <- VariableFeatures(full_analysis)
#  so_celltype <- subset(so_celltype, features = hvgs)
  
  
  celltype2 <- str_replace_all(celltype,"/","_")
  celltype2 <- str_replace_all(celltype2,"-","_")
  
  
  
  print(paste0(celltype2, ' is running.'))
  #FOR LOOP OVER VARIABLES 
  
  counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round
  count_gene <- counts_path
  meta_gene <- subset(so_celltype)@meta.data
  
  
  
  complete_idx <- complete.cases(meta_gene$sex)
  cat("Cells with complete data:", sum(complete_idx), "\n")
  
  # Step 3: Filter all your data to only include cells with complete predictor data
  meta_gene <- meta_gene[complete_idx, ]
  count_gene <- count_gene[, complete_idx]  # Note: subsetting columns for cells
  
  
  num_cells <- nrow(meta_gene)
  num_part <- unique(meta_gene$record_id) %>% length()
  
  tmp_df <- meta_gene %>% dplyr::select(record_id, group, sex) %>% filter(!duplicated(record_id))
  
  num_male <- tmp_df %>% filter(sex == 'Male') %>% nrow()
  num_female <- tmp_df %>% filter(sex == 'Female') %>% nrow()
  
  # Step 4: Create prediction matrix from the complete data
  pred_gene <- model.matrix(~sex, data = meta_gene)
  
  
  # library <- meta_gene$library_size
  library <- meta_gene$pooled_offset
  data_g_gene <- group_cell(count = count_gene, id = meta_gene$kit_id, pred = pred_gene,offset=library)
  
  if (is.null(data_g_gene)) {
    data_g_gene <- list(count = count_gene, id = meta_gene$kit_id, pred = pred_gene, offset = library)
  }
  
  #With offset
  result <- nebula(count = data_g_gene$count, id = data_g_gene$id, 
                   pred = data_g_gene$pred, ncore = 1, reml=T,model="NBLMM",output_re = T,covariance=T,offset=data_g_gene$library)
  
  
  
  #Make dataframe of final results
  full_results <- as.data.frame(result)
  #Calculate number of genes filtered out for low expression 
  
  
  full_results$num_cells <- num_cells
  full_results$num_male <- num_male
  full_results$num_female <- num_female
  
  write.table(full_results,paste0(dir.results,"Full_NEBULA_", 
                                  celltype2, "_cells__LC_pooledoffset.csv"),
              row.names=F, quote=F, sep=',')
  
  
  print(paste0(celltype2, ' is done.'))
  
  
}


celltypes_vec <- c('All', 
                   'PT', 
                   'TAL', 
                   'EC', 
                   'DCTall', 
                   'IC')

for(i  in 1:length(celltypes_vec)){
  LC_NEBULA_Analysis_full(so_subset = so_subset, dir.results = dir.results, celltype = celltypes_vec[i], genes = gene_list)
  print(celltypes_vec[i])
}









########### Volcano plots 

library(ggplot2)
library(ggbreak)
library(dplyr)

variable_names <- c('All', 'PT', 'TAL', 'EC', 'IC', 'DCTall')

for(i in c(1:length(variable_names))){
  
  sig_markers <- data.table::fread(paste0('/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/LeanControl_Only/Full_NEBULA_', 
                                          variable_names[i], '_cells__LC_pooledoffset.csv'))
  
  sig_markers <- sig_markers %>% dplyr::select(Gene = summary.gene,
                                               LogFC = summary.logFC_sexMale, 
                                               Pvalue = summary.p_sexMale)
  
  tmp_df <- sig_markers
  tmp_df$diffexp <- 'No'
  tmp_df$diffexp[tmp_df$Pvalue < 0.05 & tmp_df$LogFC > 0] <- 'Up'
  tmp_df$diffexp[tmp_df$Pvalue < 0.05 & tmp_df$LogFC < 0] <- 'Down'
  
  tmp_df <- tmp_df %>% arrange(Pvalue)
  tmp_df$label <- NA
  tmp_df$label[1:10] <- tmp_df$Gene[1:10]
  
  tmp_df <- tmp_df %>% filter(abs(LogFC) < 10)
  
  # Automatically detect if axis break is needed
  logfc_range <- range(tmp_df$LogFC, na.rm = TRUE)
  logfc_q95 <- quantile(abs(tmp_df$LogFC), 0.95, na.rm = TRUE)
  logfc_max <- max(abs(tmp_df$LogFC), na.rm = TRUE)
  
  # Determine if we need a break (outliers beyond 95th percentile create large gaps)
  needs_break <- logfc_max > (logfc_q95 * 2)
  
  # Calculate break points if needed
  if(needs_break){
    # Find gaps in the data
    sorted_logfc <- sort(abs(tmp_df$LogFC))
    # Calculate gaps between consecutive values
    gaps <- diff(sorted_logfc)
    # Find the largest gap in the outer 10% of data
    gap_threshold <- quantile(abs(tmp_df$LogFC), 0.90, na.rm = TRUE)
    outer_data_idx <- which(sorted_logfc > gap_threshold)
    
    if(length(outer_data_idx) > 1){
      outer_gaps <- gaps[outer_data_idx[-length(outer_data_idx)]]
      max_gap_idx <- which.max(outer_gaps)
      
      # Set break points around the largest gap
      break_start <- sorted_logfc[outer_data_idx[max_gap_idx]] + 0.1
      break_end <- sorted_logfc[outer_data_idx[max_gap_idx + 1]] - 0.1
      
      # Ensure breaks are reasonable
      if(break_end - break_start > 0.5){
        use_break <- TRUE
      } else {
        use_break <- FALSE
      }
    } else {
      use_break <- FALSE
    }
  } else {
    use_break <- FALSE
  }
  
  # Making graph
  if(length(unique(tmp_df$diffexp)) > 1){
    tmp_graph <- ggplot(tmp_df, aes(x= LogFC, y=-log10(Pvalue), col = diffexp, label=label))+
      geom_point()+
      geom_text(size=2, vjust = 2, color='black')+
      scale_color_manual(values = c('orange', 'grey', 'purple'),
                         labels = c('Downregulated', 'Not significant', 'Upregulated'))+
      geom_hline(yintercept = -log10(0.05), col='blue', linetype='dashed')+
      geom_vline(xintercept = c(0), col='black', linetype ='dashed')+
      theme_classic()+
      labs(x='LogFC', y='-log10 pvalue', col ='Differential Expression', 
           title = paste0('Sex Differences in ',  variable_names[i],' Cells'))
    
    # Add break if needed
    if(use_break){
      tmp_graph <- tmp_graph + scale_x_break(breaks = c(break_start, break_end), scales = 0.3)
    }
    
  } else {
    tmp_graph <- ggplot(tmp_df, aes(x= LogFC, y=-log10(Pvalue), col = diffexp, label=label))+
      geom_point()+
      geom_text(size=2, vjust = 2, color='black')+
      scale_color_manual(values = c('grey'),
                         labels = c('Not significant'))+
      geom_hline(yintercept = -log10(0.05), col='blue', linetype='dashed')+
      geom_vline(xintercept = c(0), col='black', linetype ='dashed')+
      theme_classic()+
      labs(x='LogFC', y='-log10 P-value', col ='Differential Expression', 
           title = paste0('Sex Differences in ',  variable_names[i],' Cells'))
    
    # Add break if needed
    if(use_break){
      tmp_graph <- tmp_graph + scale_x_break(breaks = c(break_start, break_end), scales = 0.3)
    }
  }
  
  pdf(paste0('/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/LeanControl_Only/VolcanoPlots_', 
             variable_names[i], '_Cells_LCOnly.pdf'))
  print(tmp_graph)
  dev.off()
  
  png(paste0('/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/LeanControl_Only/VolcanoPlots_', 
             variable_names[i], '_Cells_LCOnly.png'), width = 3000, height = 3000, res = 300)
  print(tmp_graph)
  dev.off()
  
  print(paste0('Plot done for ', variable_names[i], 
               if(use_break) ' (with axis break)' else ''))
}








































##### GSEA








# Load required libraries
library(fgsea)
library(msigdbr)
library(ggplot2)
library(dplyr)
library(patchwork)  # Optional: for combining plots

# Define cell types
celltypes_vec <- c('All', 'PT', 'TAL', 'EC', 'DCTall', 'IC')

# Get gene sets from MSigDB for human
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(hallmark_sets$gene_symbol, hallmark_sets$gs_name)

# Get GO gene sets
go_bp_sets <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
go_bp_list <- split(go_bp_sets$gene_symbol, go_bp_sets$gs_name)

go_cc_sets <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:CC")
go_cc_list <- split(go_cc_sets$gene_symbol, go_cc_sets$gs_name)

go_mf_sets <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF")
go_mf_list <- split(go_mf_sets$gene_symbol, go_mf_sets$gs_name)

# Define gene set types
geneset_types <- list(
  Hallmark = hallmark_list,
  GO_BP = go_bp_list,
  GO_CC = go_cc_list,
  GO_MF = go_mf_list
)

# Create a list to store all results
all_gsea_results <- list()

# Loop through each cell type
for(celltype in celltypes_vec) {
  
  cat("\n=== Processing", celltype, "===\n")
  
  # Load your results file for this cell type
  de_results <- read.csv(paste0(dir.results, 'Full_NEBULA_', 
                                celltype, '_cells__LC_pooledoffset.csv'))
  
  # Create ranked gene list
  ranked_genes <- de_results %>%
    dplyr::select(logFC = summary.logFC_sexMale, 
                  gene_id = summary.gene) %>% 
    arrange(desc(logFC)) %>%
    pull(logFC, name = gene_id)
  
  # Loop through each gene set type
  for(geneset_name in names(geneset_types)) {
    
    cat("Running GSEA for", geneset_name, "...\n")
    
    # Run fGSEA
    gsea_results <- fgsea(
      pathways = geneset_types[[geneset_name]],
      stats = ranked_genes,
      minSize = 15,
      maxSize = 500,
      nperm = 10000
    )
    
    # Add cell type and gene set type columns
    gsea_results$celltype <- celltype
    gsea_results$geneset_type <- geneset_name
    
    # Store results
    result_name <- paste0(celltype, "_", geneset_name)
    all_gsea_results[[result_name]] <- gsea_results
    
    # Select top significant pathways
    top_pathways <- gsea_results %>%
      filter(pval < 0.05) %>%
      arrange(pval) %>%
      slice_head(n = 20) %>%
      mutate(pathway_clean = gsub("HALLMARK_", "", pathway),
             pathway_clean = gsub("GOBP_", "", pathway_clean),
             pathway_clean = gsub("GOCC_", "", pathway_clean),
             pathway_clean = gsub("GOMF_", "", pathway_clean),
             pathway_clean = gsub("_", " ", pathway_clean),
             pathway_clean = tools::toTitleCase(tolower(pathway_clean)))
    
    # Only create plot if there are significant pathways
    if(nrow(top_pathways) > 0) {
      
      # Create dotplot
      p <- ggplot(top_pathways, aes(x = NES, y = reorder(pathway_clean, NES))) +
        geom_point(aes(size = -log10(pval), color = NES)) +
        scale_color_gradient2(
          low = "blue", 
          mid = "white", 
          high = "red", 
          midpoint = 0, 
          name = "NES"
        ) +
        scale_size_continuous(
          name = "-log10(pval)",
          range = c(2, 10)
        ) +
        theme_bw() +
        theme(
          axis.text.y = element_text(size = 9, color = "black"),
          axis.text.x = element_text(size = 9, color = "black"),
          panel.grid.major.y = element_line(color = "grey90"),
          panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 12)
        ) +
        labs(
          x = "Normalized Enrichment Score (NES)",
          y = "",
          title = paste0(celltype, " - ", geneset_name, " (Top 20)")
        ) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5)
      
      # Save individual plot
      ggsave(paste0(dir.results, 'GSEA/', geneset_name, '_gsea_', celltype, ".pdf"), 
             plot = p, width = 10, height = 8)
      
      cat("Saved plot for", celltype, "-", geneset_name, "\n")
    } else {
      cat("No significant pathways found for", celltype, "-", geneset_name, "\n")
    }
  }
}

# Combine all GSEA results into one dataframe
# Before saving, convert list columns to character strings
combined_gsea_save <- combined_gsea %>%
  mutate(leadingEdge = sapply(leadingEdge, function(x) paste(x, collapse = ";")))

# Now save combined results
write.csv(combined_gsea_save, 
          paste0(dir.results, "GSEA/all_celltypes_all_genesets_gsea_results.csv"), 
          row.names = FALSE)

# Save separate files for each gene set type
for(geneset_name in names(geneset_types)) {
  subset_results <- combined_gsea %>% 
    filter(geneset_type == geneset_name) %>%
    mutate(leadingEdge = sapply(leadingEdge, function(x) paste(x, collapse = ";")))
  
  write.csv(subset_results, 
            paste0(dir.results, "GSEA/gsea_results_", geneset_name, ".csv"), 
            row.names = FALSE)
}
cat("\n=== All analyses complete! ===\n")
cat("Total number of analyses:", length(all_gsea_results), "\n")

# Summary of significant pathways per cell type and gene set
summary_table <- combined_gsea %>%
  mutate(leadingEdge = sapply(leadingEdge, function(x) paste(x, collapse = ";"))) %>% 
  filter(pval < 0.05) %>%
  group_by(celltype, geneset_type) %>%
  summarise(n_significant = n(), .groups = 'drop') %>%
  tidyr::pivot_wider(names_from = geneset_type, values_from = n_significant, values_fill = 0)

print(summary_table)
write.csv(summary_table, 
          paste0(dir.results, "GSEA/summary_significant_pathways.csv"), 
          row.names = FALSE)














# Install required package if not already installed
if (!require("pdftools")) {
  install.packages("pdftools")
}

library(pdftools)

# Function to convert a single PDF to PNG
convert_pdf_to_png <- function(pdf_path, output_dir = NULL, dpi = 300) {
  # If no output directory specified, use same directory as PDF
  if (is.null(output_dir)) {
    output_dir <- dirname(pdf_path)
  }
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Get the base filename without extension
  base_name <- tools::file_path_sans_ext(basename(pdf_path))
  
  # Convert PDF to PNG
  # This creates one PNG per page
  png_files <- pdf_convert(
    pdf = pdf_path,
    format = "png",
    dpi = dpi,
    filenames = file.path(output_dir, paste0(base_name, "_page_%d.png"))
  )
  
  return(png_files)
}

# Convert all PDFs in a folder
convert_all_pdfs <- function(folder_path, output_dir = NULL, dpi = 300) {
  # Get all PDF files in the folder
  pdf_files <- list.files(folder_path, pattern = "\\.pdf$", 
                          full.names = TRUE, ignore.case = TRUE)
  
  if (length(pdf_files) == 0) {
    message("No PDF files found in the specified folder.")
    return(NULL)
  }
  
  message(paste("Found", length(pdf_files), "PDF file(s)"))
  
  # Convert each PDF
  all_png_files <- list()
  for (pdf in pdf_files) {
    message(paste("Converting:", basename(pdf)))
    png_files <- convert_pdf_to_png(pdf, output_dir, dpi)
    all_png_files[[basename(pdf)]] <- png_files
  }
  
  message("Conversion complete!")
  return(all_png_files)
}

# Example usage:
# Convert all PDFs in a folder
folder_path <- "/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/LeanControl_Only/GSEA/"

# Or convert to a specific output directory
 converted_files <- convert_all_pdfs(folder_path, output_dir = "/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/LeanControl_Only/GSEA/")

# Or convert with higher resolution
# converted_files <- convert_all_pdfs(folder_path, dpi = 600)









 
 #### Proteomics and Metabolomics 
 
 
 
 
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
 #library(table1)
 library(clusterProfiler)
 library('org.Hs.eg.db')
 
 
 
 



 
 
 
 harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')
 
 dat <- harmonized_data %>% 
   dplyr::select(-dob) %>% 
   arrange(date_of_screen) %>% 
   dplyr::summarise(
     across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
     across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm=T))),
     .by = c(record_id, visit)
   )
 
 dat <- dat %>% filter(group == 'Lean Control')
 
 
 data_dictionary <- readxl::read_xlsx('/Users/netio/Downloads/data_dictionary_master.xlsx')
 
 form_names <- unique(data_dictionary$form_name)
 proteo <- form_names[str_which(form_names, pattern = 'proteom')]
 metab <- form_names[str_which(form_names, pattern = 'metab')]
 
 
 data_dictionary_small <- data_dictionary %>% 
   filter()
 
 




