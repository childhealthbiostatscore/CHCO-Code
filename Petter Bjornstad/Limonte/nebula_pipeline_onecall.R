######## NEBULA PIPLINE EDITED WITH ONE NEBULA CALL ##########


# Record start time
start_time <- Sys.time()

# Nebula all genes together
seuratdata = as.SingleCellExperiment(seurat_object_hvg)
seuratdata <- scToNeb(
  obj = seuratdata,
  assay = "RNA",
  id = "record_id",
  pred = "group",
  offset = "nCount_RNA"
)
seuratdata$count = round(seuratdata$count)
seuratdata$pred = model.matrix(~group, seuratdata$pred)

# Subset count matrix to genes in genes_list
gene_idx <- rownames(seuratdata$count) %in% genes_list
count_subset <- seuratdata$count[gene_idx, ]

re = nebula(
  count_subset,
  seuratdata$id,
  pred = seuratdata$pred,
  ncore = 8,
  offset=seuratdata$offset
)


# Record end time
end_time <- Sys.time()

# Print time taken for the analysis
print(end_time - start_time)

###8.2.3 Post-Analysis Results Processing
re_summary<-re$summary
# Calculate number of genes filtered out
nebula_genes<-nrow(re_summary)
low_exp <- nrow(count_subset) - nebula_genes

# Remove NULL entries (failed analyses) from the results
re_summary$convergecode<-re$convergence
re_summary$converged<-ifelse(re_summary$convergecode==1,1,0)

# Calculate number of nonconverged genes:
nonconverge_genes<-nrow(subset(re_summary,re_summary$converged==0))

# Create a dataframe of converged nebula results:
full_results<-subset(re_summary,re_summary$converged==1)

# Calculate non-convergence percentage
nebula_nonconverged_percent <- paste0(
  round(
    (1 -
       (nrow(count_subset) - nonconverge_genes) / nrow(count_subset)) *
      100,
    3
  ),
  "%"
)

# Adjust p-values using False Discovery Rate (FDR) method
full_results$p_unadjusted<-full_results$p_group1
full_results$p_adjusted<-p.adjust(full_results$p_unadjusted, method = "fdr")

# Calculate -log10(p-values), avoid log(0) by setting minimum p-value
full_results$p10_unadjusted <- -log10(pmax(full_results$p_unadjusted, 1e-10))
full_results$p10_adjusted <- -log10(pmax(full_results$p_adjusted, 1e-10))

# Assign colors based on log fold-change p-value significance: 
full_results$color_unadjusted <- ifelse(
  full_results$p_unadjusted < 0.05 & full_results$logFC_group1 > 0,
  "lightcoral",
  ifelse(
    full_results$p_unadjusted < 0.05 & full_results$logFC_group1 < 0,
    "lightblue",
    "gray"
  )
)

full_results$color_adjusted <- ifelse(
  full_results$p_adjusted< 0.05 & full_results$logFC_group1 > 0,
  "lightcoral",
  ifelse(
    full_results$p_adjusted < 0.05 & full_results$logFC_group1 < 0,
    "lightblue",
    "gray"
  )
)

###8.3.4 Visualize Results
# ---- Metadata for Plotting ----
significant_df_adjusted <- full_results[full_results$p_adjusted < 0.05, ]
significant_df_unadjusted <- full_results[full_results$p_unadjusted < 0.05, ]

# Number of genes after filtering
Genes <- length(unique(full_results$gene))

# Number of nuclei (cells)
Nuclei <- ncol(seurat_object_hvg)

# Non-convergence rate
Nonconvergence_Rate <- nebula_nonconverged_percent

# Define x-axis limits for plotting
max <- max(full_results$logFC_group1)
min <- min(full_results$logFC_group1)

# ---- Volcano Plot ----
volcano_plot_adjusted <- ggplot(full_results, aes(x = logFC_group1, y = p10_adjusted, color = color_adjusted)) +
  geom_point(alpha = 0.7) +  # Add points with transparency
  scale_color_identity() +  # Use color directly from 'color' column
  theme_minimal() +  # Apply minimalistic theme
  labs(
    title = "FDR Adjusted Pvalues",  # Main title
    subtitle = "",  # Subtitle
    x = "logFC",  # X-axis label
    y = "-log10(P-Value)",  # Y-axis label
    color = "LogFC Direction",
    caption = paste0(
      "Genes = ", Genes, 
      ", Nuclei = ", Nuclei, 
      " (Genes Filtered out for Low Expression: ", low_exp,
      ", Non-Convergence Rate: ", Nonconvergence_Rate, 
      
      ")"
    )  # Plot caption
  ) +
  theme(
    plot.title = element_text(hjust = 0),
    axis.text.x = element_text(angle = 0, hjust = 1)
  ) +
  xlim(min, max) +  # Set x-axis limits
  # Add gene labels for significant points
  geom_text(data = significant_df_adjusted, aes(label = gene),
            vjust = 1, hjust = 1, size = 3, check_overlap = TRUE, color = "black")


volcano_plot_unadjusted <- ggplot(full_results, aes(x = logFC_group1, y = p10_unadjusted, color = color_unadjusted)) +
  geom_point(alpha = 0.7) +  # Add points with transparency
  scale_color_identity() +  # Use color directly from 'color' column
  theme_minimal() +  # Apply minimalistic theme
  labs(
    title = "Unadjusted PValues",  # Main title
    subtitle = "",  # Subtitle
    x = "logFC",  # X-axis label
    y = "-log10(P-Value)",  # Y-axis label
    color = "LogFC Direction",
    caption = paste0(
      "Genes = ", Genes, 
      ", Nuclei = ", Nuclei, 
      " (Genes Filtered out for Low Expression: ", low_exp,
      ", Non-Convergence Rate: ", Nonconvergence_Rate, 
      
      ")"
    )  # Plot caption
  ) +
  theme(
    plot.title = element_text(hjust = 0),
    axis.text.x = element_text(angle = 0, hjust = 1)
  ) +
  xlim(min, max) +  # Set x-axis limits
  # Add gene labels for significant points
  geom_text(data = significant_df_unadjusted, aes(label = gene),
            vjust = 1, hjust = 1, size = 3, check_overlap = TRUE, color = "black")
