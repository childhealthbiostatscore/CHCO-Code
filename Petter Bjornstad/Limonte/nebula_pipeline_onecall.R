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
# Remove NULL entries (failed analyses) from the results
re_summary$convergecode<-re$convergence
re_summary$converged<-ifelse(re_summary$convergecode==1,1,0)
nonconverge_genes<-subset(re_summary,re_summary$converged==0)
# Create a full results dataframe
full_results <- as.data.frame(re_summary)

# Calculate number of genes filtered out
low_exp <- length(genes_list) - length(full_results$gene)

# Filter out non-converging genes from the results
full_results<-subset(full_results,full_results$converged==1)

# Calculate non-convergence percentage
nebula_nonconverged_percent <- paste0(
  round(
    (1 -
       (length(genes_list) - nrow(nonconverge_genes)) / length(genes_list)) *
      100,
    3
  ),
  "%"
)

# Adjust p-values using False Discovery Rate (FDR) method
full_results <- full_results %>%
  mutate(fdr = p.adjust(`p_group1`, method = "fdr"))

# Calculate -log10(p-values), avoid log(0) by setting minimum p-value
full_results$PValue10 <- -log10(pmax(full_results$`p_group1`, 1e-10))

# Assign colors based on log fold-change and FDR significance
full_results$color <- ifelse(
  full_results$fdr < 0.05 & full_results$`logFC_group1` > 0,
  "lightcoral",
  ifelse(
    full_results$fdr < 0.05 & full_results$`logFC_group1` < 0,
    "lightblue",
    "gray"
  )
)

# Identify significantly differentially expressed genes (FDR < 0.05)
significant_df <- full_results[full_results$fdr < 0.05, ]

###8.3.4 Visualize Results
# ---- Metadata for Plotting ----

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

# Create the volcano plot using ggplot
volcano_plot <- ggplot(
  full_results,
  aes(x = logFC_group1, y = PValue10, color = color)
) +
  geom_point(alpha = 0.7) + # Add points with transparency
  scale_color_identity() + # Use color directly from 'color' column
  theme_minimal() + # Apply minimalistic theme
  labs(
    title = "Main title", # Main title
    subtitle = "subtitle", # Subtitle
    x = "logFC", # X-axis label
    y = "-log10(P-Value)", # Y-axis label
    color = "LogFC Direction",
    caption = paste0(
      "FDR < 0.05, Genes = ",
      Genes,
      ", Nuclei = ",
      Nuclei,
      ", Non-Convergence Rate: ",
      Nonconvergence_Rate,
      ", Genes Filtered out for Low Expression: ",
      low_exp
    ) # Plot caption
  ) +
  theme(
    plot.title = element_text(hjust = 0),
    axis.text.x = element_text(angle = 0, hjust = 1)
  ) +
  xlim(min, max) + # Set x-axis limits
  # Add gene labels for significant points
  geom_text(
    data = significant_df,
    aes(label = gene),
    vjust = 1,
    hjust = 1,
    size = 3,
    check_overlap = TRUE,
    color = "black"
  )

# Display the volcano plot
volcano_plot
