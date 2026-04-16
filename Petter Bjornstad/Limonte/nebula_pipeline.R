######## NEBULA PIPLINE FROM NOTION ##########


# ---- Parallelized Analysis ----

# Set up parallel cluster with 10 cores
cl <- makeCluster(10)
registerDoParallel(cl)

###8.2.2 Run Analysis###
# Record start time
start_time <- Sys.time()

# Perform analysis across genes in parallel
#g<-"VAMP3"
nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
  tryCatch({
    # Subset counts for the current gene
    count_gene <- counts_hvg[g, , drop = FALSE]
    
    # Subset metadata corresponding to the current gene
    meta_gene <- subset(seurat_object_hvg, features = g)@meta.data
    
    # Create predictor matrix using the factor variable
    pred_gene <- model.matrix(~group, data = meta_gene)
    
    # Group data by ID (batch information)
    data_g_gene <- group_cell(count = count_gene, id = meta_gene$record_id, pred = pred_gene)
    
    # Run the nebula model (with random effects)
    result <- nebula(count = count_gene, id = meta_gene$record_id, pred = pred_gene, ncore = 1, output_re = TRUE)
    
    # Return list containing gene name and its result
    list(gene = g, result = result)
    
  }, error = function(e) {
    # In case of error, return NULL
    NULL
  })
}

# Stop the parallel cluster
stopCluster(cl)

# Record end time
end_time <- Sys.time()

# Print time taken for the analysis
print(end_time - start_time)

###8.2.3 Post-Analysis Results Processing
# Remove NULL entries (failed analyses) from the results
nebula_results_list <- Filter(Negate(is.null), nebula_results_list)

# Assign gene names as list names
names(nebula_results_list) <- sapply(nebula_results_list, function(x) x$gene)

# Simplify the list to only contain the result objects
nebula_results_list <- lapply(nebula_results_list, function(x) x$result)

# Extract convergence status into a dataframe
hep_nebula_converged <- map_dfr(
  names(nebula_results_list),
  function(gene_name) {
    converged <- nebula_results_list[[gene_name]]$convergence
    df <- data.frame(Gene = gene_name, Convergence_Code = converged)
    return(df)
  }
)

# Extract model summaries into a dataframe
nebula_summaries <- map_dfr(
  names(nebula_results_list),
  function(gene_name) {
    df <- nebula_results_list[[gene_name]]$summary
    df <- df %>% mutate(Gene = gene_name)
    return(df)
  }
)

# Identify genes that did not converge (Convergence_Code == -40)
nonconverge_genes <- unique(hep_nebula_converged$Gene[which(hep_nebula_converged$Convergence_Code !=1)])

# ---- Final Cleanup ----

# Create a full results dataframe
full_results <- as.data.frame(nebula_summaries)

# Calculate number of genes filtered out due to low expression
low_exp <- length(genes_list) - length(full_results$gene)

# Filter out non-converging genes from the results
full_results <- full_results %>% 
  filter(!gene %in% nonconverge_genes)

# Calculate non-convergence percentage
nebula_nonconverged_percent <- paste0(round((1 - (length(genes_list)  - length(nonconverge_genes)) / length(genes_list) ) * 100, 3), "%")

# Adjust p-values using False Discovery Rate (FDR) method
full_results <- full_results %>%
  mutate(fdr = p.adjust(`p_group1`, method = "fdr"))

# Calculate -log10(p-values), avoid log(0) by setting minimum p-value
full_results$PValue10 <- -log10(pmax(full_results$`p_group1`, 1e-10))

# Assign colors based on log fold-change and FDR significance
full_results$color <- ifelse(full_results$fdr < 0.05 & full_results$`logFC_group1` > 0, "lightcoral",
                             ifelse(full_results$fdr < 0.05 & full_results$`logFC_group1` < 0, "lightblue", "gray"))

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
volcano_plot <- ggplot(full_results, aes(x = logFC_group1, y = PValue10, color = color)) +
  geom_point(alpha = 0.7) +  # Add points with transparency
  scale_color_identity() +  # Use color directly from 'color' column
  theme_minimal() +  # Apply minimalistic theme
  labs(
    title = "Main title",  # Main title
    subtitle = "subtitle",  # Subtitle
    x = "logFC",  # X-axis label
    y = "-log10(P-Value)",  # Y-axis label
    color = "LogFC Direction",
    caption = paste0(
      "FDR < 0.05, Genes = ", Genes, 
      ", Nuclei = ", Nuclei, 
      ", Non-Convergence Rate: ", Nonconvergence_Rate, 
      ", Genes Filtered out for Low Expression: ", low_exp
    )  # Plot caption
  ) +
  theme(
    plot.title = element_text(hjust = 0),
    axis.text.x = element_text(angle = 0, hjust = 1)
  ) +
  xlim(min, max) +  # Set x-axis limits
  # Add gene labels for significant points
  geom_text(data = significant_df, aes(label = gene),
            vjust = 1, hjust = 1, size = 3, check_overlap = TRUE, color = "black")

# Display the volcano plot
volcano_plot

