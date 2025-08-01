library(enrichplot)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)  # Human annotation, change for other species
library(DOSE)

library(ggplot2)
library(dplyr)

# Alternative libraries for IPA-style analysis
library(msigdbr)  # For MSigDB gene sets
library(fgsea)    # Fast GSEA implementation

library(dplyr)
library(stringr)


# ====================================================================
# METHOD 1: Using ReactomePA for Reactome pathway analysis
# ====================================================================

# Example: Prepare your gene list (differential expression results)
# Your data should have gene symbols and fold changes/statistics

results <- data.table::fread('C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/')



# Convert gene symbols to Entrez IDs (required for ReactomePA)
gene_symbols <- deg_data$gene_symbol
entrez_ids <- bitr(gene_symbols, 
                   fromType = "SYMBOL", 
                   toType = "ENTREZID", 
                   OrgDb = org.Hs.eg.db)

# Merge back with original data
deg_with_entrez <- deg_data %>%
  inner_join(entrez_ids, by = c("gene_symbol" = "SYMBOL"))

# Create named vector for GSEA (ranked gene list)
gene_list <- deg_with_entrez$log2FC
names(gene_list) <- deg_with_entrez$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

# Reactome Over-Representation Analysis (ORA)
reactome_ora <- enrichPathway(gene = names(gene_list)[abs(gene_list) > 1], # significant genes
                              pvalueCutoff = 0.05,
                              readable = TRUE,
                              organism = "human")

# Reactome Gene Set Enrichment Analysis (GSEA)
reactome_gsea <- gsePathway(geneList = gene_list,
                            pvalueCutoff = 0.05,
                            pAdjustMethod = "BH",
                            organism = "human")

# ====================================================================
# METHOD 2: Using MSigDB for IPA-style canonical pathways
# ====================================================================

# Get MSigDB gene sets (similar to IPA canonical pathways)
msigdb_sets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
# Other useful categories:
# C2:CP:KEGG - KEGG pathways
# C2:CP:BIOCARTA - BioCarta pathways  
# C2:CGP - Chemical and genetic perturbations
# C5:BP - GO Biological Process

# Convert to format for clusterProfiler
msig_t2g <- msigdb_sets %>% 
  select(gs_name, entrez_gene) %>%
  as.data.frame()

# Run enrichment analysis
msigdb_ora <- enricher(gene = names(gene_list)[abs(gene_list) > 1],
                       TERM2GENE = msig_t2g,
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

# GSEA with MSigDB
msigdb_gsea <- GSEA(geneList = gene_list,
                    TERM2GENE = msig_t2g,
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH")

# ====================================================================
# METHOD 3: Using fgsea for fast GSEA analysis
# ====================================================================

# Prepare gene sets for fgsea
fgsea_pathways <- split(msigdb_sets$entrez_gene, msigdb_sets$gs_name)

# Run fgsea
fgsea_results <- fgsea(pathways = fgsea_pathways,
                       stats = gene_list,
                       minSize = 15,
                       maxSize = 500,
                       nperm = 10000)

# ====================================================================
# VISUALIZATION: Creating dotplots
# ====================================================================

# 1. Basic dotplot for ORA results
dotplot(reactome_ora, showCategory = 20) +
  ggtitle("Reactome Pathway Over-Representation Analysis") +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 9))

# 2. Enhanced dotplot for GSEA results
dotplot(reactome_gsea, showCategory = 20, split = ".sign") +
  facet_grid(.~.sign) +
  ggtitle("Reactome Pathway GSEA") +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 9))

# 3. Custom dotplot with ggplot2
# Convert enrichment results to dataframe
if (!is.null(reactome_ora)) {
  ora_df <- as.data.frame(reactome_ora)
  
  # Create custom dotplot
  custom_dotplot_ora <- ora_df %>%
    slice_head(n = 20) %>%  # Top 20 pathways
    mutate(Description = forcats::fct_reorder(Description, Count)) %>%
    ggplot(aes(x = Count, y = Description)) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_gradient(low = "red", high = "blue", name = "Adj. P-value") +
    scale_size_continuous(name = "Gene Count", range = c(3, 8)) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 10),
      legend.position = "right"
    ) +
    labs(
      title = "Top Enriched Reactome Pathways",
      x = "Gene Count",
      y = "Pathway"
    )
  
  print(custom_dotplot_ora)
}

# 4. Custom dotplot for GSEA results
if (!is.null(reactome_gsea)) {
  gsea_df <- as.data.frame(reactome_gsea)
  
  custom_dotplot_gsea <- gsea_df %>%
    slice_head(n = 20) %>%
    mutate(Description = forcats::fct_reorder(Description, NES),
           Regulation = ifelse(NES > 0, "Upregulated", "Downregulated")) %>%
    ggplot(aes(x = NES, y = Description)) +
    geom_point(aes(size = setSize, color = p.adjust)) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.7) +
    scale_color_gradient(low = "red", high = "blue", name = "Adj. P-value") +
    scale_size_continuous(name = "Set Size", range = c(3, 8)) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 10),
      legend.position = "right"
    ) +
    labs(
      title = "Reactome Pathway GSEA Results",
      x = "Normalized Enrichment Score (NES)",
      y = "Pathway"
    )
  
  print(custom_dotplot_gsea)
}

# 5. Multi-panel comparison dotplot
create_comparison_dotplot <- function(ora_result, gsea_result, title = "Pathway Analysis") {
  
  # Prepare ORA data
  ora_df <- as.data.frame(ora_result) %>%
    slice_head(n = 15) %>%
    mutate(Analysis = "ORA",
           Score = Count,
           Description = stringr::str_wrap(Description, 40))
  
  # Prepare GSEA data  
  gsea_df <- as.data.frame(gsea_result) %>%
    slice_head(n = 15) %>%
    mutate(Analysis = "GSEA",
           Score = NES,
           Description = stringr::str_wrap(Description, 40))
  
  # Combine and plot
  combined_df <- bind_rows(
    ora_df %>% select(Description, Score, p.adjust, Analysis),
    gsea_df %>% select(Description, Score, p.adjust, Analysis)
  )
  
  ggplot(combined_df, aes(x = Score, y = Description)) +
    geom_point(aes(size = abs(Score), color = p.adjust)) +
    scale_color_gradient(low = "red", high = "blue", name = "Adj. P-value") +
    scale_size_continuous(name = "Score Magnitude", range = c(2, 6)) +
    facet_wrap(~Analysis, scales = "free_x") +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 8),
      strip.text = element_text(size = 12, face = "bold")
    ) +
    labs(title = title, x = "Score", y = "Pathway")
}

# ====================================================================
# ADDITIONAL VISUALIZATIONS
# ====================================================================

# 6. Enrichment map (network visualization)
if (!is.null(reactome_ora)) {
  # Create enrichment map
  emapplot(reactome_ora, showCategory = 30, cex_label_category = 0.6)
}

# 7. Gene-Concept Network
if (!is.null(reactome_ora)) {
  # Show relationships between genes and pathways
  cnetplot(reactome_ora, categorySize = "pvalue", foldChange = gene_list)
}

# 8. Heatmap of enriched pathways
if (!is.null(reactome_ora)) {
  heatplot(reactome_ora, showCategory = 20, foldChange = gene_list)
}