library(aws.s3)
library(jsonlite)
library(biomaRt)
library(Seurat)
library(dplyr)
library(CellChat)
library(ComplexHeatmap)
library(circlize)

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

pb90_subset <- s3readRDS(object = "data_clean/subset/pb90_ckd_analysis_subset.rds", 
                         bucket = "scrna", region = "")

# Cell chat
# Set options for CellChat
options(stringsAsFactors = FALSE)

# Initialize CellChat object from Seurat object
# Extract expression data and metadata
data.input <- GetAssayData(pb90_subset, assay = "RNA", layer = "data")
meta <- pb90_subset@meta.data

# Create CellChat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "KPMP_celltype_general")

# Set the ligand-receptor interaction database
# Use human database (change to mouse if needed with "Mouse")
CellChatDB <- CellChatDB.human

# Show available categories in the database
showDatabaseCategory(CellChatDB)

# IMPORTANT: Use all categories to include metabolic interactions
# The "Metabolic" category contains metabolic signaling
CellChatDB.use <- CellChatDB  # Use full database

# Alternatively, if you want to focus on specific categories:
# CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling", "ECM-Receptor", 
#                                                    "Cell-Cell Contact", "Metabolic"))

# Set the database in CellChat object
cellchat@DB <- CellChatDB.use

# Preprocessing the expression data
# Identify overexpressed genes and interactions
cellchat <- subsetData(cellchat) # This step is essential
options(future.globals.maxSize = 8000 * 1024^2) # Increase to 8GB
future::plan("multisession", workers = 4) # Parallelize if needed

# Identify overexpressed ligands and receptors
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Project gene expression to protein-protein interaction network (optional but recommended)
cellchat <- projectData(cellchat, PPI.human)

# Compute communication probability and infer cellular communication network
# This analyzes both paracrine (between different cell types) and autocrine (same cell type) signaling
cellchat <- computeCommunProb(cellchat, 
                              type = "triMean",  # Use triMean for more robust results
                              population.size = TRUE)  # Consider cell population size

# Filter out communications if there are only few cells in certain groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Compute communication probability at signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

# Calculate aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

# Analyze metabolic signaling specifically
# Extract metabolic interactions
metabolic_pathways <- searchPair(signaling = "metabolic", 
                                 pairLR.use = cellchat@LR$LRsig, 
                                 key = "pathway_name", 
                                 matching.exact = FALSE, 
                                 pair.only = FALSE)

# Get all available signaling pathways
all_pathways <- cellchat@netP$pathways

# Identify which pathways are metabolic
metabolic_pathway_names <- grep("metabolic|glycolysis|lipid|amino acid|nucleotide|vitamin", 
                                all_pathways, 
                                ignore.case = TRUE, 
                                value = TRUE)

print(paste("Metabolic pathways identified:", paste(metabolic_pathway_names, collapse = ", ")))

# Visualization Section
# 1. Overall communication network
groupSize <- as.numeric(table(cellchat@idents))

# Visualize number of interactions
pdf("cellchat_interaction_counts.pdf", width = 8, height = 8)
netVisual_circle(cellchat@net$count, 
                 vertex.weight = groupSize, 
                 weight.scale = TRUE, 
                 label.edge = FALSE, 
                 title.name = "Number of interactions")
dev.off()

# Visualize interaction strength/weight
pdf("cellchat_interaction_strength.pdf", width = 8, height = 8)
netVisual_circle(cellchat@net$weight, 
                 vertex.weight = groupSize, 
                 weight.scale = TRUE, 
                 label.edge = FALSE, 
                 title.name = "Interaction strength")
dev.off()

# 2. Analyze autocrine vs paracrine signaling
# Create matrix to identify autocrine (diagonal) vs paracrine (off-diagonal)
interaction_matrix <- cellchat@net$count

# Extract autocrine interactions (diagonal elements)
autocrine_interactions <- diag(interaction_matrix)
names(autocrine_interactions) <- rownames(interaction_matrix)

# Extract paracrine interactions (off-diagonal elements)
paracrine_interactions <- interaction_matrix
diag(paracrine_interactions) <- 0

# Calculate total autocrine and paracrine interactions
total_autocrine <- sum(autocrine_interactions)
total_paracrine <- sum(paracrine_interactions)

print(paste("Total autocrine interactions:", total_autocrine))
print(paste("Total paracrine interactions:", total_paracrine))

# Plot autocrine vs paracrine comparison
autocrine_df <- data.frame(
  cell_type = names(autocrine_interactions),
  autocrine = autocrine_interactions,
  paracrine_sent = rowSums(paracrine_interactions),
  paracrine_received = colSums(paracrine_interactions)
)

# Create visualization for autocrine vs paracrine
library(reshape2)
autocrine_df_long <- melt(autocrine_df, id.vars = "cell_type")

pdf("autocrine_vs_paracrine_barplot.pdf", width = 10, height = 6)
ggplot(autocrine_df_long, aes(x = cell_type, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Autocrine vs Paracrine Signaling by Cell Type",
       x = "Cell Type", y = "Number of Interactions") +
  scale_fill_manual(values = c("autocrine" = "#FF6B6B", 
                               "paracrine_sent" = "#4ECDC4", 
                               "paracrine_received" = "#45B7D1"))
dev.off()

# 3. Focus on metabolic signaling pathways
if(length(metabolic_pathway_names) > 0) {
  pdf("metabolic_signaling_heatmap.pdf", width = 10, height = 8)
  for(pathway in metabolic_pathway_names) {
    tryCatch({
      netVisual_heatmap(cellchat, 
                        signaling = pathway, 
                        color.heatmap = "Reds",
                        title.name = paste("Metabolic pathway:", pathway))
    }, error = function(e) {
      print(paste("Could not plot", pathway, ":", e$message))
    })
  }
  dev.off()
}

# 4. Identify major signaling sources and targets
# Calculate centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# Visualize dominant senders and receivers
pdf("signaling_roles_heatmap.pdf", width = 10, height = 6)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", 
                                         title = "Outgoing signaling patterns")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", 
                                         title = "Incoming signaling patterns")
ht1 + ht2
dev.off()

# 5. Analyze specific metabolic interactions in detail
# Get detailed information about metabolic signaling
metabolic_communications <- subsetCommunication(cellchat, 
                                                grep("metabolic|glycolysis|lipid", 
                                                     cellchat@LR$LRsig$pathway_name, 
                                                     ignore.case = TRUE))

if(nrow(metabolic_communications) > 0) {
  write.csv(metabolic_communications, "metabolic_communications_detailed.csv", row.names = FALSE)
  
  # Visualize top metabolic interactions
  pdf("top_metabolic_interactions.pdf", width = 12, height = 8)
  top_metabolic <- metabolic_communications %>%
    group_by(source, target) %>%
    summarise(total_prob = sum(prob, na.rm = TRUE)) %>%
    arrange(desc(total_prob)) %>%
    head(20)
  
  ggplot(top_metabolic, aes(x = reorder(paste(source, "->", target), total_prob), 
                            y = total_prob, fill = source)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    theme_minimal() +
    labs(title = "Top 20 Metabolic Signaling Interactions",
         x = "Cell Type Interaction", 
         y = "Total Communication Probability")
  dev.off()
}

# 6. Identify communication patterns
# Identify global communication patterns
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = 3)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = 3)

# Visualize communication patterns
pdf("communication_patterns.pdf", width = 12, height = 5)
netAnalysis_river(cellchat, pattern = "outgoing")
netAnalysis_river(cellchat, pattern = "incoming")
dev.off()

# 7. Compare metabolic vs non-metabolic signaling
all_communications <- subsetCommunication(cellchat)
all_communications$is_metabolic <- grepl("metabolic|glycolysis|lipid|amino acid|nucleotide", 
                                         all_communications$pathway_name, 
                                         ignore.case = TRUE)

signaling_summary <- all_communications %>%
  group_by(source, target, is_metabolic) %>%
  summarise(
    n_interactions = n(),
    mean_prob = mean(prob, na.rm = TRUE),
    total_prob = sum(prob, na.rm = TRUE)
  )

# Plot comparison
pdf("metabolic_vs_canonical_signaling.pdf", width = 10, height = 6)
ggplot(signaling_summary, aes(x = source, y = target, size = total_prob)) +
  geom_point(aes(color = is_metabolic), alpha = 0.6) +
  scale_size_continuous(range = c(1, 10)) +
  scale_color_manual(values = c("FALSE" = "#3498db", "TRUE" = "#e74c3c"),
                     labels = c("Canonical", "Metabolic")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Metabolic vs Canonical Cell-Cell Communication",
       x = "Source Cell Type", y = "Target Cell Type",
       size = "Total Probability", color = "Signaling Type")
dev.off()

# 8. Save CellChat object for future use
saveRDS(cellchat, file = "cellchat_pb90_subset_analysis.rds")

# 9. Generate comprehensive report
sink("cellchat_analysis_summary.txt")
cat("=== CellChat Analysis Summary ===\n\n")
cat("Total number of cells:", nrow(meta), "\n")
cat("Cell types analyzed:", paste(unique(cellchat@idents), collapse = ", "), "\n\n")

cat("=== Interaction Statistics ===\n")
cat("Total interactions detected:", sum(cellchat@net$count), "\n")
cat("Total autocrine interactions:", total_autocrine, "\n")
cat("Total paracrine interactions:", total_paracrine, "\n")
cat("Ratio autocrine/paracrine:", round(total_autocrine/total_paracrine, 3), "\n\n")

cat("=== Metabolic Signaling ===\n")
cat("Metabolic pathways identified:", length(metabolic_pathway_names), "\n")
cat("Metabolic pathway names:", paste(metabolic_pathway_names, collapse = ", "), "\n")
if(nrow(metabolic_communications) > 0) {
  cat("Total metabolic interactions:", nrow(metabolic_communications), "\n")
  cat("Proportion of metabolic interactions:", 
      round(nrow(metabolic_communications)/nrow(all_communications) * 100, 2), "%\n")
}

cat("\n=== Top Signaling Pathways ===\n")
pathway_counts <- table(all_communications$pathway_name)
top_pathways <- sort(pathway_counts, decreasing = TRUE)[1:10]
for(i in 1:length(top_pathways)) {
  cat(names(top_pathways)[i], ":", top_pathways[i], "interactions\n")
}
sink()

print("CellChat analysis complete! Check the generated PDF files and summary text file.")

# Optional: Interactive visualization with CellChat's Shiny app
# cellchat <- runCellChatApp(cellchat)