library(aws.s3)
library(jsonlite)
library(biomaRt)
library(Seurat)
library(dplyr)
library(CellChat)
library(ComplexHeatmap)
library(circlize)
library(reshape2)
library(ggalluvial)

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
cellchat <- createCellChat(object = data.input, meta = meta, 
                           group.by = "KPMP_celltype_general")
groupSize <- as.numeric(table(cellchat@idents))

# Set the ligand-receptor interaction database
# Use human database
# Set the database in CellChat object
cellchat@DB <- CellChatDB.human

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

# visualize the aggregated cell-cell communication network
groupSize <- as.numeric(table(cellchat@idents))


# import cellchat saved below
cellchat <- s3readRDS(bucket = "scrna",
                      object = "Projects/CKD/RH_RH2/data_clean/cellchat_pb90_subset.rds",
                      region = "")

# circle plots showing the total number of interactions between celltypes
temp_file <- tempfile(fileext = ".png")
png(temp_file, 
    width = 3000, height = 3000, res = 600, bg = "transparent")
par(mar = c(0, 0, 1.5, 0), xpd = TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge = F, title.name = "Number of interactions")
dev.off()
put_object(
  file = temp_file,
  object = "Projects/CKD/RH_RH2/Results/Figures/CellChat/cellchat_circle_plots_n_interactions.png",
  bucket = "scrna",
  region = ""
)

# circle plots showing which interactions are strongest
temp_file <- tempfile(fileext = ".png")
png(temp_file, 
    width = 3000, height = 3000, res = 600, bg = "transparent")
par(mar = c(0, 0, 1.5, 0), xpd = TRUE)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge = F, title.name = "Interaction weights/strength")
dev.off()
put_object(
  file = temp_file,
  object = "Projects/CKD/RH_RH2/Results/Figures/CellChat/cellchat_circle_plots_interactions_weights.png",
  bucket = "scrna",
  region = ""
)

# circle plots of strength of interactions by cell type
mat <- cellchat@net$weight
for (i in 1:nrow(mat)) {
  cell_name <- rownames(mat)[i]
  
  # Signals sent from cell type i
  mat_sent <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat_sent[i, ] <- mat[i, ]
  
  temp_file <- tempfile(fileext = ".png")
  png(temp_file, 
      width = 3000, height = 3000, res = 600, bg = "transparent")
  par(mar = c(0, 0, 1.5, 0), xpd = TRUE)
  netVisual_circle(mat_sent, vertex.weight = groupSize, weight.scale = T, 
                   edge.weight.max = max(mat), 
                   title.name = paste0("Signals sent from ", cell_name))
  dev.off()
  put_object(
    file = temp_file,
    object = paste0("Projects/CKD/RH_RH2/Results/Figures/CellChat/cellchat_circle_plots_sent_wt_", cell_name, ".png"),
    bucket = "scrna",
    region = ""
  )
  
  # Signals received by cell type i
  mat_recv <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat_recv[, i] <- mat[, i]
  
  temp_file <- tempfile(fileext = ".png")
  png(temp_file, 
      width = 3000, height = 3000, res = 600, bg = "transparent")
  par(mar = c(0, 0, 1.5, 0), xpd = TRUE)
  netVisual_circle(mat_recv, vertex.weight = groupSize, weight.scale = T, 
                   edge.weight.max = max(mat), 
                   title.name = paste0("Signals received by ", cell_name))
  dev.off()
  put_object(
    file = temp_file,
    object = paste0("Projects/CKD/RH_RH2/Results/Figures/CellChat/cellchat_circle_plots_rcvd_wt_", cell_name, ".png"),
    bucket = "scrna",
    region = ""
  )
}

# Analyze autocrine vs paracrine signaling
# Create matrix to identify autocrine (diagonal) vs paracrine (off-diagonal)
interaction_matrix <- cellchat@net$weight

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
autocrine_df_long <- autocrine_df %>%
  dplyr::mutate(cell_type = gsub("_", "/", cell_type),
                                     n = autocrine + paracrine_sent + paracrine_received) %>%
  arrange(desc(n)) %>%
  dplyr::mutate(cell_type = factor(cell_type, levels = unique(cell_type))) %>%
  dplyr::select(-n) %>%
  melt(id.vars = "cell_type") %>%
  dplyr::mutate(variable = case_when(variable == "autocrine" ~ "Autocrine",
                                     variable == "paracrine_sent" ~ "Paracrine Sent",
                                     variable == "paracrine_received" ~ "Paracrine Rcvd"))
  

p <- ggplot(autocrine_df_long, aes(x = cell_type, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        text = element_text(size = 15),
        legend.position = c(0.8, 0.7)) +
  labs(x = NULL, y = "Weight of Interactions",
       fill = "Interaction type") +
  scale_fill_manual(values = c("Autocrine" = "#e63946", 
                               "Paracrine Sent" = "#a8dadc", 
                               "Paracrine Rcvd" = "#457b9d"))

s3write_using_region(p, FUN = ggsave,
                     object = "Projects/CKD/RH_RH2/Results/Figures/CellChat/bar_interaction_strength.png", 
                     bucket = "scrna",
                     region = "",
                     width = 8, height = 5)

# Identify major signaling sources and targets
# Calculate centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# Visualize dominant senders and receivers
netAnalysis_signalingRole_ggplot(cellchat, pattern = "outgoing",
                                 color.heatmap = "#e76f51",
                                 filepath = "Projects/CKD/RH_RH2/Results/Figures/CellChat/outgoing_signal_patterns.png", 
                                 bucket = "scrna",
                                 height = 18, width = 10)

netAnalysis_signalingRole_ggplot(cellchat, pattern = "incoming",
                                 color.heatmap = "#4f772d",
                                 filepath = "Projects/CKD/RH_RH2/Results/Figures/CellChat/incoming_signal_patterns.png", 
                                 bucket = "scrna",
                                 height = 18, width = 10)


# Identify communication patterns
# Identify global communication patterns
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = 3)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = 3)

# Visualize communication patterns
netAnalysis_river(cellchat, pattern = "outgoing")
netAnalysis_river(cellchat, pattern = "incoming")

# Save CellChat object for future use
s3saveRDS(cellchat,
          bucket = "scrna",
          object = "Projects/CKD/RH_RH2/data_clean/cellchat_pb90_subset.rds",
          region = "",
          multipart = T)
