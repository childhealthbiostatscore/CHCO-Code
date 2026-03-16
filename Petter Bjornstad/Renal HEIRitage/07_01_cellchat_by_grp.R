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

# subset pb90 into own seurat objects (need to run CellChat for each group separately)
pb90_subset_t2d_glpn <- subset(pb90_subset, glp_t2dob == "GLP_N" & group == "Type_2_Diabetes")
pb90_subset_t2d_glpy <- subset(pb90_subset, glp_t2dob == "GLP_Y" & group == "Type_2_Diabetes")
pb90_subset_hc<- subset(pb90_subset, group == "Lean_Control")

# Cell chat
# Set options for CellChat
options(stringsAsFactors = FALSE)

# Define subsets in a named list
subsets <- list(
  t2d_glpn = pb90_subset_t2d_glpn,
  t2d_glpy = pb90_subset_t2d_glpy,
  hc = pb90_subset_hc
)

# Store results
cellchat_list <- list()

for (name in names(subsets)) {
  
  message("Processing: ", name)
  
  obj <- subsets[[name]]
  
  data.input <- GetAssayData(obj, assay = "RNA", layer = "data")
  meta <- obj@meta.data
  
  cellchat <- createCellChat(object = data.input, meta = meta,
                              group.by = "KPMP_celltype_general")
  
  cellchat@DB <- CellChatDB.human
  cellchat <- subsetData(cellchat)
  
  options(future.globals.maxSize = 8000 * 1024^2)
  future::plan("multisession", workers = 4)
  
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  
  cellchat <- computeCommunProb(cellchat,
                                 type = "triMean",
                                 population.size = TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  cellchat_list[[name]] <- cellchat
  
  message("Done: ", name)
}

cellchat_list[["hc"]] <- netAnalysis_computeCentrality(cellchat_list[["hc"]], slot.name = "netP")
cellchat_list[["t2d_glpy"]] <- netAnalysis_computeCentrality(cellchat_list[["t2d_glpy"]], slot.name = "netP")
cellchat_list[["t2d_glpn"]] <- netAnalysis_computeCentrality(cellchat_list[["t2d_glpn"]], slot.name = "netP")

s3saveRDS(cellchat_list,
          object = "Projects/CKD/RH_RH2/data_clean/cellchat_list_t2dglpyn_hc.rds",
          bucket = "scrna", region = "",
          multipart = T)

# merge cell chat objects
cellchat_merged <- mergeCellChat(cellchat_list, 
                                 add.names = names(cellchat_list))
s3saveRDS(cellchat_merged,
          object = "Projects/CKD/RH_RH2/data_clean/cellchat_merged.rds",
          bucket = "scrna", region = "",
          multipart = T)

cellchat_merged_ct <- as.data.frame(sapply(cellchat_merged@net, function(x) sum(x$count)))
colnames(cellchat_merged_ct) <- "count"
cellchat_merged_ct$dataset <- names(cellchat_merged@net)
cellchat_merged_ct$dataset = factor(cellchat_merged_ct$dataset, 
                                    labels = c("hc" = "HC",
                                               "t2d_glpn" = "T2D GLP-",
                                               "t2d_glpy" = "T2D GLP+"))
cellchat_merged_wt <- as.data.frame(sapply(cellchat_merged@net, function(x) sum(x$weight)))
cellchat_merged_wt[, 1] <- round(cellchat_merged_wt[, 1], 3)
colnames(cellchat_merged_wt) <- "weight"
cellchat_merged_wt$dataset <- names(cellchat_merged@net)
cellchat_merged_wt$dataset = factor(cellchat_merged_wt$dataset, 
                                    labels = c("hc" = "HC",
                                               "t2d_glpn" = "T2D GLP-",
                                               "t2d_glpy" = "T2D GLP+"))

# recreate compareInteractions() plot
cellchat_merged_ct_bar <- cellchat_merged_ct %>%
  ggplot(aes(x = dataset, y = count, fill = dataset)) +
  geom_col() + 
  geom_text(aes(label = count), 
            vjust = -0.3, 
            fontface = "bold",
            color = "#343a40") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 15),
        legend.position = "none") +
  labs(x = NULL, y = "N of interactions") +
  scale_fill_manual(values = c("HC" = "#588157",
                               "T2D GLP-" = "#e76f51",
                               "T2D GLP+" = "#f4a261")) 
s3write_using_region(cellchat_merged_ct_bar, 
                     FUN = ggsave,
                     object = "Projects/CKD/RH_RH2/Results/Figures/CellChat/grp_analysis/cellchat_count_bar.png",
                     bucket = "scrna",
                     region = "",
                     width = 7, height = 4)

cellchat_merged_wt_bar <- cellchat_merged_wt %>%
  ggplot(aes(x = dataset, y = weight, fill = dataset)) +
  geom_col() + 
  geom_text(aes(label = weight), 
            vjust = -0.3, 
            fontface = "bold",
            color = "#343a40") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 15),
        legend.position = "none") +
  labs(x = NULL, y = "Interaction strength") +
  scale_fill_manual(values = c("HC" = "#588157",
                               "T2D GLP-" = "#e76f51",
                               "T2D GLP+" = "#f4a261")) 
s3write_using_region(cellchat_merged_wt_bar, 
                     FUN = ggsave,
                     object = "Projects/CKD/RH_RH2/Results/Figures/CellChat/grp_analysis/cellchat_weight_bar.png",
                     bucket = "scrna",
                     region = "",
                     width = 7, height = 4)

# Circle plots
# one circle plot per grouping
for (i in 1:length(cellchat_list)) {
  obj <- cellchat_list[[i]]
  grp_name <- names(cellchat_list)[i]
  
  # --- Overall weight plot ---
  mat_wt <- obj@net$weight
  temp_file <- tempfile(fileext = ".png")
  png(temp_file, 
      width = 3000, height = 3000, res = 600, bg = "transparent")
  par(mar = c(0, 0, 1.5, 0), xpd = TRUE)
  netVisual_circle(mat_wt, vertex.weight = groupSize, weight.scale = T, 
                   label.edge = F, title.name = "Interaction weights/strength")
  dev.off()
  put_object(
    file = temp_file,
    object = paste0("Projects/CKD/RH_RH2/Results/Figures/CellChat/grp_analysis/cellchat_circle_plots_wt_", grp_name, ".png"),
    bucket = "scrna",
    region = ""
  )
  
  # --- Overall count plot ---
  mat_ct <- obj@net$count
  temp_file <- tempfile(fileext = ".png")
  png(temp_file, 
      width = 3000, height = 3000, res = 600, bg = "transparent")
  par(mar = c(0, 0, 1.5, 0), xpd = TRUE)
  netVisual_circle(mat_ct, vertex.weight = groupSize, weight.scale = T, 
                   label.edge = F, title.name = "Number of interactions")
  dev.off()
  put_object(
    file = temp_file,
    object = paste0("Projects/CKD/RH_RH2/Results/Figures/CellChat/grp_analysis/cellchat_circle_plots_ct_", grp_name, ".png"),
    bucket = "scrna",
    region = ""
  )
  # --- Per-cell-type weight plots (sent and received) ---
  for (j in 1:nrow(mat_wt)) {
    cell_name <- rownames(mat_wt)[j]
    
    # Signals sent from cell type j
    mat_sent <- matrix(0, nrow = nrow(mat_wt), ncol = ncol(mat_wt), dimnames = dimnames(mat_wt))
    mat_sent[j, ] <- mat_wt[j, ]
    
    temp_file <- tempfile(fileext = ".png")
    png(temp_file, 
        width = 3000, height = 3000, res = 600, bg = "transparent")
    par(mar = c(0, 0, 1.5, 0), xpd = TRUE)
    netVisual_circle(mat_sent, vertex.weight = groupSize, weight.scale = T, 
                     edge.weight.max = max(mat_wt), 
                     title.name = paste0("Signals sent from ", cell_name))
    dev.off()
    put_object(
      file = temp_file,
      object = paste0("Projects/CKD/RH_RH2/Results/Figures/CellChat/grp_analysis/cellchat_circle_plots_sent_wt_", grp_name, "_", cell_name, ".png"),
      bucket = "scrna",
      region = ""
    )
    
    # Signals received by cell type j
    mat_recv <- matrix(0, nrow = nrow(mat_wt), ncol = ncol(mat_wt), dimnames = dimnames(mat_wt))
    mat_recv[, j] <- mat_wt[, j]
    
    temp_file <- tempfile(fileext = ".png")
    png(temp_file, 
        width = 3000, height = 3000, res = 600, bg = "transparent")
    par(mar = c(0, 0, 1.5, 0), xpd = TRUE)
    netVisual_circle(mat_recv, vertex.weight = groupSize, weight.scale = T, 
                     edge.weight.max = max(mat_wt), 
                     title.name = paste0("Signals received by ", cell_name))
    dev.off()
    put_object(
      file = temp_file,
      object = paste0("Projects/CKD/RH_RH2/Results/Figures/CellChat/grp_analysis/cellchat_circle_plots_rcvd_wt_", grp_name, "_", cell_name, ".png"),
      bucket = "scrna",
      region = ""
    )
  }
}

# Differential interactions circle plot
temp_file <- tempfile(fileext = ".png")
png(temp_file, 
    width = 3000, height = 3000, res = 600, bg = "transparent")
par(mar = c(0, 0, 2, 0), xpd = TRUE)
netVisual_diffInteraction(cellchat_merged, weight.scale = T, comparison = c(1,2), measure = "weight",
                          title.name = "Differential interaction strength\nT2D GLP+ vs. T2D GLP-") # 1: t2d_glpn vs. 2: t2d_glpy
dev.off()
put_object(
  file = temp_file,
  object = "Projects/CKD/RH_RH2/Results/Figures/CellChat/grp_analysis/cellchat_circle_diff_t2dglpyn.png",
  bucket = "scrna",
  region = ""
)

temp_file <- tempfile(fileext = ".png")
png(temp_file, 
    width = 3000, height = 3000, res = 600, bg = "transparent")
par(mar = c(0, 0, 2, 0), xpd = TRUE)
netVisual_diffInteraction(cellchat_merged, weight.scale = T, comparison = c(3,1), measure = "weight",
                          title.name = "Differential interaction strength\nT2D GLP- vs. HC") # 3: hc vs. 1: t2d_glpn
dev.off()
put_object(
  file = temp_file,
  object = "Projects/CKD/RH_RH2/Results/Figures/CellChat/grp_analysis/cellchat_circle_diff_t2dglpn_hc.png",
  bucket = "scrna",
  region = ""
)

# Compare major sources/targets in 2D space
num.link <- sapply(cellchat_list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(cellchat_list)) {
  title = gsub("_GLPY", " GLP+", toupper(names(cellchat_list[i])))
  title = gsub("_GLPN", " GLP-", title)

  gg[[i]] <- netAnalysis_signalingRole_scatter(cellchat_list[[i]], 
                                               title = title, 
                                               weight.MinMax = weight.MinMax,
                                               font.size = 15, font.size.title = 18, 
                                               label.size = 5) +
    geom_abline(slope = 1, linetype = "dashed", color = "grey", alpha = 0.5) +
    scale_colour_manual(
      values = scales::alpha(scPalette(length(levels(cellchat@idents))), 0.8),
      drop = FALSE
    ) +
    theme(plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent", color = NA),
          legend.background = element_rect(fill = "transparent", colour = NA),
          legend.box.background = element_rect(fill = "transparent", colour = NA),
          legend.position = c(0.1,0.6)) 
  p_built <- ggplot_build(gg[[i]])
  max_x <- max(p_built$layout$panel_params[[1]]$x.range)
  max_y <- max(p_built$layout$panel_params[[1]]$y.range)
  gg[[i]] <- gg[[i]] +
    annotate("text", x = max_x, y = 0, label = "Major source", hjust = 1, vjust = -0.5,
             fontface = "italic", color = "grey") +
    annotate("text", x = 0, y = max_y, label = "Major target", hjust = -0.1, vjust = -0.5,
             fontface = "italic", color = "grey") +
    annotate("text", x = max_x, y = max_y, label = "Major source & target", hjust = 1, vjust = -0.5,
             fontface = "italic", color = "grey")
  
  s3write_using_region(gg[[i]],
                       FUN = ggsave,
                       object = paste0("Projects/CKD/RH_RH2/Results/Figures/CellChat/grp_analysis/scatter_", names(cellchat_list[i]), ".png"),
                       bucket = "scrna",
                       region = "",
                       bg = "transparent",
                       width = 7, height = 5)
}

# Identify signaling groups based on their functionality similarity
cellchat_merged <- computeNetSimilarityPairwise(cellchat_merged, type = "functional")
cellchat_merged <- netEmbedding(cellchat_merged, type = "functional")
cellchat_merged <- netClustering(cellchat_merged, type = "functional", do.parallel = F)
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat_merged, type = "functional", label.size = 3.5) +
  theme(plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", colour = NA),
        legend.box.background = element_rect(fill = "transparent", colour = NA))
s3write_using_region(FUN = ggsave,
                     object = "Projects/CKD/RH_RH2/Results/Figures/CellChat/grp_analysis/signaling_grp_functional.png",
                     bucket = "scrna",
                     region = "",
                     bg = "transparent",
                     width = 10, height = 7)

# Identify signaling groups based on structure similarity
cellchat_merged <- computeNetSimilarityPairwise(cellchat_merged, type = "structural")
cellchat_merged <- netEmbedding(cellchat_merged, type = "structural")
cellchat_merged <- netClustering(cellchat_merged, type = "structural", do.parallel = F)
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat_merged, type = "structural", label.size = 3.5)
netVisual_embeddingPairwiseZoomIn(cellchat_merged, type = "structural", nCol = 2)


# Compute and visualize the pathway distnace in the learned joint manifold
rankSimilarity(cellchat_merged, type = "functional", comparison2 = c(1,2))
rankSimilarity(cellchat_merged, type = "functional", comparison2 = c(3,1))

s3saveRDS(cellchat_merged,
          object = "Projects/CKD/RH_RH2/data_clean/cellchat_merged.rds",
          bucket = "scrna", region = "",
          multipart = T)





