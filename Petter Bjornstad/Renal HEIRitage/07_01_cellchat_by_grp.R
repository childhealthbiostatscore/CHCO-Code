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

source(file.path(git_path, "Renal HEIRitage/RH_RH2_IMPROVE_scRNA_functions.R"))
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

# cellchat_list <- s3readRDS(
#           object = "Projects/CKD/RH_RH2/data_clean/cellchat_list_t2dglpyn_hc.rds",
#           bucket = "scrna", region = "")

# merge cell chat objects
cellchat_merged <- mergeCellChat(cellchat_list, 
                                 add.names = names(cellchat_list))

# Identify signaling groups based on their functionality similarity
cellchat_merged <- computeNetSimilarityPairwise(cellchat_merged, type = "functional")
cellchat_merged <- netEmbedding(cellchat_merged, type = "functional")
cellchat_merged <- netClustering(cellchat_merged, type = "functional", do.parallel = F)

# Identify signaling groups based on structure similarity
cellchat_merged <- computeNetSimilarityPairwise(cellchat_merged, type = "structural")
cellchat_merged <- netEmbedding(cellchat_merged, type = "structural")
cellchat_merged <- netClustering(cellchat_merged, type = "structural", do.parallel = F)

s3saveRDS(cellchat_merged,
          object = "Projects/CKD/RH_RH2/data_clean/cellchat_merged.rds",
          bucket = "scrna", region = "",
          multipart = T)

# cellchat_merged <- s3readRDS(
#   object = "Projects/CKD/RH_RH2/data_clean/cellchat_merged.rds",
#   bucket = "scrna", region = "")

cellchat_merged_ct <- as.data.frame(sapply(cellchat_merged@net, function(x) sum(x$count)))
colnames(cellchat_merged_ct) <- "count"
cellchat_merged_ct$dataset <- names(cellchat_merged@net)
cellchat_merged_ct$dataset = factor(cellchat_merged_ct$dataset, 
                                    levels = c("t2d_glpy", "t2d_glpn", "hc"),
                                    labels = c("t2d_glpy" = "T2D GLP+",
                                               "t2d_glpn" = "T2D GLP-",
                                               "hc" = "HC"))
cellchat_merged_wt <- as.data.frame(sapply(cellchat_merged@net, function(x) sum(x$weight)))
cellchat_merged_wt[, 1] <- round(cellchat_merged_wt[, 1], 3)
colnames(cellchat_merged_wt) <- "weight"
cellchat_merged_wt$dataset <- names(cellchat_merged@net)
cellchat_merged_wt$dataset = factor(cellchat_merged_wt$dataset, 
                                    levels = c("t2d_glpy", "t2d_glpn", "hc"),
                                    labels = c("t2d_glpy" = "T2D GLP+",
                                               "t2d_glpn" = "T2D GLP-",
                                               "hc" = "HC"))

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

# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat_merged, type = "structural", label.size = 3.5)
netVisual_embeddingPairwiseZoomIn(cellchat_merged, type = "structural", nCol = 2)

# Compute and visualize the pathway distance in the learned joint manifold
# Top pathways have the largest function distance
rankSimilarity(cellchat_merged, type = "functional",
               comparison2 = c(1,2)) + # T2D GLP+ vs. T2D GLP-
  theme(plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", colour = NA),
        legend.box.background = element_rect(fill = "transparent", colour = NA))
s3write_using_region(FUN = ggsave,
                     object = "Projects/CKD/RH_RH2/Results/Figures/CellChat/grp_analysis/functional_dist_t2dglpyn.png",
                     bucket = "scrna",
                     region = "",
                     bg = "transparent",
                     width = 7, height = 10)

rankSimilarity(cellchat_merged, type = "functional", comparison2 = c(3,1)) + # T2D GLP- vs. HC
  theme(plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", colour = NA),
        legend.box.background = element_rect(fill = "transparent", colour = NA))
s3write_using_region(FUN = ggsave,
                     object = "Projects/CKD/RH_RH2/Results/Figures/CellChat/grp_analysis/functional_dist_t2dglpn_hc.png",
                     bucket = "scrna",
                     region = "",
                     bg = "transparent",
                     width = 7, height = 10)

# Identify altered signaling with distinct interaction strength
rankNet(cellchat_merged, mode = "comparison", measure = "weight",
        comparison = c(1,2),
        sources.use = NULL, targets.use = NULL, stacked = T, do.stat = T)

rankNet(cellchat_merged, mode = "comparison", measure = "weight",
        comparison = c(3,1),
        sources.use = NULL, targets.use = NULL, stacked = T, do.stat = T)








# ============================================================================
# CellChat: Cell-Type Level Pathway Visualizations
# Saves to S3: scrna bucket
# Base path: Projects/CKD/RH_RH2/Results/Figures/CellChat/grp_analysis/
# ============================================================================

# S3 base path
s3_base <- "Projects/CKD/RH_RH2/Results/Figures/CellChat/grp_analysis/"
s3_bucket <- "scrna"
s3_region <- ""

# Pathway lists
treatment_pathways <- c("NRG", "MPZ", "JAM", "BMP", "GDF", "ADGRE5",
                        "SEMA4", "ESAM", "VCAM", "CDH")

disease_pathways <- c("MIF", "SEMA4", "AGRN", "HSPG", "JAM", "PTPRM",
                      "FGF", "PDGF", "FN1", "CD99")

all_key_pathways <- unique(c(treatment_pathways, disease_pathways))


# ============================================================================
# HELPER: Save base R plot to S3
# ============================================================================
save_base_plot_s3 <- function(expr, s3_object,
                              width = 3000, height = 3000, res = 600,
                              mar = c(0,0,1,0)) {
  temp_file <- tempfile(fileext = ".png")
  png(temp_file, width = width, height = height, res = res, bg = "transparent")
  par(mar = mar, xpd = F)
  tryCatch({
    eval(expr)
    dev.off()
    put_object(
      file = temp_file,
      object = paste0(s3_base, s3_object),
      bucket = s3_bucket,
      region = s3_region
    )
    cat("Saved:", s3_object, "\n")
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    cat("Skipped:", s3_object, "-", e$message, "\n")
  })
  unlink(temp_file)
}

# ============================================================================
# HELPER: Save ggplot to S3
# ============================================================================
save_ggplot_s3 <- function(p, s3_object, width = 12, height = 8) {
  tryCatch({
    s3write_using_region(p,
                         FUN = ggsave,
                         object = paste0(s3_base, s3_object),
                         bucket = s3_bucket,
                         region = s3_region,
                         bg = "transparent",
                         width = width, height = height)
    cat("Saved:", s3_object, "\n")
  }, error = function(e) {
    cat("Skipped:", s3_object, "-", e$message, "\n")
  })
}


# ============================================================================
# SECTION 1: CHORD DIAGRAMS — Treatment Effect (GLP-N vs GLP-Y)
# ============================================================================
# Each condition gets its own separate plot

cat("\n========== CHORD: TREATMENT EFFECT ==========\n\n")

for (pathway in treatment_pathways) {
  
  # GLP-N
  save_base_plot_s3(
    quote(netVisual_aggregate(cellchat_glpn, signaling = pathway,
                              layout = "chord",
                              title.space = 0,
                              top = 0)),
    s3_object = paste0("chord_treatment_glpn_", pathway, ".png")
  )
  
  # GLP-Y
  save_base_plot_s3(
    quote(netVisual_aggregate(cellchat_glpy, signaling = pathway,
                              layout = "chord",
                              title.space = 0,
                              top = 0)),
    s3_object = paste0("chord_treatment_glpy_", pathway, ".png")
  )
}


# ============================================================================
# SECTION 2: CHORD DIAGRAMS — Disease Signature (HC vs GLP-N)
# ============================================================================

cat("\n========== CHORD: DISEASE SIGNATURE ==========\n\n")

for (pathway in disease_pathways) {
  
  # HC
  save_base_plot_s3(
    quote(netVisual_aggregate(cellchat_hc, signaling = pathway,
                              layout = "chord",
                              title.name = paste0(pathway, " - HC"))),
    s3_object = paste0("chord_disease_hc_", pathway, ".png")
  )
  
  # T2D GLP-N
  save_base_plot_s3(
    quote(netVisual_aggregate(cellchat_glpn, signaling = pathway,
                              layout = "chord",
                              title.name = paste0(pathway, " - T2D GLP-N"))),
    s3_object = paste0("chord_disease_glpn_", pathway, ".png")
  )
}


# ============================================================================
# SECTION 3: CHORD DIAGRAMS — Recovery (GLP-Y vs HC)
# ============================================================================

cat("\n========== CHORD: RECOVERY ==========\n\n")

for (pathway in all_key_pathways) {
  
  # GLP-Y
  save_base_plot_s3(
    quote(netVisual_aggregate(cellchat_glpy, signaling = pathway,
                              layout = "chord",
                              title.name = paste0(pathway, " - T2D GLP-Y"))),
    s3_object = paste0("chord_recovery_glpy_", pathway, ".png")
  )
  
  # HC
  save_base_plot_s3(
    quote(netVisual_aggregate(cellchat_hc, signaling = pathway,
                              layout = "chord",
                              title.name = paste0(pathway, " - HC"))),
    s3_object = paste0("chord_recovery_hc_", pathway, ".png")
  )
}


# ============================================================================
# SECTION 4: BUBBLE PLOTS — Treatment Effect (GLP-N vs GLP-Y)
# ============================================================================

cat("\n========== BUBBLE: TREATMENT EFFECT ==========\n\n")

for (pathway in treatment_pathways) {
  tryCatch({
    p <- netVisual_bubble(cellchat_merged,
                          sources.use = NULL,
                          targets.use = NULL,
                          signaling = pathway,
                          comparison = c(1, 2),
                          angle.x = 45,
                          title.name = paste0(pathway, " - GLP-N vs GLP-Y"))
    
    save_ggplot_s3(p, paste0("bubble_treatment_", pathway, ".png"))
    
    p_bar_data <- p$data %>%
      filter(source != "Other" & target != "Other") %>%
      group_by(source, target, ligand, receptor) %>%
      summarise(
        prob_glpy = sum(prob[dataset == "t2d_glpy"], na.rm = FALSE),
        prob_glpn = sum(prob[dataset == "t2d_glpn"], na.rm = FALSE),
        .groups = "drop"
      ) %>%
      mutate(
        prob_glpy = tidyr::replace_na(prob_glpy, 0),
        prob_glpn = tidyr::replace_na(prob_glpn, 0),
        diff = prob_glpy - prob_glpn,
        lr_pair = paste0(ligand, " - ", receptor),
        st_pair = paste0(source, " > ", target)
      )
    
    n_sources <- length(unique(p_bar_data$source))
    vline_positions <- seq(1.5, n_sources - 0.5, by = 1)
    
    p_bar <- p_bar_data %>%
      ggplot(aes(x = source, y = diff, fill = lr_pair)) +
      geom_vline(xintercept = vline_positions, linetype = "dashed", 
                 color = "grey60", linewidth = 0.4) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
      geom_col(position = position_dodge(width = 0.9, preserve = "single"), width = 0.8) +
      theme_minimal() +
      theme(panel.grid = element_blank(),
            text = element_text(size = 15),
            panel.border = element_rect(color = "#343a40", fill = NA, 
                                        linetype = "dashed", linewidth = 0.5),
            axis.text.x = element_text(angle = 40, hjust = 1)) +
      labs(fill = "Ligand-Receptor",
           x = "Source",
           y = "Difference in Probability\n(T2D GLP+ - T2D GLP-)") +
      facet_wrap(~ paste0("Target: ",target)) +
      scale_fill_manual(values = c(
        "#2E4F5E",
        "#3AAD9A", 
        "#E8C76A", 
        "#F0A868", 
        "#E8634A",
        "#ae2012",
        "#bb3e03",
        "#ee9b00",
        "#e9d8a6",
        "#94d2bd",
        "#0a9396",
        "#005f73"
      ))
    
    save_ggplot_s3(p_bar, paste0("bubble_bar_treatment_", pathway, ".png"))
  }, error = function(e) {
    cat("Skipped bubble treatment", pathway, "-", e$message, "\n")
  })
}

# Combined
tryCatch({
  p_all <- netVisual_bubble(cellchat_merged,
                            sources.use = NULL,
                            targets.use = NULL,
                            signaling = treatment_pathways,
                            comparison = c(1, 2),
                            angle.x = 45,
                            title.name = "All treatment pathways: GLP-N vs GLP-Y")
  
  save_ggplot_s3(p_all, "bubble_treatment_ALL.png", width = 16, height = 20)
}, error = function(e) {
  cat("Combined treatment bubble error:", e$message, "\n")
})


# ============================================================================
# SECTION 5: BUBBLE PLOTS — Disease Signature (GLP-N vs HC)
# ============================================================================

cat("\n========== BUBBLE: DISEASE SIGNATURE ==========\n\n")

for (pathway in disease_pathways) {
  tryCatch({
    p <- netVisual_bubble(cellchat_merged,
                          sources.use = NULL,
                          targets.use = NULL,
                          signaling = pathway,
                          comparison = c(1, 3),
                          angle.x = 45,
                          title.name = paste0(pathway, " - T2D GLP-N vs HC"))
    
    save_ggplot_s3(p, paste0("bubble_disease_", pathway, ".png"))
  }, error = function(e) {
    cat("Skipped bubble disease", pathway, "-", e$message, "\n")
  })
}

# Combined
tryCatch({
  p_all_disease <- netVisual_bubble(cellchat_merged,
                                    sources.use = NULL,
                                    targets.use = NULL,
                                    signaling = disease_pathways,
                                    comparison = c(1, 3),
                                    angle.x = 45,
                                    title.name = "All disease pathways: T2D GLP-N vs HC")
  
  save_ggplot_s3(p_all_disease, "bubble_disease_ALL.png", width = 16, height = 20)
}, error = function(e) {
  cat("Combined disease bubble error:", e$message, "\n")
})


# ============================================================================
# SECTION 6: BUBBLE PLOTS — Recovery (GLP-Y vs HC)
# ============================================================================

cat("\n========== BUBBLE: RECOVERY ==========\n\n")

for (pathway in all_key_pathways) {
  tryCatch({
    p <- netVisual_bubble(cellchat_merged,
                          sources.use = NULL,
                          targets.use = NULL,
                          signaling = pathway,
                          comparison = c(2, 3),
                          angle.x = 45,
                          title.name = paste0(pathway, " - GLP-Y vs HC"))
    
    save_ggplot_s3(p, paste0("bubble_recovery_", pathway, ".png"))
  }, error = function(e) {
    cat("Skipped bubble recovery", pathway, "-", e$message, "\n")
  })
}


# ============================================================================
# SECTION 7: BUBBLE PLOTS — Three-Way Comparison
# ============================================================================

cat("\n========== BUBBLE: THREE-WAY ==========\n\n")

for (pathway in all_key_pathways) {
  tryCatch({
    p <- netVisual_bubble(cellchat_merged,
                          sources.use = NULL,
                          targets.use = NULL,
                          signaling = pathway,
                          comparison = c(1, 2, 3),
                          angle.x = 45,
                          title.name = paste0(pathway, " - All conditions"))
    
    save_ggplot_s3(p, paste0("bubble_3way_", pathway, ".png"), width = 16, height = 10)
  }, error = function(e) {
    cat("Skipped 3-way bubble", pathway, "-", e$message, "\n")
  })
}


# ============================================================================
# SECTION 8: DIFFERENTIAL HEATMAPS
# ============================================================================

cat("\n========== DIFFERENTIAL HEATMAPS ==========\n\n")

# --- Treatment (GLP-N vs GLP-Y) ---
for (pathway in treatment_pathways) {
  temp_file <- tempfile(fileext = ".png")
  png(temp_file, width = 3000, height = 2400, res = 600, bg = "transparent")
  tryCatch({
    netVisual_heatmap(cellchat_merged,
                      signaling = pathway,
                      comparison = c(1, 2),
                      title.name = paste0(pathway, ": GLP-N vs GLP-Y"),
                      color.heatmap = "Spectral")
    dev.off()
    put_object(
      file = temp_file,
      object = paste0(s3_base, "heatmap_treatment_", pathway, ".png"),
      bucket = s3_bucket,
      region = s3_region
    )
    cat("Heatmap saved: treatment", pathway, "\n")
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    cat("Skipped heatmap treatment", pathway, "-", e$message, "\n")
  })
  unlink(temp_file)
}

# Overall treatment
temp_file <- tempfile(fileext = ".png")
png(temp_file, width = 3000, height = 2400, res = 600, bg = "transparent")
tryCatch({
  netVisual_heatmap(cellchat_merged,
                    comparison = c(1, 2),
                    title.name = "Overall: GLP-N vs GLP-Y",
                    color.heatmap = "Spectral")
  dev.off()
  put_object(
    file = temp_file,
    object = paste0(s3_base, "heatmap_treatment_OVERALL.png"),
    bucket = s3_bucket,
    region = s3_region
  )
  cat("Overall treatment heatmap saved\n")
}, error = function(e) {
  try(dev.off(), silent = TRUE)
  cat("Overall treatment heatmap error:", e$message, "\n")
})
unlink(temp_file)

# --- Disease (GLP-N vs HC) ---
for (pathway in disease_pathways) {
  temp_file <- tempfile(fileext = ".png")
  png(temp_file, width = 3000, height = 2400, res = 600, bg = "transparent")
  tryCatch({
    netVisual_heatmap(cellchat_merged,
                      signaling = pathway,
                      comparison = c(1, 3),
                      title.name = paste0(pathway, ": T2D GLP-N vs HC"),
                      color.heatmap = "Spectral")
    dev.off()
    put_object(
      file = temp_file,
      object = paste0(s3_base, "heatmap_disease_", pathway, ".png"),
      bucket = s3_bucket,
      region = s3_region
    )
    cat("Heatmap saved: disease", pathway, "\n")
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    cat("Skipped heatmap disease", pathway, "-", e$message, "\n")
  })
  unlink(temp_file)
}

# Overall disease
temp_file <- tempfile(fileext = ".png")
png(temp_file, width = 3000, height = 2400, res = 600, bg = "transparent")
tryCatch({
  netVisual_heatmap(cellchat_merged,
                    comparison = c(1, 3),
                    title.name = "Overall: T2D GLP-N vs HC",
                    color.heatmap = "Spectral")
  dev.off()
  put_object(
    file = temp_file,
    object = paste0(s3_base, "heatmap_disease_OVERALL.png"),
    bucket = s3_bucket,
    region = s3_region
  )
  cat("Overall disease heatmap saved\n")
}, error = function(e) {
  try(dev.off(), silent = TRUE)
  cat("Overall disease heatmap error:", e$message, "\n")
})
unlink(temp_file)


# ============================================================================
# SECTION 9: L-R PAIR CONTRIBUTIONS
# ============================================================================

cat("\n========== L-R PAIR CONTRIBUTIONS ==========\n\n")

for (pathway in all_key_pathways) {
  
  # HC
  tryCatch({
    p <- netAnalysis_contribution(cellchat_hc, signaling = pathway,
                                  title = paste0(pathway, " L-R pairs (HC)"))
    save_ggplot_s3(p, paste0("LR_hc_", pathway, ".png"), width = 8, height = 5)
  }, error = function(e) {
    cat("  No", pathway, "in HC\n")
  })
  
  # T2D GLP-N
  tryCatch({
    p <- netAnalysis_contribution(cellchat_glpn, signaling = pathway,
                                  title = paste0(pathway, " L-R pairs (T2D GLP-N)"))
    save_ggplot_s3(p, paste0("LR_glpn_", pathway, ".png"), width = 8, height = 5)
  }, error = function(e) {
    cat("  No", pathway, "in GLP-N\n")
  })
  
  # T2D GLP-Y
  tryCatch({
    p <- netAnalysis_contribution(cellchat_glpy, signaling = pathway,
                                  title = paste0(pathway, " L-R pairs (T2D GLP-Y)"))
    save_ggplot_s3(p, paste0("LR_glpy_", pathway, ".png"), width = 8, height = 5)
  }, error = function(e) {
    cat("  No", pathway, "in GLP-Y\n")
  })
}


# ============================================================================
# SECTION 10: SIGNALING ROLE SCATTER PLOTS
# ============================================================================

cat("\n========== SIGNALING ROLE SCATTER PLOTS ==========\n\n")

# These are ggplot-based in newer CellChat versions, base R in older ones.
# Try ggplot first, fall back to base R.

# GLP-N vs GLP-Y
tryCatch({
  p <- netAnalysis_signalingRole_scatter(cellchat_merged,
                                         title = "Signaling roles: GLP-N vs GLP-Y",
                                         comparison = c(1, 2))
  save_ggplot_s3(p, "signaling_roles_glpn_vs_glpy.png")
}, error = function(e) {
  # Fall back to base R saving
  save_base_plot_s3(
    quote(netAnalysis_signalingRole_scatter(cellchat_merged,
                                            title = "Signaling roles: GLP-N vs GLP-Y",
                                            comparison = c(1, 2))),
    s3_object = "signaling_roles_glpn_vs_glpy.png",
    width = 3600, height = 3000
  )
})

# GLP-N vs HC
tryCatch({
  p <- netAnalysis_signalingRole_scatter(cellchat_merged,
                                         title = "Signaling roles: T2D GLP-N vs HC",
                                         comparison = c(1, 3))
  save_ggplot_s3(p, "signaling_roles_glpn_vs_hc.png")
}, error = function(e) {
  save_base_plot_s3(
    quote(netAnalysis_signalingRole_scatter(cellchat_merged,
                                            title = "Signaling roles: T2D GLP-N vs HC",
                                            comparison = c(1, 3))),
    s3_object = "signaling_roles_glpn_vs_hc.png",
    width = 3600, height = 3000
  )
})

# GLP-Y vs HC
tryCatch({
  p <- netAnalysis_signalingRole_scatter(cellchat_merged,
                                         title = "Signaling roles: GLP-Y vs HC",
                                         comparison = c(2, 3))
  save_ggplot_s3(p, "signaling_roles_glpy_vs_hc.png")
}, error = function(e) {
  save_base_plot_s3(
    quote(netAnalysis_signalingRole_scatter(cellchat_merged,
                                            title = "Signaling roles: GLP-Y vs HC",
                                            comparison = c(2, 3))),
    s3_object = "signaling_roles_glpy_vs_hc.png",
    width = 3600, height = 3000
  )
})


cat("\n========== ALL DONE ==========\n")