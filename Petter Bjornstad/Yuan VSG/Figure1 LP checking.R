#----------------------·•●  🧽  ●•·-------------------------
#                        Figure 1
#----------------------·•●──────●•·-------------------------
# Author: Long Yuan
# Affiliation: Johns Hopkins | University of Washington
# Email: lyuan13@jhmi.edu
#-----------------------------------------------------------

library(readr)
library(dplyr)
library(gridExtra)
library(Seurat)
library(lme4)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(stats) 
#
library(lmerTest)
library(broom.mixed)
library(emmeans)
#
library(UpSetR)
library(tidyr)
#
library(ComplexHeatmap)
library(circlize)
library(scales)
library(ggridges)
library(tidytext)
library(msigdbr)
library(fgsea)
library(stringr)
library(purrr)
library(rlang)
#
library(KEGGREST)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(Matrix)
library(edgeR)
library(limma)
library(ggrepel)
#
scrna <- readRDS("/Users/pylell/Library/CloudStorage/OneDrive-UW/Bjornstad/Biostatistics Core Shared Drive/scRNA/data_raw/PB_90_RPCAFix_Metadata_LR_RM_kpmpV1labelled.rds")
scrna <- subset(scrna, subset = percent.mt < 50 & nFeature_RNA < 5000 & nFeature_RNA > 500) 

#
clin <- read_csv("/Users/pylell/Library/CloudStorage/OneDrive-UW/Bjornstad/Biostatistics Core Shared Drive/scRNA/data_clean/pb90_rhrh2improve_clinical_subset.csv", show_col_types = FALSE) %>%
  dplyr::select(record_id, visit,acr_u, gfr_bsa_plasma, gfr_raw_plasma, bmi, eGFR_CKiD_U25_avg, eGFR_fas_cr_cysc) %>%
  dplyr::mutate(visit = dplyr::recode(visit,"baseline" = "pre","12_months_post_surgery" = "post"))

#############
## IMPROVE ##
#############
improve <- subset(
  scrna, 
  subset = cohort == "IMPROVE" | 
    record_id %in% c("RH-59-T", "RH-60-T", "RH-65-T", "RH-66-T"))
improve@meta.data[improve@meta.data$record_id == "RH-65-T", ]$pre_post <- 'pre'
improve@meta.data <- improve@meta.data %>%
  tibble::rownames_to_column("cell") %>%   
  left_join(clin, by = c("record_id" = "record_id", "pre_post" = "visit")) %>%
  tibble::column_to_rownames("cell") 
improve@meta.data[improve@meta.data$record_id == "RH-59-T", ]$record_id <- 'IT_07'
improve@meta.data[improve@meta.data$record_id == "RH-60-T", ]$record_id <- 'IT_08'
improve@meta.data[improve@meta.data$record_id == "RH-65-T", ]$record_id <- 'IT_09'
improve@meta.data[improve@meta.data$record_id == "RH-66-T", ]$record_id <- 'IT_10'

########
## RH ##
########
rh <- subset(
  scrna,
  subset = cohort %in% c('RENAL HEIR', 'RENAL HEIRITAGE') &
    !scrna$record_id %in% c("RH-59-T", "RH-60-T", "RH-65-T", "RH-66-T", "RH2-07-O")
)
rh$pre_post <- ifelse(rh$cohort == "RENAL HEIR", "pre", "post")
rh@meta.data <- rh@meta.data %>%
  tibble::rownames_to_column("cell") %>%   
  left_join(clin, by = c("record_id" = "record_id")) %>%
  tibble::column_to_rownames("cell") 
#
rh@meta.data[rh@meta.data$record_id == "RH2-14-T", ]$record_id <- 'RH-23-T'
rh@meta.data[rh@meta.data$record_id == "RH2-19-T", ]$record_id <- 'RH-67-T'

#
improve@meta.data$treatment <- "VSG"
rh@meta.data$treatment <- "Standard"

healthy <- subset(scrna, 
                  subset = record_id %in% c("CRC-03","CRC-14","CRC-02","CRC-13",
                                            "CRC-10","CRC-39","CRC-40","CRC-11",
                                            "CRC-46","CRC-58","CRC-56","CRC-54"))
healthy@meta.data$treatment <- "Healthy"
healthy@meta.data$pre_post <- "healthy"
#
combined <- merge(improve, rh)
combined <- merge(combined, healthy)
#

# total ppts by group
a <- combined@meta.data[,c("record_id", "treatment", "group", "pre_post")]
a <- unique(a)
table(a$treatment)
a_vsg <- a[a$treatment == "VSG",]
a_smt <- a[a$treatment == "Standard",]
length(unique(a_vsg$record_id))
length(unique(a_smt$record_id))
length(unique(a$record_id))

combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
#
combined@meta.data <- combined@meta.data %>%
  mutate(lineage = case_when(
    KPMP_celltype %in% c("PT-S1/S2", "PT-S3", "aPT") ~ "PT",
    KPMP_celltype %in% c("aTAL", "C-TAL-1", "C-TAL-2", "dTAL") ~ "TAL",
    KPMP_celltype %in% c("CCD-PC", "CNT-PC", "dCCD-PC", "M-PC", "tPC-IC") ~ "PC",
    KPMP_celltype %in% c("IC-A", "IC-B", "aIC") ~ "IC",
    KPMP_celltype %in% c("EC-AEA", "EC-AVR", "EC-GC", "EC-PTC", "EC-LYM", "EC/VSMC") ~ "EC",
    KPMP_celltype %in% c("CD4+ T", "CD8+ T", "NK", "B", "cycT") ~ "Lymphoid",
    KPMP_celltype == "POD" ~ "POD",
    KPMP_celltype %in% c("MAC", "MON", "cDC", "pDC") ~ "Myeloid",
    KPMP_celltype %in% c("FIB", "VSMC/P") ~ "FIB/VSMC/P",
    KPMP_celltype == "SchwannCells" ~ "Schwann",
    TRUE ~ "Other"
  ))

# HVGs per layer
hvg1 <- VariableFeatures(combined, assay = "RNA", layer = "counts.1.1")
hvg2 <- VariableFeatures(combined, assay = "RNA", layer = "counts.2.1")
hvg3 <- VariableFeatures(combined, assay = "RNA", layer = "counts.2")

# union HVGs
hvg_union <- union(union(hvg1, hvg2), hvg3)

# remove mitochondrial & ribosomal genes
hvg_union <- hvg_union[!grepl("^MT-|^RPL|^RPS", hvg_union)]

data1 <- GetAssayData(improve, slot="data", layer="data")
data2 <- GetAssayData(rh, slot="data", layer="data")
data3 <- GetAssayData(healthy, slot="data", layer="data")
expr_matrix <- cbind(data1, data2, data3)
expr_matrix <- expr_matrix[hvg_union, ]

meta_improve <- improve@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  select(cell, record_id, pre_post, treatment)

meta_rh <- rh@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  select(cell, record_id, pre_post, treatment)

meta_healthy <- healthy@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  select(cell, record_id, pre_post, treatment)

meta <- bind_rows(meta_improve, meta_rh, meta_healthy)

meta$group_id <- paste(meta$record_id, meta$pre_post, meta$treatment, sep=".")

pb_list <- split(meta$cell, meta$group_id)

# compute pseudobulk means: genes × samples
pb_mat <- do.call(cbind, lapply(pb_list, function(cells) {
  Matrix::rowMeans(expr_matrix[, cells, drop=FALSE])
}))

pb_mat <- as.matrix(pb_mat)
colnames(pb_mat) <- names(pb_list)

pb_meta <- data.frame(sample = colnames(pb_mat)) %>%
  mutate(
    record_id = sapply(strsplit(sample, "\\."), `[`, 1),
    pre_post  = sapply(strsplit(sample, "\\."), `[`, 2),
    treatment = sapply(strsplit(sample, "\\."), `[`, 3)
  )

# samples × genes for PCA
pb_scaled <- t(scale(t(pb_mat)))  

pca <- prcomp(t(pb_scaled), center=TRUE, scale.=FALSE)

scores <- as.data.frame(pca$x[, 1:2])
scores$sample <- pb_meta$sample
scores$record_id <- pb_meta$record_id
pb_meta$treatment <- dplyr::recode(pb_meta$treatment, 
                                   "Standard" = "SMT")
scores$treatment <- pb_meta$treatment
scores$pre_post  <- pb_meta$pre_post
scores$pre_post <- dplyr::recode(scores$pre_post,
                                 "pre" = "Pre",
                                 "post" = "Post",
                                 "healthy" = "Healthy")
scores$pre_post <- factor(scores$pre_post, 
                          levels = c("Pre", "Post", "Healthy"))

paired_scores <- scores %>%
  filter(treatment %in% c("VSG", "SMT"),
         pre_post %in% c("Pre", "Post")) %>%
  group_by(record_id, treatment) %>%
  filter(all(c("Pre", "Post") %in% pre_post)) %>%
  ungroup()

scores_plot <- bind_rows(
  paired_scores,
  scores %>% filter(treatment == "Healthy")
)

arrow_df <- paired_scores %>%
  select(record_id, treatment, pre_post, PC1, PC2) %>%
  pivot_wider(names_from = pre_post, values_from = c(PC1, PC2)) %>%
  drop_na(PC1_Pre, PC1_Post, PC2_Pre, PC2_Post)

# Treatment colors
treatment_colors <- c(
  "Healthy" = "#7da0c8",
  "VSG"     = "#3d9f3c",
  "SMT"     = "#555d50"
)

# Timepoint shapes: healthy square; pre triangle; post circle
timepoint_shapes <- c(
  "Healthy" = 18,   # solid square
  "Pre"     = 17,   # triangle
  "Post"    = 19    # circle
)

# ----- PLOT -----
p <- ggplot(scores_plot, aes(x = PC1, y = PC2)) +
  
  # Healthy ellipse 
  stat_ellipse(
    data = subset(scores_plot, treatment == "Healthy"),
    geom = "polygon",
    aes(fill = "Healthy"),
    type = "norm", linetype = "solid",
    alpha = 0.15, color = NA
  ) +
  # VSG ellipse
  stat_ellipse(
    data = subset(scores_plot, treatment == "VSG"),
    geom = "polygon",
    aes(fill = "VSG"),
    type = "norm", linetype = "solid",
    alpha = 0.10, color = NA
  ) +
  # SMT ellipse
  stat_ellipse(
    data = subset(scores_plot, treatment == "SMT"),
    geom = "polygon",
    aes(fill = "SMT"),
    type = "norm", linetype = "solid",
    alpha = 0.10, color = NA
  ) +
  
  # ellipses use same fill colors
  scale_fill_manual(values = treatment_colors, guide = "none") +
  
  # Points
  geom_point(aes(color = treatment, shape = pre_post), 
             size = 3., alpha = 0.5) +
  
  # Arrows pre → post
  geom_segment(
    data = arrow_df,
    aes(
      x = PC1_Pre, y = PC2_Pre,
      xend = PC1_Post, yend = PC2_Post,
      color = treatment
    ),
    arrow = arrow(length = unit(0.25, "cm")),
    linewidth = 0.8, alpha = 0.7,
    show.legend = FALSE
  ) +
  
  # Labels 
  ggrepel::geom_text_repel(
    data = subset(scores_plot, treatment %in% c("VSG","SMT")),
    aes(label = record_id, color = treatment),
    size = 3.5,
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = Inf,
    min.segment.length = 0,
    show.legend = FALSE
  ) +
  
  # manual scales
  scale_color_manual(values = treatment_colors, name = "Treatment") +
  scale_shape_manual(values = timepoint_shapes, name = "Timepoint") +
  
  theme_bw() +
  theme(
    legend.position = "right",
    #legend.title = element_blank(),
    panel.grid = element_blank()
  ) +
  
  labs(
    title = "Patient-Level PCA Trajectories (Pre → Post)",
    x = "PC1",
    y = "PC2",
    color = "Treatment",
    shape = "Timepoint"
  )

p

ggsave("patient_pca_trajectory.pdf",
       p, device = cairo_pdf,
       width = 7, height = 6, units = "in")
##
scores_plot <- scores %>%
  mutate(
    group3 = case_when(
      treatment == "Healthy" ~ "HC",
      treatment == "VSG" & pre_post == "Post" ~ "Post-VSG",
      TRUE ~ "T2D"
    )
  )
scores_plot$group3 <- factor(scores_plot$group3, 
                             levels=c("HC",
                                      "T2D",
                                      "Post-VSG"))


ellipse_df_for_group <- function(df, group_label, level = 0.95) {
  library(ellipse)
  
  cov_mat <- cov(df[, c("PC1", "PC2")])
  center  <- colMeans(df[, c("PC1", "PC2")])
  
  el <- as.data.frame(ellipse(cov_mat, centre = center, level = level))
  colnames(el) <- c("x", "y")
  el$group3 <- group_label
  
  return(el)
}

hc_el      <- ellipse_df_for_group(subset(scores_plot, group3 == "HC"), "HC")
t2d_el     <- ellipse_df_for_group(subset(scores_plot, group3 == "T2D"), "T2D")
postvsg_el <- ellipse_df_for_group(subset(scores_plot, group3 == "Post-VSG"), "Post-VSG")

ellipse3_df <- dplyr::bind_rows(hc_el, t2d_el, postvsg_el)


hc_el$group3      <- "HC"
t2d_el$group3     <- "T2D"
postvsg_el$group3 <- "Post-VSG"

ellipse3_df <- dplyr::bind_rows(hc_el, t2d_el, postvsg_el)

group3_colors <- c(
  "HC"       = "#002147",
  "T2D"      = "#b87333",
  "Post-VSG" = "#c9a0dc"
)

p3 <- ggplot(scores_plot, aes(PC1, PC2)) +
  
  # filled ellipses
  geom_polygon(data=ellipse3_df,
               aes(x, y, fill=group3),
               alpha=0.15, color=NA) +
  
  # points
  geom_point(aes(color=group3), size=2, alpha=0.7) +
  
  scale_fill_manual(values=group3_colors, name = "Ellipse") +
  scale_color_manual(values=group3_colors, name = "Status") +
  
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    #legend.title = element_blank()
  ) +
  
  labs(
    title = "Sensitivity PCA: Healthy vs T2D vs Post-VSG",
    x = "PC1",
    y = "PC2"
  )

p3

ggsave("PCA_sensitivity_HC_T2D_PostVSG.pdf",
       p3, device=cairo_pdf,
       width=7, height=4.5, units="in")

##
meta_lineage <- combined@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  filter(cell %in% colnames(expr_matrix)) %>% 
  select(cell, record_id, pre_post, treatment, lineage)

# Align order to expr_matrix columns
meta_lineage <- meta_lineage[match(colnames(expr_matrix), meta_lineage$cell), ]
stopifnot(all(meta_lineage$cell == colnames(expr_matrix)))

meta_lineage$group_id <- paste(meta_lineage$record_id,
                               meta_lineage$pre_post,
                               meta_lineage$treatment,
                               meta_lineage$lineage,
                               sep = ".")

pb_lineage_list <- split(meta_lineage$cell, meta_lineage$group_id)

# genes × lineage-pseudobulk-samples
pb_lineage_mat <- do.call(cbind, lapply(pb_lineage_list, function(cells) {
  Matrix::rowMeans(expr_matrix[, cells, drop = FALSE])
}))
pb_lineage_mat <- as.matrix(pb_lineage_mat)
colnames(pb_lineage_mat) <- names(pb_lineage_list)

# metadata per lineage-pseudobulk sample
pb_lineage_meta <- data.frame(sample = colnames(pb_lineage_mat)) %>%
  mutate(
    record_id = sapply(strsplit(sample, "\\."), `[`, 1),
    pre_post  = sapply(strsplit(sample, "\\."), `[`, 2),
    treatment = sapply(strsplit(sample, "\\."), `[`, 3),
    lineage   = sapply(strsplit(sample, "\\."), `[`, 4)
  )

compute_mahalanobis_pcs <- function(expr_samples_by_genes, healthy_mask, n_pcs = 10) {
  
  # remove zero-variance genes
  vars <- apply(expr_samples_by_genes, 2, var)
  keep_genes <- vars > 0
  X <- expr_samples_by_genes[, keep_genes, drop = FALSE]
  
  X_scaled <- scale(X)
  
  # PCA
  pca <- prcomp(X_scaled, center = TRUE, scale. = FALSE)
  k <- min(n_pcs, ncol(pca$x))
  pcs <- pca$x[, 1:k, drop = FALSE]
  
  healthy_pcs <- pcs[healthy_mask, , drop = FALSE]
  
  # covariance + strong shrinkage
  covmat <- cov(healthy_pcs)
  alpha <- 0.3  
  
  covmat_reg <- covmat * (1 - alpha) + alpha * diag(diag(covmat))
  
  d <- mahalanobis(pcs, center = colMeans(healthy_pcs), cov = covmat_reg)
  return(d)
}

lineages <- setdiff(unique(pb_lineage_meta$lineage), 'Other')

dist_list <- lapply(lineages, function(lin) {
  idx <- pb_lineage_meta$lineage == lin
  mat_sub  <- t(pb_lineage_mat[, idx, drop = FALSE])       # samples × genes
  meta_sub <- pb_lineage_meta[idx, ]
  
  healthy_mask <- meta_sub$treatment == "Healthy"
  
  # skip lineages without healthy controls
  if (sum(healthy_mask) < 3) {
    warning(paste("Not enough healthy samples for lineage", lin, "- skipping"))
    return(NULL)
  }
  
  d <- compute_mahalanobis_pcs(mat_sub, healthy_mask, n_pcs = 30)
  
  out <- meta_sub
  out$distance <- d
  out
})

dist_df <- do.call(rbind, dist_list)

paired_df <- dist_df %>%
  mutate(treatment = dplyr::recode(treatment, "Standard" = "SMT")) %>%
  filter(treatment %in% c("VSG","SMT"),
         pre_post %in% c("pre","post")) %>%
  
  mutate(distance_sqrt = sqrt(distance)) %>%
  
  group_by(record_id, treatment, lineage, pre_post) %>% 
  summarise(distance_sqrt = mean(distance_sqrt), .groups = "drop") %>%
  
  pivot_wider(names_from = pre_post, values_from = distance_sqrt) %>%
  filter(!is.na(pre) & !is.na(post)) %>%
  
  mutate(
    delta = post - pre,
    patient_label = paste0(record_id, " (", treatment, ")")
  )

heat_df <- paired_df %>%
  select(patient_label, lineage, delta) %>%
  tidyr::pivot_wider(
    names_from = lineage,
    values_from = delta
  ) %>%
  as.data.frame() 

# remove patient_label column and set rownames
rownames(heat_df) <- heat_df$patient_label
heat_df$patient_label <- NULL

# remove lineages that contain any NA
heat_df <- heat_df[, colSums(is.na(heat_df)) == 0, drop = FALSE]

# final numeric matrix
heat_matrix <- as.matrix(heat_df)

col_fun <- colorRamp2(
  breaks = c(min(heat_matrix), 0, max(heat_matrix)),
  colors = c("#0c567d", "white", "#edb79c")
)

lgd <- Legend(
  title = "Δ Mahalanobis (Post - Pre)",
  col_fun = col_fun,
  at = c(min(heat_matrix), 0, max(heat_matrix)),
  labels = c("Closer to Healthy  ←", "0", "→ Further from Healthy"),
  direction = "horizontal",
  legend_width = unit(4.5, "cm"),
  title_position = "topcenter"
)

ht <- Heatmap(
  heat_matrix,
  name = "ΔD",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_names_rot = 45,
  show_heatmap_legend = FALSE,
  cell_fun = function(j,i,x,y,width,height,fill) {
    grid.text(sprintf("%.1f", heat_matrix[i,j]), x, y, gp=gpar(fontsize=9))
  }
)

cairo_pdf("lineage_delta_mahalanobis_pc30.pdf", width = 5, height = 5)
draw(ht, heatmap_legend_list = list(lgd), heatmap_legend_side = "bottom")
dev.off()

#
combined <- merge(improve, rh)
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 5000)
#
hvg_layer1 <- VariableFeatures(combined, assay = "RNA", layer = "counts.1")
hvg_layer2 <- VariableFeatures(combined, assay = "RNA", layer = "counts.2")
# Take the union which yields 6327 genes
top5000_union <- union(hvg_layer1, hvg_layer2)
#
combined@meta.data <- combined@meta.data %>%
  mutate(lineage = case_when(
    KPMP_celltype %in% c("PT-S1/S2", "PT-S3", "aPT") ~ "PT",
    KPMP_celltype %in% c("aTAL", "C-TAL-1", "C-TAL-2", "dTAL") ~ "TAL",
    KPMP_celltype %in% c("CCD-PC", "CNT-PC", "dCCD-PC", "M-PC", "tPC-IC") ~ "PC",
    KPMP_celltype %in% c("IC-A", "IC-B", "aIC") ~ "IC",
    KPMP_celltype %in% c("EC-AEA", "EC-AVR", "EC-GC", "EC-PTC", "EC-LYM", "EC/VSMC") ~ "EC",
    KPMP_celltype %in% c("CD4+ T", "CD8+ T", "NK", "B", "cycT") ~ "Lymphoid",
    KPMP_celltype %in% c("POD") ~ "POD",
    KPMP_celltype %in% c("MAC", "MON",  "cDC", "pDC") ~ "Myeloid",
    KPMP_celltype %in% c("FIB", "VSMC/P") ~ "FIB/VSMC/P",
    KPMP_celltype %in% c("SchwannCells") ~ "Schwann",
    TRUE ~ "Other"
  ))
#
meta <- combined@meta.data
meta$cell <- rownames(meta)
meta$pre_post <- factor(meta$pre_post, levels = c("pre", "post"))
meta$treatment <- factor(meta$treatment, levels = c("Standard", "VSG"))
meta$record_id <- as.factor(meta$record_id)
#
data1 <- GetAssayData(improve, slot = "data", layer = "data") #dim 31332 33170
data2 <- GetAssayData(rh, slot = "data", layer = "data") #dim 31332 72727
expr_matrix <- cbind(data1, data2) #31332 109457
expr_matrix <- expr_matrix[top5000_union, ]

#################
## Mixed Model ##
#################
# 
results <- list()
valid_lineages <- c("PT", "TAL", "PC", "EC", "IC", 
                    "POD", "Lymphoid", "Myeloid", 
                    "FIB/VSMC/P", "Schwann")

# loop through 
for (lineage in valid_lineages) {
  message("\n----------------------------")
  message("Processing lineage: ", lineage)
  message("----------------------------")
  
  meta_sub <- meta %>% filter(lineage == !!lineage)
  expr_sub <- expr_matrix[, meta_sub$cell, drop = FALSE]
  
  # loop through all genes
  for (gene in top5000_union) {
    
    if (!gene %in% rownames(expr_sub)) next
    sub_expr_all <- expr_sub[gene, , drop = FALSE]
    
    # skip if gene not expressed
    if (all(sub_expr_all == 0)) next
    
    # create DF for model
    df <- meta_sub %>% mutate(expr = as.numeric(sub_expr_all))
    
    # skip small groups
    if (nrow(df) < 50) next
    
    try({
      model <- lmer(expr ~ pre_post * treatment + (1 | record_id), data = df)
      
      tidy_mod <- tidy(model, effects = "fixed", conf.int = TRUE)
      interaction <- tidy_mod %>% filter(term == "pre_postpost:treatmentVSG")
      
      suppressMessages({
        em <- emmeans(model, ~ pre_post | treatment)
      })
      contrasts <- contrast(em, method = "pairwise")
      contrast_df <- as.data.frame(contrasts)
      
      # contrast is "pre - post"
      std_contrast <- contrast_df %>%
        filter(treatment == "Standard" & contrast == "pre - post")
      vsg_contrast <- contrast_df %>%
        filter(treatment == "VSG" & contrast == "pre - post")
      
      results[[length(results) + 1]] <- tibble(
        gene = gene,
        lineage = lineage,
        interaction_estimate = interaction$estimate,
        interaction_p = interaction$p.value,
        interaction_conf.low = interaction$conf.low,
        interaction_conf.high = interaction$conf.high,
        standard_estimate = std_contrast$estimate,
        standard_p = std_contrast$p.value,
        standard_conf.low = std_contrast$lower.CL,
        standard_conf.high = std_contrast$upper.CL,
        vsg_estimate = vsg_contrast$estimate,
        vsg_p = vsg_contrast$p.value,
        vsg_conf.low = vsg_contrast$lower.CL,
        vsg_conf.high = vsg_contrast$upper.CL,
        converged = TRUE,
        is_singular = isSingular(model, tol = 1e-4)
      )
    }, silent = TRUE)
  }
}

# combine
final_df <- bind_rows(results)
message("Total results stored: ", nrow(final_df))

# FDR correction 
if (nrow(final_df) > 0) {
  final_df <- final_df %>%
    group_by(lineage) %>%
    mutate(
      interaction_FDR = p.adjust(interaction_p, method = "BH"),
      standard_FDR = p.adjust(standard_p, method = "BH"),
      vsg_FDR = p.adjust(vsg_p, method = "BH")
    ) %>%
    ungroup()
  
  write.csv(final_df, "mixed_model_all_lineages_top5k.csv", row.names = FALSE)
}
#
# Thresholds
EST_CUTOFF <- 0.5
FDR_CUTOFF <- 0.05

# Build lineage–gene sets per contrast
get_sig_sets <- function(df, est_col, fdr_col) {
  df %>%
    filter(abs(.data[[est_col]]) > EST_CUTOFF, .data[[fdr_col]] < FDR_CUTOFF) %>%
    group_by(lineage) %>%
    summarise(genes = list(unique(gene)), .groups = "drop")
}

# Build gene sets for each comparison
vsg_sets  <- get_sig_sets(final_df, "vsg_estimate",       "vsg_FDR")
std_sets  <- get_sig_sets(final_df, "standard_estimate",  "standard_FDR")
int_sets  <- get_sig_sets(final_df, "interaction_estimate","interaction_FDR")

# Prepare list for UpSetR
sets_to_list <- function(sig_df) {
  out <- sig_df$genes
  names(out) <- sig_df$lineage
  out
}

vsg_list <- sets_to_list(vsg_sets)
std_list <- sets_to_list(std_sets)
int_list <- sets_to_list(int_sets)

# Plot
plot_upset_filtered <- function(gene_list, title) {
  if (length(gene_list) < 2) return(NULL)  # need at least 2 sets
  upset(
    fromList(gene_list),
    nsets = length(gene_list),
    nintersects = 40,
    mb.ratio = c(0.6, 0.4),
    order.by = "freq",
    mainbar.y.label = paste("Intersection size"),
    main.bar.color = "#56765E",
    matrix.color = "#CBDA99",
    sets.x.label = "Significant genes per lineage",
    text.scale = c(1.2, 1.2, 0.9, 1, 1.2, 1),
    query.legend = "bottom",
    keep.order = TRUE
  )
}

#FF770F and #000026
pdf("VSG_PrePost_UpSet.pdf", width = 6, height = 5.5)
plot_upset_filtered(vsg_list, "VSG Pre–Post")
dev.off()
#28517F and #C7E1FA
pdf("std_PrePost_UpSet.pdf", width = 6, height = 5.5)
plot_upset_filtered(std_list, "Standard Pre–Post")
dev.off()
#56765E and #CBDA99
pdf("DiD_PrePost_UpSet.pdf", width = 6, height = 5.5)
plot_upset_filtered(int_list, "Treatment × Time Interaction (DiD)")
dev.off()
#

get_top_bottom_genes <- function(df, est_col, fdr_col) {
  df %>%
    filter(abs(.data[[est_col]]) >= EST_CUTOFF,
           .data[[fdr_col]] < FDR_CUTOFF) %>%
    group_by(lineage) %>%
    arrange(desc(.data[[est_col]])) %>%
    mutate(rank_pos = row_number()) %>%
    arrange(.data[[est_col]]) %>%
    mutate(rank_neg = row_number()) %>%
    ungroup() %>%
    group_by(lineage) %>%
    filter(rank_pos <= 10 | rank_neg <= 10) %>%
    ungroup()
}

vsg_df <- get_top_bottom_genes(final_df, "vsg_estimate", "vsg_FDR")
std_df <- get_top_bottom_genes(final_df, "standard_estimate", "standard_FDR")
int_df <- get_top_bottom_genes(final_df, "interaction_estimate", "interaction_FDR")

# Combine 
all_vsg <- unique(c(vsg_df$gene, std_df$gene, int_df$gene))

# Annotate 
immune_inflammatory <- c(
  "CXCL3", "DEFB1", "SPINK1", "CD24", "CD74", "CD9",
  "HLA-C", "HLA-B", "HLA-DPA1", "HLA-DPB1", "HLA-DQB1", "HLA-DRA", "HLA-E",
  "C1QA", "C1QB", "TYROBP", "HCST", "IFITM3","NFKBIA", "JUN", "JUNB", "JUND", "FOS",
  "SERPINA1", "A2M", "GAS6", "TFPI2", "ANPEP","CLU"
)


stress_response <- c(
  "CRYAB", "HSPA2", "S100A6", "S100A10", "S100A11", "S100A13", "S100A4",
  "HIF1A", "BTG2", "DUSP1","LDHA", "ALDOB", "AKR1C1", "AKR1B1",
  "GPX3", "SOD3","FTH1", "FTL"
)

fibrosis_remodeling <- c(
  "MMP7","IGFBP5","DCN","TAGLN","ITIH5","SERPINF2","PKHD1","UTRN",
  "TM4SF1","SYNE1","PLAC9"
)


metabolic_mitochondrial <- c(
  "NQO1","ESRRG","CHCHD10","COX5B","MRPS6","ATP6V0C","ATP6V0B","ATP1B1",
  "RALBP1","DCXR","IMMP2L","SUCLG1","DPYD","FABP4","FMO3","APOE"
)


transporters_ion <- c(
  "SLC12A1","SLC26A7","SLC14A1","SLC5A3","AQP2","UMOD","FXYD2","FXYD4"
)


signaling_transcription <- c(
  "PRKG1", "PDE4D", "PDE3A", "PDE1A","ERBB4","PLCL1","MECOM","ZBTB20",
  "SOX5","ZNF521","EBF1","ID1", "ID3","NR2F2-AS1", "RORA", "NRARP",
  "PTH2R","OXR1","TFAP2A","MAML2","SKAP1", "PTPRG",
  "CALCA","RGS5","PTGER3","IRX1","HDAC9"     
)

polarity_cytoskeleton <- c(
  "DOCK4", "DOCK2", "DOCK8","MAGI1","PARD3B","TMEM158","CAPG","LRMDA",
  "LDB2","ZEB2","RBMS3","BICC1","KRT7","VIM","DIAPH2","INPP4B","PLXDC2",
  "RABGAP1L","PTPRM"
)


row_order <- c(
  immune_inflammatory,
  stress_response,
  fibrosis_remodeling,
  signaling_transcription,
  metabolic_mitochondrial,
  polarity_cytoskeleton,
  transporters_ion
)

valid_lineages <- unique(final_df$lineage)
df_full <- expand_grid(
  gene = row_order,
  lineage = valid_lineages
)
#
df_plot <- df_full %>%
  left_join(
    final_df %>%
      dplyr::select(gene, lineage,
             vsg_estimate, standard_estimate, interaction_estimate),
    by = c("gene", "lineage")
  )
#
to_matrix <- function(df, value_col) {
  df_wide <- df %>%
    dplyr::select(gene, lineage, all_of(value_col)) %>%
    tidyr::pivot_wider(names_from = lineage, values_from = all_of(value_col))
  
  mat <- as.data.frame(df_wide)
  rownames(mat) <- mat$gene
  mat <- mat[, -1, drop = FALSE]
  as.matrix(mat)
}
#
mat_vsg <- to_matrix(df_plot, "vsg_estimate")
mat_std <- to_matrix(df_plot, "standard_estimate")
mat_int <- to_matrix(df_plot, "interaction_estimate")

#
groups_list <- list(
  Inflammatory_Immune   = immune_inflammatory,
  Stress_Response       = stress_response,
  Fibrosis_Remodeling   = fibrosis_remodeling,
  Signaling_Transcription = signaling_transcription,
  Metabolic_Mitochondrial = metabolic_mitochondrial,
  Polarity_Cytoskeleton = polarity_cytoskeleton,
  Transporters_Ion      = transporters_ion
)

shared_row_order <- c()
for (grp in names(groups_list)) {
  genes <- groups_list[[grp]]
  present_genes <- intersect(genes, rownames(mat_vsg))
  missing_genes <- setdiff(genes, present_genes)
  
  if (length(present_genes) > 1) {
    submat <- mat_vsg[present_genes, , drop = FALSE]
    
    # Replace NAs with 0 (or column mean) so dist() works but structure preserved
    submat[is.na(submat)] <- 0
    
    cl <- hclust(dist(submat))
    ordered <- rownames(submat)[cl$order]
  } else {
    ordered <- present_genes
  }
  
  # Append missing genes (unclustered) at the end of this group
  shared_row_order <- c(shared_row_order, ordered, missing_genes)
}

mat_vsg <- mat_vsg[shared_row_order, , drop = FALSE]
mat_std <- mat_std[shared_row_order, , drop = FALSE]
mat_int <- mat_int[shared_row_order, , drop = FALSE]

# Build a category mapping for annotation
category_colors <- c(
  Inflammatory_Immune   = "#b2182b",
  Stress_Response       = "#ef8a62",
  Fibrosis_Remodeling   = "#fddbc7",
  Signaling_Transcription = "#d1e5f0",
  Metabolic_Mitochondrial = "#67a9cf",
  Polarity_Cytoskeleton      = "#2166ac",
  Transporters_Ion      = "#053061"
)

gene_to_group <- rep(names(groups_list), lengths(groups_list))
names(gene_to_group) <- unlist(groups_list)
row_group <- gene_to_group[shared_row_order]
row_group <- factor(row_group, levels = names(groups_list))

ha_row <- rowAnnotation(
  Function = row_group,
  col = list(Function = category_colors),
  annotation_legend_param = list(
    title = "Functional Group",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 9)
  )
)
#
col_fun <- colorRamp2(c(-1.5, 0, 1.5), c("#6f6182", "white", "#6d3c38"))

# Plot heatmap
ht_vsg <- Heatmap(
  -mat_vsg, # Flip sign so that positive = upregulated post
  name = "Estimate", 
  col=col_fun,
  cluster_rows=FALSE, cluster_columns=FALSE,
  row_order=shared_row_order,
  show_row_names=TRUE, show_column_names=TRUE,
  row_names_gp=gpar(fontsize=8),
  column_names_gp=gpar(fontsize=8),
  column_title="VSG \nΔ(post − pre)",
  heatmap_legend_param = list(
    title = "Estimated Change"
  )
)

ht_std <- Heatmap(
  -mat_std, # Flip sign so that positive = upregulated post
  name = "Estimate", 
  col=col_fun,
  cluster_rows=FALSE, cluster_columns=FALSE,
  row_order=shared_row_order,
  show_row_names=TRUE, show_column_names=TRUE,
  column_names_gp=gpar(fontsize=8),
  column_title="Standard \nΔ(post − pre)",
  show_heatmap_legend = FALSE 
)

ht_int <- Heatmap(
  mat_int,
  name = "Estimate", 
  col=col_fun,
  cluster_rows=FALSE, cluster_columns=FALSE,
  row_order=shared_row_order,
  show_row_names=TRUE, show_column_names=TRUE,
  column_names_gp=gpar(fontsize=8),
  column_title="DiD\nΔ(VSG − Standard)",
  show_heatmap_legend = FALSE 
)

# Combine horizontally
ht_list <-  ha_row + ht_vsg + ht_std + ht_int
draw(ht_list,
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom")

cairo_pdf("MixedModel_DEG_Heatmaps.pdf", width = 8, height = 19)
draw(ht_list,
     #main_heatmap = "VSG pre–post",
     heatmap_legend_side = "right",
     annotation_legend_side = "right")
dev.off()

#
species <- "Homo sapiens"

msig_all <- bind_rows(
  msigdbr(species = species, category = "H"),
  msigdbr(species = species, category = "C2", subcategory = "CP:REACTOME"),
  msigdbr(species = species, category = "C2", subcategory = "CP:KEGG"),
  msigdbr(species = species, category = "C5", subcategory = "BP")
) %>%
  mutate(gs_name = str_replace_all(gs_name, "_", " "))

# Keep immune/metabolism-relevant sets
keep_patterns <- "(IMMUNE|INFLAM|CYTOKINE|INTERFERON|COMPLEMENT|ANTIGEN|LEUKOCYTE|T CELL|B CELL|MACROPHAGE|NK|CHEMOKINE|PHAGOCYT|MHC|ADAPTIVE|INNATE|OXIDATIVE PHOSPHORYLATION|FATTY ACID|LIPID|GLYCOLYSIS|TCA|TRICARBOXYLIC|MITOCHONDR|BETA OXIDATION|OXPHOS|GLUCONEOGEN|AMINO ACID|UREA CYCLE|BILE ACID|PENTOSE|PPP|PYRUVATE|GLUTATHIONE|REDOX|CARBOHYDRATE|METABOL)"
msig_imm_met <- msig_all %>%
  filter(str_detect(toupper(gs_name), keep_patterns)) %>%
  dplyr::select(gs_name, gene_symbol)

# gene sets as a list
gmt_list <- msig_imm_met %>% group_by(gs_name) %>% summarise(genes = list(unique(gene_symbol))) %>% 
  tibble::deframe()

eps <- 1e-300

mk_ranks <- function(df, est_col, fdr_col, use = c("signed_logFDR","estimate")) {
  use <- match.arg(use)
  df %>%
    mutate(rank_score = if (use == "signed_logFDR") {
      sign(.data[[est_col]]) * (-log10(pmax(.data[[fdr_col]], eps)))
    } else {
      .data[[est_col]]
    }) %>%
    filter(!is.na(rank_score), !is.na(gene), !is.na(lineage)) %>%
    group_by(lineage) %>%
    # in case of duplicated symbols within a lineage, keep the strongest
    arrange(desc(abs(rank_score))) %>%
    distinct(gene, .keep_all = TRUE) %>%
    summarise(ranks = list(setNames(rank_score, gene)), .groups = "drop")
}

ranks_vsg <- mk_ranks(final_df, "vsg_estimate", "vsg_FDR", use = "signed_logFDR")
ranks_std <- mk_ranks(final_df, "standard_estimate", "standard_FDR", use = "signed_logFDR")
ranks_int <- mk_ranks(final_df, "interaction_estimate", "interaction_FDR", use = "signed_logFDR")

min_universe <- 30     
minSize <- 10
maxSize <- 500

drop_small_universe <- function(tab) {
  tab %>% mutate(n = lengths(ranks)) %>% filter(n >= min_universe) %>% dplyr::select(-n)
}

ranks_vsg <- drop_small_universe(ranks_vsg)
ranks_std <- drop_small_universe(ranks_std)
ranks_int <- drop_small_universe(ranks_int)

run_fgsea_list <- function(ranks_tbl, contrast_label) {
  map2_df(ranks_tbl$lineage, ranks_tbl$ranks, function(lin, rv) {
    fg <- fgsea(pathways = gmt_list, stats = rv, minSize = minSize, maxSize = maxSize, nperm = 20000)
    fg %>%
      as_tibble() %>%
      mutate(lineage = lin, contrast = contrast_label) %>%
      dplyr::select(contrast, lineage, pathway, size, NES, pval, padj, leadingEdge)
  })
}

gsea_vsg <- run_fgsea_list(ranks_vsg, "VSG")
gsea_std <- run_fgsea_list(ranks_std, "Standard")
gsea_int <- run_fgsea_list(ranks_int, "Interaction")

gsea_all <- bind_rows(gsea_vsg, gsea_std, gsea_int)

sig_cut <- 0.25
gsea_sig <- gsea_all %>%
  filter(!is.na(NES), padj < sig_cut) %>%   
  group_by(contrast) %>%
  mutate(rank_pos = dplyr::dense_rank(dplyr::desc(NES)),
         rank_neg = dplyr::dense_rank(NES)) %>%
  filter(rank_pos <= 10 | rank_neg <= 10) %>%
  ungroup() %>%
  mutate(padj_capped = pmin(padj, sig_cut))

plot_gsea_gene_heat <- function(
    gsea_sig,
    vsg_list,
    contrast_name,
    padj_cutoff = 0.25,
    palette_list = list(
      VSG = c("#6f6182", "#f7f7f7", "#6d3c38"),
      Standard = c("#4575b4", "#f7f7f7", "#d73027"),
      Interaction = c("#605f74", "#f7f7f7", "#36303c")
    ),
    title_label = NULL
) {
  # Unpack gene list
  vsg_genes <- unique(unlist(vsg_list))
  message("Total ", length(vsg_genes), " genes in vsg_list")

  # Extract pathways for selected contrast
  df <- gsea_sig %>%
    filter(contrast == contrast_name, padj <= padj_cutoff)
  
  # Ensure leadingEdge is a list-column
  if (is.character(df$leadingEdge)) {
    df <- df %>%
      mutate(leadingEdge = str_split(leadingEdge, "[,/]", simplify = FALSE))
  }
  
  # Unnest the list-column of leadingEdge genes
  df <- df %>%
    unnest(leadingEdge) %>%
    mutate(leadingEdge = str_trim(as.character(leadingEdge)))
  
  message("Found ", nrow(df), " gene–pathway pairs after unnesting")
  

  # Subset to only genes in vsg_genes
  plot_df <- df %>%
    filter(leadingEdge %in% vsg_genes) %>%
    dplyr::select(pathway, leadingEdge, NES)
  
  if (nrow(plot_df) == 0) {
    warning(paste0("⚠️ No matching genes from vsg_list found in ", contrast_name, " GSEA results."))
    return(invisible(NULL))
  }
  
  # Create matrix
  mat <- plot_df %>%
    group_by(pathway, leadingEdge) %>%
    summarise(NES = mean(NES, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = leadingEdge, values_from = NES, values_fill = 0)
  
  mat_matrix <- as.matrix(mat[,-1])
  rownames(mat_matrix) <- mat$pathway
  
  # Determine NES range dynamically for balanced color scaling
  nes_range <- range(mat_matrix, na.rm = TRUE)
  limit <- max(abs(nes_range))
  

  # Convert to color scale 
  colors <- palette_list[[contrast_name]]
  if (is.null(colors)) {
    warning(paste0("⚠️ No color palette found for contrast: ", contrast_name, 
                   ". Using default palette."))
    colors <- c("#6f6182", "#f7f7f7", "#6d3c38")
  }
  col_fun <- colorRamp2(c(-limit, 0, limit), colors)
  
  if (is.null(title_label)) {
    title_label <- contrast_name
  }
  
  # Draw Heatmap
  ht <- Heatmap(
    mat_matrix,
    name = "NES",
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    row_names_gp = gpar(fontsize = 8),
    row_names_max_width = unit(15, "cm"),
    column_names_gp = gpar(fontsize = 9, fontface = "bold"),
    column_title = paste(title_label),
    show_heatmap_legend = FALSE
  )
  
  lgd <- Legend(
    col_fun = col_fun,
    title = "NES",
    direction = "horizontal",
    title_position = "lefttop",
    at = c(-limit, 0, limit),
    labels = round(c(-limit, 0, limit), 2),
    legend_width = unit(1.5, "cm"),
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  )
  
  # Pack the legend and title together horizontally
  legend_pack <- packLegend(lgd, direction = "horizontal")
  draw(ht,
       heatmap_legend_side = "top",
       annotation_legend_side = "bottom",
       annotation_legend_list = list(legend_pack),
       padding = unit(c(10, 20, 10, 40), "mm"))
  invisible(ht)
}
pdf("GSEA_Gene_Heatmaps_VSG.pdf", width = 8, height = 3)
plot_gsea_gene_heat(gsea_sig, vsg_list, "VSG", title_label = "VSG pre vs. post")
dev.off()

pdf("GSEA_Gene_Heatmaps_STD.pdf", width = 13, height = 4.5)
plot_gsea_gene_heat(gsea_sig, std_list, "Standard", title_label = "Standard pre vs. post")
dev.off()

pdf("GSEA_Gene_Heatmaps_DiD.pdf", width = 11, height = 4.5)
plot_gsea_gene_heat(gsea_sig, int_list, "Interaction", title_label = "DiD")
dev.off()

#
chunk_vec <- function(x, size = 10) split(x, ceiling(seq_along(x) / size))

get_gene_to_compounds <- function(gene_symbols, organism = "hsa", sleep_sec = 0.2) {
  gene_symbols <- unique(gene_symbols)
  
  # Map SYMBOL → ENTREZ
  map_tbl <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = gene_symbols,
    columns = c("ENTREZID", "SYMBOL"),
    keytype = "SYMBOL"
  ) %>%
    dplyr::distinct() %>%
    dplyr::filter(!is.na(ENTREZID))
  
  if (nrow(map_tbl) == 0) {
    warning("No symbols mapped to Entrez IDs.")
    return(tibble())
  }
  
  map_tbl <- map_tbl %>%
    dplyr::mutate(KEGG_GENE = paste0(organism, ":", ENTREZID))
  
  message("Mapping via KO → Reaction → Compound")
  
  # GENE → KO 
  ko_links <- map_tbl$KEGG_GENE %>%
    purrr::map(function(g) {
      Sys.sleep(sleep_sec)
      tryCatch(keggLink("ko", g), error = function(e) character(0))
    }) %>%
    unlist()
  
  if (length(ko_links) == 0) {
    warning("No KO links found for any genes.")
    return(tibble())
  }
  
  gene_ko_df <- tibble(
    KEGG_GENE = names(ko_links),
    KO = unname(ko_links)
  )
  
  # KO → Reaction 
  ko_rxn_links <- gene_ko_df$KO %>%
    unique() %>%
    purrr::map(function(k) {
      Sys.sleep(sleep_sec)
      tryCatch(keggLink("reaction", k), error = function(e) character(0))
    }) %>%
    unlist()
  
  if (length(ko_rxn_links) == 0) {
    warning("No reactions found for any KO IDs.")
    return(tibble())
  }
  
  ko_rxn_df <- tibble(
    KO = names(ko_rxn_links),
    RN = unname(ko_rxn_links)
  )
  
  # Reaction → Compound 
  rn_cpd_links <- ko_rxn_df$RN %>%
    unique() %>%
    purrr::map(function(r) {
      Sys.sleep(sleep_sec)
      tryCatch(keggLink("cpd", r), error = function(e) character(0))
    }) %>%
    unlist()
  
  if (length(rn_cpd_links) == 0) {
    warning("No compound links found for any reactions.")
    return(tibble())
  }
  
  rn_cpd_df <- tibble(
    RN = names(rn_cpd_links),
    CPD = unname(rn_cpd_links)
  )
  
  # Merge all mappings
  gene_cpd_all <- gene_ko_df %>%
    dplyr::inner_join(ko_rxn_df, by = "KO") %>%
    dplyr::inner_join(rn_cpd_df, by = "RN") %>%
    dplyr::left_join(map_tbl, by = "KEGG_GENE") %>%
    dplyr::select(SYMBOL, ENTREZID, KO, RN, CPD) %>%
    dplyr::distinct()
  
  # Fetch compound names
  cpd_ids <- unique(gene_cpd_all$CPD)
  cpd_tbl <- tibble()
  
  if (length(cpd_ids) > 0) {
    for (cid in cpd_ids) {
      Sys.sleep(sleep_sec)
      result <- tryCatch({
        entry <- keggGet(cid)[[1]]
        tibble(
          CPD = cid,
          CompoundName = if (!is.null(entry$NAME)) entry$NAME[1] else NA_character_
        )
      }, error = function(e) {
        # fallback: just keep the ID if KEGGREST fails to parse
        tibble(CPD = cid, CompoundName = NA_character_)
      })
      cpd_tbl <- dplyr::bind_rows(cpd_tbl, result)
    }
  }
  
  # Combine 
  out <- gene_cpd_all %>%
    dplyr::left_join(cpd_tbl, by = "CPD") %>%
    dplyr::arrange(SYMBOL, CPD)
  
  return(out)
}
#
result_vsg <- get_gene_to_compounds(unique(unlist(vsg_list)))
result_std <- get_gene_to_compounds(unique(unlist(std_list)))
result_int <- get_gene_to_compounds(unique(unlist(int_list)))
# 
annotate_direction <- function(result_tbl, final_df, est_col, fdr_col, label_prefix) {
  # Build lookup table: gene → regulation category
  df_dir <- final_df %>%
    dplyr::select(gene, !!est_col := !!sym(est_col), !!fdr_col := !!sym(fdr_col)) %>%
    dplyr::mutate(
      regulation = dplyr::case_when(
        .data[[fdr_col]] >= FDR_CUTOFF ~ NA_character_,
        .data[[est_col]] < -EST_CUTOFF ~ paste("Upregulated in", label_prefix),
        .data[[est_col]] >  EST_CUTOFF ~ paste("Downregulated in", label_prefix),
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(regulation)) %>%
    dplyr::distinct(gene, regulation)
  
  # Join by SYMBOL (gene name)
  annotated <- result_tbl %>%
    dplyr::left_join(df_dir, by = c("SYMBOL" = "gene"))
  
  return(annotated)
}

result_vsg <- annotate_direction(result_vsg, final_df, "vsg_estimate", "vsg_FDR", "VSG (post vs pre)")
result_std <- annotate_direction(result_std, final_df, "standard_estimate", "standard_FDR", "Standard (post vs pre)")
result_int <- annotate_direction(result_int, final_df, "interaction_estimate", "interaction_FDR", "VSG vs Standard DiD")

write.csv(result_vsg, "gene_to_metabolites_vsg.csv", row.names = FALSE)
write.csv(result_std, "gene_to_metabolites_standard.csv", row.names = FALSE)
write.csv(result_int, "gene_to_metabolites_DiD.csv", row.names = FALSE)

## ---- Pseudobulk at Lineage level  ----
cts <- as.matrix(expr_matrix)

# One column per lineage|record_id|pre_post pseudobulk
meta$pb_key_lineage <- with(meta, paste(lineage, record_id, pre_post, sep="|"))

split_idx <- split(seq_len(nrow(meta)), meta$pb_key_lineage)
pb_cts_lineage <- do.call(cbind, lapply(split_idx, function(v) Matrix::rowSums(cts[, v, drop=FALSE])))
colnames(pb_cts_lineage) <- names(split_idx)

# sample metadata for those columns
pb_meta_lineage <- do.call(rbind, strsplit(colnames(pb_cts_lineage), "\\|")) |> as.data.frame()
colnames(pb_meta_lineage) <- c("lineage","record_id","pre_post")

# add treatment and BMI (matched by patient & time)
pb_meta_lineage <- pb_meta_lineage %>%
  left_join(meta %>% distinct(record_id, treatment), by = "record_id") %>%
  left_join(meta %>% distinct(record_id, pre_post, bmi), by = c("record_id","pre_post"))

# keep HVGs
keep_genes <- intersect(rownames(pb_cts_lineage), top5000_union)
pb_cts_lineage <- pb_cts_lineage[keep_genes, ]

#
fit_lineage_list <- list()
for (L in sort(unique(pb_meta_lineage$lineage))) {
  message("Processing lineage: ", L)
  
  sel   <- pb_meta_lineage$lineage == L
  metaL <- pb_meta_lineage[sel, , drop = FALSE]
  if (nrow(metaL) < 3) next  # too small
  
  # Dense conversion 
  cts_sub <- pb_cts_lineage[, sel, drop = FALSE]
  cts_sub <- as(cts_sub, "matrix")           
  cts_sub[is.na(cts_sub)] <- 0
  cts_sub <- cts_sub[rowSums(cts_sub) > 0, , drop = FALSE]
  if (ncol(cts_sub) < 2 || nrow(cts_sub) < 10) next
  
  # DGEList + norm factors
  dge <- edgeR::DGEList(counts = cts_sub)
  dge <- edgeR::calcNormFactors(dge)
  
  # covariates
  metaL <- metaL %>% mutate(
    BMI_c    = as.numeric(scale(bmi, scale = FALSE)),
    pre_post = factor(pre_post, levels = c("pre","post")),
    treatment = factor(treatment, levels = c("Standard","VSG"))
  )
  
  # design
  design <- model.matrix(~ BMI_c * treatment + pre_post, data = metaL)
  stopifnot(ncol(dge$counts) == nrow(design))
  

  v <- limma::voom(dge, design, plot = FALSE)
  dup <- limma::duplicateCorrelation(v, design, block = metaL$record_id)
  
  fit <- limma::lmFit(v, design, block = metaL$record_id, correlation = dup$consensus)
  fit <- limma::eBayes(fit)
  
  coef_names <- colnames(fit$coefficients)
  
  # BMI in Standard arm 
  tabs <- list()
  if ("BMI_c" %in% coef_names) {
    tab_bmi <- limma::topTable(fit, coef = "BMI_c", number = Inf, sort.by = "P") |>
      tibble::rownames_to_column("gene") |>
      dplyr::mutate(lineage = L, term = "BMI")
    tabs$BMI <- tab_bmi
  }
  
  # Interaction 
  has_int <- "BMI_c:treatmentVSG" %in% coef_names
  if (has_int) {
    tab_int <- limma::topTable(fit, coef = "BMI_c:treatmentVSG", number = Inf, sort.by = "P") |>
      tibble::rownames_to_column("gene") |>
      dplyr::mutate(lineage = L, term = "BMI:treatmentVSG")
    tabs$INT <- tab_int
    
    # Per-treatment slope in VSG = BMI_c + interaction
    colnames(design) <- make.names(colnames(design))
    colnames(fit$coefficients) <- make.names(colnames(fit$coefficients))
    
    # now the interaction term will be BMI_c.treatmentVSG instead of BMI_c:treatmentVSG
    contr.mat <- limma::makeContrasts(
      BMI_in_VSG = BMI_c + BMI_c.treatmentVSG,
      levels = design
    )
    fit2 <- limma::contrasts.fit(fit, contr.mat) |> limma::eBayes()
    tab_bmi_vsg <- limma::topTable(fit2, coef = "BMI_in_VSG", number = Inf, sort.by = "P") |>
      tibble::rownames_to_column("gene") |>
      dplyr::mutate(lineage = L, term = "BMI_in_VSG")
    tabs$BMI_VSG <- tab_bmi_vsg
  } else {
    # If only one treatment present in this lineage subset, report BMI slope for that treatment only
    present_trt <- as.character(unique(metaL$treatment))
    note <- paste0("BMI (only ", paste(present_trt, collapse="/"), " present)")
    # The BMI_c column is the slope for whatever treatment is the reference in the design
    if ("BMI_c" %in% coef_names) {
      tabs$BMI_only <- limma::topTable(fit, coef = "BMI_c", number = Inf, sort.by = "P") |>
        tibble::rownames_to_column("gene") |>
        dplyr::mutate(lineage = L, term = note)
    }
  }
  
  if (length(tabs)) {
    fit_lineage_list[[L]] <- dplyr::bind_rows(tabs)
  }
}

res_lineage <- dplyr::bind_rows(fit_lineage_list)%>%
  filter(lineage != "Other")
write_csv(res_lineage, "BMI_associated_genes_lineage.csv")

plot_df <- res_lineage %>%
  dplyr::filter(term %in% c("BMI", "BMI_in_VSG"), AveExpr > 1) %>%
  dplyr::mutate(
    Treatment = dplyr::recode(term,
                              "BMI" = "Standard",
                              "BMI_in_VSG" = "VSG"),
    label = paste0(gene, " (", lineage, ")")
  ) %>%
  dplyr::arrange(Treatment, adj.P.Val, dplyr::desc(abs(logFC))) %>%
  dplyr::group_by(Treatment) %>%
  dplyr::mutate(rank = dplyr::row_number()) %>%
  dplyr::ungroup()

pdf("BMI_lineage_genes.pdf", width = 12, height = 4)
ggplot(plot_df %>% dplyr::filter(rank <= 30),
       aes(x = rank, y = -log10(adj.P.Val),
           color = logFC, size = abs(logFC))) +
  geom_point(alpha = 0.8) +
  geom_text_repel(
    aes(label = ifelse(rank <= 10, label, "")),
    size = 3.5, max.overlaps = Inf, box.padding =0.8, segment.color = "grey70",
    force = 8, point.padding = 0.4
  ) +
  scale_color_gradient2(
    low = "#b87333", mid = "#f0ead2", high = "#242e16", midpoint = 0,
    name = "logFC"
  ) +
  facet_wrap(~Treatment, scales = "free_x") +
  theme_minimal(base_size = 12) +
  labs(
    y = "-log10(FDR)",
    x = "Ranked genes",
  ) +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
dev.off()

#
df_std <- res_lineage %>% filter(term == "BMI") %>% dplyr::select(gene, lineage, std_FDR = adj.P.Val)
df_vsg <- res_lineage %>% filter(term == "BMI_in_VSG") %>% dplyr::select(gene, lineage, vsg_FDR = adj.P.Val)

compare_df <- inner_join(df_std, df_vsg, by = c("gene", "lineage"))

lineage_colors <- c(
  "Myeloid"     = "#202127",  # deep red
  "Lymphoid"    = "#a98744",  # deep blue
  "PT"          = "#BFD3E6",
  "TAL"         = "#D1E5F0",
  "IC"          = "#568897",
  "PC"          = "#bad2e1",
  "EC"          = "#76bbbb",
  "FIB/VSMC/P"  = "#FEE090",
  "POD"         = "#c0c9b6",
  "Schwann"     = "#E5E5E5"
)

pdf("BMI_lineage_genes_scatter.pdf", width = 9, height = 6)
ggplot(compare_df, aes(x = -log10(std_FDR), y = -log10(vsg_FDR), color = lineage)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_text_repel(
    data = subset(compare_df, std_FDR < 0.05 | vsg_FDR < 0.05),
    aes(label = gene),
    size = 3, max.overlaps = 20, box.padding = 0.5, segment.color = "grey70"
  ) +
  scale_color_manual(
    name = "Lineage",
    values = lineage_colors,
    guide = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      override.aes = list(size = 4, alpha = 1)   # larger legend dots, full opacity
    )
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    legend.key.size = unit(1.6, "lines"),
    legend.key.width = unit(1.4, "lines"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13, face = "bold"),
    panel.grid.minor = element_blank(),
    strip.text = element_blank()
  ) +
  labs(
    x = expression(-log[10]*"(FDR)[Standard]"),
    y = expression(-log[10]*"(FDR)[VSG]")
  )
dev.off()
