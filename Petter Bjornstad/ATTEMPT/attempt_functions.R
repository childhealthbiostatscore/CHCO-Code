
# ATTEMPT analysis related functions


# ---- Volcano Plot Functions ----

# ===========================================================================
# Function: plot_volcano
# ===========================================================================
plot_volcano <- function(data, fc, p_col, title_suffix, x_axis, y_axis, file_suffix, p_thresh = 0.05) {
  set.seed(1)
  top_pos <- data %>%
    dplyr::filter(!!sym(fc) > 0 & !!sym(p_col) < p_thresh) %>%
    dplyr::arrange(!!sym(p_col)) %>%
    slice_head(n=20)
  
  top_neg <- data %>%
    dplyr::filter(!!sym(fc) < 0 & !!sym(p_col) < p_thresh) %>%
    dplyr::arrange(!!sym(p_col)) %>%
    slice_head(n=20)
  
  data <- data %>%
    dplyr::mutate(top_color = case_when(Gene %in% top_pos$Gene ~ "#f28482",
                                        Gene %in% top_neg$Gene ~ "#457b9d",
                                        TRUE ~ "#ced4da"),
                  top_size = if_else(Gene %in% c(top_pos$Gene, top_neg$Gene), 1.3, 1),
                  top_lab  = if_else(Gene %in% c(top_pos$Gene, top_neg$Gene), Gene, ""))
  
  p <- ggplot(data, aes(x = !!sym(fc), y = -log10(!!sym(p_col)))) +
    geom_hline(yintercept = -log10(p_thresh), linetype = "dashed", color = "darkgrey") +
    geom_point(alpha = 0.5, aes(color = top_color, size = top_size)) +
    geom_text_repel(aes(label = top_lab, color = top_color),
                    size = 3, max.overlaps = Inf,
                    force = 6, segment.alpha = 0.3, segment.size = 0.3) +
    labs(title = paste(title_suffix),
         x = paste(x_axis),
         y = paste(y_axis)) +
    scale_size_continuous(range = c(1, 1.3)) + 
    scale_color_manual(values = c("#457b9d"="#457b9d", "#ced4da"="#ced4da", "#f28482"="#f28482")) +
    theme_minimal() +
    theme(legend.title = element_blank(),
          panel.grid = element_blank(),
          text = element_text(size = 15),
          title = element_text(size = 9)) +
    guides(color = "none", size = "none")
  
  ggsave(paste0("/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/ATTEMPT/Results/Figures/Volcano Plots/", file_suffix, ".jpeg"), plot = p, width = 7, height = 5)
  return(p)
}

# ===========================================================================
# Function: plot_volcano_associations
# ===========================================================================
plot_volcano_associations <- function(clin_results, fc, p_col, title_suffix, 
                                      x_axis, y_axis, file_suffix, p_thresh = 0.05,  treatment_results, treatment_p_col, 
                                      positive_text = "Positive with Dapagliflozin", 
                                      negative_text = "Negative with Dapagliflozin",
                                      formula = "\u0394 treatment", color_by = "treatment_direction",
                                      legend_position = c(0.3,0.8)) {
  
  set.seed(1)
  
  # Add treatment effect directions
  treatment_results <- treatment_results %>%
    mutate(treatment_direction = case_when(`p_treatmentDapagliflozin:visitPOST` < 0.05 &
                                             `logFC_treatmentDapagliflozin:visitPOST` > 0 ~ "Up w/ Dapagliflozin",
                                           `p_treatmentDapagliflozin:visitPOST` < 0.05 &
                                             `logFC_treatmentDapagliflozin:visitPOST` < 0 ~ "Down w/ Dapagliflozin",
                                           TRUE ~ "NS"))
  
  # Genes significant in both models
  sig_clin_genes <- clin_results %>%
    filter(!!sym(p_col) < p_thresh) %>%
    pull(Gene)
  
  sig_treat_genes <- treatment_results %>%
    filter(!!sym(treatment_p_col) < p_thresh) %>%
    pull(Gene)
  
  sig_both_genes <- intersect(sig_clin_genes, sig_treat_genes)
  
  # Top 5 left pos / neg and right pos / neg
  top_5_left_pos <- clin_results %>%
    filter(Gene %in% sig_both_genes, !!sym(fc) < 0) %>%
    left_join(treatment_results %>% 
                dplyr::select(Gene, `logFC_treatmentDapagliflozin:visitPOST`), by = "Gene") %>%
    filter(`logFC_treatmentDapagliflozin:visitPOST` > 0) %>%
    arrange(!!sym(p_col)) %>%
    slice_head(n = 5) %>%
    pull(Gene)
  
  top_5_left_neg <- clin_results %>%
    filter(Gene %in% sig_both_genes, !!sym(fc) < 0) %>%
    left_join(treatment_results %>% 
                dplyr::select(Gene, `logFC_treatmentDapagliflozin:visitPOST`), by = "Gene") %>%
    filter(`logFC_treatmentDapagliflozin:visitPOST` < 0) %>%
    arrange(!!sym(p_col)) %>%
    slice_head(n = 5) %>%
    pull(Gene)
  
  top_5_right_pos <- clin_results %>%
    filter(Gene %in% sig_both_genes, !!sym(fc) > 0) %>%
    left_join(treatment_results %>% 
                dplyr::select(Gene, `logFC_treatmentDapagliflozin:visitPOST`), by = "Gene") %>%
    filter(`logFC_treatmentDapagliflozin:visitPOST` > 0) %>%
    arrange(!!sym(p_col)) %>%
    slice_head(n = 5) %>%
    pull(Gene)
  
  top_5_right_neg <- clin_results %>%
    filter(Gene %in% sig_both_genes, !!sym(fc) > 0) %>%
    left_join(treatment_results %>% 
                dplyr::select(Gene, `logFC_treatmentDapagliflozin:visitPOST`), by = "Gene") %>%
    filter(`logFC_treatmentDapagliflozin:visitPOST` < 0) %>%
    arrange(!!sym(p_col)) %>%
    slice_head(n = 5) %>%
    pull(Gene)
  
  # Combine for plotting
  top_label_genes <- c(top_5_left_pos, top_5_left_neg, top_5_right_pos, top_5_right_neg)
  
  message("Number of genes significant in both: ", length(sig_both_genes))
  message("Number of labeled genes: ", length(top_label_genes))
  
  clin_results <- clin_results %>%
    mutate(shape_var_plot = case_when(
      Gene %in% sig_both_genes & !!sym(fc) < 0 & treatment_results$treatment_direction[match(Gene, treatment_results$Gene)] == "Up w/ Dapagliflozin" ~ "Left_Pos",
      Gene %in% sig_both_genes & !!sym(fc) < 0 & treatment_results$treatment_direction[match(Gene, treatment_results$Gene)] == "Down w/ Dapagliflozin" ~ "Left_Neg",
      Gene %in% sig_both_genes & !!sym(fc) > 0 & treatment_results$treatment_direction[match(Gene, treatment_results$Gene)] == "Up w/ Dapagliflozin" ~ "Right_Pos",
      Gene %in% sig_both_genes & !!sym(fc) > 0 & treatment_results$treatment_direction[match(Gene, treatment_results$Gene)] == "Down w/ Dapagliflozin" ~ "Right_Neg",
      TRUE ~ "NS"
    ),
    color_var_plot = case_when(
      Gene %in% sig_both_genes ~ treatment_results$treatment_direction[match(Gene, treatment_results$Gene)],
      TRUE ~ "NS"
    ),
    top_lab = if_else(Gene %in% top_label_genes, Gene, "")
    )
  
  clin_results$shape_var_plot <- factor(clin_results$shape_var_plot, 
                                        levels = c("Left_Pos", "Left_Neg", "Right_Pos", "Right_Neg", "NS"))
  
  clin_results$color_var_plot <- factor(clin_results$color_var_plot, 
                                        levels = c("Up w/ Dapagliflozin", "Down w/ Dapagliflozin", "NS"))
  
  # Max and min for annotation arrows
  max_fc <- max(clin_results[[fc]], na.rm = TRUE)
  min_fc <- min(clin_results[[fc]], na.rm = TRUE)
  
  # Plot
  p <- ggplot(clin_results, aes(x = !!sym(fc), y = -log10(!!sym(p_col)))) +
    geom_hline(yintercept = -log10(p_thresh), linetype = "dashed", color = "darkgrey") +
    geom_point(alpha = 0.3, aes(fill = color_var_plot, 
                                color = color_var_plot,
                                shape = shape_var_plot), size = 2) +
    scale_shape_manual(values = c("Left_Pos" = 23, "Left_Neg" = 23, "Right_Pos" = 22, "Right_Neg" = 22, "NS" = 21),
                       labels = c("Left_Pos" = "Left Pos", "Left_Neg" = "Left Neg", 
                                  "Right_Pos" = "Right Pos", "Right_Neg" = "Right Neg", "NS" = "NS")) +
    scale_color_manual(values = c("Up w/ Dapagliflozin" = "#f28482", "Down w/ Dapagliflozin" = "#457b9d", "NS" = "#ced4da")) +
    scale_fill_manual(values = c("Up w/ Dapagliflozin" = "#f28482", "Down w/ Dapagliflozin" = "#457b9d", "NS" = "#ced4da")) +
    guides(color = guide_legend(title = NULL), fill = "none") +
    geom_text_repel(aes(label = top_lab, color = color_var_plot),
                    size = 3, max.overlaps = Inf, force = 10,
                    segment.alpha = 0.5, segment.size = 0.4,
                    min.segment.length = 0, box.padding = 0.6, point.padding = 0.4,
                    segment.color = "#ced4da") +
    labs(title = paste(title_suffix),
         x = paste(x_axis),
         y = paste(y_axis)) +
    theme_minimal() +
    annotate("segment", 
             x=max_fc/8, 
             xend=(max_fc*7)/8, 
             y=-1.5,
             col="darkgrey", arrow=arrow(length=unit(0.2, "cm"))) +
    annotate("text", 
             x=mean(c(max_fc/8, (max_fc*7)/8)), 
             y=-2, 
             label=positive_text,
             size=3, color="#343a40") +
    annotate("segment", 
             x=min_fc/8, 
             xend=(min_fc*7)/8, 
             y=-1.5,
             col="darkgrey", arrow=arrow(length=unit(0.2, "cm"))) +
    annotate("text", 
             x=mean(c(min_fc/8, (min_fc*7)/8)), 
             y=-2, 
             label=negative_text,
             size=3, color="#343a40") +
    annotate("text", 
             x=max_fc * 0.95,
             y=-3.5, 
             hjust = 1,
             label=paste0("Formula: ~ ", formula, " + (1|subject)"),
             size=3, color="#343a40") +
    scale_y_continuous(expand=c(0,0)) +
    coord_cartesian(ylim = c(0, 15), clip="off") +
    theme(legend.title = element_blank(),
          panel.grid = element_blank(),
          text = element_text(size = 9),
          legend.position = legend_position,
          legend.justification = c(1, 1),
          plot.margin = margin(t = 10, r = 20, b = 28, l = 20),
          axis.title.x = element_text(margin = margin(t = 32))) + 
    guides(shape = "none")
  
  # Save
  ggsave(paste0("/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/ATTEMPT/Results/Figures/Volcano Plots/", file_suffix, ".jpeg"), plot = p, width = 7, height = 5)
  
  return(p)
}

# ===========================================================================
# Function: plot_volcano_concordance
# ===========================================================================

plot_volcano_concordance <- function(clin_results, fc, p_col, 
                                     x_axis, y_axis, file_suffix, p_thresh = 0.05,  treatment_results, treatment_p_col, 
                                     positive_text = "Positive with Dapagliflozin", 
                                     negative_text = "Negative with Dapagliflozin",
                                     formula = "\u0394 treatment", color_by = "treatment_direction",
                                     legend_position = "top",
                                     clinical_direction = NULL,
                                     caption_text = "Point position reflects association with clinical variable; point color indicates treatment effect direction.\nPoints are colored if concordant with clinical variable direction after treatment. \nUp to top 50 from each direction are labeled.",
                                     cell_type = "") {
  
  set.seed(1)
  
  # Add treatment effect directions
  treatment_results <- treatment_results %>%
    mutate(treatment_direction = case_when(`p_treatmentDapagliflozin:visitPOST` < 0.05 &
                                             `logFC_treatmentDapagliflozin:visitPOST` > 0 ~ "Up w/ Dapagliflozin",
                                           `p_treatmentDapagliflozin:visitPOST` < 0.05 &
                                             `logFC_treatmentDapagliflozin:visitPOST` < 0 ~ "Down w/ Dapagliflozin",
                                           TRUE ~ "NS/NC"))
  
  # Calculate concordant genes based on clinical direction
  concordant_genes <- c()
  if (!is.null(clinical_direction)) {
    if (clinical_direction == "-") {
      # For negative clinical direction: clinical - and treatment +
      concordant_genes <- clin_results %>%
        left_join(treatment_results %>% 
                    dplyr::select(Gene, `logFC_treatmentDapagliflozin:visitPOST`, `p_treatmentDapagliflozin:visitPOST`), by = "Gene") %>%
        filter((!!sym(fc) > 0 & !!sym(p_col) < 0.05) & (`logFC_treatmentDapagliflozin:visitPOST` < 0 & `p_treatmentDapagliflozin:visitPOST` < 0.05)|
                 (!!sym(fc) < 0 & !!sym(p_col) < 0.05) & (`logFC_treatmentDapagliflozin:visitPOST` > 0 & `p_treatmentDapagliflozin:visitPOST` < 0.05)) %>%
        pull(Gene)
    } else if (clinical_direction == "+") {
      # For positive clinical direction: clinical + and treatment -
      concordant_genes <- clin_results %>%
        left_join(treatment_results %>% 
                    dplyr::select(Gene, `logFC_treatmentDapagliflozin:visitPOST`, `p_treatmentDapagliflozin:visitPOST`), by = "Gene") %>%
        filter((!!sym(fc) > 0 & !!sym(p_col) < 0.05) & (`logFC_treatmentDapagliflozin:visitPOST` > 0 & `p_treatmentDapagliflozin:visitPOST` < 0.05)|
                 (!!sym(fc) < 0 & !!sym(p_col) < 0.05) & (`logFC_treatmentDapagliflozin:visitPOST` < 0 & `p_treatmentDapagliflozin:visitPOST` < 0.05)) %>%
        pull(Gene)
    }
  }
  
  n_concordant <- length(concordant_genes)
  sig_clin_genes <- clin_results %>%
    filter(!!sym(p_col) < p_thresh) %>%
    pull(Gene)
  
  sig_treat_genes <- treatment_results %>%
    filter(!!sym(treatment_p_col) < p_thresh) %>%
    pull(Gene)
  
  sig_both_genes <- intersect(sig_clin_genes, sig_treat_genes)
  
  # Label top 50 concordant genes on each side
  top_50_left <- clin_results %>%
    filter(Gene %in% concordant_genes, !!sym(fc) < 0) %>%
    arrange(!!sym(p_col)) %>%
    slice_head(n = 50) %>%
    pull(Gene)
  
  top_50_right <- clin_results %>%
    filter(Gene %in% concordant_genes, !!sym(fc) > 0) %>%
    arrange(!!sym(p_col)) %>%
    slice_head(n = 50) %>%
    pull(Gene)
  
  # Combine for plotting
  top_label_genes <- c(top_50_left, top_50_right)
  
  message("Number of genes significant in both: ", length(sig_both_genes))
  message("Number of concordant genes: ", n_concordant)
  message("Number of labeled genes: ", length(top_label_genes))
  
  clin_results <- clin_results %>%
    mutate(
      # Color only concordant genes
      color_var_plot = case_when(
        Gene %in% concordant_genes ~ treatment_results$treatment_direction[match(Gene, treatment_results$Gene)],
        TRUE ~ "NS/NC"
      ),
      # Shape based on significance and direction
      shape_var_plot = case_when(
        Gene %in% sig_both_genes & !!sym(fc) < 0 & treatment_results$treatment_direction[match(Gene, treatment_results$Gene)] == "Up w/ Dapagliflozin" ~ "Left_Pos",
        Gene %in% sig_both_genes & !!sym(fc) < 0 & treatment_results$treatment_direction[match(Gene, treatment_results$Gene)] == "Down w/ Dapagliflozin" ~ "Left_Neg",
        Gene %in% sig_both_genes & !!sym(fc) > 0 & treatment_results$treatment_direction[match(Gene, treatment_results$Gene)] == "Up w/ Dapagliflozin" ~ "Right_Pos",
        Gene %in% sig_both_genes & !!sym(fc) > 0 & treatment_results$treatment_direction[match(Gene, treatment_results$Gene)] == "Down w/ Dapagliflozin" ~ "Right_Neg",
        TRUE ~ "NS/NC"
      ),
      # Label only top concordant genes
      top_lab = if_else(Gene %in% top_label_genes, Gene, "")
    )
  
  clin_results$shape_var_plot <- factor(clin_results$shape_var_plot, 
                                        levels = c("Left_Pos", "Left_Neg", "Right_Pos", "Right_Neg", "NS/NC"))
  
  clin_results$color_var_plot <- factor(clin_results$color_var_plot, 
                                        levels = c("Up w/ Dapagliflozin", "Down w/ Dapagliflozin", "NS/NC"))
  
  # Max and min for annotation arrows
  max_fc <- max(clin_results[[fc]], na.rm = TRUE)
  min_fc <- min(clin_results[[fc]], na.rm = TRUE)
  
  # Get y-axis max for dynamic scaling
  y_max <- max(-log10(clin_results[[p_col]]), na.rm = TRUE) * 1.1
  
  # Plot (no background shading)
  p <- ggplot(clin_results, aes(x = !!sym(fc), y = -log10(!!sym(p_col)))) +
    geom_hline(yintercept = -log10(p_thresh), linetype = "dashed", color = "darkgrey") +
    geom_point(alpha = 0.3, aes(fill = color_var_plot, 
                                color = color_var_plot,
                                shape = shape_var_plot), size = 2) +
    scale_shape_manual(values = c("Left_Pos" = 22, "Left_Neg" = 22, "Right_Pos" = 22, "Right_Neg" = 22, "NS/NC" = 22),
                       labels = c("Left_Pos" = "Left Pos", "Left_Neg" = "Left Neg", 
                                  "Right_Pos" = "Right Pos", "Right_Neg" = "Right Neg", "NS/NC" = "NS/NC")) +
    scale_color_manual(values = c("Up w/ Dapagliflozin" = "#f28482", "Down w/ Dapagliflozin" = "#457b9d", "NS/NC" = "#ced4da")) +
    scale_fill_manual(values = c("Up w/ Dapagliflozin" = "#f28482", "Down w/ Dapagliflozin" = "#457b9d", "NS/NC" = "#ced4da")) +
    guides(color = guide_legend(title = NULL, 
                                override.aes = list(label = ""), 
                                nrow = 1), 
           fill = "none",
           shape = "none") +
    geom_text_repel(aes(label = top_lab, color = color_var_plot),
                    size = 3, max.overlaps = Inf, force = 10,
                    segment.alpha = 0.5, segment.size = 0.4,
                    min.segment.length = 0, box.padding = 0.6, point.padding = 0.4,
                    segment.color = "#ced4da") +
    labs(title = NULL,
         x = paste(x_axis),
         y = paste(y_axis),
         caption = if (!is.null(clinical_direction)) {
           paste0("Formula: ~ ", formula, " + (1|subject)", 
                  "\n\n Cell type: ", 
                  cell_type, if(cell_type != "") " | " else "", 
                  "Concordant genes: n = ", n_concordant,
                  "\n\n", caption_text)
         } else NULL) +
    theme_minimal() +
    annotate("segment", 
             x=max_fc/8, 
             xend=(max_fc*7)/8, 
             y=-y_max * 0.09,
             col="darkgrey", arrow=arrow(length=unit(0.2, "cm"))) +
    annotate("text", 
             x=mean(c(max_fc/8, (max_fc*7)/8)), 
             y=-y_max * 0.14, 
             label=positive_text,
             size=3, color="#343a40") +
    annotate("segment", 
             x=min_fc/8, 
             xend=(min_fc*7)/8, 
             y=-y_max * 0.09,
             col="darkgrey", arrow=arrow(length=unit(0.2, "cm"))) +
    annotate("text", 
             x=mean(c(min_fc/8, (min_fc*7)/8)), 
             y=-y_max * 0.14, 
             label=negative_text,
             size=3, color="#343a40") +
    scale_y_continuous(expand=c(0,0)) +
    coord_cartesian(ylim = c(0, y_max), clip="off") +
    theme(legend.title = element_blank(),
          panel.grid = element_blank(),
          text = element_text(size = 9, family = "Arial"),
          legend.position = legend_position,
          legend.justification = if(is.character(legend_position) && legend_position == "bottom") c(0.5, 0) else c(1, 1),
          legend.direction = "horizontal",
          legend.spacing.x = unit(0.3, 'cm'),
          plot.margin = margin(t = 10, r = 20, b = 25, l = 20),
          axis.title.x = element_text(margin = margin(t = 28)),
          plot.caption = element_text(size = 8.5, hjust = 0.5, margin = margin(t = 8), family = "Arial"),
          legend.margin = margin(t = 5, b = 5)) + 
    guides(shape = "none")
  
  # Save
  ggsave(paste0("/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/ATTEMPT/Results/Figures/Volcano Plots/", file_suffix, "_concordance.jpeg"), plot = p, width = 7, height = 5)
  
  return(p)
}


# ---- Pathway Plot Functions ----


# ===========================================================================
# Function: matrix_to_list
# ===========================================================================
matrix_to_list <- function(pws){
  pws.l <- list()
  for (pw in colnames(pws)) {
    pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
  }
  return(pws.l)
}
# ===========================================================================
# Function: prepare_gmt
# ===========================================================================
prepare_gmt <- function(gmt_file, genes_in_data, savefile = FALSE){
  # for debug
  #file <- gmt_files[1]
  #genes_in_data <- df$gene_symbol
  
  # Read in gmt file
  gmt <- gmtPathways(gmt_file)
  hidden <- unique(unlist(gmt))
  
  # Convert gmt file to a matrix with the genes as rows and for each go annotation (columns) the values are 0 or 1
  mat <- matrix(NA, dimnames = list(hidden, names(gmt)),
                nrow = length(hidden), ncol = length(gmt))
  for (i in 1:dim(mat)[2]){
    mat[,i] <- as.numeric(hidden %in% gmt[[i]])
  }
  
  #Subset to the genes that are present in our data to avoid bias
  hidden1 <- intersect(genes_in_data, hidden)
  mat <- mat[hidden1, colnames(mat)[which(colSums(mat[hidden1,])>5)]] # filter for gene sets with more than 5 genes annotated
  # And get the list again
  final_list <- matrix_to_list(mat) # for this we use the function we previously defined
  
  if(savefile){
    saveRDS(final_list, file = paste0(gsub('.gmt', '', gmt_file), '_subset_', format(Sys.time(), '%d%m'), '.RData'))
  }
  
  print('Wohoo! .gmt conversion successfull!:)')
  return(final_list)
}

# List of acronyms to uppercase
acronyms <- c("rna", "dna", "mtor", "foxo", "ppar", "nmd", "fgfr", "robo", 
              "bhl", "cov", "jak", "stat", "wnt", "hiv", "bcl", "mapk",
              "pt", "tal", "pc", "ic", "ec", "fibvsmcp")
special_mixed <- c("rrna", "mrna", "trna", "gtpase", "atpase", "robos", "slits", "fibvsmcp")
special_replacements <- c("rRNA", "mRNA", "tRNA", "GTPase", "ATPase", "ROBOs", "SLITs", "FIB/VSMC/P")


# ===========================================================================
# Function: replace_mixed_case
# ===========================================================================

replace_mixed_case <- function(text, from, to) {
  for (i in seq_along(from)) {
    pattern <- paste0("\\b", from[i], "\\b")
    text <- str_replace_all(text, regex(pattern, ignore_case = TRUE), to[i])
  }
  return(text)
}


# ===========================================================================
# Function: capitalize_acronyms
# ===========================================================================

capitalize_acronyms <- function(text, terms) {
  for (term in terms) {
    pattern <- paste0("\\b", term, "\\b")
    replacement <- toupper(term)
    text <- str_replace_all(text, regex(pattern, ignore_case = TRUE), replacement)
  }
  return(text)
}


# ===========================================================================
# Function: filter_redundant_pathways
# ===========================================================================

filter_redundant_pathways <- function(gsea_result, overlap_pct = 0.3) {
  # Check required columns
  if (!all(c("pathway", "leadingEdge", "padj") %in% colnames(gsea_result))) {
    stop("Input data frame must have 'pathway', 'leadingEdge', and 'padj' columns.")
  }
  
  # Extract leading edge and name them
  leading_edges <- gsea_result$leadingEdge
  names(leading_edges) <- gsea_result$pathway
  
  # Compute Jaccard overlap matrix
  overlap_matrix <- sapply(leading_edges, function(x) 
    sapply(leading_edges, function(y) 
      length(intersect(x, y)) / length(union(x, y))
    )
  )
  
  # Identify pairs with high overlap
  overlap_pairs <- which(overlap_matrix > overlap_pct & lower.tri(overlap_matrix), arr.ind = TRUE)
  
  # If no overlaps found, return original
  if (nrow(overlap_pairs) == 0) {
    message("No redundant pathways found with overlap above threshold.")
    return(gsea_result)
  }
  
  # Create overlap graph
  library(igraph)
  edges <- data.frame(
    from = rownames(overlap_matrix)[overlap_pairs[,1]],
    to   = colnames(overlap_matrix)[overlap_pairs[,2]]
  )
  g <- graph_from_data_frame(edges, directed = FALSE)
  components <- components(g)
  
  # Cluster terms and select representative from each
  redundant_clusters <- split(names(components$membership), components$membership)
  representative_terms <- sapply(redundant_clusters, function(cluster_terms) {
    sub <- gsea_result[gsea_result$pathway %in% cluster_terms, ]
    sub[which.min(sub$padj), "pathway"]
  })
  
  # Keep only representative + non-overlapping terms
  all_overlapping <- unique(unlist(redundant_clusters))
  non_overlapping <- setdiff(gsea_result$pathway, all_overlapping)
  final_terms <- c(non_overlapping, representative_terms)
  
  # Return filtered GSEA result
  gsea_result[gsea_result$pathway %in% final_terms, ]
}

# ===========================================================================
# Function: plot_fgsea_transpose
# ===========================================================================

plot_fgsea_transpose <- function(fgsea_res,
                                 top_n = 30,
                                 title = "Top Enriched Pathways",
                                 xmin = 0,
                                 xmax = 3,
                                 xnudge = (xmax - xmin)/100,
                                 text1 = 6.5,
                                 text2 = 18,
                                 text3 = 20) {
  
  fgsea_res <- fgsea_res %>%
    arrange(pval) %>%
    head(top_n) %>%
    mutate(
      direction = case_when((NES < 0 & pval <= 0.05 ~ "Negative"), 
                            (NES > 0 & pval <= 0.05 ~ "Positive"),
                            (NES < 0 & pval > 0.05 ~ "Negative p > 0.05"), 
                            (NES > 0 & pval > 0.05 ~ "Positive p > 0.05")),
      face = case_when((NES < 0 & pval <= 0.05 ~ "bold"), 
                       (NES > 0 & pval <= 0.05 ~ "bold"),
                       (NES < 0 & pval > 0.05 ~ "plain"), 
                       (NES > 0 & pval > 0.05 ~ "plain")),
      pathway_clean = str_remove(pathway, "^KEGG_"), 
      pathway_clean = str_remove(pathway_clean, "^REACTOME_"), 
      pathway_clean = str_remove(pathway_clean, "^GOBP_"), 
      pathway_clean = str_remove(pathway_clean, "^GOMF_"), 
      pathway_clean = str_replace_all(pathway_clean, "_", " "),
      pathway_clean = str_to_sentence(pathway_clean),
      pathway_clean = str_replace_all(pathway_clean, "\\bi\\b", "I"),
      pathway_clean = str_replace_all(pathway_clean, "\\bii\\b", "II"),
      pathway_clean = str_replace_all(pathway_clean, "\\biii\\b", "III"),
      pathway_clean = str_replace_all(pathway_clean, "\\biv\\b", "IV"),
      pathway_clean = str_replace_all(pathway_clean, "\\bv\\b", "V"),
      pathway_clean = str_replace_all(pathway_clean, regex("\\(immune\\)", ignore_case = TRUE), "(IMMUNE)"),
      pathway_clean = capitalize_acronyms(pathway_clean, acronyms),
      pathway_clean = replace_mixed_case(pathway_clean, special_mixed, special_replacements),
      pathway_clean = paste0(pathway_clean, " (", size, ")")
    ) %>%
    arrange(pval)
  
  fgsea_res$pathway_clean <- reorder(fgsea_res$pathway_clean, -abs(fgsea_res$NES))
  
  fgsea_res %>%
    ggplot(aes(x = abs(NES), y = fct_rev(pathway_clean), label = pathway_clean)) +
    geom_point(aes(size = -log10(pval), color = direction, alpha = 0.8)) +
    # geom_vline(xintercept = 2, linetype = "dashed") +
    geom_text(aes(group = pathway_clean, color = direction, fontface = face), 
              hjust = 0, size = text1, nudge_x = xnudge) +
    scale_size_binned() +
    scale_color_manual(values = c("Positive" = "#c75146", "Negative" = "#2c7da0", 
                                  "Positive p > 0.05" = "#e18c80", "Negative p > 0.05" = "#7ab6d1")) +
    scale_x_continuous(limits = c(xmin, xmax), expand = expansion(mult = c(0, 0))) +
    scale_y_discrete(expand = expansion(add = 1)) +
    labs(
      x = "NES",
      y = "Pathways",
      color = "Direction",
      size = "-log(p-value)",
      title = title
    ) +
    guides(alpha = "none") +
    theme_bw() +
    theme(
      axis.text.y = element_blank(),
      panel.grid = element_blank(),
      axis.text.x = element_text(size = text3),
      axis.title = element_text(size = text3),
      axis.ticks.y = element_blank(), 
      legend.position = c(0.9, 0.2),
      legend.background = element_blank(),
      legend.box.background = element_rect(color = "black"),
      legend.title = element_text(size = text2),
      legend.text = element_text(size = text2),
      title = element_text(size = text3)
    )
}

# ---- Summary Functions ----

# ===========================================================================
# Function: run_cell_type_analysis
# ===========================================================================

# Function to run complete analysis for a given cell type
run_cell_type_analysis <- function(cell_type, 
                                   input_path, 
                                   input_suffix,
                                   output_base_path,
                                   bg_path,
                                   bucket = "attempt",
                                   region = "") {
  
  # Print status
  cat("Starting analysis for cell type:", cell_type, "\n")
  
  # Create cell type specific paths
  cell_type_lower <- tolower(cell_type)
  
  # Read in nebula results
  input_file <- file.path(input_path, cell_type, "nebula", 
                          paste0(cell_type_lower, input_suffix))
  
  croc_res <- s3readRDS(input_file, bucket = bucket, region = region)
  
  # Check convergence
  croc_res_convergence <- map_dfr(
    names(croc_res),
    function(gene_name) {
      converged <- croc_res[[gene_name]]$convergence
      df <- data.frame(Gene = gene_name,
                       Convergence_Code = converged)
      return(df)
    }
  )
  
  croc_res_converged <- croc_res_convergence %>%
    filter(Convergence_Code >= -10)
  
  # Combine results
  croc_res_combined <- map_dfr(
    names(croc_res),
    function(gene_name) {
      df <- croc_res[[gene_name]]$summary
      df <- df %>% 
        filter(gene_name %in% croc_res_converged$Gene) %>%
        mutate(Gene = gene_name)
      return(df)
    }
  )
  
  # Calculate FDR
  croc_res_combined <- croc_res_combined %>%
    ungroup() %>%
    mutate(fdr_interaction = p.adjust(p_groupType_1_Diabetes, method = "fdr"))
  
  # Save results
  output_csv <- file.path(output_base_path, "Results", "NEBULA", 
                          paste0("CROC_", cell_type_lower, "_nebula_res.csv"))
  write.csv(croc_res_combined, output_csv, row.names = FALSE)
  
  # Add direction columns for visualization
  croc_res_combined <- croc_res_combined %>%
    mutate(group_direction = case_when(
      p_groupType_1_Diabetes < 0.05 & logFC_groupType_1_Diabetes > 0 ~ "Positive",
      p_groupType_1_Diabetes < 0.05 & logFC_groupType_1_Diabetes < 0 ~ "Negative",
      TRUE ~ "NS"
    ),
    group_direction_fdr = case_when(
      fdr_interaction < 0.05 & logFC_groupType_1_Diabetes > 0 ~ "Positive",
      fdr_interaction < 0.05 & logFC_groupType_1_Diabetes < 0 ~ "Negative",
      TRUE ~ "NS"
    ))
  
  # Create volcano plots
  plot_volcano(croc_res_combined, 
               "logFC_groupType_1_Diabetes", 
               "p_groupType_1_Diabetes",
               paste(cell_type, "Volcano Plot"), 
               "logFC group", 
               "-log10(p-value)",
               file.path("nebula", paste0("croc_", cell_type_lower, "_deg")),
               formula = "group",
               positive_text = "Upregulated in T1D compared to LC",
               negative_text = "Downregulated in T1D compared to LC")
  
  plot_volcano(croc_res_combined, 
               "logFC_groupType_1_Diabetes", 
               "fdr_interaction",
               paste(cell_type, "Volcano Plot (FDR)"), 
               "logFC group", 
               "-log10(FDR adjusted p-value)",
               file.path("nebula", paste0("croc_", cell_type_lower, "_deg_fdr")),
               formula = "group",
               positive_text = "Upregulated in T1D compared to LC",
               negative_text = "Downregulated in T1D compared to LC")
  
  # GSEA Analysis
  # Load pathway files
  gmt_files <- list.files(path = bg_path, pattern = '.gmt', full.names = TRUE)
  kegg_legacy <- prepare_gmt(gmt_files[1], unique(croc_res_combined$Gene), savefile = FALSE)
  reactome <- prepare_gmt(gmt_files[3], unique(croc_res_combined$Gene), savefile = FALSE)
  go <- prepare_gmt(gmt_files[4], unique(croc_res_combined$Gene), savefile = FALSE)
  
  # Rank genes by logFC
  rankings <- croc_res_combined$logFC_groupType_1_Diabetes
  names(rankings) <- croc_res_combined$Gene
  rankings <- sort(rankings, decreasing = TRUE)
  
  # Run GSEA
  set.seed(1234)
  
  kegg_legacy_res <- fgsea(pathways = kegg_legacy,
                           stats = rankings,
                           scoreType = 'std', 
                           minSize = 3,
                           maxSize = 500,
                           nproc = 1)
  
  reactome_res <- fgsea(pathways = reactome,
                        stats = rankings,
                        scoreType = 'std', 
                        minSize = 3,
                        maxSize = 500,
                        nproc = 1)
  
  go_res <- fgsea(pathways = go,
                  stats = rankings,
                  scoreType = 'std', 
                  minSize = 5,
                  maxSize = 500,
                  nproc = 1)
  
  # Summary table
  fgsea_summary <- data.frame(
    "KEGG_Legacy" = c(sum(kegg_legacy_res[, padj < 0.05], na.rm = TRUE), 
                      sum(kegg_legacy_res[, pval < 0.05], na.rm = TRUE)),
    "REACTOME" = c(sum(reactome_res[, padj < 0.05], na.rm = TRUE), 
                   sum(reactome_res[, pval < 0.05], na.rm = TRUE)),
    "GO" = c(sum(go_res[, padj < 0.05], na.rm = TRUE), 
             sum(go_res[, pval < 0.05], na.rm = TRUE))
  )
  rownames(fgsea_summary) <- c("adj.pval", "p.val")
  
  # Plot pathways
  # KEGG
  plot_fgsea_transpose(kegg_legacy_res, 
                       title = paste(cell_type, "Top 30 KEGG Pathways"), 
                       xmin = 1, xmax = 3)
  
  ggsave(file.path(output_base_path, "Results", "Figures", "Pathways", "nebula",
                   paste0("croc_", cell_type_lower, "_res_top30_kegg_pathways.jpeg")),
         width = 27.5, height = 14, scale = 1)
  
  # REACTOME
  plot_fgsea_transpose(reactome_res, 
                       title = paste(cell_type, "Top 30 REACTOME Pathways"), 
                       xmin = 1.45, xmax = 2.6)
  
  ggsave(file.path(output_base_path, "Results", "Figures", "Pathways", "nebula",
                   paste0("croc_", cell_type_lower, "_res_top30_reactome_pathways.jpeg")),
         width = 27.5, height = 14, scale = 1)
  
  # GO
  plot_fgsea_transpose(go_res, 
                       title = paste(cell_type, "Top 30 GO Pathways"), 
                       xmin = 1.4, xmax = 5)
  
  ggsave(file.path(output_base_path, "Results", "Figures", "Pathways", "nebula",
                   paste0("croc_", cell_type_lower, "_res_top30_go_pathways.jpeg")),
         width = 27.5, height = 14, scale = 1)
  
  # Return results
  return(list(
    nebula_results = croc_res_combined,
    fgsea_summary = fgsea_summary,
    kegg_results = kegg_legacy_res,
    reactome_results = reactome_res,
    go_results = go_res,
    rankings = rankings
  ))
}
