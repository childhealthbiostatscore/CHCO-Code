
# ATTEMPT analysis related functions

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
                                           TRUE ~ "NS"))
  
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
        TRUE ~ "NS"
      ),
      # Shape based on significance and direction
      shape_var_plot = case_when(
        Gene %in% sig_both_genes & !!sym(fc) < 0 & treatment_results$treatment_direction[match(Gene, treatment_results$Gene)] == "Up w/ Dapagliflozin" ~ "Left_Pos",
        Gene %in% sig_both_genes & !!sym(fc) < 0 & treatment_results$treatment_direction[match(Gene, treatment_results$Gene)] == "Down w/ Dapagliflozin" ~ "Left_Neg",
        Gene %in% sig_both_genes & !!sym(fc) > 0 & treatment_results$treatment_direction[match(Gene, treatment_results$Gene)] == "Up w/ Dapagliflozin" ~ "Right_Pos",
        Gene %in% sig_both_genes & !!sym(fc) > 0 & treatment_results$treatment_direction[match(Gene, treatment_results$Gene)] == "Down w/ Dapagliflozin" ~ "Right_Neg",
        TRUE ~ "NS"
      ),
      # Label only top concordant genes
      top_lab = if_else(Gene %in% top_label_genes, Gene, "")
    )
  
  clin_results$shape_var_plot <- factor(clin_results$shape_var_plot, 
                                        levels = c("Left_Pos", "Left_Neg", "Right_Pos", "Right_Neg", "NS"))
  
  clin_results$color_var_plot <- factor(clin_results$color_var_plot, 
                                        levels = c("Up w/ Dapagliflozin", "Down w/ Dapagliflozin", "NS"))
  
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
    scale_shape_manual(values = c("Left_Pos" = 22, "Left_Neg" = 22, "Right_Pos" = 22, "Right_Neg" = 22, "NS" = 22),
                       labels = c("Left_Pos" = "Left Pos", "Left_Neg" = "Left Neg", 
                                  "Right_Pos" = "Right Pos", "Right_Neg" = "Right Neg", "NS" = "NS")) +
    scale_color_manual(values = c("Up w/ Dapagliflozin" = "#f28482", "Down w/ Dapagliflozin" = "#457b9d", "NS" = "#ced4da")) +
    scale_fill_manual(values = c("Up w/ Dapagliflozin" = "#f28482", "Down w/ Dapagliflozin" = "#457b9d", "NS" = "#ced4da")) +
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
             y=-y_max * 0.08,
             col="darkgrey", arrow=arrow(length=unit(0.2, "cm"))) +
    annotate("text", 
             x=mean(c(max_fc/8, (max_fc*7)/8)), 
             y=-y_max * 0.14, 
             label=positive_text,
             size=3, color="#343a40") +
    annotate("segment", 
             x=min_fc/8, 
             xend=(min_fc*7)/8, 
             y=-y_max * 0.08,
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
  ggsave(paste0("/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/ATTEMPT/Results/Figures/Volcano Plots/", file_suffix, ".jpeg"), plot = p, width = 7, height = 5)
  
  return(p)
}
