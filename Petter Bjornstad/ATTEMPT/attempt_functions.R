

color_5 <- c("#264653", "#2a9d8f", "#e9c46a", "#f4a261", "#e76f51")

# ATTEMPT analysis related functions

# ---- Analysis Functions ----


#----------------------------------------------------------#
# run_doubletfinder_custom
#----------------------------------------------------------#
# run_doubletfinder_custom runs Doublet_Finder() and returns a dataframe with the cell IDs and a column with either 'Singlet' or 'Doublet'
run_doubletfinder_custom <- function(seu_sample_subset, multiplet_rate = NULL){
  # for debug
  #seu_sample_subset <- samp_split[[1]]
  # Print sample number
  print(paste0("Sample ", unique(seu_sample_subset[['SampleID']]), '...........')) 
  
  if(is.null(multiplet_rate)){
    print('multiplet_rate not provided....... estimating multiplet rate from cells in dataset')
    
    # 10X multiplet rates table
    #https://rpubs.com/kenneditodd/doublet_finder_example
    multiplet_rates_10x <- data.frame('Multiplet_rate'= c(0.004, 0.008, 0.0160, 0.023, 0.031, 0.039, 0.046, 0.054, 0.061, 0.069, 0.076),
                                      'Loaded_cells' = c(800, 1600, 3200, 4800, 6400, 8000, 9600, 11200, 12800, 14400, 16000),
                                      'Recovered_cells' = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000))
    
    print(multiplet_rates_10x)
    
    multiplet_rate <- multiplet_rates_10x %>% dplyr::filter(Recovered_cells < nrow(seu_sample_subset@meta.data)) %>% 
      dplyr::slice(which.max(Recovered_cells)) %>% # select the min threshold depending on your number of samples
      dplyr::select(Multiplet_rate) %>% as.numeric(as.character()) # get the expected multiplet rate for that number of recovered cells
    
    print(paste('Setting multiplet rate to', multiplet_rate))
  }
  
  # Pre-process seurat object with standard seurat workflow --- 
  sample <- NormalizeData(seu_sample_subset)
  sample <- FindVariableFeatures(sample)
  sample <- ScaleData(sample)
  sample <- RunPCA(sample, nfeatures.print = 10)
  
  # Find significant PCs
  stdv <- sample[["pca"]]@stdev
  percent_stdv <- (stdv/sum(stdv)) * 100
  cumulative <- cumsum(percent_stdv)
  co1 <- which(cumulative > 90 & percent_stdv < 5)[1] 
  co2 <- sort(which((percent_stdv[1:length(percent_stdv) - 1] - 
                       percent_stdv[2:length(percent_stdv)]) > 0.1), 
              decreasing = T)[1] + 1
  min_pc <- min(co1, co2)
  
  # Finish pre-processing with min_pc
  sample <- RunUMAP(sample, dims = 1:min_pc)
  sample <- FindNeighbors(object = sample, dims = 1:min_pc)              
  sample <- FindClusters(object = sample, resolution = 0.1)
  
  # pK identification (no ground-truth) 
  #introduces artificial doublets in varying props, merges with real data set and 
  # preprocesses the data + calculates the prop of artficial neighrest neighbours, 
  # provides a list of the proportion of artificial nearest neighbours for varying
  # combinations of the pN and pK
  sweep_list <- paramSweep(sample, PCs = 1:min_pc, sct = FALSE)   
  sweep_stats <- summarizeSweep(sweep_list)
  bcmvn <- find.pK(sweep_stats) # computes a metric to find the optimal pK value (max mean variance normalised by modality coefficient)
  # Optimal pK is the max of the bimodality coefficient (BCmvn) distribution
  optimal.pk <- bcmvn %>% 
    dplyr::filter(BCmetric == max(BCmetric)) %>%
    dplyr::select(pK)
  optimal.pk <- as.numeric(as.character(optimal.pk[[1]]))
  
  ## Homotypic doublet proportion estimate
  annotations <- sample@meta.data$seurat_clusters # use the clusters as the user-defined cell types
  homotypic.prop <- modelHomotypic(annotations) # get proportions of homotypic doublets
  
  nExp.poi <- round(multiplet_rate * nrow(sample@meta.data)) # multiply by number of cells to get the number of expected multiplets
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop)) # expected number of doublets
  
  # run DoubletFinder
  sample <- doubletFinder(seu = sample, 
                          PCs = 1:min_pc, 
                          pK = optimal.pk, # the neighborhood size used to compute the number of artificial nearest neighbours
                          nExp = nExp.poi.adj) # number of expected real doublets
  # change name of metadata column with Singlet/Doublet information
  colnames(sample@meta.data)[grepl('DF.classifications.*', colnames(sample@meta.data))] <- "doublet_finder"
  
  # Subset and save
  # head(sample@meta.data['doublet_finder'])
  # singlets <- subset(sample, doublet_finder == "Singlet") # extract only singlets
  # singlets$ident
  double_finder_res <- sample@meta.data['doublet_finder'] # get the metadata column with singlet, doublet info
  double_finder_res <- rownames_to_column(double_finder_res, "row_names") # add the cell IDs as new column to be able to merge correctly
  return(double_finder_res)
}



# ===========================================================================
# Function: process_nebula_results
# ===========================================================================

process_nebula_results <- function(nebula_list, 
                                   pval_col = "p_treatmentDapagliflozin:visitPOST", 
                                   convergence_cut = -10) {
  # Extract convergence codes
  convergence_df <- purrr::map_dfr(names(nebula_list), function(gene_name) {
    convergence_code <- nebula_list[[gene_name]]$convergence
    data.frame(Gene = gene_name, Convergence_Code = convergence_code)
  })
  
  # Filter to converged models
  converged_genes <- convergence_df %>%
    filter(Convergence_Code >= convergence_cut) %>%
    pull(Gene)
  
  # Combine model summary results
  summary_df <- purrr::map_dfr(converged_genes, function(gene_name) {
    nebula_list[[gene_name]]$summary %>%
      mutate(Gene = gene_name)
  })
  
  # Add FDR adjustment
  if (pval_col %in% names(summary_df)) {
    summary_df <- summary_df %>%
      mutate(fdr = p.adjust(.data[[pval_col]], method = "fdr"))
  } else {
    warning(paste("Column", pval_col, "not found in summary data. FDR not computed."))
    summary_df$fdr <- NA
  }
  
  # Extract overdispersion estimates
  overdisp_df <- purrr::map_dfr(names(nebula_list), function(gene_name) {
    od <- nebula_list[[gene_name]]$overdispersion
    od$Gene <- gene_name
    od
  })
  
  return(list(
    convergence     = convergence_df,
    results         = summary_df,
    overdispersion  = overdisp_df
  ))
}


# ===========================================================================
# Function: run_nebula_attempt
# ===========================================================================

run_nebula_attempt <- function(so,
                               trait            = "",
                               extra_covars     = "treatment",
                               subject_var      = "subject_id",
                               offset_var       = "pooled_offset",
                               assay_layer      = "counts",
                               n_cores          = max(parallel::detectCores() - 1, 1),
                               aws_s3           = NULL,
                               s3_bucket        = NULL,
                               s3_key           = NULL) {
  
  stopifnot(trait %in% colnames(so@meta.data))
  
  # ── 1. Keep only cells with non-missing trait ──────────────────────────────
  md            <- so@meta.data
  keep_cells    <- rownames(md)[!is.na(md[[trait]])]
  so_subset     <- so[, keep_cells]
  
  # ── 2. Pull counts and gene list ──────────────────────────────────────────
  counts_mat    <- round(GetAssayData(so_subset, layer = assay_layer))
  genes_list    <- rownames(counts_mat)
  
  # ── 3. Spin up parallel backend ───────────────────────────────────────────
  cl            <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  
  on.exit({          # make *sure* we clean up
    parallel::stopCluster(cl)
  }, add = TRUE)
  
  start_time <- Sys.time()
  
  # ── 4. Per-gene nebula fits ───────────────────────────────────────────────
  nebula_res <- foreach::foreach(
    g = genes_list,
    .packages      = c("nebula", "Matrix"),
    .errorhandling = "pass"
  ) %dopar% {
    
    warn <- err <- NULL
    res  <- NULL
    
    tryCatch({
      count_gene <- counts_mat[g, , drop = FALSE]
      meta_gene  <- subset(so_subset, features = g)@meta.data
      
      pred_formula <- reformulate(c(trait, extra_covars), response = NULL) # ~ trait  extras
      pred_gene    <- model.matrix(pred_formula, data = meta_gene)
      
      data_g       <- list(count = count_gene,
                           id    = meta_gene[[subject_var]],
                           pred  = pred_gene)
      
      res <- withCallingHandlers(
        nebula::nebula(count      = data_g$count,
                       id         = data_g$id,
                       pred       = data_g$pred,
                       model      = "NBLMM",
                       output_re  = TRUE,
                       covariance = TRUE,
                       reml       = TRUE,
                       offset     = if (!is.null(offset_var)) meta_gene[[offset_var]] else NULL,
                       ncore      = 1),
        warning = function(w) { warn <<- conditionMessage(w); invokeRestart("muffleWarning") }
      )
      
    }, error = function(e) {
      err <<- conditionMessage(e)
    })
    
    list(gene = g, result = res, warning = warn, error = err)
  }
  
  # ── 5. Collate warnings / errors ──────────────────────────────────────────
  for (x in nebula_res) {
    if (!is.null(x$warning)) message(sprintf("⚠️  Warning for %s: %s", x$gene, x$warning))
    if (!is.null(x$error))   message(sprintf("⛔ Error   for %s: %s", x$gene, x$error))
  }
  
  # ── 6. Keep successful fits only & name the list ──────────────────────────
  fits <- lapply(nebula_res, `[[`, "result")
  names(fits) <- vapply(nebula_res, `[[`, "", "gene")
  fits        <- Filter(Negate(is.null), fits)
  
  # ── 7. Report & (optionally) persist to S3 ────────────────────────────────
  drop_pct <- 100 * (1 - length(fits) / length(genes_list))
  message(sprintf("%0.2f%% of genes were dropped (low expression / errors).", drop_pct))
  
  if (!is.null(aws_s3) && !is.null(s3_bucket) && !is.null(s3_key)) {
    tmp <- tempfile(fileext = ".rds")
    saveRDS(fits, tmp)
    aws_s3$upload_file(tmp, Bucket = s3_bucket, Key = s3_key)
    unlink(tmp)
    message(sprintf("⬆️  Results uploaded to s3://%s/%s", s3_bucket, s3_key))
  }
  
  end_time <- Sys.time()
  message(sprintf("Finished in %.1f minutes.", as.numeric(difftime(end_time, start_time, units = "mins"))))
  
  invisible(fits)
}

# ---- Plot Functions ----
# ===========================================================================
# Function: make_comp_plot
# ===========================================================================

make_comp_plot <- function(attempt_df,
                           croc_df,
                           attempt_p_cut = 0.05,
                           croc_p_cut   = 0.05,
                           save_path    = NULL,
                           csv_path     = NULL,  
                           width        = 10,
                           height       = 5,
                           caption      = NULL,
                           FC = "logFC_treatmentDapagliflozin:visitPOST") {
  
  ## 1. Harmonize & flag significance ------------------------------------------
  df_attempt <- attempt_df %>%
    dplyr::rename(p_fdr = fdr,
                  FC    = FC) %>%
    mutate(direction = case_when(p_fdr < attempt_p_cut & FC < 0 ~ "-",
                                 p_fdr < attempt_p_cut & FC > 0 ~ "+"),
           source = "ATTEMPT") %>%
    dplyr::select(Gene, direction, p_fdr, source, FC)
  
  df_croc <- croc_df %>%
    dplyr::rename(p_fdr = p_groupType_1_Diabetes,
                  FC    = logFC_groupType_1_Diabetes) %>%
    mutate(direction = case_when(p_fdr < croc_p_cut & FC < 0 ~ "-",
                                 p_fdr < croc_p_cut & FC > 0 ~ "+"),
           source = "CROCODILE") %>%
    dplyr::select(Gene, direction, p_fdr, source, FC)
  
  pt_comp_df <- bind_rows(df_attempt, df_croc) %>%
    filter(!is.na(direction))
  
  ## 2. Wide table → effect category ------------------------------------------
  pt_comp_df_wide <- pt_comp_df %>%
    dplyr::select(Gene, source, direction) %>%
    pivot_wider(names_from  = source,
                values_from = direction) %>%
    mutate(effect = case_when(
      (ATTEMPT == "+" & CROCODILE == "+") | 
        (ATTEMPT == "-" & CROCODILE == "-") ~ "Inconsistent",
      ATTEMPT == "+" & CROCODILE == "-" ~ "Reversed towards +\n(T1D depleted)",
      ATTEMPT == "-" & CROCODILE == "+" ~ "Reversed towards -\n(T1D elevated)"
    )) %>%
    mutate(effect = factor(effect,
                           levels = c("Inconsistent",
                                      "Reversed towards +\n(T1D depleted)",
                                      "Reversed towards -\n(T1D elevated)")))
  
  ## Count inconsistent genes for caption
  n_inconsistent <- sum(pt_comp_df_wide$effect == "Inconsistent", na.rm = TRUE)
  
  ## Save inconsistent genes to CSV if requested
  if (!is.null(csv_path)) {
    inconsistent_genes <- pt_comp_df_wide %>%
      filter(effect == "Inconsistent") %>%
      mutate(consistency_type = case_when(
        ATTEMPT == "+" & CROCODILE == "+" ~ "Both_upregulated",
        ATTEMPT == "-" & CROCODILE == "-" ~ "Both_downregulated"
      )) %>%
      arrange(Gene)
    
    write.csv(inconsistent_genes, csv_path, row.names = FALSE)
  }
  
  ## 3. Bar plot ---------------------------------------------------------------
  effect_cols <- c("Reversed towards +\n(T1D depleted)" = "#3a5a40",
                   "Reversed towards -\n(T1D elevated)" = "#a3b18a")
  
  # Update caption to include inconsistent count
  if (is.null(caption)) {
    caption <- paste0("Number of inconsistent genes: ", n_inconsistent)
  } else {
    caption <- paste0(caption, "\nNumber of inconsistent genes: ", n_inconsistent)
  }
  
  pt_comp_bar <- pt_comp_df_wide %>%
    filter(!is.na(effect) & effect != "Inconsistent") %>%  # Filter out inconsistent genes
    ggplot(aes(x = effect, fill = effect)) +
    geom_bar_rounded() +
    geom_text(stat  = "count",
              aes(label = after_stat(count)),
              vjust = 1.5, hjust = 0.5, color = "white") +
    theme_bw() +
    theme(panel.grid  = element_blank(),
          panel.border = element_blank(),
          text         = element_text(size = 15),
          legend.position = "none",
          axis.text.x  = element_text(angle = 50, hjust = 1),
          axis.ticks   = element_blank(),
          plot.caption = element_text(hjust = 0.5, size = 12)) +
    labs(y = "Count", x = NULL,
         caption = caption) +
    scale_fill_manual(values = effect_cols)
  
  ## 4. Label colouring for reversed genes ------------------------------------
  reversed_cols <- effect_cols[c(
    "Reversed towards +\n(T1D depleted)",
    "Reversed towards -\n(T1D elevated)"
  )]
  
  pt_gene_color_df <- pt_comp_df_wide %>%
    filter(str_detect(effect, "Reversed")) %>%
    mutate(Gene_label = ifelse(
      str_detect(effect, "depleted"),
      paste0("<span style='color:", reversed_cols[1], "'>", Gene, "</span>"),
      paste0("<span style='color:", reversed_cols[2], "'>", Gene, "</span>")
    )) %>%
    dplyr::select(Gene, Gene_label)
  
  ## 5. Dot plot ---------------------------------------------------------------
  pt_comp_dot <- pt_comp_df %>%
    left_join(pt_gene_color_df, by = "Gene") %>%
    filter(Gene %in% pt_gene_color_df$Gene) %>%
    mutate(group = ifelse(source == "ATTEMPT",
                          "Dapa vs. Placebo", "T1D vs. HC"),
           Gene_label = dplyr::coalesce(Gene_label, Gene)) %>%
    ggplot(aes(x = Gene_label, y = group,
               color = direction, size = abs(FC))) +
    geom_point() +
    theme_bw() +
    scale_color_manual(values = c("+" = "#f28482", "-" = "#457b9d")) +
    theme(panel.border = element_blank(),
          panel.grid   = element_blank(),
          legend.position = "top",
          legend.box      = "vertical",
          axis.text.x  = element_markdown(angle = 50, hjust = 1),
          axis.ticks   = element_blank(),
          text = element_text(size = 15)) +
    scale_y_discrete(position = "right") +
    labs(x = NULL, y = NULL, color = "Direction")
  
  ## 6. Combine & optionally save ---------------------------------------------
  combined <- pt_comp_bar + pt_comp_dot
  
  if (!is.null(save_path)) {
    ggsave(save_path, combined, width = width, height = height)
  }
  
  return(invisible(combined))
}
# ===========================================================================
# Function: make_plot_df
# ===========================================================================

make_plot_df <- function(res_df, celltype_label, pos_genes, neg_genes) {
  res_df %>%
    dplyr::select(Gene, 
                  logFC = `logFC_treatmentDapagliflozin:visitPOST`, 
                  pval = `p_treatmentDapagliflozin:visitPOST`, 
                  fdr = fdr) %>%
    mutate(celltype = celltype_label,
           direction = case_when(
             logFC > 0 ~ "+",
             logFC < 0 ~ "-"
           )) %>%
    filter(Gene %in% c(pos_genes, neg_genes))
}

# ===========================================================================
# Function: plot_volcano
# ===========================================================================
plot_volcano <- function(data, fc, p_col, title = NULL, x_axis, y_axis, file_suffix, p_thresh = 0.05,
                         positive_text = "Positive with Dapagliflozin", 
                         negative_text = "Negative with Dapagliflozin",
                         formula = "group", legend_position = c(0.8, 0.9),
                         cell_type = "",
                         output_base_path = "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/ATTEMPT/Results/Figures/Volcano Plots/") {
  set.seed(1)
  top_pos <- data %>%
    dplyr::filter(!!sym(fc) > 0 & !!sym(p_col) < p_thresh) %>%
    dplyr::arrange(!!sym(p_col))
  
  n_pos <- nrow(top_pos)
    
  top_pos <- top_pos %>%
    slice_head(n=20)
  
  top_neg <- data %>%
    dplyr::filter(!!sym(fc) < 0 & !!sym(p_col) < p_thresh) %>%
    dplyr::arrange(!!sym(p_col))
  
  n_neg <- nrow(top_neg)
    
  top_neg <- top_neg %>%
    slice_head(n=20)
  
  data <- data %>%
    dplyr::mutate(top_color = case_when(Gene %in% top_pos$Gene ~ "#f28482",
                                        Gene %in% top_neg$Gene ~ "#457b9d",
                                        TRUE ~ "#ced4da"),
                  top_size = if_else(Gene %in% c(top_pos$Gene, top_neg$Gene), 1.3, 1),
                  top_lab  = if_else(Gene %in% c(top_pos$Gene, top_neg$Gene), Gene, ""))
  
  # Max and min for annotation arrows
  max_fc <- max(data[[fc]], na.rm = TRUE)
  min_fc <- min(data[[fc]], na.rm = TRUE)
  
  # Get y-axis max for dynamic scaling
  y_max <- max(-log10(data[[p_col]]), na.rm = TRUE) * 1.1
  
  
  p <- ggplot(data, aes(x = !!sym(fc), y = -log10(!!sym(p_col)))) +
    geom_hline(yintercept = -log10(p_thresh), linetype = "dashed", color = "darkgrey") +
    geom_point(alpha = 0.5, aes(color = top_color, size = top_size)) +
    geom_text_repel(aes(label = top_lab, color = top_color),
                    size = 3, max.overlaps = Inf,
                    force = 6, segment.alpha = 0.3, segment.size = 0.3) +
    labs(title = paste(title),
         x = paste(x_axis),
         y = paste(y_axis),
         caption = if (!is.null(cell_type)) {
           paste0("Formula: ~ ", formula, " + (1|subject)", 
                  "\n\n Cell type: ", 
                  cell_type, if(cell_type != "") " | " else "", 
                  "Positive n = ", n_pos, " | Negative n = ", n_neg)
         } else NULL) +
    scale_size_continuous(range = c(1, 1.3)) + 
    scale_color_manual(values = c("#457b9d"="#457b9d", "#ced4da"="#ced4da", "#f28482"="#f28482")) +
    theme_minimal() +
    theme(legend.title = element_blank(),
          panel.grid = element_blank(),
          text = element_text(size = 15),
          title = element_text(size = 9)) +
    guides(color = "none", size = "none")  +
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
          axis.title.x = element_text(margin = margin(t = 32)),
          plot.caption = element_text(size = 8.5, hjust = 0.5, margin = margin(t = 8), family = "Arial"),
          legend.margin = margin(t = 5, b = 5))
    
  
  ggsave(paste0(output_base_path, file_suffix, ".jpeg"), plot = p, width = 7, height = 5)
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
  
  # Get y-axis max for dynamic scaling
  y_max <- max(-log10(clin_results[[p_col]]), na.rm = TRUE) * 1.1
  
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
                                     arrow_label = "",
                                     formula = "\u0394 treatment", color_by = "treatment_direction",
                                     legend_position = "top",
                                     clinical_direction = NULL,
                                     caption_text = paste0("Point position reflects association with clinical variable; point color indicates treatment effect direction.\nPoints are colored if concordant with clinical variable direction after treatment. \nUp to top ", top_n, " from each direction are labeled."),
                                     cell_type = "", top_n = 50) {
  
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
    slice_head(n = top_n) %>%
    pull(Gene)
  
  top_50_right <- clin_results %>%
    filter(Gene %in% concordant_genes, !!sym(fc) > 0) %>%
    arrange(!!sym(p_col)) %>%
    slice_head(n = top_n) %>%
    pull(Gene)
  
  # Combine for plotting
  top_label_genes <- c(top_50_left, top_50_right)
  
  message("Number of genes significant in both: ", length(sig_both_genes))
  message("Number of concordant genes: ", n_concordant)
  message("Number of labeled genes: ", length(top_label_genes))
  
  clin_results <- clin_results %>%
    mutate(
      color_var_plot = case_when(
        Gene %in% concordant_genes ~ treatment_results$treatment_direction[match(Gene, treatment_results$Gene)],
        TRUE ~ "NS/NC"
      ),
      shape_var_plot = case_when(
        Gene %in% sig_both_genes & !!sym(fc) < 0 & treatment_results$treatment_direction[match(Gene, treatment_results$Gene)] == "Up w/ Dapagliflozin" ~ "Left_Pos",
        Gene %in% sig_both_genes & !!sym(fc) < 0 & treatment_results$treatment_direction[match(Gene, treatment_results$Gene)] == "Down w/ Dapagliflozin" ~ "Left_Neg",
        Gene %in% sig_both_genes & !!sym(fc) > 0 & treatment_results$treatment_direction[match(Gene, treatment_results$Gene)] == "Up w/ Dapagliflozin" ~ "Right_Pos",
        Gene %in% sig_both_genes & !!sym(fc) > 0 & treatment_results$treatment_direction[match(Gene, treatment_results$Gene)] == "Down w/ Dapagliflozin" ~ "Right_Neg",
        TRUE ~ "NS/NC"
      ),
      top_lab = if_else(Gene %in% top_label_genes, Gene, ""),
      # Updated arrow logic based on clinical direction - only for concordant genes
      arrow_symbol = case_when(
        # Only show arrows for concordant genes
        Gene %in% concordant_genes & !is.null(clinical_direction) & clinical_direction == "+" & !!sym(fc) < 0 ~ "\u2193",
        Gene %in% concordant_genes & !is.null(clinical_direction) & clinical_direction == "+" & !!sym(fc) > 0 ~ "\u2191",
        Gene %in% concordant_genes & !is.null(clinical_direction) & clinical_direction == "-" & !!sym(fc) < 0 ~ "\u2191",
        Gene %in% concordant_genes & !is.null(clinical_direction) & clinical_direction == "-" & !!sym(fc) > 0 ~ "\u2193",
        TRUE ~ ""
      )
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
    geom_text(aes(label = arrow_symbol, color = color_var_plot,
                  vjust = ifelse(arrow_symbol == "\u2193", 1, -0.5)),
              size = 6, family = "Arial") +
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
    annotate("segment",
             x = 0, xend = 0,
             y = if(!is.null(clinical_direction) && clinical_direction == "+") y_max * 0.4 else y_max * 0.9, 
             yend = if(!is.null(clinical_direction) && clinical_direction == "+") y_max * 0.9 else y_max * 0.4,
             size = 10, linejoin = "mitre",
             color = if(!is.null(clinical_direction) && clinical_direction == "+") "#f7c1bf" else "#a9c9dd", 
             arrow = arrow(type = "closed")) +
    annotate("text", x = 0, y = ((y_max*0.9 + y_max*0.4)/2),
             label = arrow_label,
             fontface = "bold",
             color = if(!is.null(clinical_direction) && clinical_direction == "+") "#f28482" else "#457b9d") +
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

# ===========================================================================
# Function: create_gene_expression_plots
# ===========================================================================


create_gene_expression_plots <- function(main_results,
                                         subtype_results_list,
                                         cell_type_labels = NULL,
                                         cell_type_order = NULL,
                                         cell_type_prefix = "PT",
                                         n_top_genes = 20,
                                         output_dir = ".",
                                         save_plots = TRUE,
                                         output_prefix = cell_type_prefix,
                                         logfc_col = "logFC_treatmentDapagliflozin:visitPOST",
                                         fdr_col = "fdr",
                                         volcano_pval_col = "fdr") {
  
  # Load required libraries
  require(tidyverse)
  require(ggtext)
  require(patchwork)
  
  # Determine cell type labels
  if (is.null(cell_type_labels)) {
    # Use automatic labeling with prefix
    cell_type_labels <- paste0(cell_type_prefix, "-", names(subtype_results_list))
  } else {
    # Validate that the number of labels matches the number of subtypes
    if (length(cell_type_labels) != length(subtype_results_list)) {
      stop("Length of cell_type_labels must match length of subtype_results_list")
    }
  }
  
  # Create a mapping between list names and labels
  label_mapping <- setNames(cell_type_labels, names(subtype_results_list))
  
  # Determine cell type order
  if (is.null(cell_type_order)) {
    # Use reverse order of labels (for bottom to top display)
    cell_type_order <- rev(cell_type_labels)
  } else {
    # Validate that all labels are included in the order
    if (!all(cell_type_labels %in% cell_type_order) || !all(cell_type_order %in% cell_type_labels)) {
      stop("cell_type_order must contain exactly the same labels as cell_type_labels")
    }
  }
  
  # Determine output prefix
  if (is.null(output_prefix)) {
    output_prefix <- ifelse(is.null(cell_type_labels), cell_type_prefix, "celltypes")
  }
  
  # Select top positive genes
  top_pos_genes <- main_results %>%
    arrange(!!sym(fdr_col)) %>%
    filter(!!sym(logfc_col) > 0) %>%
    head(n_top_genes) %>%
    pull(Gene)
  
  # Select top negative genes
  top_neg_genes <- main_results %>%
    arrange(!!sym(fdr_col)) %>%
    filter(!!sym(logfc_col) < 0) %>%
    head(n_top_genes) %>%
    pull(Gene)
  
  # Create combined plot data
  combined_plot_dat <- purrr::imap_dfr(subtype_results_list, function(res_df, subtype_key) {
    label <- label_mapping[subtype_key]
    make_plot_df(res_df, label, top_pos_genes, top_neg_genes)
  })
  
  # Set factor levels for ordering
  combined_plot_dat$celltype <- factor(combined_plot_dat$celltype, levels = cell_type_order)
  combined_plot_dat$Gene <- factor(combined_plot_dat$Gene, levels = c(top_neg_genes, top_pos_genes))
  
  # Create color vector for x-axis labels
  gene_levels <- levels(combined_plot_dat$Gene)
  x_colors <- setNames(
    c(rep("#457b9d", n_top_genes), rep("#f28482", n_top_genes)),
    gene_levels
  )
  
  # Create dot plot
  dot_plot <- create_dot_plot(combined_plot_dat, x_colors, n_top_genes)
  
  # Get main logFC values for the selected genes
  main_logfc_values <- main_results %>%
    filter(Gene %in% c(top_neg_genes, top_pos_genes)) %>%
    dplyr::select(Gene, main_logFC = !!sym(logfc_col))
  
  # Create consistency score data with weighted scores
  consistency_scores <- calculate_weighted_consistency_scores(
    combined_plot_dat, 
    main_logfc_values,
    top_neg_genes, 
    top_pos_genes, 
    cell_type_labels
  )
  
  # Create bar plot
  bar_plot <- create_bar_plot(consistency_scores)
  
  # Create volcano plot if requested
  volcano_plot <- NULL
  if (!is.null(volcano_pval_col)) {
    volcano_plot <- plot_volcano(
      main_results, 
      logfc_col, 
      volcano_pval_col,
      title = NULL,
      paste0("logFC ", gsub("logFC_", "", logfc_col)), 
      "-log10(FDR adjusted p-value)",
      paste0("nebula/", tolower(output_prefix), "_placebo_pvalue_nebula_reml_pooled"),
      cell_type = output_prefix
    )
  }
  
  # Combine plots
  if (!is.null(volcano_plot)) {
    layout <- c(
      area(1, 1, 1, 4),
      area(2, 1, 5, 4),
      area(1, 5, 5, 7)
    )
    combined_plot <- bar_plot + dot_plot + volcano_plot + plot_layout(design = layout)
    
    if (save_plots) {
      output_path <- file.path(output_dir, paste0(output_prefix, "_subtypes_nebula_scores_volcano.jpeg"))
      ggsave(output_path, width = 12, height = 7, plot = combined_plot)
    }
  }
  
  if (save_plots) {
    layout <- c(
      area(1, 1, 1, 1),
      area(2, 1, 5, 1)
    )
    combined_plot <- bar_plot + dot_plot + plot_layout(design = layout)
    
    if (save_plots) {
      output_path <- file.path(output_dir, paste0(output_prefix, "_subtypes_nebula_scores.jpeg"))
      ggsave(output_path, width = 7, height = 5, plot = combined_plot)
    }
  }
  
  # Save individual plots if requested
  if (save_plots) {
    ggsave(file.path(output_dir, paste0(output_prefix, "_subtypes_nebula.jpeg")), 
           width = 7, height = 5, plot = dot_plot)
  }
  
  # Return all components
  return(list(
    dot_plot = dot_plot,
    bar_plot = bar_plot,
    volcano_plot = volcano_plot,
    combined_plot = combined_plot,
    top_pos_genes = top_pos_genes,
    top_neg_genes = top_neg_genes,
    consistency_scores = consistency_scores,
    combined_plot_data = combined_plot_dat
  ))
}

#' Create dot plot for gene expression
create_dot_plot <- function(plot_data, x_colors, n_genes_per_direction) {
  plot_data %>%
    ggplot(aes(x = Gene, y = celltype, size = abs(logFC), color = direction)) +
    annotate("rect", 
             xmin = -Inf, 
             xmax = 0.5 + length(unique(plot_data$Gene)) / 2, 
             ymin = -Inf, 
             ymax = Inf, 
             fill = "#dceef5", 
             alpha = 0.5) +  
    annotate("rect", 
             xmin = 0.5 + length(unique(plot_data$Gene)) / 2, 
             xmax = Inf, 
             ymin = -Inf, 
             ymax = Inf, 
             fill = "#f9dada", 
             alpha = 0.5) +
    geom_point() +
    scale_color_manual(values = c("+" = "#f28482", "-" = "#457b9d")) +
    scale_x_discrete(
      guide = guide_axis(angle = 70), 
      labels = function(x) {
        mapply(function(label, col) {
          glue::glue("<span style='color:{col}'>{label}</span>")
        }, x, x_colors[x], SIMPLIFY = TRUE)
      }
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_markdown(angle = 70, hjust = 1),
      panel.border = element_blank(),
      text = element_text(size = 10),
      legend.key = element_rect(fill = NA),
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.position = c(0.5, .95),
      legend.background = element_rect(fill = NA, color = NA),
      legend.title = element_text(size = 8),
      plot.margin = margin(t = 0, r = 5, b = 0, l = 5)
    ) +
    labs(x = NULL, y = NULL, color = "Direction")
}

#' Calculate weighted consistency scores for genes based on logFC differences
calculate_weighted_consistency_scores <- function(plot_data, main_logfc_values, neg_genes, pos_genes, cell_type_labels) {
  # Join with main logFC values
  plot_data_with_main <- plot_data %>%
    left_join(main_logfc_values, by = "Gene")
  
  # Calculate weighted scores
  weighted_scores <- plot_data_with_main %>%
    group_by(Gene) %>%
    dplyr::summarize(
      # Calculate weighted score based on logFC ratio
      # For genes with same direction as main: positive contribution
      # For genes with opposite direction: stronger negative penalty
      # Normalize by number of cell types
      score = {
        main_fc <- dplyr::first(main_logFC)
        if (abs(main_fc) < 0.001) {
          # Handle near-zero main logFC
          0
        } else {
          # Calculate raw weights
          raw_weights <- logFC / main_fc
          
          # Apply penalty factor for opposite direction
          # Same direction: use raw weight
          # Opposite direction: apply penalty multiplier (e.g., 2x)
          penalty_factor <- 2  # Adjust this to control penalty strength
          
          weights <- ifelse(
            sign(logFC) == sign(main_fc),
            raw_weights,  # Same direction: positive contribution
            raw_weights * penalty_factor  # Opposite direction: amplified negative contribution
          )
          
          # Cap absolute weights at 2 to avoid extreme values
          weights <- pmax(pmin(weights, 2), -2)
          
          # Calculate mean
          mean_weight <- mean(weights, na.rm = TRUE)
          
          # Additional penalty: if any subtype is opposite direction, apply reduction factor
          has_opposite <- any(sign(logFC) != sign(main_fc), na.rm = TRUE)
          if (has_opposite) {
            # Count proportion of opposite direction subtypes
            prop_opposite <- sum(sign(logFC) != sign(main_fc), na.rm = TRUE) / length(logFC)
            # Apply additional reduction based on proportion of opposite effects
            mean_weight <- mean_weight * (1 - 0.5 * prop_opposite)
          }
          
          mean_weight
        }
      },
      main_logFC = dplyr::first(main_logFC),
      .groups = "drop"
    ) %>%
    # Ensure genes are in the correct order
    mutate(
      direction = case_when(main_logFC < 0 ~ "-", main_logFC > 0 ~ "+"),
      Gene = factor(Gene, levels = c(neg_genes, pos_genes))
    ) %>%
    arrange(Gene)
  
  return(weighted_scores)
}

#' Create bar plot for consistency scores
create_bar_plot <- function(scores_data) {
  scores_data %>%
    ggplot(aes(x = Gene, y = score, fill = direction)) +
    geom_hline(yintercept = 1, linetype = "dashed", size = 0.4, color = "#a7c957") +
    geom_hline(yintercept = 0.5, linetype = "dashed", size = 0.4, color = "#f4a261") +
    geom_col_rounded(width = 0.8) +
    geom_text(
      aes(label = Gene),
      vjust = 0.5, 
      hjust = ifelse(scores_data$score >= 0, 1.3, -0.3),
      size = 2,
      angle = 90,
      color = "white"
    ) +
    scale_fill_manual(values = c("+" = "#f28482", "-" = "#457b9d")) +
    scale_y_continuous(
      breaks = c(-2, -1, 0, 1, 2),
      labels = c("-2", "-1", "0", "1", "2"),
      limits = c(min(-0.1, min(scores_data$score, na.rm = TRUE) * 1.1), 
                 max(0.1, max(scores_data$score, na.rm = TRUE) * 1.1))
    ) +
    theme_bw() + 
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_blank(),
      panel.border = element_blank(),
      text = element_text(size = 10),
      axis.text.y = element_text(size = 6),
      axis.title.y = element_text(size = 8),
      legend.position = "none",
      legend.background = element_rect(fill = NA, color = NA),
      plot.margin = margin(t = 5, r = 5, b = 0, l = 5)
    ) +
    labs(x = NULL, y = "Weighted \nconsistency score")
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
# Function: t1dhc_run_cell_type_analysis
# ===========================================================================

# Function to run complete analysis for a given cell type
trt_run_cell_type_analysis <- function(cell_type, 
                                         input_path, 
                                         input_suffix,
                                         output_base_path,
                                         output_prefix,
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
  
  res <- s3readRDS(input_file, bucket = bucket, region = region)
  
  # Check convergence
  res_convergence <- map_dfr(
    names(res),
    function(gene_name) {
      converged <- res[[gene_name]]$convergence
      df <- data.frame(Gene = gene_name,
                       Convergence_Code = converged)
      return(df)
    }
  )
  
  res_converged <- res_convergence %>%
    filter(Convergence_Code >= -10)
  
  # Combine results
  res_combined <- map_dfr(
    names(res),
    function(gene_name) {
      df <- res[[gene_name]]$summary
      df <- df %>% 
        filter(gene_name %in% res_converged$Gene) %>%
        mutate(Gene = gene_name)
      return(df)
    }
  )
  
  # Calculate FDR
  res_combined <- res_combined %>%
    ungroup() %>%
    mutate(fdr = p.adjust(`p_treatmentDapagliflozin:visitPOST`, method = "fdr"))
  
  res_combined <- subset(res_combined, abs(`logFC_treatmentDapagliflozin:visitPOST`) < 10)
  
  # Save results
  output_csv <- file.path(output_base_path, "Results", "NEBULA", 
                          paste0(output_prefix, cell_type_lower, "_nebula_res.csv"))
  write.csv(res_combined, output_csv, row.names = FALSE)
  
}
# ===========================================================================
# Function: t1dhc_run_cell_type_analysis
# ===========================================================================

# Function to run complete analysis for a given cell type
t1dhc_run_cell_type_analysis <- function(cell_type, 
                                   input_path, 
                                   input_suffix,
                                   output_base_path,
                                   output_prefix,
                                   bg_path,
                                   plot_title = "Volcano Plot",
                                   bucket = "attempt",
                                   region = "") {
  
  # Print status
  cat("Starting analysis for cell type:", cell_type, "\n")
  
  # Create cell type specific paths
  cell_type_lower <- tolower(cell_type)
  
  # Read in nebula results
  input_file <- file.path(input_path, cell_type, "nebula", 
                          paste0(cell_type_lower, input_suffix))
  
  res <- s3readRDS(input_file, bucket = bucket, region = region)
  
  # Check convergence
  res_convergence <- map_dfr(
    names(res),
    function(gene_name) {
      converged <- res[[gene_name]]$convergence
      df <- data.frame(Gene = gene_name,
                       Convergence_Code = converged)
      return(df)
    }
  )
  
  res_converged <- res_convergence %>%
    filter(Convergence_Code >= -10)
  
  # Combine results
  res_combined <- map_dfr(
    names(res),
    function(gene_name) {
      df <- res[[gene_name]]$summary
      df <- df %>% 
        filter(gene_name %in% res_converged$Gene) %>%
        mutate(Gene = gene_name)
      return(df)
    }
  )
  
  # Calculate FDR
  res_combined <- res_combined %>%
    ungroup() %>%
    mutate(fdr = p.adjust(p_groupType_1_Diabetes, method = "fdr"))
  
  res_combined <- subset(res_combined, abs(logFC_groupType_1_Diabetes) < 10)
  # Save results
  
  output_csv <- file.path(output_base_path, "Results", "NEBULA", 
                          paste0(output_prefix, cell_type_lower, "_nebula_res.csv"))
  write.csv(res_combined, output_csv, row.names = FALSE)
  
  # Add direction columns for visualization
  res_combined <- res_combined %>%
    mutate(group_direction = case_when(
      p_groupType_1_Diabetes < 0.05 & logFC_groupType_1_Diabetes > 0 ~ "Positive",
      p_groupType_1_Diabetes < 0.05 & logFC_groupType_1_Diabetes < 0 ~ "Negative",
      TRUE ~ "NS"
    ),
    group_direction_fdr = case_when(
      fdr < 0.05 & logFC_groupType_1_Diabetes > 0 ~ "Positive",
      fdr < 0.05 & logFC_groupType_1_Diabetes < 0 ~ "Negative",
      TRUE ~ "NS"
    ))
  
  # Create volcano plots
  plot_volcano(res_combined, 
               "logFC_groupType_1_Diabetes", 
               "p_groupType_1_Diabetes",
               NULL,
               "logFC group", 
               "-log10(p-value)",
               "_deg",
               formula = "group",
               positive_text = "Upregulated in T1D compared to LC",
               negative_text = "Downregulated in T1D compared to LC",
               cell_type = cell_type,
               output_base_path = file.path(output_base_path, "Results", "Figures", "Volcano", "nebula",
                                            paste0(output_prefix, cell_type_lower)))
  
  plot_volcano(res_combined, 
               "logFC_groupType_1_Diabetes", 
               "fdr",
               NULL,
               "logFC group", 
               "-log10(FDR adjusted p-value)",
               "_deg_fdr",
               formula = "group",
               positive_text = "Upregulated in T1D compared to LC",
               negative_text = "Downregulated in T1D compared to LC",
               cell_type = cell_type,
               output_base_path = file.path(output_base_path, "Results", "Figures", "Volcano", "nebula",
                                            paste0(output_prefix, cell_type_lower)))
  
  # GSEA Analysis
  # Load pathway files
  gmt_files <- list.files(path = bg_path, pattern = '.gmt', full.names = TRUE)
  kegg_legacy <- prepare_gmt(gmt_files[1], unique(res_combined$Gene), savefile = FALSE)
  reactome <- prepare_gmt(gmt_files[3], unique(res_combined$Gene), savefile = FALSE)
  go <- prepare_gmt(gmt_files[4], unique(res_combined$Gene), savefile = FALSE)
  
  # Rank genes by logFC
  rankings <- res_combined$logFC_groupType_1_Diabetes
  names(rankings) <- res_combined$Gene
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
                   paste0(output_prefix, cell_type_lower, "_res_top30_kegg_pathways.jpeg")),
         width = 27.5, height = 14, scale = 1)
  
  # REACTOME
  plot_fgsea_transpose(reactome_res, 
                       title = paste(cell_type, "Top 30 REACTOME Pathways"), 
                       xmin = 1.45, xmax = 2.6)
  
  ggsave(file.path(output_base_path, "Results", "Figures", "Pathways", "nebula",
                   paste0(output_prefix, cell_type_lower, "_res_top30_reactome_pathways.jpeg")),
         width = 27.5, height = 14, scale = 1)
  
  # GO
  plot_fgsea_transpose(go_res, 
                       title = paste(cell_type, "Top 30 GO Pathways"), 
                       xmin = 1.4, xmax = 5)
  
  ggsave(file.path(output_base_path, "Results", "Figures", "Pathways", "nebula",
                   paste0(output_prefix, cell_type_lower, "_res_top30_go_pathways.jpeg")),
         width = 27.5, height = 14, scale = 1)
  
  # Return results
  return(list(
    nebula_results = res_combined,
    fgsea_summary = fgsea_summary,
    kegg_results = kegg_legacy_res,
    reactome_results = reactome_res,
    go_results = go_res,
    rankings = rankings
  ))
}

# slingshot
# ===========================================================================
# Function: slingshot_setup
# ===========================================================================

slingshot_setup <- function(seurat_obj, celltype_prefix) {
  library(SingleCellExperiment)
  library(slingshot)
  library(uwot)
  library(mclust)
  library(ggplot2)
  library(RColorBrewer)
  
  # Extract counts and convert to SCE
  counts <- GetAssayData(seurat_obj, layer = "counts")
  sce <- SingleCellExperiment(assays = List(counts = counts), colData = seurat_obj@meta.data)
  
  # Gene filtering
  gene_filter <- apply(assays(sce)$counts, 1, function(x) sum(x >= 3) >= 10)
  sce <- sce[gene_filter, ]
  
  # Normalization using quantile normalization (FQ)
  FQnorm <- function(counts){
    rk <- apply(counts, 2, rank, ties.method = 'min')
    counts.sort <- apply(counts, 2, sort)
    refdist <- apply(counts.sort, 1, median)
    norm <- apply(rk, 2, function(r) { refdist[r] })
    rownames(norm) <- rownames(counts)
    return(norm)
  }
  assays(sce)$norm <- FQnorm(assays(sce)$counts)
  
  # PCA
  pca <- prcomp(t(log1p(assays(sce)$norm)), scale. = FALSE)
  
  # Variance explained
  var_explained <- pca$sdev^2 / sum(pca$sdev^2)
  
  # Elbow plot as ggplot
  elbow_df <- data.frame(PC = 1:50, VarianceExplained = var_explained[1:50])
  elbow_plot <- ggplot(elbow_df, aes(x = PC, y = VarianceExplained)) +
    geom_point(size = 2) +
    geom_line() +
    labs(title = paste("Elbow Plot:", celltype_prefix),
         x = "Principal Component",
         y = "Proportion of Variance Explained") +
    theme_bw()
  
  # Return necessary outputs
  return(list(
    sce = sce,
    pca = pca,
    var_explained = var_explained,
    elbow_plot = elbow_plot
  ))
}

# ===========================================================================
# Function: run_slingshot
# ===========================================================================

run_slingshot <- function(sce, pca_obj, n_pcs = 6, start_cluster = NULL, end_cluster = NULL, cluster_label = "celltype") {
  library(slingshot)
  library(uwot)
  
  # Get top PCs
  rd1 <- pca_obj$x[, 1:n_pcs]
  
  # Visualize PCA (optional)
  plot(rd1, col = rgb(0, 0, 0, 0.5), pch = 16, asp = 1, main = "PCA")
  
  # Run UMAP
  umap_mat <- uwot::umap(t(log1p(assays(sce)$norm)))
  colnames(umap_mat) <- c("UMAP1", "UMAP2")
  
  # Visualize UMAP (optional)
  plot(umap_mat, col = rgb(0, 0, 0, 0.5), pch = 16, asp = 1, main = "UMAP")
  
  # Add to reducedDims
  reducedDims(sce) <- SimpleList(PCA = rd1, UMAP = umap_mat)
  
  # Run Slingshot
  sce_sl <- slingshot(sce, clusterLabels = cluster_label, reducedDim = "PCA", 
                      start.clus = start_cluster,
                      end.clus = end_cluster)
  
  return(sce_sl)
}

# ===========================================================================
# Function: plot_slingshot_trajectory
# ===========================================================================

plot_slingshot_trajectory <- function(sce_sl, 
                                      celltype_levels, 
                                      custom_colors = color_5, 
                                      cluster_label = "celltype",
                                      bucket = "attempt",
                                      celltype_suffix = NULL,
                                      title = "Slingshot trajectory",
                                      lineage = 1) {
  library(ggplot2)
  library(RColorBrewer)
  library(slingshot)
  
  # Extract PCA and pseudotime
  pca_df <- as.data.frame(reducedDims(sce_sl)$PCA)
  
  # Dynamically extract cluster labels
  pca_df$cluster_label <- factor(colData(sce_sl)[[cluster_label]], levels = celltype_levels)
  pca_df$pseudotime <- slingPseudotime(sce_sl)[, lineage]
  
  # Get curve
  curves <- slingCurves(sce_sl)
  curve_df <- as.data.frame(curves[[lineage]]$s)
  lineage_obj <- slingshot::slingLineages(sce_sl)
  
  # Generate ggplot
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = cluster_label)) +
    geom_point(alpha = 0.3, size = 1) +
    geom_path(data = curve_df, aes(x = PC1, y = PC2), 
              color = "black", size = 1.2, inherit.aes = FALSE) +
    theme_bw() + 
    theme(panel.grid = element_blank()) + 
    labs(color = NULL, title = title) +
    scale_color_manual(values = custom_colors)
  
  # Save if bucket is provided
  if (!is.null(bucket)) {
    temp_file <- tempfile(fileext = ".jpeg")
    ggsave(filename = temp_file, plot = p, width = 7, height = 5)
    s3$upload_file(temp_file, bucket, paste0("slingshot/attempt_pca_", tolower(celltype_suffix), "_slingshot.jpeg"))
  }
  
  return(list(
    pca_plot = p,
    curves = curves,
    pca_df = pca_df,
    curve_df = curve_df,
    lineage = lineage_obj
  ))
}
# ===========================================================================
# Function: create_pseudotime_df
# ===========================================================================

create_pseudotime_df <- function(sce_obj, lineage = 1) {
  # Check if slingPseudotime exists
  if (!"slingPseudotime_1" %in% names(colData(sce_obj))) {
    stop("The input object does not contain slingPseudotime. Run Slingshot first.")
  }
  
  # Extract relevant metadata and pseudotime
  pseudotime_df <- data.frame(
    pseudotime = slingPseudotime(sce_obj)[, lineage],
    treatment = sce_obj$treatment,
    visit = sce_obj$visit,
    celltype = sce_obj$celltype,
    visit_treatment = paste(sce_obj$visit, sce_obj$treatment)
  )
  
  # Reorder factor levels
  pseudotime_df$visit_treatment <- factor(pseudotime_df$visit_treatment, 
                                          levels = c("PRE Placebo", "POST Placebo", 
                                                     "PRE Dapagliflozin", "POST Dapagliflozin"))
  return(pseudotime_df)
}


# ===========================================================================
# Function: plot_pseudotime_violin
# ===========================================================================

plot_pseudotime_violin <- function(df, s3_folder = "slingshot", 
                                   celltype_suffix = "celltype") {
  pseudotime_df <- df
  # Calculate medians
  median_df <- pseudotime_df %>%
    group_by(visit_treatment) %>%
    summarise(median_pt = median(pseudotime, na.rm = TRUE), .groups = "drop")
  
  # Define colors
  fill_colors <- c("PRE Placebo" = "#fbc4ab", "POST Placebo" = "#f4978e",
                   "PRE Dapagliflozin" = "#ccd5ae", "POST Dapagliflozin" = "#828e82")
  
  # Create plot
  p <- ggplot(pseudotime_df, aes(x = visit_treatment, y = pseudotime)) +
    geom_violin(aes(fill = visit_treatment, color = visit_treatment), trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.1, outlier.shape = NA, aes(color = visit_treatment)) +
    geom_line(data = subset(median_df, grepl("Placebo", visit_treatment)), 
              aes(x = visit_treatment, y = median_pt, group = 1), 
              color = "#343a40", linewidth = 0.5, linetype = "dashed") +
    geom_line(data = subset(median_df, grepl("Dapagliflozin", visit_treatment)), 
              aes(x = visit_treatment, y = median_pt, group = 1), 
              color = "#343a40", linewidth = 0.5, linetype = "dashed") +
    theme_classic() +
    labs(x = NULL, y = "Pseudotime", fill = NULL, color = NULL) +
    scale_fill_manual(values = fill_colors) +
    scale_color_manual(values = fill_colors) +
    theme(legend.position = "none",
          text = element_text(size = 15),
          axis.text.x = element_text(angle = 30, hjust = 1))
  
  # Save and upload to S3
  temp_file <- tempfile(fileext = ".jpeg")
  ggsave(filename = temp_file, plot = p, width = 7, height = 5)
  s3$upload_file(temp_file, "attempt", file.path(s3_folder, paste0(tolower(celltype_suffix), "_attempt_slingshot_violin.jpeg")))
  return(p)
}

# ===========================================================================
# Function: plot_slingshot_3d
# ===========================================================================

library(plotly)
library(htmlwidgets)

plot_slingshot_3d <- function(pca_df, curve_df, custom_colors = color_5,
                              s3_folder = "slingshot", 
                              cluster_label = "celltype",
                              celltype_suffix = "celltype") {
  
  # Create plot
  pt_plotly <- plot_ly() %>%
    add_trace(
      data = pca_df,
      x = ~PC1, y = ~PC2, z = ~PC3,
      type = "scatter3d",
      mode = "markers",
      color = ~cluster_label,
      colors = custom_colors,
      marker = list(size = 2),
      name = "Cells") %>%
    add_trace(
      data = curve_df,
      x = ~PC1, y = ~PC2, z = ~PC3,
      type = "scatter3d",
      mode = "lines",
      line = list(color = "black", width = 4),
      name = "Trajectory")
  
  # Save and upload to S3
  temp_file <- tempfile(fileext = ".html")
  saveWidget(pt_plotly, file = temp_file, selfcontained = TRUE)
  s3$upload_file(temp_file, "attempt", file.path(s3_folder, paste0("attempt_", tolower(celltype_suffix), "_slingshot_plotly.html")))
  return(pt_plotly)
}

# ===========================================================================
# Function: plot_and_test_pseudotime_distribution
# ===========================================================================

plot_and_test_pseudotime_distribution <- function(df,
                                                  sce_object,
                                                  pseudotime_var = "pseudotime",
                                                  visit_treatment_var = "visit_treatment",
                                                  s3_folder = "slingshot",
                                                  filename_suffix = "celltype") {
  library(condiments)
  BiocParallel::register(BiocParallel::SerialParam())
  # Set colors
  group_colors <- c("PRE Placebo" = "#fbc4ab",
                    "POST Placebo" = "#f4978e",
                    "PRE Dapagliflozin" = "#ccd5ae",
                    "POST Dapagliflozin" = "#828e82")
  
  # Create plot
  p <- ggplot(df, aes(x = .data[[pseudotime_var]],
                      fill = .data[[visit_treatment_var]],
                      color = .data[[visit_treatment_var]])) +
    geom_density(alpha = 0.3) +
    theme_minimal() +
    labs(x = "Pseudotime", y = "Density", fill = NULL, color = NULL) +
    theme(
      legend.position = c(0.5, 0.85),
      legend.direction = "horizontal",
      legend.text = element_text(size = 10),
      panel.grid = element_blank(),
      text = element_text(size = 15)
    ) +
    scale_fill_manual(values = group_colors) +
    scale_color_manual(values = group_colors)
  
  # Save and upload
  temp_file <- tempfile(fileext = ".jpeg")
  ggsave(filename = temp_file, plot = p, width = 7, height = 5)
  s3$upload_file(temp_file, "attempt", file.path(s3_folder, paste0("attempt_", filename_suffix, "_slingshot_density_trtvisit.jpeg")))
  
  # Run progressionTest
  test_result <- progressionTest(sce_object, conditions = df[[visit_treatment_var]])
  print(test_result)
  return(p)
  invisible(list(plot = p, test_result = test_result))
}

# ===========================================================================
# Function: plot_pseudotime_density_faceted_by_treatment
# ===========================================================================

plot_pseudotime_density_faceted_by_treatment <- function(df,
                                                         pseudotime_var = "slingPseudotime_1",
                                                         visit_col = "visit",
                                                         treatment_col = "treatment",
                                                         visit_treatment_col = "visit_treatment",
                                                         s3_folder = "slingshot",
                                                         filename_suffix = "celltype") {
  # Quantile summary
  raw_summary <- df %>%
    as.data.frame() %>%
    group_by(.data[[treatment_col]], .data[[visit_col]]) %>%
    summarise(
      p25 = quantile(.data[[pseudotime_var]], probs = 0.25, na.rm = TRUE),
      p65 = quantile(.data[[pseudotime_var]], probs = 0.65, na.rm = TRUE),
      p85 = quantile(.data[[pseudotime_var]], probs = 0.85, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(visit_treatment = paste(.data[[visit_col]], .data[[treatment_col]]))
  
  # Color palette
  fill_colors <- c("PRE Placebo" = "#fbc4ab", 
                   "POST Placebo" = "#f4978e",
                   "PRE Dapagliflozin" = "#ccd5ae", 
                   "POST Dapagliflozin" = "#828e82")
  
  # Plot
  p <- ggplot(df, aes(x = .data[[pseudotime_var]], fill = .data[[visit_treatment_col]])) +
    geom_density(alpha = 0.5, aes(color = .data[[visit_treatment_col]])) +
    facet_wrap(vars(.data[[treatment_col]]), strip.position = "bottom") +
    geom_vline(data = raw_summary, aes(xintercept = p25, color = visit_treatment), linetype = "dashed") +
    geom_vline(data = raw_summary, aes(xintercept = p65, color = visit_treatment), linetype = "dashed") +
    geom_vline(data = raw_summary, aes(xintercept = p85, color = visit_treatment), linetype = "dashed") +
    theme_minimal() +
    labs(x = "Pseudotime", y = "Density", color = NULL, fill = NULL) +
    theme(panel.grid = element_blank(),
          text = element_text(size = 15),
          legend.position = c(0.45, 0.95),
          legend.direction = "horizontal",
          legend.text = element_text(size = 12)) +
    scale_fill_manual(values = fill_colors) +
    scale_color_manual(values = fill_colors)
  
  # Save and upload
  temp_file <- tempfile(fileext = ".jpeg")
  ggsave(filename = temp_file, plot = p, width = 7, height = 5)
  s3$upload_file(temp_file, "attempt", file.path(s3_folder, paste0("attempt_", filename_suffix, "_slingshot_density_trt.jpeg")))
  
  return(p)
}

# ===========================================================================
# Function: plot_delta_percentile_heatmap
# ===========================================================================

plot_delta_percentile_heatmap <- function(df,
                                          percentile_prefix = "rq",
                                          visit_col = "visit",
                                          treatment_col = "treatment",
                                          s3_folder = "slingshot",
                                          filename_suffix = "celltype",
                                          width = 5,
                                          height = 5) {
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  # Prepare data for heatmap
  heatmap_df <- df %>%
    pivot_longer(cols = starts_with(percentile_prefix),
                 names_to = "q",
                 names_prefix = percentile_prefix,
                 values_to = "value") %>%
    pivot_wider(names_from = .data[[visit_col]], values_from = value) %>%
    mutate(delta = POST - PRE)
  
  # Plot
  p <- ggplot(heatmap_df, aes(x = .data[[treatment_col]], y = q, fill = delta)) +
    geom_tile() +
    scale_fill_gradient2(low = "#89c2d9", mid = "white", high = "#ee7674", midpoint = 0) +
    labs(x = NULL, y = "Percentile", fill = "\u0394\n(POST-PRE)") +
    theme_minimal() +
    theme(
      legend.title.position = "top",
      legend.title = element_text(hjust = 0.5, size = 10),
      legend.position = c(1.08, 0.5),
      panel.grid = element_blank(),
      text = element_text(size = 15),
      legend.text = element_text(size = 10),
      plot.margin = unit(c(0, 1.8, 0, 0), "cm")
    )
  
  # Save and upload
  temp_file <- tempfile(fileext = ".jpeg")
  ggsave(filename = temp_file, plot = p, width = width, height = height)
  s3$upload_file(temp_file, "attempt", file.path(s3_folder, paste0("attempt_", filename_suffix, "_slingshot_percentile_heatmap.jpeg")))
  
  return(p)
}

# ===========================================================================
# Function: analyze_pseudotime_by_clinvar
# ===========================================================================

library(dplyr)
library(ggplot2)

analyze_pseudotime_by_clinvar <- function(df,
                                          clinical_var,          # unquoted column name of clinical variable
                                          pseudotime_var,        # unquoted column name of pseudotime
                                          subject_id = "subject_id",
                                          visit_col = "visit",
                                          pre_label = "PRE",
                                          post_label = "POST",
                                          bin_probs = 2,
                                          bin_levels_to_compare = c(1, 2),
                                          caption_clinical_var = "",
                                          bucket = "attempt",
                                          celltype_suffix = "",
                                          filesuffix = "") {
  
  clinical_var <- ensym(clinical_var)
  pseudotime_var <- ensym(pseudotime_var)
  clinical_var_chr <- rlang::as_string(clinical_var)
  
  # Step 1: Calculate subject-level delta of the clinical variable
  delta_df <- df %>%
    group_by(.data[[subject_id]]) %>%
    summarise(
      delta_value = mean(.data[[clinical_var_chr]][.data[[visit_col]] == post_label], na.rm = TRUE) -
        mean(.data[[clinical_var_chr]][.data[[visit_col]] == pre_label], na.rm = TRUE),
      .groups = "drop"
    )
  
  # Step 2: Join delta back to the original dataframe
  df <- df %>%
    left_join(delta_df, by = subject_id)
  
  # Step 3: Bin the clinical variable by quantiles
  df_binned <- if (identical(bin_probs, "direction")) {
    df %>%
      mutate(clinical_bin = case_when(
        .data[[clinical_var_chr]] > 0 ~ "+",
        .data[[clinical_var_chr]] < 0 ~ "-",
        .data[[clinical_var_chr]] == 0 ~  "No Change"
      )) %>%
      mutate(clinical_bin = factor(clinical_bin, levels = c("-", "No Change", "+")))
  } else {
    df %>%
      mutate(clinical_bin = cut(.data[[clinical_var_chr]],
                                breaks = quantile(.data[[clinical_var_chr]], 
                                                  probs = seq(0, 1, 1/bin_probs), na.rm = TRUE),
                                include.lowest = TRUE)) %>%
      filter(!is.na(clinical_bin))
  }
  
  # Step 4: Plot density of pseudotime by clinical bins
  n_bins <- nlevels(droplevels(df_binned$clinical_bin))
  custom_colors <- NULL
  if (n_bins == 2) {
    custom_colors <- c("#457b9d", "#f28482")
    names(custom_colors) <- levels(droplevels(df_binned$clinical_bin))
  }
  
  p <- df_binned %>%
    filter(!is.na(clinical_bin)) %>%
    ggplot(aes(x = !!pseudotime_var, fill = clinical_bin)) +
    geom_density(alpha = 0.5) +
    theme_minimal() +
    labs(
      y = "Density", x = "Pseudotime", fill = NULL, color = NULL,
      caption = paste0(caption_clinical_var, " across pseudotime")) +
    theme(panel.grid = element_blank(), 
          text = element_text(size = 15),
          legend.position = c(0.5,0.95),
          legend.direction = "horizontal") +
    {
      if (!is.null(custom_colors)) {
        list(
          scale_fill_manual(values = custom_colors)
        )
      } else {
        NULL
      }
    }
  
  if (!is.null(bucket)) {
    temp_file <- tempfile(fileext = ".jpeg") # need to create a temporary file
    ggsave(filename = temp_file, width = 7, height = 5)
    s3$upload_file(temp_file, bucket, paste0("slingshot/attempt_density_", 
                                             tolower(celltype_suffix), "_", clinical_var, 
                                             filesuffix,
                                             "_slingshot.jpeg"))
  }
  
  print(p)
  
  # Step 5: KS test between selected bin levels
  bin_levels <- levels(droplevels(df_binned$clinical_bin))
  if (length(bin_levels) >= max(bin_levels_to_compare)) {
    pt1 <- df_binned %>%
      filter(clinical_bin == bin_levels[bin_levels_to_compare[1]]) %>%
      pull(!!pseudotime_var)
    
    pt2 <- df_binned %>%
      filter(clinical_bin == bin_levels[bin_levels_to_compare[2]]) %>%
      pull(!!pseudotime_var)
    
    ks <- ks.test(pt1, pt2)
    print(ks)
  } else {
    warning("Not enough bins to compare the requested levels.")
    ks <- NULL
  }
  
  return(invisible(list(plot = p, ks_test = ks, delta_df = delta_df)))
}

# ===========================================================================
# Function: plot_clinvar_pseudotime_arrows
# ===========================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(rlang)
library(grid)

plot_clinvar_pseudotime_arrows <- function(df,
                                           pseudotime_var = "slingPseudotime_1",
                                           subject_col = "subject_id",
                                           visit_col = "visit",
                                           treatment_col = "treatment",
                                           visit_treatment_col = "visit_treatment",
                                           clinical_var,                # unquoted, e.g. hba1c
                                           clinical_var_label = NULL,   # for axis label
                                           percentile_filter = 50,
                                           percentiles = c(25, 50, 65, 85),
                                           shape_pre = 16,
                                           shape_post = 17,
                                           celltype_suffix = "",
                                           bucket = "attempt",
                                           color_palette = c("PRE Placebo" = "#fbc4ab",
                                                             "POST Placebo" = "#f4978e",
                                                             "PRE Dapagliflozin" = "#ccd5ae",
                                                             "POST Dapagliflozin" = "#828e82",
                                                             "Consistent" = "#160f29",
                                                             "Inconsistent" = "#e5e5e5")) {
  
  clinical_var <- rlang::ensym(clinical_var)
  pseudotime_var_chr <- rlang::as_string(rlang::ensym(pseudotime_var))
  clinical_var_chr <- rlang::as_string(clinical_var)
  
  # Step 1: Pseudotime percentiles
  med <- df %>%
    group_by(.data[[subject_col]], .data[[visit_col]]) %>%
    summarise(across(
      .data[[pseudotime_var_chr]],
      list("25" = ~quantile(., 0.25, na.rm = TRUE),
           "50" = ~quantile(., 0.50, na.rm = TRUE),
           "65" = ~quantile(., 0.65, na.rm = TRUE),
           "85" = ~quantile(., 0.85, na.rm = TRUE)),
      .names = "pseudotime_{fn}"
    ), .groups = "drop")
  
  # Step 2: Construct plot_df
  plot_df <- df %>%
    dplyr::select(all_of(c(subject_col, visit_col, treatment_col, clinical_var_chr, visit_treatment_col))) %>%
    distinct(.data[[subject_col]], .data[[visit_col]], .keep_all = TRUE) %>%
    left_join(med, by = c(subject_col, visit_col)) %>%
    pivot_longer(
      cols = starts_with("pseudotime_"),
      names_to = "percentile",
      names_prefix = "pseudotime_",
      values_to = "pseudotime"
    ) %>%
    mutate(percentile = as.numeric(percentile))
  
  # Step 3: Format wide for arrow plot
  arrow_df <- plot_df %>%
    dplyr::select(all_of(c(subject_col, treatment_col, clinical_var_chr, "percentile", visit_col, "pseudotime", visit_treatment_col))) %>%
    pivot_wider(names_from = .data[[visit_col]],
                values_from = c(pseudotime, !!clinical_var, !!visit_treatment_col)) %>%
    mutate(expectations = case_when(
      pseudotime_POST - pseudotime_PRE > 0 & !!sym(paste0(clinical_var_chr, "_POST")) - !!sym(paste0(clinical_var_chr, "_PRE")) > 0 ~ "Consistent",
      pseudotime_POST - pseudotime_PRE < 0 & !!sym(paste0(clinical_var_chr, "_POST")) - !!sym(paste0(clinical_var_chr, "_PRE")) < 0 ~ "Consistent",
      TRUE ~ "Inconsistent"
    ))
  
  # Step 4: Plot
  p <- arrow_df %>%
    filter(percentile == percentile_filter) %>%
    ggplot() +
    geom_segment(aes(x = pseudotime_PRE,
                     y = !!sym(paste0(clinical_var_chr, "_PRE")),
                     xend = pseudotime_POST,
                     yend = !!sym(paste0(clinical_var_chr, "_POST")),
                     color = expectations),
                 arrow = arrow(length = unit(0.5, "cm")),
                 alpha = 1) +
    geom_point(aes(x = pseudotime_PRE,
                   y = !!sym(paste0(clinical_var_chr, "_PRE")),
                   color = .data[[paste0(visit_treatment_col, "_PRE")]]),
               shape = shape_pre, size = 4, alpha = 0.7) +
    geom_point(aes(x = pseudotime_POST,
                   y = !!sym(paste0(clinical_var_chr, "_POST")),
                   color = .data[[paste0(visit_treatment_col, "_POST")]]),
               shape = shape_post, size = 4, alpha = 0.7) +
    scale_color_manual(values = color_palette) +
    labs(x = "Pseudotime",
         y = clinical_var_label %||% clinical_var_chr,
         color = NULL) +
    theme_minimal(base_size = 15) +
    theme(panel.grid = element_blank())
  
  if (!is.null(bucket)) {
    temp_file <- tempfile(fileext = ".jpeg") # need to create a temporary file
    ggsave(filename = temp_file, width = 7, height = 5)
    s3$upload_file(temp_file, bucket, paste0("slingshot/attempt_arrow_scatter", 
                                             tolower(celltype_suffix), "_", clinical_var, "_slingshot.jpeg"))
  }
  return(p)
}

# ===========================================================================
# Function: plot_clinvar_pseudotime_arrows
# ===========================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(rlang)
library(grid)

plot_clinvar_pseudotime_arrows_delta <- function(df,
                                           pseudotime_var = "slingPseudotime_1",
                                           subject_col = "subject_id",
                                           visit_col = "visit",
                                           treatment_col = "treatment",
                                           visit_treatment_col = "visit_treatment",
                                           clinical_var,                # unquoted, e.g. hba1c
                                           clinical_var_label = NULL,   # for axis label
                                           percentile_filter = 50,
                                           percentiles = c(25, 50, 65, 85),
                                           shape_pre = 16,
                                           shape_post = 17,
                                           bucket = "attempt",
                                           celltype_suffix = "",
                                           color_palette = c("PRE Placebo" = "#fbc4ab",
                                                             "POST Placebo" = "#f4978e",
                                                             "PRE Dapagliflozin" = "#ccd5ae",
                                                             "POST Dapagliflozin" = "#828e82",
                                                             "Healthy state" = "#8d99ae",
                                                             "Injured state" = "#e5e5e5")) {
  
  clinical_var <- rlang::ensym(clinical_var)
  pseudotime_var_chr <- rlang::as_string(rlang::ensym(pseudotime_var))
  clinical_var_chr <- rlang::as_string(clinical_var)
  
  # Step 1: Pseudotime percentiles
  med <- df %>%
    group_by(.data[[subject_col]], .data[[visit_col]]) %>%
    summarise(across(
      .data[[pseudotime_var_chr]],
      list("25" = ~quantile(., 0.25, na.rm = TRUE),
           "50" = ~quantile(., 0.50, na.rm = TRUE),
           "65" = ~quantile(., 0.65, na.rm = TRUE),
           "85" = ~quantile(., 0.85, na.rm = TRUE)),
      .names = "pseudotime_{fn}"
    ), .groups = "drop")
  
  # Step 2: Construct plot_df
  plot_df <- df %>%
    dplyr::select(all_of(c(subject_col, visit_col, treatment_col, clinical_var_chr, visit_treatment_col))) %>%
    distinct(.data[[subject_col]], .data[[visit_col]], .keep_all = TRUE) %>%
    left_join(med, by = c(subject_col, visit_col)) %>%
    pivot_longer(
      cols = starts_with("pseudotime_"),
      names_to = "percentile",
      names_prefix = "pseudotime_",
      values_to = "pseudotime"
    ) %>%
    mutate(percentile = as.numeric(percentile))
  
  # Step 3: Format wide for arrow plot
  arrow_df <- plot_df %>%
    dplyr::select(all_of(c(subject_col, treatment_col, clinical_var_chr, "percentile", visit_col, "pseudotime", visit_treatment_col))) %>%
    pivot_wider(names_from = .data[[visit_col]],
                values_from = c(pseudotime, !!clinical_var, !!visit_treatment_col)) %>%
    mutate(expectations = case_when(
      pseudotime_POST - pseudotime_PRE < 0 ~ "Healthy state",
      TRUE ~ "Injured state"
    ))
  
  # Step 4: Plot
  p <- arrow_df %>%
    filter(percentile == percentile_filter & 
             !is.na(expectations) & 
             (!is.na(.data[[paste0(visit_treatment_col, "_PRE")]])| !is.na(.data[[paste0(visit_treatment_col, "_POST")]]))) %>%
    ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "#2b2d42") +
    geom_segment(aes(x = pseudotime_PRE,
                     y = !!sym(paste0(clinical_var_chr, "_PRE")),
                     xend = pseudotime_POST,
                     yend = !!sym(paste0(clinical_var_chr, "_POST")),
                     color = expectations),
                 arrow = arrow(length = unit(0.5, "cm")),
                 alpha = 1) +
    geom_point(aes(x = pseudotime_PRE,
                   y = !!sym(paste0(clinical_var_chr, "_PRE")),
                   color = .data[[paste0(visit_treatment_col, "_PRE")]]),
               shape = shape_pre, size = 4, alpha = 0.7) +
    geom_point(aes(x = pseudotime_POST,
                   y = !!sym(paste0(clinical_var_chr, "_POST")),
                   color = .data[[paste0(visit_treatment_col, "_POST")]]),
               shape = shape_post, size = 4, alpha = 0.7) +
    scale_color_manual(values = color_palette) +
    labs(x = "Pseudotime",
         y = clinical_var_label %||% clinical_var_chr,
         color = NULL) +
    theme_minimal(base_size = 15) +
    theme(panel.grid = element_blank())
  
  if (!is.null(bucket)) {
    temp_file <- tempfile(fileext = ".jpeg") # need to create a temporary file
    ggsave(filename = temp_file, width = 7, height = 5)
    s3$upload_file(temp_file, bucket, paste0("slingshot/attempt_arrow_scatter", 
                                             tolower(celltype_suffix), "_", clinical_var, "_slingshot.jpeg"))
  }

  return(p)
}