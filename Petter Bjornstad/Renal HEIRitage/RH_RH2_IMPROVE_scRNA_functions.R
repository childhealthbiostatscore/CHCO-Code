library(Seurat)
library(nebula)
library(Matrix)
library(foreach)
library(doParallel)
library(parallel)
library(dplyr)
library(purrr)
library(aws.s3)
library(jsonlite)

# pb90_subset_clean <- s3readRDS(object = "data_clean/subset/pb90_ckd_analysis_subset.rds", bucket = "scrna", region = "")

# Set up environment for Kopah
user <- Sys.info()[["user"]]

if (user == "choiyej") { # local version
  root_path <- "/Users/choiyej/Library/CloudStorage/OneDrive-UW/Bjornstad/Biostatistics Core Shared Drive"
  git_path <- "/Users/choiyej/GitHub/CHCO-Code/Petter Bjornstad"
  keys <- fromJSON("/Users/choiyej/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Bjornstad Pyle Lab/keys.json")
} else if (user == "rameshsh") { # hyak version
  root_path <- ""
  git_path <- "/mmfs1/gscratch/togo/rameshsh/CHCO-Code/Petter Bjornstad"
} else if (user == "yejichoi") { # hyak version
  root_path <- "/mmfs1/gscratch/togo/yejichoi/"
  git_path <- "/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad"
  keys <- fromJSON("/mmfs1/home/yejichoi/keys.json")
} else if (user == "pylell") {
  root_path <- "/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive"
  git_path <- "/Users/pylell/Documents/GitHub/CHCO-Code/Petter Bjornstad"
  keys <- fromJSON("/mmfs1/home/pylell/keys.json")
} else {
  stop("Unknown user: please specify root path for this user.")
}

Sys.setenv(
  "AWS_ACCESS_KEY_ID" = keys$MY_ACCESS_KEY,
  "AWS_SECRET_ACCESS_KEY" = keys$MY_SECRET_KEY,
  "AWS_DEFAULT_REGION" = "",
  "AWS_REGION" = "",
  "AWS_S3_ENDPOINT" = "s3.kopah.uw.edu"
)

# ===========================================================================
# Function: run_nebula_parallel
# ===========================================================================

run_nebula_parallel <- function(seurat_obj,
                                n_cores = 10,
                                layer = "counts",
                                subject_id_col = "record_id",
                                offset_col = "pooled_offset",
                                formula = ~ group,
                                model = "NBLMM",
                                reml = 1,
                                output_re = TRUE,
                                covariance = TRUE,
                                s3_bucket = "scrna",
                                s3_key = NULL,
                                verbose = TRUE, 
                                group = F) {
  
  # Extract counts and gene list
  counts_mat <- round(GetAssayData(seurat_obj, layer = layer))
  genes_list <- rownames(counts_mat)
  
  # Set up parallel backend
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  # Ensure cleanup on exit
  on.exit({
    stopCluster(cl)
  }, add = TRUE)
  
  start_time <- Sys.time()
  
  # Run nebula in parallel
  nebula_results_list <- foreach(g = genes_list, 
                                 .packages = c("nebula", "Matrix"),
                                 .errorhandling = "pass") %dopar% {
                                   warn <- NULL
                                   err <- NULL
                                   res <- NULL
                                   
                                   tryCatch({
                                     # Subset data for single gene
                                     count_gene <- counts_mat[g, , drop = FALSE]
                                     meta_gene <- subset(seurat_obj, features = g)@meta.data
                                     
                                     # Create model matrix
                                     pred_gene <- model.matrix(formula, data = meta_gene)
                                     
                                     # Group cells
                                     if (group) {
                                       data_g_gene <- group_cell(count = count_gene, 
                                                                 id = meta_gene[[subject_id_col]], 
                                                                 pred = pred_gene)
                                     } else {
                                       data_g_gene <- list(count = count_gene, 
                                                           id = meta_gene[[subject_id_col]], 
                                                           pred = pred_gene)
                                     }
                                     
                                     # Run nebula with warning handling
                                     res <- withCallingHandlers({
                                       nebula(count = data_g_gene$count, 
                                              id = data_g_gene$id, 
                                              pred = data_g_gene$pred,
                                              ncore = 1, 
                                              output_re = output_re, 
                                              covariance = covariance, 
                                              reml = reml, 
                                              model = model, 
                                              offset = meta_gene[[offset_col]])
                                     }, warning = function(w) {
                                       warn <<- conditionMessage(w)
                                       invokeRestart("muffleWarning")
                                     })
                                   }, error = function(e) {
                                     err <<- conditionMessage(e)
                                   })
                                   
                                   list(gene = g, result = res, warning = warn, error = err)
                                 }
  
  # Report warnings and errors if verbose
  if (verbose) {
    for (res in nebula_results_list) {
      if (!is.null(res$warning)) {
        cat(sprintf("Warning for gene %s: %s\n", res$gene, res$warning))
      }
      if (!is.null(res$error)) {
        cat(sprintf("Error for gene %s: %s\n", res$gene, res$error))
      }
    }
  }
  
  # Clean up results
  names(nebula_results_list) <- sapply(nebula_results_list, function(x) x$gene)
  nebula_results_list <- lapply(nebula_results_list, function(x) x$result)
  nebula_results_list <- Filter(Negate(is.null), nebula_results_list)
  
  # Calculate runtime and non-convergence rate
  end_time <- Sys.time()
  runtime <- difftime(end_time, start_time, units = "mins")
  
  nebula_nonconverged_percent <- (length(genes_list) - length(nebula_results_list)) / length(genes_list)
  
  if (verbose) {
    cat(sprintf("\nRuntime: %.2f minutes\n", runtime))
    cat(sprintf("%.2f%% of genes filtered due to low expression or convergence issues\n", 
                nebula_nonconverged_percent * 100))
  }
  
  # Save to S3 if requested
  if (!is.null(s3_key) && !is.null(s3_bucket)) {
    s3saveRDS(nebula_results_list, bucket = s3_bucket, object = s3_key, region = "")
    
    if (verbose) {
      cat(sprintf("Results uploaded to s3://%s/%s\n", s3_bucket, s3_key))
    }
  }
  
  # Return results and metadata
  return(list(
    results = nebula_results_list,
    n_genes_tested = length(genes_list),
    n_genes_converged = length(nebula_results_list),
    nonconverged_percent = nebula_nonconverged_percent,
    runtime_minutes = as.numeric(runtime)
  ))
}

# ===========================================================================
# Function: process_nebula_results
# ===========================================================================

process_nebula_results <- function(nebula_list, 
                                   pval_col = "p_groupType_1_Diabetes", 
                                   logfc_col = "logFC_groupType_1_Diabetes",
                                   convergence_cut = -10,
                                   logfc_cut = 10) {
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
      dplyr::mutate(Gene = gene_name)
  }) %>%
    filter(abs(.data[[logfc_col]]) < logfc_cut)
  
  # Add FDR adjustment
  if (pval_col %in% names(summary_df)) {
    summary_df <- summary_df %>%
      dplyr::mutate(fdr = p.adjust(.data[[pval_col]], method = "fdr")) 
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
# Function: plot_volcano (for transcripts)
# ===========================================================================
plot_volcano <- function(data, fc, p_col, title = NULL, x_axis, y_axis, file_suffix, p_thresh = 0.05,
                         positive_text = "Positive with DKD", 
                         negative_text = "Negative with DKD",
                         formula = "group", legend_position = c(0.8, 0.9),
                         full_formula = F,
                         text_size = 20, caption_size = 8.5,
                         cell_type = "",
                         output_base_path = file.path(root_path, "Renal HERITAGE/Results/Figures/Volcano Plots/"),
                         geom_text_size = 4,
                         arrow_padding = 0.13,
                         arrow_text_padding = 0.18,
                         legend_text_size = 10,
                         caption_padding = 15,
                         x_title_padding_t = 32,
                         genes_to_label= NULL,
                         volcano_force = 6,
                         volcano_box_padding = 0,
                         off_chart_threshold = 0.95,  # New parameter: proportion of y_max to consider "off chart"
                         off_chart_y_position = 0.85,  # New parameter: where to place off-chart labels
                         off_chart_arrow_length = 0.02) {  # New parameter: length of off-chart arrows
  
  set.seed(1)
  
  # Add epsilon for log transformation
  epsilon <- 1e-300
  
  # Calculate -log10(p) with epsilon
  data <- data %>%
    mutate(neg_log_p = -log10(!!sym(p_col) + epsilon))
  
  # Get y-axis max for dynamic scaling
  y_max <- max(data$neg_log_p, na.rm = TRUE) * 1.1
  y_cutoff <- y_max * off_chart_threshold  # Points above this are "off chart"
  
  top_pos <- data %>%
    dplyr::filter(!!sym(fc) > 0 & !!sym(p_col) < p_thresh) %>%
    dplyr::arrange(!!sym(p_col))
  
  n_pos <- nrow(top_pos)
  
  top_pos_n <- top_pos %>%
    filter(if (!is.null(genes_to_label)) Gene %in% genes_to_label else TRUE) %>%
    dplyr::slice_head(n=20)
  
  top_neg <- data %>%
    dplyr::filter(!!sym(fc) < 0 & !!sym(p_col) < p_thresh) %>%
    dplyr::arrange(!!sym(p_col))
  
  n_neg <- nrow(top_neg)
  
  top_neg_n <- top_neg %>%
    filter(if (!is.null(genes_to_label)) Gene %in% genes_to_label else TRUE) %>%
    dplyr::slice_head(n=20)
  
  # Identify off-chart genes
  off_chart_genes <- data %>%
    filter(Gene %in% c(top_pos_n$Gene, top_neg_n$Gene) & neg_log_p > y_cutoff) %>%
    mutate(
      is_positive = !!sym(fc) > 0,
      # Spread out x positions for off-chart labels
      x_position = if_else(is_positive,
                           !!sym(fc) + seq(from = 0.1, by = 0.2, length.out = n()),
                           !!sym(fc) - seq(from = 0.1, by = 0.2, length.out = n())),
      y_position = y_max * off_chart_y_position
    )
  
  # Separate on-chart and off-chart genes for labeling
  on_chart_genes <- c(top_pos_n$Gene, top_neg_n$Gene)[!c(top_pos_n$Gene, top_neg_n$Gene) %in% off_chart_genes$Gene]
  
  data <- data %>%
    dplyr::mutate(
      top_color = case_when(
        Gene %in% top_pos$Gene ~ "#f28482",
        Gene %in% top_neg$Gene ~ "#457b9d",
        TRUE ~ "#ced4da"
      ),
      top_size = if_else(Gene %in% c(top_pos$Gene, top_neg$Gene), 1.3, 1),
      # Only label on-chart genes normally
      top_lab  = if_else(Gene %in% on_chart_genes, Gene, ""),
      # Cap display values at y_cutoff for plotting
      display_neg_log_p = pmin(neg_log_p, y_cutoff)
    ) %>%
    filter(abs(!!sym(fc)) < 10)
  
  # Max and min for annotation arrows
  max_fc <- max(data[[fc]], na.rm = TRUE)
  min_fc <- min(data[[fc]], na.rm = TRUE)
  
  p <- ggplot(data, aes(x = !!sym(fc), y = display_neg_log_p)) +
    geom_hline(yintercept = -log10(p_thresh), linetype = "dashed", color = "darkgrey") +
    geom_point(alpha = 0.5, aes(color = top_color), size = 3) +
    # Regular labels for on-chart genes
    geom_label_repel(seed = 1, 
                     data = filter(data, top_lab != ""),
                     aes(label = top_lab, color = top_color),
                     fontface = "bold",
                     size = geom_text_size, max.overlaps = Inf,
                     force = volcano_force, segment.alpha = 0.3, segment.size = 0.3,
                     box.padding = volcano_box_padding,
                     fill = fill_alpha("white", 0.7),     # background color of the label box
                     label.size = 0
                     # bg.color = "black",
                     # bg.r = .15
    ) +
    # Add arrows for off-chart genes
    {if(nrow(off_chart_genes) > 0) {
      list(
        geom_segment(
          data = off_chart_genes,
          aes(x = !!sym(fc), y = y_cutoff * 0.98,
              xend = !!sym(fc), yend = y_cutoff - (y_max * off_chart_arrow_length)),
          arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
          color = "black",
          size = 0.6
        ),
        geom_text(
          data = off_chart_genes,
          aes(x = x_position, y = y_position, label = Gene),
          size = geom_text_size * 0.9,
          hjust = if_else(off_chart_genes$is_positive, 0, 1),
          color = "black",
          fontface = "italic"
        ),
        # Add p-value annotation for off-chart genes
        geom_text(
          data = off_chart_genes,
          aes(x = x_position, y = y_position - (y_max * 0.03), 
              label = paste0("p=", format(!!sym(p_col), scientific = TRUE, digits = 2))),
          size = geom_text_size * 0.7,
          hjust = if_else(off_chart_genes$is_positive, 0, 1),
          color = "darkgrey"
        )
      )
    }} +
    labs(title = paste(title),
         x = paste(x_axis),
         y = paste(y_axis),
         caption = if (!is.null(cell_type) & !full_formula & !is.null(formula)) {
           paste0("Formula: ~ ", formula, " + (1|subject)", 
                  "\n\nCell type: ", 
                  cell_type, if(cell_type != "") " | " else "", 
                  "Positive n = ", n_pos, " | Negative n = ", n_neg,
                  if(nrow(off_chart_genes) > 0) paste0("\n", nrow(off_chart_genes), " gene(s) with p-values near zero (arrows indicate off-scale values)") else "")
         } else if(!is.null(cell_type) & full_formula) {
           bquote(
             atop(
               Expression[ij] == (beta[0] + b[0*i]) +
                 beta[1]*visit[ij] +
                 beta[2]*treatment[i] +
                 beta[3]*(visit[ij]*treatment[i]) +
                 epsilon[ij],
               "Cell type:" ~ .(cell_type) ~ "|" ~ 
                 "Negative n =" ~ .(n_neg) ~ "|" ~ 
                 "Positive n =" ~ .(n_pos) ~
                 .(if(nrow(off_chart_genes) > 0) paste0(" | ", nrow(off_chart_genes), " off-scale") else "")
             )
           )
         } else {
           paste0("\nCell type: ", 
                  cell_type, if(cell_type != "") " | " else "", 
                  "Positive n = ", n_pos, " | Negative n = ", n_neg,
                  if(nrow(off_chart_genes) > 0) paste0("\n", nrow(off_chart_genes), " gene(s) with p-values near zero (arrows indicate off-scale values)") else "")
         }) +
    scale_size_continuous(range = c(1, 1.3)) + 
    scale_color_manual(values = c("#457b9d"="#457b9d", "#ced4da"="#ced4da", "#f28482"="#f28482")) +
    # theme_minimal() +
    guides(color = "none", size = "none")  +
    annotate("segment", 
             x=max_fc/8, 
             xend=(max_fc*7)/8, 
             y=-y_max * arrow_padding,
             col="darkgrey", arrow=arrow(length=unit(0.2, "cm"))) +
    annotate("text", 
             x=mean(c(max_fc/8, (max_fc*7)/8)), 
             y=-y_max * arrow_text_padding, 
             label=positive_text,
             size=geom_text_size, color="#343a40") +
    annotate("segment", 
             x=min_fc/8, 
             xend=(min_fc*7)/8, 
             y=-y_max * arrow_padding,
             col="darkgrey", arrow=arrow(length=unit(0.2, "cm"))) +
    annotate("text", 
             x=mean(c(min_fc/8, (min_fc*7)/8)), 
             y=-y_max * arrow_text_padding, 
             label=negative_text,
             size=geom_text_size, color="#343a40") +
    scale_y_continuous(expand=c(0,0)) +
    coord_cartesian(ylim = c(0, y_max), clip="off") +
    theme(legend.title = element_blank(),
          panel.grid = element_blank(),
          text = element_text(size = text_size, family = "Arial"),
          title = element_text(size = legend_text_size, family = "Arial"),
          legend.position = legend_position,
          legend.justification = if(is.character(legend_position) && legend_position == "bottom") c(0.5, 0) else c(1, 1),
          legend.direction = "horizontal",
          legend.spacing.x = unit(0.05, 'cm'),
          plot.margin = margin(t = 10, r = 20, b = 25, l = 20),
          axis.title.x = element_text(margin = margin(t = x_title_padding_t)),
          plot.caption = element_text(size = caption_size, hjust = 0.5, margin = margin(t = caption_padding)),
          legend.margin = margin(t = 5, b = 5),
          panel.background = element_blank(),
          legend.background = element_blank())
  
  if (!is.null(output_base_path)) {
    ggsave(paste0(output_base_path, file_suffix, ".png"), bg = "transparent", plot = p, width = 7, height = 5)
  }
  return(p)
}

make_overlap_table <- function(
    modA, modB, modC, modD, modE,
    sig = c("fdr", "p"), alpha = 0.05
) {
  sig <- match.arg(sig)
  
  # map each model to its code + label; allow NULL (e.g., no C)
  model_map <- list(
    A = list(df = modA, label = "DKDn_HC"),
    B = list(df = modB, label = "DKD_HC"),
    C = list(df = modC, label = "DKDy_DKDn"),
    D = list(df = modD, label = "GLPn_HC"),
    E = list(df = modE, label = "GLPy_GLPn")
  )
  
  # bind to long format with code + label so we can control wide names
  long <- imap_dfr(model_map, function(m, code) {
    if (is.null(m$df)) return(NULL)
    m$df %>%
      dplyr::select(Gene, logFC, p.value, fdr) %>%
      # keep the best row per Gene within each model (ties: lowest FDR, then largest |logFC|)
      group_by(Gene) %>%
      arrange(fdr, desc(abs(logFC)), p.value, .by_group = TRUE) %>%
      dplyr::slice(1) %>%
      ungroup() %>%
      mutate(code = code, label = m$label)
  })
  
  # apply significance filter
  long_filt <- long %>%
    filter(if (sig == "fdr") fdr < alpha else p.value < alpha)
  
  # pivot to wide with desired name pattern: A_logFC_DKDn_HC, A_p_DKDn_HC, A_fdr_DKDn_HC
  wide <- long_filt %>%
    pivot_wider(
      id_cols = Gene,
      names_from = c(code, label),
      values_from = c(logFC, p.value, fdr),
      names_glue = "{code}_{.value}_{label}",
      values_fill = NA
    )
  
  # figure out which logFC columns exist (some models may be absent)
  logfc_cols <- grep("^[A-E]_logFC_", names(wide), value = TRUE)
  # count how many models each gene shows up in (non-missing logFC)
  wide <- wide %>%
    mutate(
      n_overlap = if (length(logfc_cols)) rowSums(!is.na(pick(all_of(logfc_cols)))) else 0
    )
  
  # build overlaps string like "A+B-D-..." in A→E order, skipping missing columns
  # helper to get sign tag per code if present
  code_order <- c("A","B","C","D","E")
  tag_df <- map_dfc(code_order, function(code) {
    col <- grep(paste0("^", code, "_logFC_"), names(wide), value = TRUE)
    if (length(col) == 0) {
      rep("", nrow(wide)) |> as.data.frame() |> setNames(code)
    } else {
      x <- wide[[col]]
      tag <- ifelse(is.na(x), "", ifelse(x > 0, paste0(code, "+"), paste0(code, "-")))
      as.data.frame(tag) |> setNames(code)
    }
  })
  
  overlaps_chr <- do.call(paste0, as.list(tag_df))
  wide$overlaps <- overlaps_chr
  
  # optional: put overlaps & n_overlap up front
  cols_front <- c("Gene", "n_overlap", "overlaps")
  wide %>%
    relocate(any_of(cols_front))
}


# ===========================================================================
# Function: prepare_umap_metadata
# ===========================================================================
prepare_umap_metadata <- function(
    seurat_obj,
    genes,
    celltype_col = "KPMP_celltype",
    umap_reduction = "umap.harmony",
    color_palette = "Set3",        # RColorBrewer palette name
    custom_colors = NULL           # named vector to override palette entirely
) {
  
  # Determine cell type levels
  celltypes <- levels(seurat_obj@meta.data[[celltype_col]])
  n_celltypes <- length(celltypes)
  
  # Build color mapping
  if (!is.null(custom_colors)) {
    # Use provided named color vector as-is
    celltype_colors <- custom_colors
  } else {
    # Generate from RColorBrewer palette, expanding as needed
    umap_colors <- colorRampPalette(brewer.pal(min(brewer.pal.info[color_palette, "maxcolors"], n_celltypes), color_palette))(n_celltypes)
    celltype_colors <- setNames(umap_colors, celltypes)
  }
  
  # Fetch gene expression
  expr_df <- Seurat::FetchData(seurat_obj, vars = genes) %>%
    as.data.frame()
  colnames(expr_df) <- genes  # ensure clean names
  
  # Build metadata with UMAP embeddings and expression
  umap_embed_cols <- colnames(seurat_obj@reductions[[umap_reduction]]@cell.embeddings)
  metadata <- seurat_obj@meta.data %>%
    cbind(seurat_obj@reductions[[umap_reduction]]@cell.embeddings) %>%
    cbind(expr_df)
  colnames(metadata)[colnames(metadata) %in% umap_embed_cols] <- c("umapharmony_1", "umapharmony_2")
  
  # Compute cell type centers
  centers <- metadata %>%
    group_by(.data[[celltype_col]]) %>%
    summarise(x = median(umapharmony_1),
              y = median(umapharmony_2), .groups = "drop")
  
  list(
    metadata = metadata,
    centers = centers,
    celltype_colors = celltype_colors
  )
}

# ===========================================================================
# Function: plot_feature_umap
# ===========================================================================
plot_feature_umap <- function(
    genes,
    metadata,          # pb90_subset_meta equivalent
    centers,           # UMAP label centers df with x, y, KPMP_celltype
    celltype_colors,   # named vector of colors
    pct_threshold = 10,  # minimum % expressed to highlight
    save_to_s3 = FALSE,
    s3_prefix = "Projects/CKD/RH_RH2/Results/Figures/UMAP/REMODEL genes/",
    s3_suffix = "",
    s3_bucket = "scrna",
    plot_height = 10,
    plot_width = 10
) {
  
  plots <- list()
  
  for (gene in genes) {
    
    # Calculate expression stats per cell type
    expr_stats <- metadata %>%
      group_by(KPMP_celltype) %>%
      summarise(
        pct_expressed = round(mean(.data[[gene]] > 0) * 100, 1),
        n_total = n(),
        .groups = "drop"
      ) %>%
      filter(pct_expressed > pct_threshold)
    
    # Tag centers by whether they pass expression threshold
    centers_filtered <- centers %>%
      mutate(pct_expressed_5 = case_when(
        KPMP_celltype %in% expr_stats$KPMP_celltype ~ "Y",
        TRUE ~ "N"
      ))
    
    # Build caption
    caption_text <- expr_stats %>%
      arrange(desc(pct_expressed)) %>%
      mutate(
        color = celltype_colors[KPMP_celltype],
        label = paste0("<span style='color:", color, "'>&#9679; </span>**", KPMP_celltype, "**: ", pct_expressed, "%"),
        group = ceiling(row_number() / 5)
      ) %>%
      group_by(group) %>%
      summarise(line = paste(label, collapse = " | "), .groups = "drop") %>%
      rbind(data.frame(group = 0, line = paste0("**% Expressed (>", pct_threshold, "%):**"))) %>%
      arrange(group) %>%
      pull(line) %>%
      paste(collapse = "<br><br>")
    
    # Build plot
    feature_p <- metadata %>%
      dplyr::mutate(feature_color_logic = case_when(
        .data[[gene]] > 0 ~ KPMP_celltype,
        TRUE ~ "No"
      )) %>%
      ggplot(aes(x = umapharmony_1, y = umapharmony_2, color = feature_color_logic)) +
      geom_point(alpha = 0.4) +
      geom_text_repel(
        data = centers_filtered,
        aes(x = x, y = y, label = KPMP_celltype),
        inherit.aes = FALSE,
        size = 4,
        fontface = ifelse(centers_filtered$pct_expressed_5 == "Y", "bold", "plain"),
        color = ifelse(centers_filtered$pct_expressed_5 == "Y", "black", "#495057"),
        max.overlaps = Inf
      ) +
      scale_color_manual(values = c(celltype_colors, "No" = "#edede9")) +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        text = element_text(size = 15),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.caption = element_markdown(size = 12, hjust = 0.5),
        axis.title.x = element_text(hjust = 0),
        axis.title.y = element_text(hjust = 0)
      ) +
      labs(x = "UMAP1", y = "UMAP2", title = gene, caption = caption_text)
    
    print(feature_p)
    plots[[gene]] <- feature_p
    
    if (save_to_s3) {
      s3write_using_region(
        feature_p,
        FUN = ggsave,
        object = paste0(s3_prefix, gene, "_FeaturePlot", s3_suffix, ".png"),
        bucket = s3_bucket,
        region = "",
        height = plot_height,
        width = plot_width
      )
    }
  }
  
  invisible(plots)  # Return list of plots silently
}
