


########################
# specify user for paths
########################

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
} else if (user == "leidholt") {
    root_path <- c("/mmfs1/gscratch/togo/leidholt/")
    git_path <- "/mmfs1/gscratch/togo/leidholt/CHCO-Code/"
    keys <- fromJSON("/mmfs1/home/leidholt/keys.json")
} else if (user == "sleidholt") {
    root_path <- c("/Users/Shared/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/")
    dir.results <- c("/Users/Shared/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/T1D Adiposity/Results/Compiled Results/Proteomics Results")
    git_path <- "/Users/sleidholt/Documents/GitHub/CHCO-Code/Petter Bjornstad/"
    keys <- fromJSON("/Users/Shared/OneDrive - UW/Personal_storage/keys.json")
} else if (user == "savanahleidholt") {
    root_path <- c("/Users/savanahleidholt/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/")
    dir.results <- c("/Users/savanahleidholt/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/T1D Adiposity/Results/Compiled Results/Proteomics Results")
    git_path <- "/Users/savanahleidholt/Desktop/CHCO-Code/Petter Bjornstad"
    keys <- fromJSON("/Users/savanahleidholt/OneDrive - UW/Personal_storage/keys.json")
} else {
    stop("Unknown user: please specify root path for this user.")
}


#---------------------------------------------------------
#Function for plotting clinical data of interest
#---------------------------------------------------------
plot_clinical_boxplot <- function(data, outcome_var, outcome_label) {
  plot_df <- data %>%
    dplyr::select(group_bmi, age, sex, study, dplyr::all_of(outcome_var)) %>%
    dplyr::mutate(
      y = as.numeric(.data[[outcome_var]]),
      group_bmi = factor(group_bmi),
      sex = factor(sex),
      study = factor(study)
    ) %>%
    tidyr::drop_na(y, group_bmi)
  
  ggplot2::ggplot(plot_df,
    ggplot2::aes(
      x = group_bmi,
      y = y,
      fill = group_bmi
    )
  ) +
    ggplot2::geom_boxplot(
      width = 0.65,
      alpha = 0.75,
      outlier.shape = NA
    ) +
    ggplot2::geom_jitter(
      width = 0.15,
      size = 2.5,
      alpha = 0.75
    ) +
    ggplot2::labs(
      title = outcome_label,
      x = NULL,
      y = outcome_label
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "LC_Normal" = "cyan4",
        "LC_Overweight_Obese" = "darkred",
        "T1D_Normal" = "gold",
        "T1D_Overweight" = "plum3",
        "T1D_Obese" = "black"
      )
    ) +
    ggplot2::theme_bw(base_size = 18) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      ),
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5
      ),
      panel.grid = ggplot2::element_blank()
    )
}

#-------------------------------------
#PCA plotting function for group_bmi
#-------------------------------------
make_pca <- function(mat, meta_df, title_str, add_ellipses = TRUE) {
  stopifnot(identical(rownames(mat), meta_df$record_id))
  
  pca <- prcomp(mat, center = TRUE, scale. = TRUE)
  vexp <- round(100 * (pca$sdev^2) / sum(pca$sdev^2), 1)
  
  pca_df <- tibble(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    group_bmi = meta_df$group_bmi,
    study = meta_df$study
  )
  
  p <- ggplot(
    pca_df,
    aes(
      x = PC1,
      y = PC2,
      color = group_bmi
    )
  ) +
    {
      if (add_ellipses) {
        stat_ellipse(
          aes(group = group_bmi, color = group_bmi),
          type = "norm",
          level = 0.68,
          linewidth = 1,
          linetype = "solid",
          alpha = 0.8
        )
      }
    } +
    geom_point(aes(shape = study), size = 3, alpha = 0.85) +
    scale_color_manual(
      values = c(
        LC_Normal  = "cyan4",
        LC_Overweight_Obese = "darkred",
        T1D_Normal = "gold",
        T1D_Overweight = "plum3",
        T1D_Obese = "black"
      ),
      name = "Disease-adiposity group"
    ) +
    labs(
      title = title_str,
      x = paste0("PC1 (", vexp[1], "%)"),
      y = paste0("PC2 (", vexp[2], "%)"),
      color = "Disease-adiposity group",
      shape = "Study"
    ) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  return(p)
}

#----------------------
#all T1D limma function
#----------------------
run_t1d_global_adiposity_limma <- function(mat, meta_df, design_formula, model_label) {
  
  design <- model.matrix(design_formula, data = meta_df)
  colnames(design) <- make.names(colnames(design))
  
  message("Design columns for ", model_label, ":")
  print(colnames(design))
  
  fit <- lmFit(mat, design)
  fit <- eBayes(fit)
  
  group_cols <- grep("^group_bmi", colnames(design), value = TRUE)
  
  topTable(
    fit,
    coef = group_cols,
    number = Inf,
    adjust.method = "BH",
    sort.by = "F"
  ) %>%
    rownames_to_column("protein_id") %>%
    mutate(
      model = model_label,
      contrast = "T1D_global_adiposity",
      direction = case_when(
        adj.P.Val < 0.05 ~ "Global adiposity-associated difference",
        TRUE ~ "Not significant"
      )
    )
}


#-------------------------------------------------------
#Limma function for all differences & pairwise comparisons
#--------------------------------------------------------
run_group_bmi_limma <- function(mat, meta_df, design_formula, model_label) {
  
  design <- model.matrix(design_formula, data = meta_df)
  colnames(design) <- make.names(colnames(design))
  
  message("Design columns for ", model_label, ":")
  print(colnames(design))
  
  fit <- lmFit(mat, design)
  
  
  #global f-test
  fit_global <- eBayes(fit)
  
  group_cols <- grep("^group_bmi", colnames(design))
  
  global_results <- topTable(
    fit_global,
    coef = group_cols,
    number = Inf,
    adjust.method = "BH",
    sort.by = "F"
  ) %>%
    rownames_to_column("protein_id") %>%
    mutate(
      model = model_label,
      contrast = "Global_group_bmi",
      direction = case_when(
        adj.P.Val < 0.05 ~ "Globally significant",
        TRUE ~ "Not significant"
      )
    )
  
  
  #pairwise comparison
  contrast_matrix <- makeContrasts(
    T1DNormal_vs_LCNormal =
      group_bmiT1D_Normal - group_bmiLC_Normal,
    
    T1DOverObese_vs_LCOverObese =
      ((group_bmiT1D_Overweight + group_bmiT1D_Obese) / 2) -
      group_bmiLC_Overweight_Obese,
    
    T1DNormal_vs_T1DOverweight =
      group_bmiT1D_Normal - group_bmiT1D_Overweight,
    
    T1DNormal_vs_T1DObese =
      group_bmiT1D_Normal - group_bmiT1D_Obese,
    
    T1DObese_vs_T1DOverweight =
      group_bmiT1D_Obese - group_bmiT1D_Overweight,
    
    levels = design
  )
  
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  
  pairwise_results <- lapply(colnames(contrast_matrix), function(cn) {
    topTable(
      fit2,
      coef = cn,
      number = Inf,
      adjust.method = "BH"
    ) %>%
      rownames_to_column("protein_id") %>%
      mutate(
        model = model_label,
        contrast = cn,
        direction = case_when(
          adj.P.Val < 0.05 & logFC > 0 ~ "Higher in first group",
          adj.P.Val < 0.05 & logFC < 0 ~ "Higher in second group",
          TRUE ~ "Not significant"
        )
      )
  }) %>%
    bind_rows()
  
  bind_rows(global_results, pairwise_results)
}

#-----------------------------------------------
#heatmap function for all proteomic differences
#----------------------------------------------
make_global_group_bmi_heatmap <- function(results_df,
                                          expr_mat,
                                          meta_df,
                                          model_name,
                                          sample_id_col,
                                          padj_cutoff = 0.05,
                                          top_n = 50) {
  
  sig_prots <- results_df %>%
    filter(
      model == model_name,
      contrast == "Global_group_bmi",
      adj.P.Val < padj_cutoff
    ) %>%
    arrange(adj.P.Val) %>%
    slice_head(n = top_n) %>%
    pull(protein_id)
  
  heat_mat <- expr_mat[sig_prots, , drop = FALSE]
  
  heat_mat <- t(scale(t(heat_mat)))
  heat_mat <- heat_mat[complete.cases(heat_mat), , drop = FALSE]
  
  annot_df <- meta_df %>%
    dplyr::select(all_of(sample_id_col), group_bmi) %>%
    mutate(
      group_bmi = factor(
        group_bmi,
        levels = c(
          "LC_Normal",
          "LC_Overweight_Obese",
          "T1D_Normal",
          "T1D_Overweight",
          "T1D_Obese"
        )
      )
    ) %>%
    as.data.frame()
  
  rownames(annot_df) <- annot_df[[sample_id_col]]
  annot_df[[sample_id_col]] <- NULL
  
  shared_samples <- intersect(colnames(heat_mat), rownames(annot_df))
  
  heat_mat <- heat_mat[, shared_samples, drop = FALSE]
  annot_df <- annot_df[shared_samples, , drop = FALSE]
  
  stopifnot(ncol(heat_mat) > 0)
  stopifnot(nrow(annot_df) > 0)
  stopifnot(all(colnames(heat_mat) == rownames(annot_df)))
  
  row_labels <- map_tbl %>%
    filter(protein_id %in% rownames(heat_mat)) %>%
    distinct(protein_id, .keep_all = TRUE)
  
  gene_labels <- row_labels$gene_symbol
  names(gene_labels) <- row_labels$protein_id
  
  rownames(heat_mat) <- ifelse(
    is.na(gene_labels[rownames(heat_mat)]),
    rownames(heat_mat),
    gene_labels[rownames(heat_mat)]
  )
  
  group_colors <- setNames(
    viridisLite::viridis(5),
    c(
      "LC_Normal",
      "LC_Overweight_Obese",
      "T1D_Normal",
      "T1D_Overweight",
      "T1D_Obese"
    )
  )
  
  ann_colors <- list(
    group_bmi = group_colors
  )
  
  pheatmap::pheatmap(
    heat_mat,
    annotation_col = annot_df,
    annotation_colors = ann_colors,
    show_colnames = FALSE,
    fontsize_row = 9,
    fontsize_col = 10,
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    border_color = NA,
    main = paste0(
      "Top ",
      nrow(heat_mat),
      " globally different proteins across group-BMI strata"
    )
  )
}

#-------------------------------------
#volcano plots for DEPs
#-------------------------------------
make_volcano <- function(results_df,
                         contrast_name,
                         title_prefix = NULL,
                         p_cutoff = 0.05,
                         fc_cutoff = 0.5,
                         label_n = 15) {
  
  df <- results_df %>%
    left_join(map_tbl, by = "protein_id") %>%
    
    mutate(
      gene_symbol = ifelse(
        is.na(gene_symbol) | gene_symbol == "",
        protein_id,
        gene_symbol
      )
    ) %>%
    dplyr::filter(contrast == contrast_name) %>%
    dplyr::mutate(
      neg_log10_adj_p = -log10(adj.P.Val),
      sig_group = dplyr::case_when(
        adj.P.Val < p_cutoff & logFC > fc_cutoff  ~ "Up",
        adj.P.Val < p_cutoff & logFC < -fc_cutoff ~ "Down",
        TRUE ~ "Not significant"
      )
    )
  
  n_up <- sum(df$sig_group == "Up", na.rm = TRUE)
  n_down <- sum(df$sig_group == "Down", na.rm = TRUE)
  n_ns <- sum(df$sig_group == "Not significant", na.rm = TRUE)
  
  legend_labs <- c(
    "Up" = paste0("Upregulated (n = ", n_up, ")"),
    "Down" = paste0("Downregulated (n = ", n_down, ")"),
    "Not significant" = paste0("Not significant (n = ", n_ns, ")")
  )
  
  top_labels <- df %>%
    dplyr::filter(sig_group %in% c("Up", "Down")) %>%
    dplyr::arrange(dplyr::desc(abs(t))) %>%
    dplyr::distinct(gene_symbol, .keep_all = TRUE) %>%
    dplyr::slice_head(n = label_n)
  
  ggplot(df, aes(x = logFC, y = neg_log10_adj_p)) +
    
    geom_point(aes(color = sig_group), alpha = 0.75, size = 2) +
    
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "solid", alpha = 0.5) +
    
    geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed") +
    
    ggrepel::geom_text_repel(
      data = top_labels,
      aes(label = gene_symbol),
      size = 3.5,
      max.overlaps = Inf
    ) +
    
    scale_color_manual(
      values = c(
        "Up" = "orange3",
        "Down" = "purple4",
        "Not significant" = "grey70"
      ),
      labels = legend_labs
    ) +
    
    coord_cartesian(xlim = c(-2, 3)) +  
    
    scale_x_continuous(breaks = seq(-2, 3, by = 1)) +
    
    labs(
      title = paste0(title_prefix, "\n", contrast_name),
      x = "log2 Fold Change",
      y = "-log10 adjusted P-value",
      color = NULL
    ) +
    
    theme_minimal(base_size = 22) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", linewidth = 1),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(
        face = "bold",
        size = 12
      ),
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12)
    )
}

#---------------------------------------
#making a list of volcano plots/contrast
#---------------------------------------

make_volcano_list <- function(results_df, model_name, label_n = 15) {
  
  contrasts <- unique(results_df$contrast)
  
  plot_list <- lapply(contrasts, function(ct) {
    make_volcano(
      results_df = results_df,
      contrast_name = ct,
      title_prefix = model_name,
      label_n = label_n
    )
  })
  
  names(plot_list) <- contrasts
  
  plot_list
}

#--------------
#fgsea funciton
#--------------
run_fgsea_model <- function(res_df,
                            model_name,
                            pathways,
                            pathway_type,
                            protein_gene_map,
                            min_size = 10,
                            max_size = 500,
                            eps = 0,
                            nPermSimple = 10000) {
  
  res_annot <- res_df %>%
    dplyr::left_join(protein_gene_map, by = "protein_id") %>%
    dplyr::filter(
      !is.na(gene_symbol),
      gene_symbol != "",
      !is.na(t)
    )
  
  fgsea_out <- res_annot %>%
    dplyr::group_by(contrast) %>%
    dplyr::group_split() %>%
    purrr::map_dfr(function(df) {
      
      contrast_name <- unique(df$contrast)
      
      ranks <- df %>%
        dplyr::group_by(gene_symbol) %>%
        dplyr::slice_max(order_by = abs(t), n = 1, with_ties = FALSE) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(desc(t))
      
      ranks_vec <- ranks$t
      names(ranks_vec) <- ranks$gene_symbol
      
      fgsea::fgsea(
        pathways = pathways,
        stats = ranks_vec,
        minSize = min_size,
        maxSize = max_size
      ) %>%
        tibble::as_tibble() %>%
        dplyr::mutate(
          model = model_name,
          contrast = contrast_name,
          pathway_type = pathway_type
        )
    })
  
  fgsea_out
}


#----------------------------------
#plotting function for GSEA pathways
#----------------------------------
plot_fgsea_top20 <- function(fgsea_df,
                             model_name,
                             contrast_name,
                             pathway_type_name = NULL,
                             top_n = 20,
                             padj_cutoff = 0.05) {
  
  df <- fgsea_df %>%
    filter(
      model == model_name,
      contrast == contrast_name
    )
  
  if (!is.null(pathway_type_name)) {
    df <- df %>%
      filter(pathway_type == pathway_type_name)
  }
  
  df_top <- df %>%
    filter(!is.na(padj), !is.na(NES)) %>%
    arrange(padj) %>%
    slice_head(n = top_n) %>%
    mutate(
      pathway_clean = pathway %>%
        str_remove("^REACTOME_") %>%
        str_replace_all("_", " ") %>%
        str_replace("^GO ", "") %>%
        str_to_sentence() %>%
        stringr::str_wrap(width = 60),
      
      pathway_clean = forcats::fct_reorder(pathway_clean, NES),
      
      direction = ifelse(NES > 0, "Up", "Down"),
      
      # transform FDR so smaller = bigger
      fdr_size = -log10(padj)
    )
  
  ggplot(df_top, aes(y = pathway_clean)) +
    geom_segment(
      aes(x = 0, xend = NES, yend = pathway_clean),
      color = "grey70",
      linewidth = 1
    ) +
    
    geom_point(
      aes(x = NES, size = fdr_size, color = direction)
    ) +
    
    geom_text(
      aes(x = 0, label = pathway_clean),
      hjust = ifelse(df_top$NES > 0, 1, 0),
      size = 5
    ) +
    
    geom_vline(xintercept = 0, linetype = "dashed") +
    
    scale_color_manual(
      values = c(
        "Up" = "orange3",
        "Down" = "purple4"
      ),
      name = "Direction"
    ) +
    
    scale_size_continuous(
      name = "-log10(FDR)",
      range = c(6, 14)
    ) +
    
    scale_x_continuous(
      limits = c(-4, 4),
      breaks = seq(-4, 4, by = 2)
    ) +
    
    labs(
      title = paste(model_name, contrast_name, pathway_type_name, sep = "\n"),
      x = "Normalized enrichment score",
      y = NULL
    ) +
    
    theme_minimal(base_size = 22) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", linewidth = 1),
      
      axis.text.x = element_text(size = 14),
      axis.title.x = element_text(size = 14, face = "bold"),
      
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      
      plot.title = element_text(
        face = "bold",
        size = 11
      ),
      
      legend.position = "right",
      
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14)
    )
}

#----------------------------------------
#ElasticNet/LASSO Regression
#----------------------------------------
run_elastic_net <- function(
    contrast_name,
    limma_df,
    meta_df,
    X_mat,
    protein2symbol = NULL,
    group_var = "group_bmi",
    top_n = 100,
    sig_only = TRUE,
    p_cutoff = 0.05,
    fc_cutoff = 0.5,
    covars = c("age", "sex", "study"),
    alpha_grid = seq(0.1, 1, by = 0.1),
    nfolds = 10,
    seed = 1835,
    min_per_class = 5,
    maxit = 1e6
) {
  
  cat("\n============================\n")
  cat("Running contrast:", contrast_name, "\n")
  
  contrast_groups <- list(
    T1D_Normal_vs_LC_Normal = list(
      positive = "T1D_Normal",
      negative = "LC_Normal"
    ),
    T1D_OverObese_vs_LC_OverObese = list(
      positive = c("T1D_Overweight", "T1D_Obese"),
      negative = "LC_Overweight_Obese"
    ),
    T1D_Normal_vs_T1D_Overweight = list(
      positive = "T1D_Normal",
      negative = "T1D_Overweight"
    ),
    T1D_Normal_vs_T1D_Obese = list(
      positive = "T1D_Normal",
      negative = "T1D_Obese"
    ),
    T1D_Obese_vs_T1D_Overweight = list(
      positive = "T1D_Obese",
      negative = "T1D_Overweight"
    )
  )
  
  if (!contrast_name %in% names(contrast_groups)) {
    stop("Contrast not found in contrast_groups: ", contrast_name)
  }
  
  g_pos <- contrast_groups[[contrast_name]]$positive
  g_neg <- contrast_groups[[contrast_name]]$negative
  
  meta_sub <- meta_df %>%
    dplyr::filter(.data[[group_var]] %in% c(g_pos, g_neg)) %>%
    dplyr::mutate(
      record_id = as.character(record_id),
      y = ifelse(.data[[group_var]] %in% g_pos, 1, 0),
      sex = factor(sex),
      study = factor(study)
    )
  
  if (nrow(meta_sub) == 0) {
    stop("No samples found for contrast: ", contrast_name)
  }
  
  tab <- table(meta_sub$y)
  print(tab)
  
  if (any(tab < min_per_class)) {
    stop(
      "Too few samples in at least one class for contrast: ",
      contrast_name,
      " | counts = ",
      paste(names(tab), tab, collapse = ", ")
    )
  }
  
  ids <- intersect(meta_sub$record_id, rownames(X_mat))
  
  meta_sub <- meta_sub %>%
    dplyr::filter(record_id %in% ids)
  
  meta_sub <- meta_sub[match(ids, meta_sub$record_id), , drop = FALSE]
  X_mat_sub <- X_mat[ids, , drop = FALSE]
  
  stopifnot(identical(rownames(X_mat_sub), meta_sub$record_id))
  
  limma_sub <- limma_df %>%
    dplyr::filter(contrast == contrast_name)
  
  if (nrow(limma_sub) == 0) {
    stop("No limma rows found for contrast: ", contrast_name)
  }
  
  if (sig_only) {
    limma_sub <- limma_sub %>%
      dplyr::filter(
        adj.P.Val < p_cutoff,
        abs(logFC) > fc_cutoff
      )
  }
  
  if (nrow(limma_sub) == 0) {
    stop("No significant limma DEPs found for contrast: ", contrast_name)
  }
  
  if (!is.null(protein2symbol)) {
    limma_sub <- limma_sub %>%
      dplyr::left_join(protein2symbol, by = "protein_id") %>%
      dplyr::mutate(
        gene_symbol = dplyr::if_else(
          is.na(gene_symbol) | gene_symbol == "",
          protein_id,
          gene_symbol
        )
      ) %>%
      dplyr::filter(protein_id %in% colnames(X_mat_sub)) %>%
      dplyr::arrange(dplyr::desc(abs(t)), adj.P.Val, P.Value) %>%
      dplyr::group_by(gene_symbol) %>%
      dplyr::slice_head(n = 1) %>%
      dplyr::ungroup()
  } else {
    limma_sub <- limma_sub %>%
      dplyr::filter(protein_id %in% colnames(X_mat_sub)) %>%
      dplyr::arrange(dplyr::desc(abs(t)), adj.P.Val, P.Value) %>%
      dplyr::mutate(gene_symbol = protein_id)
  }
  
  top_feats_tbl <- limma_sub %>%
    dplyr::arrange(dplyr::desc(abs(t)), adj.P.Val, P.Value) %>%
    dplyr::slice_head(n = top_n)
  
  top_feats <- top_feats_tbl$protein_id
  
  cat("Input limma DEPs after redundancy filtering:", nrow(limma_sub), "\n")
  cat("Top proteins used in elastic net:", length(top_feats), "\n")
  
  if (length(top_feats) < 5) {
    stop("Fewer than 5 overlapping proteins for contrast: ", contrast_name)
  }
  
  covars_present <- covars[covars %in% colnames(meta_sub)]
  
  cat("Covariates used:", paste(covars_present, collapse = ", "), "\n")
  
  if (length(covars_present) > 0) {
    Z <- model.matrix(
      as.formula(paste("~", paste(covars_present, collapse = " + "))),
      data = meta_sub
    )
  } else {
    Z <- matrix(1, nrow = nrow(meta_sub), ncol = 1)
    colnames(Z) <- "(Intercept)"
  }
  
  Xp <- X_mat_sub[, top_feats, drop = FALSE]
  
  X <- cbind(Z, Xp)
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  
  y <- as.numeric(meta_sub$y)
  
  penalty_factor <- c(rep(0, ncol(Z)), rep(1, ncol(Xp)))
  
  sds <- apply(X, 2, sd)
  keep <- is.finite(sds) & sds > 0
  
  X <- X[, keep, drop = FALSE]
  penalty_factor <- penalty_factor[keep]
  
  cat("Final elastic net matrix:", nrow(X), "samples x", ncol(X), "features\n")
  
  nfolds_use <- min(nfolds, sum(y == 0), sum(y == 1))
  
  cat("nfolds used:", nfolds_use, "\n")
  
  if (nfolds_use < 2) {
    stop("Not enough samples per class for cross-validation.")
  }
  
  set.seed(seed)
  
  fit_list <- lapply(alpha_grid, function(a) {
    cv <- tryCatch(
      glmnet::cv.glmnet(
        x = X,
        y = y,
        family = "binomial",
        alpha = a,
        nfolds = nfolds_use,
        type.measure = "deviance",
        penalty.factor = penalty_factor,
        standardize = TRUE,
        maxit = maxit
      ),
      error = function(e) {
        message("glmnet failed at alpha = ", a, ": ", conditionMessage(e))
        NULL
      }
    )
    
    if (is.null(cv)) return(NULL)
    
    list(alpha = a, cvfit = cv)
  })
  
  fit_list <- Filter(Negate(is.null), fit_list)
  
  if (length(fit_list) == 0) {
    stop("All glmnet fits failed.")
  }
  
  best_idx <- which.min(sapply(fit_list, function(z) min(z$cvfit$cvm)))
  best <- fit_list[[best_idx]]
  
  coef_mat <- as.matrix(coef(best$cvfit, s = "lambda.1se"))
  
  coef_df <- tibble::tibble(
    protein_id = rownames(coef_mat),
    beta = as.numeric(coef_mat)
  ) %>%
    dplyr::filter(beta != 0)
  
  selected_proteins <- intersect(top_feats, coef_df$protein_id)
  
  coef_df_annot <- coef_df %>%
    dplyr::left_join(
      top_feats_tbl %>%
        dplyr::select(protein_id, gene_symbol, logFC, t, P.Value, adj.P.Val),
      by = "protein_id"
    )
  
  selected_proteins_annot <- top_feats_tbl %>%
    dplyr::filter(protein_id %in% selected_proteins) %>%
    dplyr::select(protein_id, gene_symbol, logFC, t, P.Value, adj.P.Val)
  
  p_hat <- as.numeric(
    predict(best$cvfit, newx = X, s = "lambda.1se", type = "response")
  )
  
  auc_val <- as.numeric(
    pROC::auc(pROC::roc(y, p_hat, quiet = TRUE))
  )
  
  cat("Best alpha:", best$alpha, "\n")
  cat("Training AUC:", auc_val, "\n")
  cat("Selected proteins:", length(selected_proteins), "\n")
  
  list(
    contrast = contrast_name,
    positive_group = g_pos,
    negative_group = g_neg,
    best_alpha = best$alpha,
    lambda_1se = best$cvfit$lambda.1se,
    lambda_min = best$cvfit$lambda.min,
    auc_train = auc_val,
    n_input_proteins = length(top_feats),
    top_feats = top_feats,
    top_feats_tbl = top_feats_tbl,
    selected_proteins = selected_proteins,
    selected_proteins_annot = selected_proteins_annot,
    coef_df = coef_df,
    coef_df_annot = coef_df_annot,
    cvfit = best$cvfit
  )
}

