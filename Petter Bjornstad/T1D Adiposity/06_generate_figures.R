#!/usr/bin/env Rscript
# =============================================================================
# Generate volcano plots and GSEA pathway plots for all NEBULA results
# =============================================================================
# Run this locally after all NEBULA jobs complete and metadata is compiled.
# Reads processed results from S3, generates:
#   1. Volcano plots (p-value < 0.05 and FDR < 0.05)
#   2. GSEA hallmark pathway plots (styled after plot_gsea_results)
# Saves figures locally to results_dir/Figures/
# =============================================================================

library(aws.s3)
library(jsonlite)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(fgsea)
library(msigdbr)
library(tidyverse)
library(stringr)
library(purrr)
library(Hmisc)

# --- S3 setup ---
user <- Sys.info()[["user"]]
if (user == "choiyej") {
  keys <- fromJSON("/Users/choiyej/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Bjornstad Pyle Lab/keys.json")
  results_dir <- "/Users/choiyej/Library/CloudStorage/OneDrive-UW/Bjornstad/T1D Adiposity"
} else if (user == "yejichoi") {
  keys <- fromJSON("/mmfs1/home/yejichoi/keys.json")
  results_dir <- "/mmfs1/gscratch/togo/yejichoi/T1D_Adiposity"
} else if (user == "pylell") {
  keys <- fromJSON("/mmfs1/home/pylell/keys.json")
  results_dir <- "/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/T1D Adiposity"
} else {
  stop("Unknown user: please specify keys path and results dir.")
}

Sys.setenv(
  "AWS_ACCESS_KEY_ID" = keys$MY_ACCESS_KEY,
  "AWS_SECRET_ACCESS_KEY" = keys$MY_SECRET_KEY,
  "AWS_DEFAULT_REGION" = "",
  "AWS_REGION" = "",
  "AWS_S3_ENDPOINT" = "s3.kopah.uw.edu"
)

bucket <- "t1d.adiposity"
s3_base <- "results/nebula"

# --- Create figure directories ---
fig_dirs <- c(
  file.path(results_dir, "Results/Figures/Volcano Plots/pval"),
  file.path(results_dir, "Results/Figures/Volcano Plots/fdr"),
  file.path(results_dir, "Results/Figures/Pathways"),
  file.path(results_dir, "Results/Figures/Pathways Lollipop")
)
for (d in fig_dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# --- Load Hallmark gene sets ---
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H") %>%
  split(.$gs_name) %>%
  lapply(function(x) x$gene_symbol)

# --- Read compiled metadata ---
meta_compiled <- tryCatch({
  tmp <- tempfile(fileext = ".csv")
  save_object(object = paste0(s3_base, "/csv/run_metadata_compiled.csv"),
              bucket = bucket, region = "", file = tmp)
  df <- read.csv(tmp, stringsAsFactors = FALSE)
  unlink(tmp)
  df
}, error = function(e) {
  message("Could not read compiled metadata: ", e$message)
  data.frame()
})

if (nrow(meta_compiled) == 0) {
  cat("No completed results found in metadata. Run compile_metadata.R first.\n")
  quit(save = "no", status = 0)
}

cat(sprintf("Processing %d analysis-celltype combinations...\n", nrow(meta_compiled)))

# Track success/failure counts
n_success <- 0
n_skipped <- 0
n_error <- 0

# =============================================================================
# HELPER FUNCTIONS (matching RH_RH2_IMPROVE_functions.R style)
# =============================================================================

clean_pathway_names <- function(pathways) {
  cleaned <- gsub("^REACTOME_", "", pathways)
  cleaned <- gsub("^GOBP_", "", cleaned)
  cleaned <- gsub("^KEGG_", "", cleaned)
  cleaned <- gsub("^HALLMARK_", "", cleaned)
  cleaned <- gsub("_", " ", cleaned)
  cleaned <- tools::toTitleCase(tolower(cleaned))

  uppercase_words <- c(
    "\\bI\\b", "\\bIi\\b", "\\bIii\\b", "\\bIv\\b", "\\bV\\b", "\\bVi\\b",
    "\\bVii\\b", "\\bViii\\b", "\\bIx\\b", "\\bX\\b",
    "\\bTca\\b", "\\bAtp\\b", "\\bAdp\\b", "\\bAmp\\b", "\\bGtp\\b", "\\bGdp\\b",
    "\\bNad\\b", "\\bNadh\\b", "\\bFad\\b", "\\bFadh2\\b", "\\bCoa\\b",
    "\\bDna\\b", "\\bRna\\b", "\\bMrna\\b", "\\bTrna\\b", "\\bRrna\\b",
    "\\bEr\\b", "\\bUpr\\b", "\\bNf\\b", "\\bHif\\b", "\\bMhc\\b",
    "\\bTgf\\b", "\\bEgf\\b", "\\bVegf\\b", "\\bPdgf\\b", "\\bFgf\\b",
    "\\bRos\\b", "\\bRns\\b", "\\bNo\\b", "\\bNos\\b", "\\bInos\\b", "\\bEnos\\b", "\\bNnos\\b",
    "\\Mapk\\b"
  )

  for (pattern in uppercase_words) {
    word <- gsub("\\\\b", "", pattern)
    replacement <- toupper(word)
    cleaned <- gsub(pattern, replacement, cleaned, ignore.case = TRUE)
  }

  cleaned <- gsub("\\bOf\\b", "of", cleaned)
  cleaned <- gsub("\\bBy\\b", "by", cleaned)
  cleaned <- gsub("\\bThe\\b", "the", cleaned)
  cleaned <- gsub("^the", "The", cleaned)
  cleaned <- gsub("\\bO Linked\\b", "O-Linked", cleaned, ignore.case = TRUE)

  return(cleaned)
}

shorten_pathway_names <- function(pathway_names, max_length = 40, aggressive = FALSE) {
  word_replacements <- c(
    "Hallmark " = "",
    "Role of " = "",
    "Function of " = "",
    "Eukaryotic " = "",
    "in the" = "in",
    "by the" = "by",
    "of the" = "of",
    "Phosphorylation" = "Phosph",
    "Metabolism" = "Metab",
    "Extracellular" = "EC",
    "Intracellular" = "IC",
    "Dysfunction" = "Adaptation",
    "Alternative" = "Alt",
    "Dependent" = "Dep",
    "Independent" = "Indep",
    "Associated" = "Assoc",
    "Inflammatory" = "Inflam",
    "Inflammation" = "Inflam",
    "Lymphocyte" = "Lymph",
    "Lymphocytes" = "Lymph",
    "Fibroblast" = "Fib",
    "Fibroblasts" = "Fib",
    "Endothelial" = "Endo",
    "Epithelial" = "Epi",
    "Transition" = "Trans",
    "Rheumatoid Arthritis" = "RA",
    "Multiple Sclerosis" = "MS",
    "Alzheimer's Disease" = "Alzheimer's",
    "Parkinson's Disease" = "Parkinson's",
    "Respiratory" = "Resp",
    "Gastrointestinal" = "GI",
    "Central Nervous System" = "CNS",
    "Oxidative" = "Ox",
    "between" = "b/w",
    "through" = "via",
    "and" = "&"
  )

  aggressive_replacements <- c(
    " the " = " ",
    " in " = " ",
    " of " = " ",
    " to " = "\u2192",
    " from " = "\u2190"
  )

  molecule_replacements <- c(
    "Natural Killer" = "NK",
    "Tumor Necrosis Factor" = "TNF",
    "Transforming Growth Factor" = "TGF",
    "Platelet Derived Growth Factor" = "PDGF",
    "Epidermal Growth Factor" = "EGF",
    "Vascular Endothelial Growth Factor" = "VEGF",
    "Insulin-like Growth Factor" = "IGF",
    "Interferon" = "IFN",
    "Interleukin" = "IL",
    "G-Protein Coupled" = "GPCR",
    "G Protein Coupled" = "GPCR",
    "Activator Protein" = "AP",
    "Electron Transport" = "e- Transport",
    "Nitric Oxide" = "NO",
    "Reactive Oxygen Species" = "ROS",
    "Amino Acid" = "AA",
    "Fatty Acid" = "FA"
  )

  apply_replacements <- function(text, replacements) {
    for (pattern in names(replacements)) {
      text <- gsub(pattern, replacements[pattern], text, ignore.case = FALSE)
    }
    return(text)
  }

  step6_truncated <- rep(FALSE, length(pathway_names))

  shortened <- sapply(seq_along(pathway_names), function(i) {
    pathway <- pathway_names[i]
    if (is.na(pathway)) return(NA)

    result <- apply_replacements(pathway, molecule_replacements)
    result <- apply_replacements(result, word_replacements)

    if (nchar(result) > max_length && aggressive) {
      result <- apply_replacements(result, aggressive_replacements)
    }

    if (nchar(result) > max_length) {
      result <- gsub(" \\([^)]+\\)", "", result)
    }

    if (nchar(result) > max_length) {
      step6_truncated[i] <<- TRUE
      result <- paste0(substr(result, 1, max_length - 3), "...")
    }

    return(result)
  })

  return(list(
    shortened = shortened,
    step6_truncated = step6_truncated
  ))
}

# =============================================================================
# VOLCANO PLOT FUNCTION
# =============================================================================
make_volcano <- function(data, p_col, fc, title = NULL,
                         test_name,
                         sig_type = "pval", fdr_col = NULL,
                         p_thresh = 0.05,
                         cell_type = "",
                         formula_text = NULL,
                         cohort_text = NULL,
                         positive_text = "Positive",
                         negative_text = "Negative",
                         text_size = 20,
                         geom_text_size = 4,
                         caption_size = 8.5,
                         legend_text_size = 10,
                         volcano_force = 6,
                         volcano_box_padding = 0,
                         off_chart_threshold = 0.95,
                         off_chart_y_position = 0.85,
                         off_chart_arrow_length = 0.02) {
  
  if (!p_col %in% names(data) || !fc %in% names(data)) return(NULL)
  
  data <- data %>% dplyr::filter(!is.na(.data[[p_col]]) & !is.na(.data[[fc]]))
  if (nrow(data) == 0) return(NULL)
  
  # Apply significance filter based on sig_type
  if (sig_type == "fdr") {
    if (is.null(fdr_col) || !fdr_col %in% names(data)) return(NULL)
    data <- data %>% dplyr::filter(!is.na(.data[[fdr_col]]))
    data$is_sig <- data[[fdr_col]] < p_thresh
  } else {
    data$is_sig <- data[[p_col]] < p_thresh
  }
  
  set.seed(1)
  epsilon <- 1e-300
  
  # For FDR plots, use adjusted p-values on the y-axis; otherwise use raw p-values
  y_col <- if (sig_type == "fdr") fdr_col else p_col
  data <- data %>%
    dplyr::mutate(neg_log_p = -log10(.data[[y_col]] + epsilon))
  
  y_max <- max(data$neg_log_p, na.rm = TRUE) * 1.1
  y_cutoff <- y_max * off_chart_threshold
  
  top_pos <- data %>%
    dplyr::filter(.data[[fc]] > 0 & is_sig) %>%
    dplyr::arrange(.data[[y_col]])
  n_pos <- nrow(top_pos)
  top_pos_n <- top_pos %>% dplyr::slice_head(n = 20)
  
  top_neg <- data %>%
    dplyr::filter(.data[[fc]] < 0 & is_sig) %>%
    dplyr::arrange(.data[[y_col]])
  n_neg <- nrow(top_neg)
  top_neg_n <- top_neg %>% dplyr::slice_head(n = 20)
  
  # Identify off-chart genes
  off_chart_genes <- data %>%
    dplyr::filter(Gene %in% c(top_pos_n$Gene, top_neg_n$Gene) & neg_log_p > y_cutoff) %>%
    dplyr::mutate(
      is_positive = .data[[fc]] > 0,
      x_position = ifelse(is_positive,
                          .data[[fc]] + seq(from = 0.1, by = 0.2, length.out = n()),
                          .data[[fc]] - seq(from = 0.1, by = 0.2, length.out = n())),
      y_position = y_max * off_chart_y_position
    )
  
  on_chart_genes <- c(top_pos_n$Gene, top_neg_n$Gene)[!c(top_pos_n$Gene, top_neg_n$Gene) %in% off_chart_genes$Gene]
  
  data <- data %>%
    dplyr::mutate(
      top_color = case_when(
        Gene %in% top_pos$Gene ~ "#f28482",
        Gene %in% top_neg$Gene ~ "#457b9d",
        TRUE ~ "#ced4da"
      ),
      top_size = if_else(Gene %in% c(top_pos$Gene, top_neg$Gene), 1.3, 1),
      top_lab  = if_else(Gene %in% on_chart_genes, Gene, ""),
      display_neg_log_p = pmin(neg_log_p, y_cutoff)
    ) %>%
    dplyr::filter(abs(.data[[fc]]) < 10)
  
  max_fc <- max(data[[fc]], na.rm = TRUE)
  min_fc <- min(data[[fc]], na.rm = TRUE)
  
  # Build caption
  has_formula <- !is.null(formula_text) && !is.na(formula_text)
  has_cohort  <- !is.null(cohort_text) && !is.na(cohort_text)
  
  caption_text <- paste0(
    if (has_formula && has_cohort) {
      paste0("Formula: ", formula_text, " + (1|subject) | Cohort: ", cohort_text, "\n\n")
    } else if (has_formula) {
      paste0("Formula: ", formula_text, " + (1|subject)\n\n")
    } else if (has_cohort) {
      paste0("Cohort: ", cohort_text, "\n\n")
    } else { "" },
    "Cell type: ", cell_type,
    if (cell_type != "") " | " else "",
    "Positive n = ", n_pos, " | Negative n = ", n_neg,
    if (nrow(off_chart_genes) > 0) paste0("\n", nrow(off_chart_genes), " gene(s) with p-values near zero (arrows indicate off-scale values)") else ""
  )
  
  # Y-axis label
  y_label <- ifelse(sig_type == "fdr", "-log10(adj. p-value)", "-log10(p-value)")
  
  p <- ggplot(data, aes(x = .data[[fc]], y = display_neg_log_p)) +
    geom_hline(yintercept = -log10(p_thresh), linetype = "dashed", color = "darkgrey") +
    geom_point(alpha = 0.5, aes(color = top_color), size = 3) +
    geom_label_repel(seed = 1,
                     data = dplyr::filter(data, top_lab != ""),
                     aes(label = top_lab, color = top_color),
                     fontface = "bold",
                     size = geom_text_size, max.overlaps = Inf,
                     force = volcano_force, segment.alpha = 0.3, segment.size = 0.3,
                     box.padding = volcano_box_padding,
                     fill = fill_alpha("white", 0.7),
                     label.size = 0
    )
  
  # Add arrows for off-chart genes
  if (nrow(off_chart_genes) > 0) {
    p <- p +
      geom_segment(
        data = off_chart_genes,
        aes(x = .data[[fc]], y = y_cutoff * 0.98,
            xend = .data[[fc]], yend = y_cutoff - (y_max * off_chart_arrow_length)),
        arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
        color = "black", size = 0.6
      ) +
      geom_text(
        data = off_chart_genes,
        aes(x = x_position, y = y_position, label = Gene),
        size = geom_text_size * 0.9,
        hjust = ifelse(off_chart_genes$is_positive, 0, 1),
        color = "black", fontface = "italic"
      ) +
      geom_text(
        data = off_chart_genes,
        aes(x = x_position, y = y_position - (y_max * 0.03),
            label = paste0("p=", format(.data[[y_col]], scientific = TRUE, digits = 2))),
        size = geom_text_size * 0.7,
        hjust = ifelse(off_chart_genes$is_positive, 0, 1),
        color = "darkgrey"
      )
  }
  
  p <- p +
    labs(title = title,
         x = paste0("log2 FC ", test_name),
         y = y_label,
         caption = caption_text) +
    scale_size_continuous(range = c(1, 1.3)) +
    scale_color_manual(values = c("#457b9d" = "#457b9d", "#ced4da" = "#ced4da", "#f28482" = "#f28482")) +
    guides(color = "none", size = "none") +
    annotate("segment",
             x = max_fc / 8, xend = (max_fc * 7) / 8,
             y = -y_max * 0.13,
             col = "darkgrey", arrow = arrow(length = unit(0.2, "cm"))) +
    annotate("text",
             x = mean(c(max_fc / 8, (max_fc * 7) / 8)),
             y = -y_max * 0.18,
             label = positive_text,
             size = geom_text_size, color = "#343a40") +
    annotate("segment",
             x = min_fc / 8, xend = (min_fc * 7) / 8,
             y = -y_max * 0.13,
             col = "darkgrey", arrow = arrow(length = unit(0.2, "cm"))) +
    annotate("text",
             x = mean(c(min_fc / 8, (min_fc * 7) / 8)),
             y = -y_max * 0.18,
             label = negative_text,
             size = geom_text_size, color = "#343a40") +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, y_max), clip = "off") +
    theme(legend.title = element_blank(),
          panel.grid = element_blank(),
          text = element_text(size = text_size),
          title = element_text(size = legend_text_size),
          plot.margin = margin(t = 10, r = 20, b = 25, l = 20),
          axis.title.x = element_text(margin = margin(t = 32)),
          plot.caption = element_text(size = caption_size, hjust = 0.5, margin = margin(t = 15)),
          panel.background = element_rect(fill = "transparent", color = NA),
          plot.background  = element_rect(fill = "transparent", color = NA),
          legend.background = element_rect(fill = "transparent", color = NA),
          legend.box.background = element_rect(fill = "transparent", color = NA))
  
  return(p)
}

# =============================================================================
# GSEA FUNCTION (styled after plot_gsea_results from RH_RH2_IMPROVE_functions.R)
# =============================================================================
run_gsea_and_plot <- function(df, logfc_col, pval_col, cell_name,
                              test_name, cohort_text,
                               pathways = hallmark_sets,
                               reference = "hallmark",
                               top_n = 20,
                               max_pathway_length = 45,
                               caption_width = 70,
                               p_threshold = 0.05,
                               min_x_limit = 3,
                               low_color = "#89c2d9",
                               mid_color = "white",
                               high_color = "#ee7674",
                               show_truncated_in_caption = TRUE) {

  if (!logfc_col %in% names(df) || !pval_col %in% names(df)) return(NULL)

  # Create ranked gene list (signed -log10 p-value)
  df <- df %>%
    dplyr::filter(!is.na(.data[[logfc_col]]) & !is.na(.data[[pval_col]]) & .data[[pval_col]] > 0) %>%
    dplyr::mutate(rank_stat = sign(.data[[logfc_col]]) * -log10(.data[[pval_col]]))

  gene_ranks <- setNames(df$rank_stat, df$Gene)
  gene_ranks <- sort(gene_ranks, decreasing = TRUE)
  gene_ranks <- gene_ranks[!duplicated(names(gene_ranks))]

  if (length(gene_ranks) < 10) return(NULL)

  # Run fgsea
  fgsea_res <- tryCatch(
    fgsea(pathways = pathways, stats = gene_ranks, minSize = 15, maxSize = 500),
    error = function(e) { message(e$message); NULL }
  )

  if (is.null(fgsea_res) || nrow(fgsea_res) == 0) return(NULL)

  # Top pathways by p-value
  top_pathways <- fgsea_res %>%
    dplyr::arrange(pval) %>%
    head(top_n) %>%
    dplyr::mutate(
      leadingEdge_size = lengths(leadingEdge),
      leadingEdge_fraction = leadingEdge_size / size
    )

  # Shorten and clean pathway names
  shorten_result <- shorten_pathway_names(top_pathways$pathway, max_length = max_pathway_length)

  # Append leading edge fraction percentage
  shorten_result$shortened <- paste0(shorten_result$shortened, " (", round(top_pathways$leadingEdge_fraction * 100), "%)")

  top_pathways <- top_pathways %>%
    dplyr::mutate(
      was_step6_truncated = shorten_result$step6_truncated,
      shortened_pathway = shorten_result$shortened,
      clean_pathway = clean_pathway_names(shortened_pathway),
      neg_log_p = -log10(pval)
    )

  # Build caption
  base_caption <- paste0("Cell type: ", gsub("_", " ", cell_name),
                         " | Model: ", test_name,
                         " | Cohort: ", cohort_text,
                          " | Reference: ", toupper(reference))

  if (show_truncated_in_caption) {
    truncated_pathways <- top_pathways %>%
      dplyr::filter(was_step6_truncated) %>%
      dplyr::mutate(
        full_clean = clean_pathway_names(pathway),
        full_clean_wrapped = str_wrap(full_clean, width = caption_width, indent = 0, exdent = 4)
      ) %>%
      dplyr::select(clean_pathway, full_clean_wrapped)

    if (nrow(truncated_pathways) > 0) {
      truncated_text <- paste(
        apply(truncated_pathways, 1, function(x) {
          paste0(x["full_clean_wrapped"])
        }),
        collapse = "\n"
      )
      full_caption <- paste0(base_caption, "\n\nTruncated pathways:\n", truncated_text)
    } else {
      full_caption <- base_caption
    }
  } else {
    full_caption <- base_caption
  }

  # Set factor levels for plotting order
  top_pathways$clean_pathway <- factor(top_pathways$clean_pathway,
                                        levels = rev(top_pathways$clean_pathway))

  # Determine x-axis limits
  if (!is.null(min_x_limit)) {
    actual_max <- max(top_pathways$neg_log_p, na.rm = TRUE)
    x_upper_limit <- max(actual_max, min_x_limit)
  } else {
    x_upper_limit <- NULL
  }
  
  # Create plot (matching plot_gsea_results style)
  p <- top_pathways %>%
    ggplot(aes(y = clean_pathway, x = neg_log_p, fill = NES)) +
    geom_col(width = 0.9) +
    geom_vline(xintercept = -log10(p_threshold), linetype = "dashed", color = "#aaaaaa") +
    geom_text(aes(label = clean_pathway),
              x = -log10(p_threshold) + 0.1, hjust = 0,
              fontface = "bold",
              color = "#2b2b2b") +
    scale_fill_gradient2(low = low_color, mid = mid_color, high = high_color,
                          midpoint = 0,
                          guide = guide_colorbar(barheight = 0.4, barwidth = 8)) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "top",
          axis.text.y = element_blank(),
          legend.title = element_text(vjust = 0.8),
          plot.caption = element_text(hjust = 0, size = 8, lineheight = 1.2),
          panel.background = element_rect(fill = "transparent", color = NA),
          plot.background  = element_rect(fill = "transparent", color = NA),
          legend.background = element_rect(fill = "transparent", color = NA),
          legend.box.background = element_rect(fill = "transparent", color = NA)) +
    labs(y = NULL,
         x = "-log(p-value)",
         fill = "NES",
         caption = full_caption,
         title = cell_name)

  # Apply x-axis limits
  if (!is.null(min_x_limit)) {
    p <- p + scale_x_continuous(limits = c(0, x_upper_limit), expand = c(0, 0))
  } else {
    p <- p + scale_x_continuous(expand = expansion(mult = c(0, 0)))
  }

  return(p)
}

# =============================================================================
# GSEA LOLLIPOP FUNCTION
# =============================================================================
make_gsea_lollipop <- function(df, logfc_col, pval_col, cell_name,
                               pathways = hallmark_sets,
                               test_name, cohort_text,
                               reference = "hallmark",
                               top_n = 20,
                               max_pathway_length = 45,
                               p_threshold = 0.05,
                               low_color = "#219ebc",
                               mid_color = "white",
                               high_color = "#d00000",
                               point_size_range = c(3, 10),
                               text_size_val = 5,
                               base_text_size = 15) {
  
  if (!logfc_col %in% names(df) || !pval_col %in% names(df)) return(NULL)
  
  # Create ranked gene list (signed -log10 p-value)
  df <- df %>%
    dplyr::filter(!is.na(.data[[logfc_col]]) & !is.na(.data[[pval_col]]) & .data[[pval_col]] > 0) %>%
    dplyr::mutate(rank_stat = sign(.data[[logfc_col]]) * -log10(.data[[pval_col]]))
  
  gene_ranks <- setNames(df$rank_stat, df$Gene)
  gene_ranks <- sort(gene_ranks, decreasing = TRUE)
  gene_ranks <- gene_ranks[!duplicated(names(gene_ranks))]
  
  if (length(gene_ranks) < 10) return(NULL)
  
  # Run fgsea
  fgsea_res <- tryCatch(
    fgsea(pathways = pathways, stats = gene_ranks, minSize = 15, maxSize = 500),
    error = function(e) { message(e$message); NULL }
  )
  
  if (is.null(fgsea_res) || nrow(fgsea_res) == 0) return(NULL)
  
  # Filter significant, take top_n by p-value, then arrange by NES for display
  plot_df <- fgsea_res %>%
    dplyr::filter(!is.na(pval)) %>%
    dplyr::filter(pval < p_threshold & NES != 0) %>%
    dplyr::arrange(pval) %>%
    head(top_n) %>%
    dplyr::arrange(NES)
  
  if (nrow(plot_df) == 0) return(NULL)
  
  # Clean pathway names and assign row number for y position
  plot_df <- plot_df %>%
    dplyr::mutate(
      clean_name = clean_pathway_names(pathway),
      neg_log_p = -log10(pval),
      n = row_number(),
      text_col = ifelse(NES > 0, high_color, low_color)
    )
  
  # Shorten long names
  shorten_result <- shorten_pathway_names(plot_df$clean_name, max_length = max_pathway_length)
  plot_df$display_name <- shorten_result$shortened
  
  # Caption
  base_caption <- paste0("Cell type: ", gsub("_", " ", cell_name),
                         " | Model: ", test_name,
                         " | Cohort: ", cohort_text,
                         " | Reference: ", toupper(reference))
  
  # --- Dynamic x-axis expansion based on label lengths ---
  nes_range <- diff(range(plot_df$NES))
  if (nes_range == 0) nes_range <- 1  # guard against single-direction pathways
  
  pos_labels <- plot_df$display_name[plot_df$NES > 0]
  neg_labels <- plot_df$display_name[plot_df$NES < 0]
  
  max_pos_chars <- if (length(pos_labels) > 0) max(nchar(pos_labels)) else 0
  max_neg_chars <- if (length(neg_labels) > 0) max(nchar(neg_labels)) else 0
  
  char_to_nes <- text_size_val * 0.016
  pos_label_width <- max_pos_chars * char_to_nes
  neg_label_width <- max_neg_chars * char_to_nes
  
  nes_min <- min(plot_df$NES)
  nes_max <- max(plot_df$NES)
  total_span <- max(abs(nes_min), abs(nes_max), 1)
  
  left_expand  <- max(0.02, pos_label_width / total_span * 0.5)
  right_expand <- max(0.02, neg_label_width / total_span * 0.5)
  
  left_expand  <- min(left_expand, 0.6)
  right_expand <- min(right_expand, 0.6)
  
  # --- Build unique breaks/labels for the color legend (min, 0, max) ---
  # The key fix: explicitly set `limits` to always include 0 AND all NES values.
  # Without this, scale_color_gradient2 auto-sets limits to the data range,
  
  # which drops out-of-range breaks (e.g., 0 when all NES > 0) while keeping
  # the corresponding labels, causing the breaks/labels length mismatch error.
  color_limits <- c(min(nes_min, 0), max(nes_max, 0))
  
  brk <- c(nes_min, 0, nes_max)
  lbl <- c(as.character(round(nes_min, 1)), "0", as.character(round(nes_max, 1)))
  
  # Remove duplicates (e.g., when nes_min == nes_max, or nes_min/nes_max == 0)
  dup <- !duplicated(round(brk, 10))
  brk <- brk[dup]
  lbl <- lbl[dup]
  
  # Build plot (matching panther_baseline_manuscript.qmd style)
  p <- plot_df %>%
    ggplot(aes(x = NES, y = n, size = neg_log_p, color = NES)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#495057") +
    geom_segment(aes(xend = 0), linewidth = 1) +
    geom_point() +
    scale_size(range = point_size_range) +
    scale_color_gradient2(low = low_color, mid = mid_color, high = high_color, midpoint = 0,
                          limits = color_limits,
                          breaks = brk, labels = lbl,
                          guide = guide_colorbar(barwidth = 8, barheight = 0.5)) +
    geom_text(data = subset(plot_df, NES > 0),
              aes(x = -0.02, y = n, label = display_name),
              color = subset(plot_df, NES > 0)$text_col,
              hjust = 1, inherit.aes = FALSE, size = text_size_val) +
    geom_text(data = subset(plot_df, NES < 0),
              aes(x = 0.02, y = n, label = display_name),
              color = subset(plot_df, NES < 0)$text_col,
              hjust = 0, inherit.aes = FALSE, size = text_size_val) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          text = element_text(size = base_text_size),
          axis.text.y = element_blank(),
          legend.position = "top",
          legend.title.position = "top",
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.caption = element_text(hjust = 0, size = 10, lineheight = 1.2),
          panel.background = element_rect(fill = "transparent", color = NA),
          plot.background  = element_rect(fill = "transparent", color = NA),
          legend.background = element_rect(fill = "transparent", color = NA),
          legend.box.background = element_rect(fill = "transparent", color = NA)) +
    labs(y = NULL, x = "NES",
         color = "NES", size = "-log10(p)",
         caption = base_caption,
         title = cell_name) +
    scale_x_continuous(expand = expansion(mult = c(left_expand, right_expand)))
  
  return(p)
}

# =============================================================================
# DEG BUTTERFLY PLOT FUNCTION
# =============================================================================
# sig_type: "pval" counts genes with raw p < 0.05; "fdr" counts genes with FDR < 0.05
plot_deg_butterfly <- function(sig_df_summary, analysis_type_filter, ct_names, output_file,
                               sig_type = "pval") {
  count_cols <- if (sig_type == "fdr") c("Positive_fdr", "Negative_fdr") else c("Positive", "Negative")
  
  sig_df_clean <- sig_df_summary %>%
    filter(analysis_type == analysis_type_filter) %>%
    dplyr::rename(Pos = !!count_cols[1], Neg = !!count_cols[2]) %>%
    dplyr::mutate(Total = Pos + Neg) %>%
    arrange(desc(Total)) %>%
    pivot_longer(cols = c(Pos, Neg), names_to = "Direction", values_to = "Count") %>%
    dplyr::mutate(
      Direction = ifelse(Direction == "Neg", "Negative", "Positive"),
      Count = ifelse(Direction == "Negative", -Count, Count),
      celltype = factor(celltype, levels = rev(unique(celltype)))
    )
  
  sig_label <- ifelse(sig_type == "fdr", "FDR < 0.05", "p < 0.05")
  
  p <- sig_df_clean %>%
    filter(celltype %in% ct_names) %>%
    ggplot(aes(x = Count, y = celltype, fill = Direction)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = 0, color = "white", linewidth = 0.6) +
    geom_text(aes(label = celltype, x = 0), size = 4) +
    scale_x_continuous(labels = abs, expand = expansion(mult = c(0.01, 0.01))) +
    scale_fill_manual(values = c("Negative" = "#c0d6df", "Positive" = "#f4978e")) +
    labs(x = paste0("Number of DE genes (", sig_label, ")"), y = NULL, fill = NULL,
         title = analysis_type_filter) +
    theme(legend.position = "top", panel.background = element_blank(), text = element_text(size = 15),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave(output_file, plot = p, width = 7, height = 7, bg = "transparent")
  p
}

# =============================================================================
# FILTER METADATA: exclude low cell counts and "Other" celltype
# =============================================================================
meta_compiled <- meta_compiled %>%
  dplyr::filter(total_cells >= 1000, celltype != "Other")

cat(sprintf("After filtering (total_cells >= 1000, celltype != 'Other'): %d combinations\n",
            nrow(meta_compiled)))

if (nrow(meta_compiled) == 0) {
  cat("No combinations remain after filtering. Exiting.\n")
  quit(save = "no", status = 0)
}

# =============================================================================
# PHASE 1: BUTTERFLY SCAN — run first to survey DEG counts across cell types
# =============================================================================
cat("\n=== Phase 1: Scanning DEG counts for butterfly plots ===\n")
deg_summary_list <- list()

for (i in seq_len(nrow(meta_compiled))) {
  row <- meta_compiled[i, ]
  analysis_type <- row$analysis_type
  celltype <- row$celltype
  s3_key <- row$s3_results_key
  pval_col <- row$pval_column
  logfc_col <- row$logfc_column
  
  cat(sprintf("  [%d/%d] Scanning %s - %s\n", i, nrow(meta_compiled), analysis_type, celltype))
  
  df <- tryCatch(
    s3readRDS(object = s3_key, bucket = bucket, region = ""),
    error = function(e) {
      message(paste("    Error reading results:", e$message))
      NULL
    }
  )
  if (is.null(df)) next
  
  # Determine FDR column
  fdr_candidates <- grep("p_.*_fdr$|padj|q_value|^fdr$", names(df), value = TRUE)
  fdr_col <- if (length(fdr_candidates) > 0) fdr_candidates[1] else NULL
  
  # Count DEGs at p < 0.05
  sig_pval <- df %>%
    dplyr::filter(!is.na(.data[[pval_col]]) & !is.na(.data[[logfc_col]]) & .data[[pval_col]] < 0.05)
  
  pos_pval <- sum(sig_pval[[logfc_col]] > 0, na.rm = TRUE)
  neg_pval <- sum(sig_pval[[logfc_col]] < 0, na.rm = TRUE)
  
  # Count DEGs at FDR < 0.05
  pos_fdr <- 0
  neg_fdr <- 0
  if (!is.null(fdr_col) && fdr_col %in% names(df)) {
    sig_fdr <- df %>%
      dplyr::filter(!is.na(.data[[fdr_col]]) & !is.na(.data[[logfc_col]]) & .data[[fdr_col]] < 0.05)
    pos_fdr <- sum(sig_fdr[[logfc_col]] > 0, na.rm = TRUE)
    neg_fdr <- sum(sig_fdr[[logfc_col]] < 0, na.rm = TRUE)
  }
  
  # Determine resolution group from celltype_variable column in metadata
  ct_var <- row$celltype_variable
  resolution <- dplyr::case_when(
    ct_var == "KPMP_celltype"          ~ "high_res",
    ct_var == "KPMP_celltype_general"  ~ "low_res",
    ct_var == "KPMP_celltype_general2" ~ "low_res2",
    TRUE                               ~ "other"
  )
  
  deg_summary_list[[length(deg_summary_list) + 1]] <- data.frame(
    analysis_type = analysis_type,
    celltype = celltype,
    celltype_variable = ct_var,
    resolution = resolution,
    Positive = pos_pval, Negative = neg_pval,
    Positive_fdr = pos_fdr, Negative_fdr = neg_fdr,
    stringsAsFactors = FALSE
  )
}

# Generate butterfly plots (pval and FDR) for each analysis type × resolution
if (length(deg_summary_list) > 0) {
  deg_summary_df <- bind_rows(deg_summary_list)
  
  # Create subdirectories for each resolution
  for (res in unique(deg_summary_df$resolution)) {
    dir.create(file.path(results_dir, "Results/Figures/Butterfly", res),
               recursive = TRUE, showWarnings = FALSE)
  }
  
  # Save summary CSV for inspection
  write.csv(deg_summary_df,
            file.path(results_dir, "Results/Figures/Butterfly/deg_summary.csv"),
            row.names = FALSE)
  cat(sprintf("  DEG summary saved: deg_summary.csv (%d rows)\n", nrow(deg_summary_df)))
  
  # Get unique analysis_type × resolution combinations
  combos <- deg_summary_df %>%
    dplyr::distinct(analysis_type, resolution)
  
  cat(sprintf("\n  Generating butterfly plots for %d analysis × resolution combinations...\n",
              nrow(combos)))
  
  for (j in seq_len(nrow(combos))) {
    at  <- combos$analysis_type[j]
    res <- combos$resolution[j]
    at_norm <- tolower(gsub("[/ ]+", "_", at))
    
    # Subset to this analysis_type + resolution
    sub_df <- deg_summary_df %>%
      dplyr::filter(analysis_type == at, resolution == res)
    ct_names <- unique(sub_df$celltype)
    
    # p-value butterfly
    tryCatch({
      plot_deg_butterfly(sub_df,
                         analysis_type_filter = at,
                         ct_names = ct_names,
                         output_file = file.path(results_dir, "Results/Figures/Butterfly", res,
                                                 paste0(at_norm, "_butterfly_pval.png")),
                         sig_type = "pval")
      cat(sprintf("    %s [%s] (pval)\n", at, res))
    }, error = function(e) {
      message(sprintf("    Error: %s [%s] (pval): %s", at, res, e$message))
    })
    
    # FDR butterfly
    tryCatch({
      plot_deg_butterfly(sub_df,
                         analysis_type_filter = at,
                         ct_names = ct_names,
                         output_file = file.path(results_dir, "Results/Figures/Butterfly", res,
                                                 paste0(at_norm, "_butterfly_fdr.png")),
                         sig_type = "fdr")
      cat(sprintf("    %s [%s] (fdr)\n", at, res))
    }, error = function(e) {
      message(sprintf("    Error: %s [%s] (fdr): %s", at, res, e$message))
    })
  }
}

cat("\n=== Phase 1 complete. Review butterfly plots to decide which analyses to focus on. ===\n")
cat("=== Phase 2: Generating per-celltype figures (volcano, GSEA, lollipop) ===\n\n")

# =============================================================================
# PHASE 2: MAIN LOOP — Generate per-celltype figures
# =============================================================================
for (i in seq_len(nrow(meta_compiled))) {
  row <- meta_compiled[i, ]
  analysis_type <- row$analysis_type
  celltype <- row$celltype
  s3_key <- row$s3_results_key
  pval_col <- row$pval_column
  logfc_col <- row$logfc_column
  
  cell_norm <- tolower(gsub("/", "_", celltype))
  cat(sprintf("[%d/%d] %s - %s\n", i, nrow(meta_compiled), analysis_type, celltype))
  
  # Read processed results from S3
  df <- tryCatch(
    s3readRDS(object = s3_key, bucket = bucket, region = ""),
    error = function(e) {
      message(paste("  Error reading results:", e$message))
      NULL
    }
  )
  if (is.null(df)) {
    n_error <- n_error + 1
    next
  }
  
  # Determine FDR column
  fdr_candidates <- grep("p_.*_fdr$|padj|q_value|^fdr$", names(df), value = TRUE)
  fdr_col <- if (length(fdr_candidates) > 0) fdr_candidates[1] else NULL
  
  # Determine formula text for caption (extract from analysis_type)
  formula_text <- sub("^~", "~ ", row$formula_used)  # If available in metadata; NULL otherwise
  test_group <- sub(".*?([A-Z].*)", "\\1", gsub("_", " ", row$pval_column))
  test_name <- paste0(test_group, " vs. ", gsub("_", " ", row$reference_level))
  cohort <- ifelse(
    grepl("group\\s*==\\s*'Type 1 Diabetes'", row$subset_condition) &
      grepl("group\\s*==\\s*'Lean Control", row$subset_condition),
    "T1D + HC",
    ifelse(
      grepl("group\\s*==\\s*'Type 1 Diabetes'", row$subset_condition),
      "T1D",
      ""
    )
  )
  
  # --- Volcano plot (p-value) ---
  vp_pval <- make_volcano(df, p_col = pval_col, fc = logfc_col,
                          sig_type = "pval",
                          cell_type = celltype,
                          test_name = test_name,
                          formula_text = formula_text,
                          cohort_text = cohort,
                          positive_text = paste0("Positive with ", test_group),
                          negative_text = paste0("Negative with ", test_group))
  if (!is.null(vp_pval)) {
    ggsave(file.path(results_dir, "Results/Figures/Volcano Plots/pval",
                     paste0(analysis_type, "_", cell_norm, ".png")),
           plot = vp_pval, width = 7, height = 5, dpi = 300, bg = "transparent")
  }
  
  # --- Volcano plot (FDR) — skip if no FDR-significant DEGs ---
  fdr_row <- deg_summary_df %>%
    dplyr::filter(analysis_type == !!analysis_type, celltype == !!celltype) %>%
    dplyr::slice(1)
  has_fdr_sig <- nrow(fdr_row) > 0 && (fdr_row$Positive_fdr + fdr_row$Negative_fdr) > 0
  
  if (has_fdr_sig) {
    vp_fdr <- make_volcano(df, p_col = pval_col, fc = logfc_col,
                           sig_type = "fdr", fdr_col = fdr_col,
                           cell_type = celltype,
                           test_name = test_name,
                           formula_text = formula_text,
                           cohort_text = cohort,
                           positive_text = paste0("Positive with ", test_group),
                           negative_text = paste0("Negative with ", test_group))
    if (!is.null(vp_fdr)) {
      ggsave(file.path(results_dir, "Results/Figures/Volcano Plots/fdr",
                       paste0(analysis_type, "_", cell_norm, "_fdr.png")),
             plot = vp_fdr, width = 7, height = 5, dpi = 300, bg = "transparent")
    }
  }
  
  # --- GSEA (styled after plot_gsea_results) ---
  gsea_plot <- run_gsea_and_plot(df, logfc_col, pval_col,
                                 cell_name = celltype,
                                 test_name = test_name,
                                 cohort_text = cohort,
                                 pathways = hallmark_sets,
                                 reference = "hallmark")
  if (!is.null(gsea_plot)) {
    ggsave(file.path(results_dir, "Results/Figures/Pathways",
                     paste0("hallmark_", analysis_type, "_", cell_norm, "_gsea.png")),
           plot = gsea_plot, width = 10, height = 7, dpi = 300, bg = "transparent")
  }
  
  # --- GSEA lollipop (styled after panther_baseline_manuscript.qmd) ---
  lollipop_plot <- make_gsea_lollipop(df, logfc_col, pval_col,
                                      cell_name = celltype,
                                      test_name = test_name,
                                      cohort_text = cohort,
                                      pathways = hallmark_sets,
                                      reference = "hallmark")
  if (!is.null(lollipop_plot)) {
    ggsave(file.path(results_dir, "Results/Figures/Pathways Lollipop",
                     paste0("hallmark_", analysis_type, "_", cell_norm, "_lollipop.png")),
           plot = lollipop_plot, width = 10, dpi = 300, bg = "transparent")
  }
  
  n_success <- n_success + 1
}

cat("\n=== Figure generation summary ===\n")
cat(sprintf("  Total in metadata (after filtering): %d\n", nrow(meta_compiled)))
cat(sprintf("  Successfully processed: %d\n", n_success))
cat(sprintf("  Failed to read from S3: %d\n", n_error))
cat(sprintf("  Unique analysis types: %d / 94\n", n_distinct(meta_compiled$analysis_type)))
cat(sprintf("  Unique cell types: %d\n", n_distinct(meta_compiled$celltype)))
cat(sprintf("  Butterfly plots: %d (pval) + %d (fdr)\n",
            if (exists("deg_summary_df")) n_distinct(deg_summary_df$analysis_type) else 0,
            if (exists("deg_summary_df")) n_distinct(deg_summary_df$analysis_type) else 0))
cat("\nFigure generation complete!\n")
