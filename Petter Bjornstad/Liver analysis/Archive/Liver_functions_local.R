# Define genes
# Get genes_subset of interest (TCA cycle)

#TCA Cycle genes
genes <- c(
  "ACO1", "ACO2", "IDH1", "IDH2", "IDH3A", "IDH3B", "IDH3G", "OGDH", "OGDHL",
  "SUCLA2", "SUCLG1", "SUCLG2", "SDHA", "FH", "MDH1", "MDH2",  "CS"
)
# genes = sort(genes)
# all.genes = rownames(so)

acetylcoa_genes <- c(
  "PDHA1", "PDHB", "ACADM", "ACADS", "ACADL", "BCAT1", "BCAT2", "GOT1", "GOT2", 
  "ACSS2", "ACSS1", "PDC","PC"
)

ckd_genes <- c(
  "GLS", "GLUD1", "EGLN2", "EGLN1", "EGLN3", "EPO", "PC", "PCK1", "PCK2", "FBP1", 
  "G6PC", "IDH1", "IDH2", "IDH3A", "OGDH"
)

ir_genes <- c(
  "TFAM", "ACLY", "DGAT1", "SPTLC1", "IDH1", "PRKAA2", "NFE2L2", 
  "INSR", "IRS1", "IRS2", "PDHA1", "PDHB", "ACACA", "FASN",
  "PPARGC1A", "SLC2A4", "PC", "GLUD1", "CACNA1C", 
  "IL6", "TNF", "MRPL12", "MTOR", "FOXP3"
)

sglt2_genes <- c(
  "SLC2A2", "INSR", "ACADM", "ACADS",
  "HMGCS2", "BDH1", "NOX4", "SOD1", "PPARGC1A", "NOS3"
)

oxy_phos_genes <- c(
  "NDUFS6",  "SDHB", "SDHC", "SDHD",
  "UQCRC1", "UQCRC2", "COX4I1", "COX4I2", "ATP5PF"
)

leptin_adipo_genes <- c(
  "LEPR", "LEPROT", "ADIPOR1", "ADIPOR2", "CDH13"
)

# function for de.markers
de.markers <- function(seurat_object, genes, group.by, id1, id2, celltype, extension){
  m = FindMarkers(seurat_object, features = genes,group.by = group.by,ident.1 = id1, 
                  ident.2 = id2, subset.ident = celltype, verbose = F, logfc.threshold=0.001,
                  min.pct = 0.001)
  m$p_val_adj = p.adjust(m$p_val,method = "fdr")
  m <- m %>% 
    rownames_to_column('gene') %>%
    arrange(p_val) %>%
    column_to_rownames('gene') %>%
    dplyr::select(avg_log2FC,pct.1,pct.2,p_val,p_val_adj) %>%
    filter(!is.na(p_val))
  
  genes_subset <- rownames(m)[m$p_val <= 0.05]
  
  if (length(genes_subset) > 0){
    assign(paste0("genes_subset", extension), genes_subset, envir = .GlobalEnv)
  }
  assign(paste0("m", extension), m, envir = .GlobalEnv)
  return(knitr::kable(m, digits = 3
  ))
}

GeomSplitViolin <- ggproto(
  "GeomSplitViolin", 
  GeomViolin, 
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    data <- transform(data, 
                      xminv = x - violinwidth * (x - xmin), 
                      xmaxv = x + violinwidth * (xmax - x))
    grp <- data[1,'group']
    newdata <- plyr::arrange(
      transform(data, x = if(grp%%2==1) xminv else xmaxv), 
      if(grp%%2==1) y else -y
    )
    newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
    newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
    if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
      stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
      quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
      aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
      aesthetics$alpha <- rep(1, nrow(quantiles))
      both <- cbind(quantiles, aesthetics)
      quantile_grob <- GeomPath$draw_panel(both, ...)
      ggplot2:::ggname("geom_split_violin", 
                       grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
    } else {
      ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
    }
  }
)

geom_split_violin <- function (mapping = NULL, 
                               data = NULL, 
                               stat = "ydensity", 
                               position = "identity", ..., 
                               draw_quantiles = NULL, 
                               trim = TRUE, 
                               scale = "area", 
                               na.rm = FALSE, 
                               show.legend = NA, 
                               inherit.aes = TRUE) {
  layer(data = data, 
        mapping = mapping, 
        stat = stat, 
        geom = GeomSplitViolin, 
        position = position, 
        show.legend = show.legend, 
        inherit.aes = inherit.aes, 
        params = list(trim = trim, 
                      scale = scale, 
                      draw_quantiles = draw_quantiles, 
                      na.rm = na.rm, ...)
  )
}

split.vp <- function(seurat_object, genes, filepath, color1 = "#6c9a8b", color2 = "#e8998d") {
  for (i in 1:length(genes)){
    cat("\n")
    cat("###", genes[i])
    cat("\n")  
    d = VlnPlot(seurat_object, features = genes[i], split.by = "Group", idents = "PT", split.plot = F, pt.size = 0) 
    d = d$data
    p = ggplot(d,aes(x=ident, y = !!sym(genes[i]), fill=split))+
      geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.9), 
                  size = 0.2, alpha = 0.3, show.legend = T, aes(color = split)) +
      geom_split_violin(trim = T) +
      theme_bw()+
      theme(legend.title = element_blank(),
            axis.title.x = element_blank(), 
            axis.text.x=element_blank(),
            plot.title = element_text()) +
      labs(title = (paste0(genes[i]," in  PT Cells")),
           y = "Expression") +
      scale_fill_manual(values=c(color1, color2)) + 
      scale_color_manual(values=c(color1, color2))
    print(p)
    cat("\n")
    # Save
    ggsave(filename = paste0(filepath,"Violin_",genes[i],".jpeg"),plot = p,scale = 5,
           width = 800,height = 600,units = "px")
  }
}

split.vp.combined <- function(seurat_object, genes, filepath, color1 = "#6c9a8b", color2 = "#e8998d", idents = NULL) {
  compiled_d = data.frame()
  for (i in 1:length(genes)){
    d = VlnPlot(seurat_object, features = genes[i], split.by = "Group", idents = idents, split.plot = F, pt.size = 0) 
    d = d$data
    d$genename = colnames(d)[1]
    colnames(d)[1] <- "expression"
    compiled_d = rbind(compiled_d, d)
  }
  p = 
    ggplot(compiled_d,aes(x=genename, y = expression, fill=split)) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.9), 
                size = 0.2, alpha = 0.3, show.legend = T, aes(color = split)) +    
    geom_split_violin(scale = "width", trim = T) +
    theme_bw() +
    theme(legend.title = element_blank(),
          axis.title.x = element_blank(), 
          plot.title = element_text()) +
    labs(y = "Expression") +
    scale_fill_manual(values=c(color1, color2)) + 
    scale_color_manual(values=c(color1, color2))
  print(p)
  cat("\n")
  # Save
  ggsave(filename = paste0(filepath,"Violin_combined",".jpeg"),plot = p,scale = 5,
         width = 1000,height = 600,units = "px")
  
}

# Formatted dot plot function
# colorlow = "#8ecae6", colormid = "#fcbf49", colorhigh = "#d90429"
dp.formatted <- function(seurat_object, genes, celltype, group.by, m,
                         colorlow = "#83c5be", colormid = "#f4f1bb", colorhigh = "#d90429")
{
  pt.combined <- DotPlot(seurat_object,
                         features = genes,idents = celltype, group.by = group.by,
                         scale = F, cols = "RdYlBu"
  )$data 
  
  pt.plot <- pt.combined %>% 
    ggplot(aes(x=features.plot, y = id, color = avg.exp.scaled, size = pct.exp)) + 
    geom_point() +
    theme_bw() +
    scale_color_gradient2(low = colorlow, mid = colormid, high = colorhigh, midpoint = 2,
                          guide = guide_colorbar(label.vjust = 0.8, ticks = F, draw.ulim = T, draw.llim = T),
                          limits = c(0,4)) +
    scale_size(range = c(0,4), 
               limits = c(1,80)) +
    theme(panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8, vjust = 0.5),
          legend.spacing.x = unit(.1, "cm"),
          legend.direction = "horizontal") +
    guides(size = guide_legend(label.position = "bottom",
                               title.position = "top"),
           color = guide_colorbar(label.position = "bottom",
                                  title.position = "top")) +
    labs(color = "Scaled average expression",
         size = "Expression (%) ") + 
    scale_y_discrete(limits=rev)
  
  pt.table <- m %>%
    filter(rownames(m) %in% genes) %>%
    dplyr::mutate(p_val_rounded = round(p_val, 4),
                  p_val = p_format(p_val_rounded, trailing.zero = T, accuracy = 0.001, digits = 3),
                  p_val_adj_rounded = round(p_val_adj, 4),
                  p_val_adj = p_format(p_val_adj_rounded, trailing.zero = T, accuracy = 0.001, digits = 3),
                  pct.1 = sprintf("%.3f", pct.1),
                  pct.2 = sprintf("%.3f", pct.2),
                  avg_log2FC = sprintf("%.3f", avg_log2FC)) %>% 
    dplyr::select(pct.1, pct.2, avg_log2FC, p_val, p_val_adj) 
  gg.pt.table <- ggtexttable(pt.table,
                             cols = c("T1D", "HC", "Log2FC", "p-value", "q-value"),
                             theme = ttheme("blank")) %>%
    tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 1) %>%
    tab_add_title("% Expressed in PT Cells")
  
  pt.plot_table <- ggarrange(pt.plot, NULL, gg.pt.table,
                             nrow = 1, widths = c(1,-0.1,1), common.legend = F,
                             legend = "top")
  
}

metabo_sc <- function(data, gene, transcripts, gene_name) {
  plots <- list()
  for (i in 1:length(transcripts)) {
    plot <- subset(data, apply(!is.na(data[, c(transcripts[i], gene)]), 1, all)) %>%
      ggplot(aes_string(y = gene, x = transcripts[i], color = "group")) +
      geom_smooth(method = lm, se = F, aes(color=NULL), color = "black",
                  linetype = "dashed") +
      geom_point() +
      geom_smooth(method = lm, se = F) +
      labs(color = "Group",
           y = gene_name,
           x = if(is.null(label(data[[transcripts[i]]]))) (transcripts[i]) else label(data[[transcripts[i]]])) +
      theme_bw()
    plots[[i]] <- plot
  }
  ggarrange(plotlist = plots, common.legend = T)
}

add_direction <- function(df) {
  df <- df %>%
    mutate(direction = case_when(
      (avg_log2FC > 0 & p_val < 0.05) ~ "Upregulated", 
      (avg_log2FC <= 0 & p_val < 0.05) ~ "Downregulated",
      TRUE ~ "NS"
    ))
  df$direction <- factor(df$direction, levels = c("Downregulated", "Upregulated", "NS") )
  return(df)
}

